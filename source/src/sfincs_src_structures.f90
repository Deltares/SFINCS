module sfincs_src_structures
   !
   ! Point structures that move water between two grid cells by user-specified
   ! rules rather than by momentum conservation:
   !    type 1 - pump          (fixed discharge)
   !    type 2 - culvert       (bidirectional, weir-like)
   !    type 3 - check valve   (unidirectional culvert)
   !    type 4 - controlled gate, water-level triggered
   !    type 5 - controlled gate, schedule triggered
   !
   ! These used to live in sfincs_discharges.f90 alongside the river point
   ! discharges read from src/dis/netsrcdis. They have been split out so that
   ! each module has a single responsibility.
   !
   ! Runtime handoff to the continuity module is via the cell-wise qsrc(np)
   ! array (in sfincs_data): this module accumulates qq on intake (nmindrn_in)
   ! and outfall (nmindrn_out) cells. Per-structure signed discharge is also
   ! stored in qdrain(ndrn) for his output.
   !
   ! Concurrency: qsrc updates use atomic because two structures (or a river
   ! source and a structure) can land in the same cell.
   !
   use sfincs_log
   use sfincs_error
   !
   private :: parse_action_kind, parse_rule_lhs, parse_comparator, parse_rule_rhs, parse_structure_type, to_lower, check_required
   private :: initialize_src_structures_legacy
   !
   ! ------------------------------------------------------------------
   ! Named constants for the keyword-based src structure input.
   !
   ! These are scaffolding for a future TOML/YAML reader; no runtime
   ! code consumes them yet.
   ! ------------------------------------------------------------------
   !
   ! Structure type codes
   !
   integer, parameter :: structure_pump        = 1
   integer, parameter :: structure_check_valve = 2
   integer, parameter :: structure_culvert     = 3
   integer, parameter :: structure_gate        = 4
   !
   ! Action kind codes
   !
   integer, parameter :: ACTION_OPEN     = 1
   integer, parameter :: ACTION_CLOSE    = 2
   !
   ! Rule left-hand-side kind codes
   !
   integer, parameter :: RULE_LHS_Z1     = 1
   integer, parameter :: RULE_LHS_Z2     = 2
   !
   ! Rule comparator codes
   !
   integer, parameter :: CMP_LT          = 1
   integer, parameter :: CMP_LE          = 2
   integer, parameter :: CMP_GT          = 3
   integer, parameter :: CMP_GE          = 4
   integer, parameter :: CMP_EQ          = 5
   integer, parameter :: CMP_NE          = 6
   !
   ! Rule right-hand-side kind codes
   !
   integer, parameter :: RULE_RHS_PAR1   = 1
   integer, parameter :: RULE_RHS_PAR2   = 2
   integer, parameter :: RULE_RHS_PAR3   = 3
   integer, parameter :: RULE_RHS_CONST  = 4
   !
   ! ------------------------------------------------------------------
   ! Derived types for the keyword-based src structure input.
   !
   ! Scaffolding only - not yet wired into any reader or the runtime.
   ! ------------------------------------------------------------------
   !
   type :: t_src_action
      !
      integer :: kind         ! ACTION_OPEN / ACTION_CLOSE
      real    :: value        ! payload (e.g. target state / timing), unused for now
      !
   end type t_src_action
   !
   type :: t_src_rule
      !
      integer :: lhs_kind     ! RULE_LHS_*
      integer :: comparator   ! CMP_*
      integer :: rhs_kind     ! RULE_RHS_*
      real    :: rhs_value    ! only used when rhs_kind == RULE_RHS_CONST
      !
   end type t_src_rule
   !
   type :: t_src_structure
      !
      ! Identification (populated by the TOML reader). id is required,
      ! name is a human-friendly label and optional.
      !
      character(len=:), allocatable :: id
      character(len=:), allocatable :: name
      !
      ! Structure kind (one of the structure_* codes)
      !
      integer :: structure_type
      !
      ! Geometry - single representative point (x, y), and two paired
      ! coords: src_1/src_2 (the old source/sink pair) and obs_1/obs_2.
      !
      real :: x,       y
      real :: src_1_x, src_1_y
      real :: src_2_x, src_2_y
      real :: obs_1_x, obs_1_y
      real :: obs_2_x, obs_2_y
      !
      ! State
      !
      integer :: status       ! 0/1/2/3 - meaning reserved for later
      !
      ! Parameters
      !
      ! q              - pump discharge
      ! width          - gate width
      ! sill_elevation - gate sill elevation
      ! mannings_n     - gate Manning's n
      ! zmin           - gate min water level for open
      ! zmax           - gate max water level for open
      ! t_close        - gate closing time (seconds)
      ! cd, par1, par2, par3 - generic parameters (use depends on type)
      !
      real :: q
      real :: width
      real :: sill_elevation
      real :: mannings_n
      real :: zmin
      real :: zmax
      real :: t_close
      real :: cd
      real :: par1
      real :: par2
      real :: par3
      !
      ! Actions and rules
      !
      type(t_src_action), allocatable :: actions(:)
      type(t_src_rule),   allocatable :: rules(:)
      !
   end type t_src_structure
   !
   ! ------------------------------------------------------------------
   ! Module-level storage for structures parsed from a TOML input file.
   !
   ! Populated by the dispatcher when the drn file parses as TOML.
   ! Not yet consumed by any downstream runtime code - wiring is a later
   ! step. The legacy path continues to populate the flat arrays in
   ! sfincs_data (drainage_type, drainage_params, etc.).
   ! ------------------------------------------------------------------
   !
   type(t_src_structure), allocatable :: src_structures(:)
   !

contains
   !
   subroutine initialize_src_structures()
   !
   ! Dispatcher for the src_structures / drainage input file.
   !
   ! Probes the file with toml-f. If it parses as TOML, the TOML reader
   ! populates the module-level src_structures(:) array. If toml-f rejects
   ! it, falls back to the legacy fixed-column reader, which populates the
   ! drainage_* arrays in sfincs_data.
   !
   ! If a file parses as TOML but fails semantic validation (e.g. a
   ! missing required field), that is treated as a hard error: we do NOT
   ! fall back to the legacy reader, because the file was already
   ! unambiguously TOML.
   !
   use sfincs_data
   use tomlf, only : toml_table, toml_error, toml_load
   !
   implicit none
   !
   type(toml_table), allocatable :: probe_top
   type(toml_error), allocatable :: probe_err
   integer                       :: ierr_toml
   logical                       :: ok
   !
   if (drnfile(1:4) == 'none') return
   !
   ok = check_file_exists(drnfile, 'Drainage points drn file', .true.)
   !
   ! Probe: try to parse as TOML. This is a cheap check; on success we
   ! discard the probe table and let read_toml_src_structures re-parse,
   ! which keeps the two code paths decoupled.
   !
   call toml_load(probe_top, drnfile, error=probe_err)
   !
   if (.not. allocated(probe_err)) then
      !
      ! TOML path
      !
      if (allocated(probe_top)) deallocate(probe_top)
      !
      call read_toml_src_structures(drnfile, src_structures, ierr_toml)
      !
      if (ierr_toml /= 0) then
         !
         ! File was valid TOML but failed semantic validation; do NOT
         ! fall back to legacy.
         !
         write(logstr,'(a,a,a)')' Error ! Failed to load TOML src_structures file ', trim(drnfile), ' !'
         call stop_sfincs(logstr, -1)
         !
      endif
      !
      return
      !
   else
      !
      ! Legacy path
      !
      deallocate(probe_err)
      if (allocated(probe_top)) deallocate(probe_top)
      !
      call initialize_src_structures_legacy()
      !
      return
      !
   endif
   !
   end subroutine
   !
   !
   subroutine initialize_src_structures_legacy()
   !
   ! Parse drnfile in the fixed-column legacy format and populate
   ! drainage_type/_params/_status/_distance/_fraction_open,
   ! nmindrn_in(ndrn), nmindrn_out(ndrn), and the output buffer qdrain(ndrn).
   !
   use sfincs_data
   use quadtree
   !
   implicit none
   !
   real*4, dimension(:), allocatable :: xsrc_drn, ysrc_drn
   real*4, dimension(:), allocatable :: xsnk,     ysnk
   real*4    :: dummy, xsnk_tmp, ysnk_tmp, xsrc_tmp, ysrc_tmp
   integer   :: idrn, nmq, stat, npars
   logical   :: ok
   character(len=256) :: drainage_line
   !
   ndrn = 0
   !
   if (drnfile(1:4) == 'none') return
   !
   ok = check_file_exists(drnfile, 'Drainage points drn file', .true.)
   !
   ! Count lines
   !
   open(501, file=trim(drnfile))
   !
   do while (.true.)
      !
      read(501, *, iostat=stat) dummy
      if (stat < 0) exit
      ndrn = ndrn + 1
      !
   enddo
   !
   rewind(501)
   !
   if (ndrn <= 0) then
      !
      close(501)
      return
      !
   endif
   !
   write(logstr,'(a,a,a,i0,a)')' Reading ', trim(drnfile), ' (', ndrn, ' drainage points found) ...'
   call write_log(logstr, 0)
   !
   allocate(xsrc_drn(ndrn))
   allocate(ysrc_drn(ndrn))
   allocate(xsnk(ndrn))
   allocate(ysnk(ndrn))
   !
   allocate(nmindrn_in(ndrn))
   allocate(nmindrn_out(ndrn))
   allocate(qdrain(ndrn))
   allocate(drainage_type(ndrn))
   allocate(drainage_params(ndrn, 6))
   allocate(drainage_status(ndrn))
   allocate(drainage_distance(ndrn))
   allocate(drainage_fraction_open(ndrn))
   !
   nmindrn_in             = 0
   nmindrn_out            = 0
   qdrain                 = 0.0
   drainage_params        = 0.0
   drainage_distance      = 0.0
   drainage_fraction_open = 1.0   ! initially fully open (could be refined from zmin/zmax)
   drainage_status        = 1     ! 0=closed, 1=open, 2=closing, 3=opening
   !
   do idrn = 1, ndrn
      !
      read(501, '(a)') drainage_line
      !
      ! Determine drainage type first (5th integer in the line)
      !
      read(drainage_line, *, iostat=stat) xsnk_tmp, ysnk_tmp, xsrc_tmp, ysrc_tmp, drainage_type(idrn)
      !
      npars = 0
      !
      if (drainage_type(idrn) == 1 .or. drainage_type(idrn) == 2 .or. drainage_type(idrn) == 3) then
         !
         npars = 1      ! pump, culvert, or check valve
         !
      elseif (drainage_type(idrn) == 4 .or. drainage_type(idrn) == 5) then
         !
         npars = 6      ! controlled gate (width, sill, manning, zmin/tclose, zmax/topen, closing time)
         !
      endif
      !
      if (npars == 0) then
         !
         write(logstr,'(a,i0,a)') 'Drainage type ', drainage_type(idrn), ' not recognized !'
         call stop_sfincs(logstr, -1)
         !
      endif
      !
      if (npars == 1) then
         !
         read(drainage_line, *, iostat=stat) xsnk(idrn), ysnk(idrn), xsrc_drn(idrn), ysrc_drn(idrn), &
              drainage_type(idrn), drainage_params(idrn, 1)
         !
      elseif (npars == 6) then
         !
         read(drainage_line, *, iostat=stat) xsnk(idrn), ysnk(idrn), xsrc_drn(idrn), ysrc_drn(idrn), &
              drainage_type(idrn), drainage_params(idrn, 1), drainage_params(idrn, 2), &
              drainage_params(idrn, 3), drainage_params(idrn, 4), drainage_params(idrn, 5), &
              drainage_params(idrn, 6)
         !
      endif
      !
      if (stat /= 0) then
         !
         write(logstr,'(a,i0,a,i0,a)') 'Drainage type ', drainage_type(idrn), ' requires ', npars, ' parameters !'
         call stop_sfincs(logstr, -1)
         !
      endif
      !
   enddo
   !
   close(501)
   !
   ! Map intake / outfall points to cell indices and compute centre-to-centre
   ! distance (needed by controlled-gate types 4 and 5).
   !
   do idrn = 1, ndrn
      !
      nmq = find_quadtree_cell(xsnk(idrn), ysnk(idrn))
      if (nmq > 0) nmindrn_in(idrn)  = index_sfincs_in_quadtree(nmq)
      !
      nmq = find_quadtree_cell(xsrc_drn(idrn), ysrc_drn(idrn))
      if (nmq > 0) nmindrn_out(idrn) = index_sfincs_in_quadtree(nmq)
      !
      if (nmindrn_in(idrn) > 0 .and. nmindrn_out(idrn) > 0) then
         !
         xsnk_tmp = z_xz(nmindrn_in(idrn))
         ysnk_tmp = z_yz(nmindrn_in(idrn))
         xsrc_tmp = z_xz(nmindrn_out(idrn))
         ysrc_tmp = z_yz(nmindrn_out(idrn))
         drainage_distance(idrn) = sqrt( (xsrc_tmp - xsnk_tmp)**2 + (ysrc_tmp - ysnk_tmp)**2 )
         !
      endif
      !
   enddo
   !
   deallocate(xsrc_drn)
   deallocate(ysrc_drn)
   deallocate(xsnk)
   deallocate(ysnk)
   !
   if (any(nmindrn_in == 0) .or. any(nmindrn_out == 0)) then
      !
      write(logstr,'(a)') 'Warning ! For some sink/source drainage points no matching active grid cell was found!'
      call write_log(logstr, 0)
      write(logstr,'(a)') 'Warning ! These points will be skipped, please check your input!'
      call write_log(logstr, 0)
      !
   endif
   !
   end subroutine
   !
   !
   subroutine update_src_structures(t, dt, tloop)
   !
   ! Compute discharges through each drainage structure, accumulate them
   ! into qsrc(np) (intake: -qq, outfall: +qq), and store per-structure
   ! signed discharge in qdrain(ndrn) for his output.
   !
   ! Called AFTER update_discharges, which zeros qsrc first.
   !
   ! Atomic updates on qsrc(nm) guard against two structures (or a river
   ! and a structure) writing to the same cell under parallel execution.
   !
   use sfincs_data
   !
   implicit none
   !
   real*8  :: t
   real*4  :: dt
   real    :: tloop
   !
   integer :: count0, count1, count_rate, count_max
   integer :: idrn, nmin, nmout
   real*4  :: qq, qq0
   real*4  :: dzds, frac, wdt, zsill, zmin, zmax, mng, hgate, dfrac, tcls, topen, tclose
   !
   if (ndrn <= 0) return
   !
   call system_clock(count0, count_rate, count_max)
   !
   !$acc parallel loop present( z_volume, zs, zb, qsrc, qdrain, &
   !$acc                        nmindrn_in, nmindrn_out, &
   !$acc                        drainage_type, drainage_params, &
   !$acc                        drainage_distance, drainage_status, drainage_fraction_open ) &
   !$acc              private( nmin, nmout, qq, qq0, dzds, frac, wdt, zsill, &
   !$acc                       zmin, zmax, mng, hgate, dfrac, tcls, topen, tclose )
   !$omp parallel do &
   !$omp   private( nmin, nmout, qq, qq0, dzds, frac, wdt, zsill, &
   !$omp            zmin, zmax, mng, hgate, dfrac, tcls, topen, tclose ) &
   !$omp   schedule ( static )
   do idrn = 1, ndrn
      !
      nmin  = nmindrn_in(idrn)
      nmout = nmindrn_out(idrn)
      !
      if (nmin > 0 .and. nmout > 0) then
         !
         select case(drainage_type(idrn))
            !
            case(1)
               !
               ! Pump
               !
               qq = drainage_params(idrn, 1)
               !
            case(2)
               !
               ! Culvert (bidirectional)
               !
               if (zs(nmin) > zs(nmout)) then
                  !
                  qq =  drainage_params(idrn, 1) * sqrt(zs(nmin)  - zs(nmout))
                  !
               else
                  !
                  qq = -drainage_params(idrn, 1) * sqrt(zs(nmout) - zs(nmin))
                  !
               endif
               !
            case(3)
               !
               ! Check valve (culvert, but flow only from intake to outfall)
               !
               if (zs(nmin) > zs(nmout)) then
                  !
                  qq =  drainage_params(idrn, 1) * sqrt(zs(nmin)  - zs(nmout))
                  !
               else
                  !
                  qq = -drainage_params(idrn, 1) * sqrt(zs(nmout) - zs(nmin))
                  !
               endif
               !
               qq = max(qq, 0.0)
               !
            case(4)
               !
               ! Controlled gate - opens when intake water level is between zmin and zmax.
               !
               wdt   = drainage_params(idrn, 1)                        ! width
               zsill = drainage_params(idrn, 2)                        ! sill elevation
               mng   = drainage_params(idrn, 3)                        ! Manning's n
               zmin  = drainage_params(idrn, 4)                        ! min water level for open
               zmax  = drainage_params(idrn, 5)                        ! max water level for open
               tcls  = drainage_params(idrn, 6)                        ! closing time (s)
               !
               dzds  = (zs(nmout) - zs(nmin)) / drainage_distance(idrn)
               frac  = drainage_fraction_open(idrn)
               hgate = max(max(zs(nmin), zs(nmout)) - zsill, 0.0)
               dfrac = dt / tcls
               !
               qq0 = qdrain(idrn) / (wdt * max(frac, 0.001))           ! previous discharge per unit width, ignoring fraction
               !
               if (drainage_status(idrn) == 0) then
                  !
                  if (zs(nmin) > zmin .and. zs(nmin) < zmax) drainage_status(idrn) = 3
                  !
               elseif (drainage_status(idrn) == 1) then
                  !
                  if (zs(nmin) <= zmin .or. zs(nmin) >= zmax) drainage_status(idrn) = 2
                  !
               endif
               !
               if (drainage_status(idrn) == 2) then
                  !
                  frac = frac - dfrac
                  !
                  if (frac < 0.0) then
                     !
                     frac = 0.0
                     drainage_status(idrn) = 0
                     !
                  endif
                  !
               elseif (drainage_status(idrn) == 3) then
                  !
                  frac = frac + dfrac
                  !
                  if (frac > 1.0) then
                     !
                     frac = 1.0
                     drainage_status(idrn) = 1
                     !
                  endif
                  !
               endif
               !
               drainage_fraction_open(idrn) = frac
               !
               qq = (qq0 - g * hgate * dzds * dt) / (1.0 + g * mng**2 * dt * abs(qq0) / hgate**(7.0 / 3.0))
               qq = qq * wdt * frac
               !
            case(5)
               !
               ! Controlled gate - schedule triggered (one open/close window).
               !
               wdt    = drainage_params(idrn, 1)                       ! width
               zsill  = drainage_params(idrn, 2)                       ! sill elevation
               mng    = drainage_params(idrn, 3)                       ! Manning's n
               tclose = drainage_params(idrn, 4)                       ! time wrt tref to close
               topen  = drainage_params(idrn, 5)                       ! time wrt tref to open
               tcls   = drainage_params(idrn, 6)                       ! closing time (s)
               !
               dzds  = (zs(nmout) - zs(nmin)) / drainage_distance(idrn)
               frac  = drainage_fraction_open(idrn)
               hgate = max(max(zs(nmin), zs(nmout)) - zsill, 0.0)
               dfrac = dt / tcls
               !
               qq0 = qdrain(idrn) / (wdt * max(frac, 0.001))
               !
               if (drainage_status(idrn) == 0) then
                  !
                  if (t >= topen) drainage_status(idrn) = 3
                  !
               elseif (drainage_status(idrn) == 1) then
                  !
                  if (t >= tclose .and. t < topen) drainage_status(idrn) = 2
                  !
               endif
               !
               if (drainage_status(idrn) == 2) then
                  !
                  frac = frac - dfrac
                  !
                  if (frac < 0.0) then
                     !
                     frac = 0.0
                     drainage_status(idrn) = 0
                     !
                  endif
                  !
               elseif (drainage_status(idrn) == 3) then
                  !
                  frac = frac + dfrac
                  !
                  if (frac > 1.0) then
                     !
                     frac = 1.0
                     drainage_status(idrn) = 1
                     !
                  endif
                  !
               endif
               !
               drainage_fraction_open(idrn) = frac
               !
               qq = (qq0 - g * hgate * dzds * dt) / (1.0 + g * mng**2 * dt * abs(qq0) / hgate**(7.0 / 3.0))
               qq = qq * wdt * frac
               !
         end select
         !
         ! Relaxation: blend new and previous discharge to damp oscillations.
         !
         qq = 1.0 / (structure_relax / dt) * qq + (1.0 - (1.0 / (structure_relax / dt))) * qdrain(idrn)
         !
         ! Limit discharge by available volume in the intake / outfall cell.
         !
         if (subgrid) then
            !
            if (qq > 0.0) then
               !
               qq = min(qq,  max(z_volume(nmin),  0.0) / dt)
               !
            else
               !
               qq = max(qq, -max(z_volume(nmout), 0.0) / dt)
               !
            endif
            !
         else
            !
            if (qq > 0.0) then
               !
               qq = min(qq,  max((zs(nmin)  - zb(nmin))  * cell_area(z_flags_iref(nmin)),  0.0) / dt)
               !
            else
               !
               qq = max(qq, -max((zs(nmout) - zb(nmout)) * cell_area(z_flags_iref(nmout)), 0.0) / dt)
               !
            endif
            !
         endif
         !
         qdrain(idrn) = qq
         !
         ! Accumulate into cell-wise qsrc. Atomic guards against multiple
         ! structures (or a river and a structure) in the same cell.
         !
         !$acc atomic update
         !$omp atomic
         qsrc(nmin)  = qsrc(nmin)  - qq
         !$acc atomic update
         !$omp atomic
         qsrc(nmout) = qsrc(nmout) + qq
         !
      endif
      !
   enddo
   !$omp end parallel do
   !$acc end parallel loop
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0 * (count1 - count0) / count_rate
   !
   end subroutine
   !
   !
   subroutine read_toml_src_structures(filename, structures, ierr)
   !
   ! Parse a TOML input file describing point source structures into an
   ! allocatable array of t_src_structure.
   !
   ! The TOML schema is an array of tables under the key "src_structure":
   !
   !    [[src_structure]]
   !    id   = "gate_south"          ! required, string
   !    name = "South tide gate"     ! optional, string
   !    type = "gate"                ! required, one of pump/check_valve/culvert/gate
   !    x = ... ; y = ...            ! optional single-point coord
   !    src_1_x = ... ; src_1_y = ... ; src_2_x = ... ; src_2_y = ...
   !    obs_1_x = ... ; obs_1_y = ... ; obs_2_x = ... ; obs_2_y = ...
   !    status = 0
   !    q = ...                      ! pump discharge
   !    width = ... ; sill_elevation = ... ; mannings_n = ...
   !    zmin = ... ; zmax = ... ; t_close = ...
   !    cd = ... ; par1 = ... ; par2 = ... ; par3 = ...
   !    actions = [ { kind = "open",  value = 10.0 }, ... ]
   !    rules   = [ { lhs = "z1", comparator = ">", rhs = "par1" }, ... ]
   !
   ! Per-type required keys (enforced on parse):
   !    pump        : x, y, q
   !    check_valve : src_1_x, src_1_y, src_2_x, src_2_y, par1
   !    culvert     : src_1_x, src_1_y, src_2_x, src_2_y, par1
   !    gate        : src_1_x, src_1_y, src_2_x, src_2_y,
   !                  width, sill_elevation, mannings_n, zmin, zmax, t_close
   !
   ! On success, structures is allocated to the exact number of entries
   ! (can be 0). On any I/O or parse failure, structures is left
   ! unallocated and ierr is non-zero.
   !
   ! This routine does not modify module state; it is the caller's job to
   ! decide what to do with the parsed array.
   !
   use tomlf
   !
   implicit none
   !
   character(len=*), intent(in)                    :: filename
   type(t_src_structure), allocatable, intent(out) :: structures(:)
   integer, intent(out)                            :: ierr
   !
   type(toml_table), allocatable    :: top
   type(toml_error), allocatable    :: err
   type(toml_array), pointer        :: arr_structs, arr_actions, arr_rules
   type(toml_table), pointer        :: tbl_struct, tbl_entry
   character(len=:), allocatable    :: id_str, name_str, type_str, kind_str, lhs_str, cmp_str, rhs_str
   integer                          :: n_struct, n_act, n_rule, i, j, stat
   real                             :: rval
   !
   ierr = 0
   !
   ! Parse the file. toml_load returns an allocatable table; on failure the
   ! table is not allocated and err carries the diagnostic.
   !
   call toml_load(top, filename, error=err)
   !
   if (allocated(err)) then
      !
      write(logstr,'(a,a,a,a)')' Error ! Failed to parse TOML file ', trim(filename), ': ', trim(err%message)
      call write_log(logstr, 1)
      ierr = 1
      return
      !
   endif
   !
   if (.not. allocated(top)) then
      !
      write(logstr,'(a,a)')' Error ! Could not load TOML file ', trim(filename)
      call write_log(logstr, 1)
      ierr = 1
      return
      !
   endif
   !
   ! Look for the top-level array of tables "src_structure". If it is not
   ! present at all, treat that as "zero entries" (empty but valid).
   !
   nullify(arr_structs)
   call get_value(top, 'src_structure', arr_structs, requested=.false., stat=stat)
   !
   if (.not. associated(arr_structs)) then
      !
      allocate(structures(0))
      return
      !
   endif
   !
   if (.not. is_array_of_tables(arr_structs)) then
      !
      write(logstr,'(a,a)')' Error ! Key "src_structure" must be an array of tables in ', trim(filename)
      call write_log(logstr, 1)
      ierr = 1
      return
      !
   endif
   !
   n_struct = len(arr_structs)
   allocate(structures(n_struct))
   !
   do i = 1, n_struct
      !
      nullify(tbl_struct)
      call get_value(arr_structs, i, tbl_struct, stat=stat)
      !
      if (.not. associated(tbl_struct)) then
         !
         write(logstr,'(a,i0,a)')' Error ! src_structure entry ', i, ' is not a table'
         call write_log(logstr, 1)
         call cleanup_on_error()
         return
         !
      endif
      !
      ! Required id string
      !
      if (allocated(id_str)) deallocate(id_str)
      call get_value(tbl_struct, 'id', id_str, stat=stat)
      !
      if (.not. allocated(id_str)) then
         !
         write(logstr,'(a,i0,a,a)')' Error ! Missing required "id" in src_structure entry ', i, &
              ' of ', trim(filename)
         call write_log(logstr, 1)
         ierr = 1
         call cleanup_on_error()
         return
         !
      endif
      !
      structures(i)%id = id_str
      !
      ! Optional name
      !
      if (allocated(name_str)) deallocate(name_str)
      call get_value(tbl_struct, 'name', name_str, stat=stat)
      if (allocated(name_str)) structures(i)%name = name_str
      !
      ! Required type string, mapped to structure_* code
      !
      if (allocated(type_str)) deallocate(type_str)
      call get_value(tbl_struct, 'type', type_str, stat=stat)
      !
      if (.not. allocated(type_str)) then
         !
         write(logstr,'(a,i0,a,a)')' Error ! Missing required "type" in src_structure entry ', i, &
              ' of ', trim(filename)
         call write_log(logstr, 1)
         ierr = 1
         call cleanup_on_error()
         return
         !
      endif
      !
      call parse_structure_type(type_str, structures(i)%structure_type, ierr)
      !
      if (ierr /= 0) then
         !
         write(logstr,'(a,a,a,i0)')' Error ! Unknown structure type "', trim(type_str), &
              '" in src_structure entry ', i
         call write_log(logstr, 1)
         call cleanup_on_error()
         return
         !
      endif
      !
      ! Per-type required-field validation. Checked by key presence
      ! (has_key) so that 0.0 remains a legal input value.
      !
      select case (structures(i)%structure_type)
         !
         case (structure_pump)
            !
            call check_required(tbl_struct, [ character(len=14) :: &
                 'x', 'y', 'q' ], id_str, ierr)
            !
         case (structure_check_valve)
            !
            call check_required(tbl_struct, [ character(len=14) :: &
                 'src_1_x', 'src_1_y', 'src_2_x', 'src_2_y', 'par1' ], id_str, ierr)
            !
         case (structure_culvert)
            !
            call check_required(tbl_struct, [ character(len=14) :: &
                 'src_1_x', 'src_1_y', 'src_2_x', 'src_2_y', 'par1' ], id_str, ierr)
            !
         case (structure_gate)
            !
            call check_required(tbl_struct, [ character(len=14) :: &
                 'src_1_x', 'src_1_y', 'src_2_x', 'src_2_y', &
                 'width', 'sill_elevation', 'mannings_n', &
                 'zmin', 'zmax', 't_close' ], id_str, ierr)
            !
      end select
      !
      if (ierr /= 0) then
         !
         call cleanup_on_error()
         return
         !
      endif
      !
      ! Coordinates - all default to 0.0 if missing. A structure may use only
      ! the single point (x, y), or only the paired coords.
      !
      call get_value(tbl_struct, 'x',       structures(i)%x,       0.0, stat=stat)
      call get_value(tbl_struct, 'y',       structures(i)%y,       0.0, stat=stat)
      call get_value(tbl_struct, 'src_1_x', structures(i)%src_1_x, 0.0, stat=stat)
      call get_value(tbl_struct, 'src_1_y', structures(i)%src_1_y, 0.0, stat=stat)
      call get_value(tbl_struct, 'src_2_x', structures(i)%src_2_x, 0.0, stat=stat)
      call get_value(tbl_struct, 'src_2_y', structures(i)%src_2_y, 0.0, stat=stat)
      call get_value(tbl_struct, 'obs_1_x', structures(i)%obs_1_x, 0.0, stat=stat)
      call get_value(tbl_struct, 'obs_1_y', structures(i)%obs_1_y, 0.0, stat=stat)
      call get_value(tbl_struct, 'obs_2_x', structures(i)%obs_2_x, 0.0, stat=stat)
      call get_value(tbl_struct, 'obs_2_y', structures(i)%obs_2_y, 0.0, stat=stat)
      !
      ! State
      !
      call get_value(tbl_struct, 'status',         structures(i)%status,         0, stat=stat)
      !
      ! Named physical parameters
      !
      call get_value(tbl_struct, 'q',              structures(i)%q,              0.0, stat=stat)
      call get_value(tbl_struct, 'width',          structures(i)%width,          0.0, stat=stat)
      call get_value(tbl_struct, 'sill_elevation', structures(i)%sill_elevation, 0.0, stat=stat)
      call get_value(tbl_struct, 'mannings_n',     structures(i)%mannings_n,     0.0, stat=stat)
      call get_value(tbl_struct, 'zmin',           structures(i)%zmin,           0.0, stat=stat)
      call get_value(tbl_struct, 'zmax',           structures(i)%zmax,           0.0, stat=stat)
      call get_value(tbl_struct, 't_close',        structures(i)%t_close,        0.0, stat=stat)
      !
      ! Generic parameters (kept for future use / rule rhs)
      !
      call get_value(tbl_struct, 'cd',             structures(i)%cd,             0.0, stat=stat)
      call get_value(tbl_struct, 'par1',           structures(i)%par1,           0.0, stat=stat)
      call get_value(tbl_struct, 'par2',           structures(i)%par2,           0.0, stat=stat)
      call get_value(tbl_struct, 'par3',           structures(i)%par3,           0.0, stat=stat)
      !
      ! Optional actions array
      !
      nullify(arr_actions)
      call get_value(tbl_struct, 'actions', arr_actions, requested=.false., stat=stat)
      !
      if (associated(arr_actions)) then
         !
         n_act = len(arr_actions)
         allocate(structures(i)%actions(n_act))
         !
         do j = 1, n_act
            !
            nullify(tbl_entry)
            call get_value(arr_actions, j, tbl_entry, stat=stat)
            !
            if (.not. associated(tbl_entry)) then
               !
               write(logstr,'(a,i0,a,i0,a)')' Error ! actions entry ', j, &
                    ' of src_structure ', i, ' is not a table'
               call write_log(logstr, 1)
               call cleanup_on_error()
               return
               !
            endif
            !
            if (allocated(kind_str)) deallocate(kind_str)
            call get_value(tbl_entry, 'kind', kind_str, stat=stat)
            !
            if (.not. allocated(kind_str)) then
               !
               write(logstr,'(a,i0,a,i0)')' Error ! Missing "kind" in actions entry ', j, &
                    ' of src_structure ', i
               call write_log(logstr, 1)
               call cleanup_on_error()
               return
               !
            endif
            !
            call parse_action_kind(kind_str, structures(i)%actions(j)%kind, ierr)
            !
            if (ierr /= 0) then
               !
               write(logstr,'(a,a,a,i0,a,i0)')' Error ! Unknown action kind "', trim(kind_str), &
                    '" in actions entry ', j, ' of src_structure ', i
               call write_log(logstr, 1)
               call cleanup_on_error()
               return
               !
            endif
            !
            rval = 0.0
            call get_value(tbl_entry, 'value', rval, stat=stat)
            structures(i)%actions(j)%value = rval
            !
         enddo
         !
      else
         !
         allocate(structures(i)%actions(0))
         !
      endif
      !
      ! Optional rules array
      !
      nullify(arr_rules)
      call get_value(tbl_struct, 'rules', arr_rules, requested=.false., stat=stat)
      !
      if (associated(arr_rules)) then
         !
         n_rule = len(arr_rules)
         allocate(structures(i)%rules(n_rule))
         !
         do j = 1, n_rule
            !
            nullify(tbl_entry)
            call get_value(arr_rules, j, tbl_entry, stat=stat)
            !
            if (.not. associated(tbl_entry)) then
               !
               write(logstr,'(a,i0,a,i0,a)')' Error ! rules entry ', j, &
                    ' of src_structure ', i, ' is not a table'
               call write_log(logstr, 1)
               call cleanup_on_error()
               return
               !
            endif
            !
            if (allocated(lhs_str)) deallocate(lhs_str)
            if (allocated(cmp_str)) deallocate(cmp_str)
            if (allocated(rhs_str)) deallocate(rhs_str)
            !
            call get_value(tbl_entry, 'lhs',        lhs_str, stat=stat)
            call get_value(tbl_entry, 'comparator', cmp_str, stat=stat)
            call get_value(tbl_entry, 'rhs',        rhs_str, stat=stat)
            !
            if (.not. allocated(lhs_str) .or. .not. allocated(cmp_str) .or. &
                .not. allocated(rhs_str)) then
               !
               write(logstr,'(a,i0,a,i0)')' Error ! rules entry ', j, &
                    ' needs lhs/comparator/rhs keys in src_structure ', i
               call write_log(logstr, 1)
               ierr = 1
               call cleanup_on_error()
               return
               !
            endif
            !
            call parse_rule_lhs(lhs_str, structures(i)%rules(j)%lhs_kind, ierr)
            !
            if (ierr /= 0) then
               !
               write(logstr,'(a,a,a,i0,a,i0)')' Error ! Unknown rule lhs "', trim(lhs_str), &
                    '" in rules entry ', j, ' of src_structure ', i
               call write_log(logstr, 1)
               call cleanup_on_error()
               return
               !
            endif
            !
            call parse_comparator(cmp_str, structures(i)%rules(j)%comparator, ierr)
            !
            if (ierr /= 0) then
               !
               write(logstr,'(a,a,a,i0,a,i0)')' Error ! Unknown comparator "', trim(cmp_str), &
                    '" in rules entry ', j, ' of src_structure ', i
               call write_log(logstr, 1)
               call cleanup_on_error()
               return
               !
            endif
            !
            call parse_rule_rhs(rhs_str, structures(i)%rules(j)%rhs_kind, ierr)
            !
            if (ierr /= 0) then
               !
               write(logstr,'(a,a,a,i0,a,i0)')' Error ! Unknown rule rhs "', trim(rhs_str), &
                    '" in rules entry ', j, ' of src_structure ', i
               call write_log(logstr, 1)
               call cleanup_on_error()
               return
               !
            endif
            !
            rval = 0.0
            !
            if (structures(i)%rules(j)%rhs_kind == RULE_RHS_CONST) then
               !
               call get_value(tbl_entry, 'rhs_value', rval, stat=stat)
               !
            endif
            !
            structures(i)%rules(j)%rhs_value = rval
            !
         enddo
         !
      else
         !
         allocate(structures(i)%rules(0))
         !
      endif
      !
   enddo
   !
   contains
      !
      subroutine cleanup_on_error()
      !
      if (allocated(structures)) deallocate(structures)
      !
      end subroutine
      !
   end subroutine
   !
   !
   subroutine check_required(table, keys, id_str, ierr)
   !
   ! Verify that every key in "keys" is present in the TOML table. Missing
   ! keys are reported to the log (naming the structure id and key) and
   ! ierr is set non-zero. Presence is checked via has_key so that a legal
   ! value of 0.0 is not mistaken for "missing".
   !
   use tomlf
   !
   implicit none
   !
   type(toml_table), pointer,    intent(in)    :: table
   character(len=*),             intent(in)    :: keys(:)
   character(len=*),             intent(in)    :: id_str
   integer,                      intent(inout) :: ierr
   !
   integer :: k
   !
   do k = 1, size(keys)
      !
      if (.not. table%has_key(trim(keys(k)))) then
         !
         write(logstr,'(a,a,a,a,a)')' Error ! src_structure "', trim(id_str), &
              '" is missing required key "', trim(keys(k)), '"'
         call write_log(logstr, 1)
         ierr = 1
         !
      endif
      !
   enddo
   !
   end subroutine
   !
   !
   subroutine parse_action_kind(str, code, ierr)
   !
   ! Translate a TOML action "kind" string to one of the ACTION_* codes.
   !
   implicit none
   !
   character(len=*), intent(in) :: str
   integer, intent(out)         :: code
   integer, intent(out)         :: ierr
   !
   character(len=:), allocatable :: s
   !
   ierr = 0
   code = 0
   s    = to_lower(str)
   !
   select case (s)
      !
      case ('open')
         !
         code = ACTION_OPEN
         !
      case ('close')
         !
         code = ACTION_CLOSE
         !
      case default
         !
         ierr = 1
         !
   end select
   !
   end subroutine
   !
   !
   subroutine parse_structure_type(str, code, ierr)
   !
   ! Translate a TOML "type" string to one of the structure_* codes.
   !
   implicit none
   !
   character(len=*), intent(in) :: str
   integer, intent(out)         :: code
   integer, intent(out)         :: ierr
   !
   character(len=:), allocatable :: s
   !
   ierr = 0
   code = 0
   s    = to_lower(str)
   !
   select case (s)
      !
      case ('pump')
         !
         code = structure_pump
         !
      case ('check_valve')
         !
         code = structure_check_valve
         !
      case ('culvert')
         !
         code = structure_culvert
         !
      case ('gate')
         !
         code = structure_gate
         !
      case default
         !
         ierr = 1
         !
   end select
   !
   end subroutine
   !
   !
   subroutine parse_rule_lhs(str, code, ierr)
   !
   ! Translate a TOML rule "lhs" string to one of the RULE_LHS_* codes.
   !
   implicit none
   !
   character(len=*), intent(in) :: str
   integer, intent(out)         :: code
   integer, intent(out)         :: ierr
   !
   character(len=:), allocatable :: s
   !
   ierr = 0
   code = 0
   s    = to_lower(str)
   !
   select case (s)
      !
      case ('z1')
         !
         code = RULE_LHS_Z1
         !
      case ('z2')
         !
         code = RULE_LHS_Z2
         !
      case default
         !
         ierr = 1
         !
   end select
   !
   end subroutine
   !
   !
   subroutine parse_comparator(str, code, ierr)
   !
   ! Translate a TOML "comparator" string to one of the CMP_* codes.
   !
   implicit none
   !
   character(len=*), intent(in) :: str
   integer, intent(out)         :: code
   integer, intent(out)         :: ierr
   !
   ierr = 0
   code = 0
   !
   select case (trim(str))
      !
      case ('<')
         !
         code = CMP_LT
         !
      case ('<=')
         !
         code = CMP_LE
         !
      case ('>')
         !
         code = CMP_GT
         !
      case ('>=')
         !
         code = CMP_GE
         !
      case ('==')
         !
         code = CMP_EQ
         !
      case ('!=')
         !
         code = CMP_NE
         !
      case default
         !
         ierr = 1
         !
   end select
   !
   end subroutine
   !
   !
   subroutine parse_rule_rhs(str, code, ierr)
   !
   ! Translate a TOML rule "rhs" string to one of the RULE_RHS_* codes.
   !
   implicit none
   !
   character(len=*), intent(in) :: str
   integer, intent(out)         :: code
   integer, intent(out)         :: ierr
   !
   character(len=:), allocatable :: s
   !
   ierr = 0
   code = 0
   s    = to_lower(str)
   !
   select case (s)
      !
      case ('par1')
         !
         code = RULE_RHS_PAR1
         !
      case ('par2')
         !
         code = RULE_RHS_PAR2
         !
      case ('par3')
         !
         code = RULE_RHS_PAR3
         !
      case ('const')
         !
         code = RULE_RHS_CONST
         !
      case default
         !
         ierr = 1
         !
   end select
   !
   end subroutine
   !
   !
   function to_lower(str) result(lower)
   !
   ! Return a lowercase copy of str (ASCII only).
   !
   implicit none
   !
   character(len=*), intent(in)  :: str
   character(len=:), allocatable :: lower
   !
   integer :: k, ic
   !
   lower = str
   !
   do k = 1, len(lower)
      !
      ic = iachar(lower(k:k))
      !
      if (ic >= iachar('A') .and. ic <= iachar('Z')) then
         !
         lower(k:k) = achar(ic + 32)
         !
      endif
      !
   enddo
   !
   end function

end module
