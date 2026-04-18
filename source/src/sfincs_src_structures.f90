module sfincs_src_structures
   !
   ! Point structures that move water between two grid cells by user-specified
   ! rules rather than by momentum conservation:
   !    type 1 - pump          (fixed discharge)
   !    type 2 - culvert       (bidirectional)
   !    type 3 - check valve   (unidirectional culvert)
   !    type 4 - controlled gate, water-level triggered
   !    type 5 - controlled gate, schedule triggered
   !
   ! These used to live in sfincs_discharges.f90 alongside the river point
   ! discharges read from src/dis/netsrcdis. They have been split out so that
   ! each module has a single responsibility.
   !
   ! Runtime handoff to the continuity module is via the cell-wise qsrc(np)
   ! array (in sfincs_data): this module accumulates qq on intake (struc_nm_in)
   ! and outfall (struc_nm_out) cells. Per-structure signed discharge is also
   ! stored in qstruc(nr_src_structures) for his output.
   !
   ! Concurrency: qsrc updates use atomic because two structures (or a river
   ! source and a structure) can land in the same cell.
   !
   use sfincs_log
   use sfincs_error
   use sfincs_rule_expression, only: add_rule, evaluate_rule, finalize_rule_storage
   !
   private :: parse_structure_type, to_lower, check_required
   private :: convert_legacy_to_toml
   private :: allocate_struc_flat_arrays, finalize_src_structures_state
   private :: marshal_src_structures_to_flat_arrays
   private :: initialize_gate_status, write_src_structures_log_summary
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
   ! ------------------------------------------------------------------
   ! Derived type for the TOML-based src structure input.
   !
   ! Gate open/close triggers are described by small boolean expressions
   ! in strings (e.g. "(z1<0.5 | z2-z1>0.05) & z2<1.5"). Those strings
   ! live here as raw characters on the derived type; the parser runs
   ! during marshalling and emits bytecode into the shared rule_*
   ! streams owned by the sfincs_rule_expression module.
   ! ------------------------------------------------------------------
   !
   type :: t_src_structure
      !
      ! Identification (populated by the TOML reader). name is the sole
      ! identifier and is required for every structure type.
      !
      character(len=:), allocatable :: name
      !
      ! Structure kind (one of the structure_* codes)
      !
      integer :: structure_type
      !
      ! Geometry - src_1/src_2 define the intake/outfall cell pair;
      ! obs_1/obs_2 are optional and default to src_1/src_2 in the
      ! marshal when the TOML reader did not see the keys (tracked via
      ! has_obs_1 / has_obs_2).
      !
      real :: src_1_x, src_1_y
      real :: src_2_x, src_2_y
      real :: obs_1_x, obs_1_y
      real :: obs_2_x, obs_2_y
      logical :: has_obs_1
      logical :: has_obs_2
      !
      ! Parameters
      !
      ! q                 - pump discharge
      ! qmax              - maximum discharge magnitude (safety clamp)
      ! width             - gate width
      ! sill_elevation    - gate sill elevation
      ! mannings_n        - gate Manning's n
      ! opening_duration  - time (s) to go from closed to fully open
      ! closing_duration  - time (s) to go from open to fully closed
      ! flow_coef         - culvert / check_valve flow coefficient
      !
      real :: q
      real :: qmax
      real :: width
      real :: sill_elevation
      real :: mannings_n
      real :: opening_duration
      real :: closing_duration
      real :: flow_coef
      !
      ! Gate control rule expressions (raw strings; parsed by marshal).
      ! Either or both may be unallocated, meaning "no trigger for this action".
      !
      character(len=:), allocatable :: rule_open
      character(len=:), allocatable :: rule_close
      !
   end type t_src_structure
   !
   ! ------------------------------------------------------------------
   ! Module-level storage for structures parsed from a TOML input file.
   !
   ! Populated by the dispatcher when the drn file parses as TOML.
   ! Not yet consumed by any downstream runtime code - wiring is a later
   ! step. The legacy path continues to populate the flat arrays below
   ! directly (struc_type, struc_q, struc_flow_coef, etc.).
   ! ------------------------------------------------------------------
   !
   type(t_src_structure), allocatable :: src_structures(:)   ! intermediate derived-type array; flattened + deallocated by marshal_src_structures_to_flat_arrays on the toml path (gpu deep-copy avoidance).
   !
   ! ------------------------------------------------------------------
   ! Module-level runtime state for src structures (moved from sfincs_data).
   ! Populated by the legacy reader or by marshal_src_structures_to_flat_arrays
   ! from the TOML path; consumed by update_src_structures and the his output.
   ! Public so downstream modules (sfincs_openacc, sfincs_output, sfincs_ncoutput,
   ! sfincs_lib) can reference them.
   ! ------------------------------------------------------------------
   !
   ! Meta / name
   !
   integer, parameter :: struc_name_len = 128          ! max length of struct name strings
   character(len=struc_name_len), dimension(:), allocatable, public :: struc_name
   !
   ! Kind / state
   !
   integer*1, dimension(:), allocatable, public :: struc_type
   integer*1, dimension(:), allocatable, public :: struc_status
   real*4,    dimension(:), allocatable, public :: struc_distance
   real*4,    dimension(:), allocatable, public :: struc_fraction_open
   !
   ! Cell mapping
   !
   integer, public :: nr_src_structures
   integer*4, dimension(:), allocatable, public :: struc_nm_in     ! (nr_src_structures) intake  (sink)   cell indices
   integer*4, dimension(:), allocatable, public :: struc_nm_out    ! (nr_src_structures) outfall (source) cell indices
   integer*4, dimension(:), allocatable, public :: struc_nm_obs_1  ! (nr_src_structures) obs_1 cell indices (gate rule inputs; defaults to src_1 cell)
   integer*4, dimension(:), allocatable, public :: struc_nm_obs_2  ! (nr_src_structures) obs_2 cell indices (gate rule inputs; defaults to src_2 cell)
   !
   ! Gate transition timer (simulation time at which current status was entered).
   ! Only meaningful for structure_gate; ignored for other types.
   !
   real*4,    dimension(:), allocatable, public :: struc_t_state
   !
   ! Coordinates
   !
   real*4, dimension(:), allocatable, public :: struc_src_1_x, struc_src_1_y
   real*4, dimension(:), allocatable, public :: struc_src_2_x, struc_src_2_y
   real*4, dimension(:), allocatable, public :: struc_obs_1_x, struc_obs_1_y
   real*4, dimension(:), allocatable, public :: struc_obs_2_x, struc_obs_2_y
   !
   ! Named parameters
   !
   real*4, dimension(:), allocatable, public :: struc_q                  ! pump discharge
   real*4, dimension(:), allocatable, public :: struc_qmax               ! max discharge magnitude (safety clamp)
   real*4, dimension(:), allocatable, public :: struc_flow_coef          ! culvert / check_valve flow coefficient
   real*4, dimension(:), allocatable, public :: struc_width              ! gate width
   real*4, dimension(:), allocatable, public :: struc_sill_elevation     ! gate sill elevation
   real*4, dimension(:), allocatable, public :: struc_mannings_n         ! gate Manning's n
   real*4, dimension(:), allocatable, public :: struc_opening_duration   ! gate opening duration (s)
   real*4, dimension(:), allocatable, public :: struc_closing_duration   ! gate closing duration (s)
   !
   ! Runtime state
   !
   real*4, dimension(:), allocatable, public :: qstruc                   ! (nr_src_structures) signed discharge per structure, mirrors the qsrc pattern
   !
   ! ------------------------------------------------------------------
   ! Per-structure rule ids into the registry owned by sfincs_rule_expression.
   ! A rule_id of 0 means "no rule; never fires".
   !
   ! struc_rule_open_src / struc_rule_close_src hold the raw source strings
   ! (for log emission only); these do not need to travel to GPU.
   ! ------------------------------------------------------------------
   !
   integer, dimension(:), allocatable, public :: struc_rule_open       ! (nr_src_structures) rule_id for open action, 0 = no rule
   integer, dimension(:), allocatable, public :: struc_rule_close      ! (nr_src_structures) rule_id for close action, 0 = no rule
   !
   integer, parameter :: struc_rule_src_len = 256
   character(len=struc_rule_src_len), dimension(:), allocatable, public :: struc_rule_open_src
   character(len=struc_rule_src_len), dimension(:), allocatable, public :: struc_rule_close_src
   !
contains
   !
   subroutine initialize_src_structures()
   !
   ! Dispatcher for the src_structures / drainage input file.
   !
   ! Probes the file with toml-f. If it parses as TOML, the TOML reader
   ! populates the module-level src_structures(:) array. If toml-f rejects
   ! it, the file is assumed to be in the legacy fixed-column format and
   ! is transcribed on-the-fly into a TOML sibling file, which is then
   ! read via the same TOML path. This keeps only one parser alive in the
   ! source tree.
   !
   ! If a file parses as TOML but fails semantic validation (e.g. a
   ! missing required field), that is treated as a hard error.
   !
   use sfincs_data
   use tomlf, only : toml_table, toml_error, toml_load
   !
   implicit none
   !
   type(toml_table), allocatable :: probe_top
   type(toml_error), allocatable :: probe_err
   integer                       :: ierr_toml, ierr_conv
   logical                       :: ok, is_toml
   character(len=512)            :: toml_path
   integer                       :: n, p
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
   is_toml = .not. allocated(probe_err)
   !
   if (allocated(probe_err)) deallocate(probe_err)
   if (allocated(probe_top)) deallocate(probe_top)
   !
   if (is_toml) then
      !
      ! TOML path: read drnfile directly.
      !
      toml_path = drnfile
      !
   else
      !
      ! Legacy path: transcribe to a TOML sibling file, then fall through
      ! to the TOML reader. Derived path: if drnfile ends in ".drn"
      ! (case-insensitive), insert ".toml" before the extension; else
      ! append ".toml".
      !
      n = len_trim(drnfile)
      p = 0
      !
      if (n >= 4) then
         !
         if (to_lower(drnfile(n-3:n)) == '.drn') then
            !
            p = n - 3
            !
         endif
         !
      endif
      !
      if (p > 0) then
         !
         toml_path = drnfile(1:p-1) // '.toml' // drnfile(p:n)
         !
      else
         !
         toml_path = drnfile(1:n) // '.toml'
         !
      endif
      !
      call convert_legacy_to_toml(drnfile, trim(toml_path), ierr_conv)
      !
      if (ierr_conv /= 0) then
         !
         write(logstr,'(a,a,a)')' Error ! Failed to convert legacy drn file "', trim(drnfile), &
              '" to TOML; see preceding log entries for the reason'
         call stop_sfincs(trim(logstr), -1)
         !
      endif
      !
   endif
   !
   call read_toml_src_structures(trim(toml_path), src_structures, ierr_toml)
   !
   if (ierr_toml /= 0) then
      !
      write(logstr,'(a,a,a)')' Error ! Failed to load TOML src_structures file ', trim(toml_path), ' !'
      call stop_sfincs(trim(logstr), -1)
      !
   endif
   !
   ! Flatten the parsed derived-type array into the module-level
   ! struc_* 1D arrays, then deallocate src_structures(:).
   !
   call marshal_src_structures_to_flat_arrays()
   !
   ! Dump a per-structure description block to the log file.
   !
   call write_src_structures_log_summary()
   !
   ! Seed gate status + fraction_open from the initial water-level field.
   ! zs(:) has already been populated by initialize_domain -> initialize_hydro
   ! -> set_initial_conditions by the time we get here, so obs-point lookups
   ! against zs are valid. Emitted after the summary so the per-gate init
   ! status lines trail the structure block they annotate.
   !
   call initialize_gate_status()
   !
   end subroutine
   !
   !
   subroutine allocate_struc_flat_arrays(n)
   !
   ! Allocate every struc_* flat array to size n and initialise defaults.
   ! Used by both the legacy reader and the TOML marshal helper.
   ! Defensively deallocates first so re-entry is safe.
   !
   use sfincs_data
   !
   implicit none
   !
   integer, intent(in) :: n
   !
   if (allocated(struc_nm_in)) deallocate(struc_nm_in)
   if (allocated(struc_nm_out)) deallocate(struc_nm_out)
   if (allocated(struc_nm_obs_1)) deallocate(struc_nm_obs_1)
   if (allocated(struc_nm_obs_2)) deallocate(struc_nm_obs_2)
   if (allocated(qstruc)) deallocate(qstruc)
   if (allocated(struc_type)) deallocate(struc_type)
   if (allocated(struc_distance)) deallocate(struc_distance)
   if (allocated(struc_status)) deallocate(struc_status)
   if (allocated(struc_fraction_open)) deallocate(struc_fraction_open)
   if (allocated(struc_t_state)) deallocate(struc_t_state)
   if (allocated(struc_name)) deallocate(struc_name)
   if (allocated(struc_src_1_x)) deallocate(struc_src_1_x)
   if (allocated(struc_src_1_y)) deallocate(struc_src_1_y)
   if (allocated(struc_src_2_x)) deallocate(struc_src_2_x)
   if (allocated(struc_src_2_y)) deallocate(struc_src_2_y)
   if (allocated(struc_obs_1_x)) deallocate(struc_obs_1_x)
   if (allocated(struc_obs_1_y)) deallocate(struc_obs_1_y)
   if (allocated(struc_obs_2_x)) deallocate(struc_obs_2_x)
   if (allocated(struc_obs_2_y)) deallocate(struc_obs_2_y)
   if (allocated(struc_q)) deallocate(struc_q)
   if (allocated(struc_qmax)) deallocate(struc_qmax)
   if (allocated(struc_flow_coef)) deallocate(struc_flow_coef)
   if (allocated(struc_width)) deallocate(struc_width)
   if (allocated(struc_sill_elevation)) deallocate(struc_sill_elevation)
   if (allocated(struc_mannings_n)) deallocate(struc_mannings_n)
   if (allocated(struc_opening_duration)) deallocate(struc_opening_duration)
   if (allocated(struc_closing_duration)) deallocate(struc_closing_duration)
   if (allocated(struc_rule_open))       deallocate(struc_rule_open)
   if (allocated(struc_rule_close))      deallocate(struc_rule_close)
   if (allocated(struc_rule_open_src))   deallocate(struc_rule_open_src)
   if (allocated(struc_rule_close_src))  deallocate(struc_rule_close_src)
   !
   allocate(struc_nm_in(n))
   allocate(struc_nm_out(n))
   allocate(struc_nm_obs_1(n))
   allocate(struc_nm_obs_2(n))
   allocate(qstruc(n))
   allocate(struc_type(n))
   allocate(struc_distance(n))
   allocate(struc_status(n))
   allocate(struc_fraction_open(n))
   allocate(struc_t_state(n))
   allocate(struc_name(n))
   allocate(struc_src_1_x(n))
   allocate(struc_src_1_y(n))
   allocate(struc_src_2_x(n))
   allocate(struc_src_2_y(n))
   allocate(struc_obs_1_x(n))
   allocate(struc_obs_1_y(n))
   allocate(struc_obs_2_x(n))
   allocate(struc_obs_2_y(n))
   allocate(struc_q(n))
   allocate(struc_qmax(n))
   allocate(struc_flow_coef(n))
   allocate(struc_width(n))
   allocate(struc_sill_elevation(n))
   allocate(struc_mannings_n(n))
   allocate(struc_opening_duration(n))
   allocate(struc_closing_duration(n))
   allocate(struc_rule_open(n))
   allocate(struc_rule_close(n))
   allocate(struc_rule_open_src(n))
   allocate(struc_rule_close_src(n))
   !
   struc_rule_open      = 0
   struc_rule_close     = 0
   struc_rule_open_src  = ' '
   struc_rule_close_src = ' '
   !
   struc_nm_in          = 0
   struc_nm_out         = 0
   struc_nm_obs_1       = 0
   struc_nm_obs_2       = 0
   qstruc               = 0.0
   struc_type           = 0
   struc_distance       = 0.0
   struc_fraction_open  = 1.0   ! 1.0 => no-op multiplier for non-gate types; gates get their real value from initialize_gate_status
   struc_status         = 0     ! 0=closed, 1=open, 2=opening, 3=closing
   struc_t_state        = 0.0
   struc_name           = ' '
   struc_src_1_x        = 0.0
   struc_src_1_y        = 0.0
   struc_src_2_x        = 0.0
   struc_src_2_y        = 0.0
   struc_obs_1_x        = 0.0
   struc_obs_1_y        = 0.0
   struc_obs_2_x        = 0.0
   struc_obs_2_y        = 0.0
   struc_q              = 0.0
   struc_qmax           = 1.0e30
   struc_flow_coef      = 1.0
   struc_width          = 0.0
   struc_sill_elevation = 0.0
   struc_mannings_n     = 0.024
   struc_opening_duration = 600.0
   struc_closing_duration = 600.0
   !
   end subroutine
   !
   !
   subroutine finalize_src_structures_state()
   !
   ! Shared post-processing for both the legacy and TOML paths. Looks up
   ! intake / outfall cell indices from the struc_src_1_* / _2_* coords
   ! and computes centre-to-centre distance.
   !
   use sfincs_data
   use quadtree
   !
   implicit none
   !
   integer :: istruc, nmq
   real*4  :: xsnk_tmp, ysnk_tmp, xsrc_tmp, ysrc_tmp
   !
   do istruc = 1, nr_src_structures
      !
      nmq = find_quadtree_cell(struc_src_1_x(istruc), struc_src_1_y(istruc))
      if (nmq > 0) struc_nm_in(istruc)  = index_sfincs_in_quadtree(nmq)
      !
      nmq = find_quadtree_cell(struc_src_2_x(istruc), struc_src_2_y(istruc))
      if (nmq > 0) struc_nm_out(istruc) = index_sfincs_in_quadtree(nmq)
      !
      ! obs cell indices feed the gate rule evaluator. The marshal has
      ! already defaulted obs_*_x/y to src_*_x/y when the TOML reader
      ! did not see the keys, so this lookup gives us obs_1 == src_1
      ! and obs_2 == src_2 for those cases without extra branching.
      !
      nmq = find_quadtree_cell(struc_obs_1_x(istruc), struc_obs_1_y(istruc))
      if (nmq > 0) struc_nm_obs_1(istruc) = index_sfincs_in_quadtree(nmq)
      !
      nmq = find_quadtree_cell(struc_obs_2_x(istruc), struc_obs_2_y(istruc))
      if (nmq > 0) struc_nm_obs_2(istruc) = index_sfincs_in_quadtree(nmq)
      !
      if (struc_nm_in(istruc) > 0 .and. struc_nm_out(istruc) > 0) then
         !
         xsnk_tmp = z_xz(struc_nm_in(istruc))
         ysnk_tmp = z_yz(struc_nm_in(istruc))
         xsrc_tmp = z_xz(struc_nm_out(istruc))
         ysrc_tmp = z_yz(struc_nm_out(istruc))
         struc_distance(istruc) = sqrt( (xsrc_tmp - xsnk_tmp)**2 + (ysrc_tmp - ysnk_tmp)**2 )
         !
      endif
      !
   enddo
   !
   if (any(struc_nm_in == 0) .or. any(struc_nm_out == 0)) then
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
   subroutine marshal_src_structures_to_flat_arrays()
   !
   ! Copy the module-level src_structures(:) array (populated by
   ! read_toml_src_structures) into the struc_* flat arrays, then run
   ! the shared post-processing and deallocate src_structures(:).
   !
   ! The TOML and legacy paths are mutually exclusive, so the flat arrays
   ! should not yet be allocated when this is called; allocate_struc_flat_arrays
   ! defensively deallocates any residual allocation first.
   !
   ! Rule expressions (rule_open / rule_close) are handed to
   ! sfincs_rule_expression's add_rule, which appends bytecode to its
   ! shared stream, registers a new rule, and returns an integer rule_id
   ! per structure. finalize_rule_storage is called at the end to shrink
   ! the stream and registry to fit (and to allocate zero-length arrays
   ! when no rules were seen).
   !
   ! ------------------------------------------------------------------
   ! why does this marshal exist?
   !
   ! the runtime reads all src-structure state from flat per-struct
   ! arrays (the struc_* family: struc_type, struc_q, struc_flow_coef, ...).
   ! the toml reader, however, naturally produces a derived-type array
   ! src_structures(:) of t_src_structure, which carries allocatable
   ! components: character(len=:), allocatable :: name, plus the
   ! nested actions(:) and rules(:) arrays.
   !
   ! nvfortran's openacc implicit deep-copy of derived types that
   ! contain allocatable components has been unreliable in practice:
   ! pushing a type(...), allocatable :: arr(:) with nested allocatables
   ! to the device tends to produce runtime issues. flat arrays of
   ! primitive types (real, integer, fixed-length character) copy
   ! cleanly across !$acc enter data copyin(...), so we keep the live
   ! runtime state in those.
   !
   ! this routine is the one-shot bridge: toml -> src_structures(:)
   ! -> struc_* flat arrays -> deallocate(src_structures). after it
   ! runs, nothing of the derived-type array survives, so no gpu
   ! region ever sees a problematic allocatable-in-derived-type.
   !
   ! the legacy fixed-column reader populates the same struc_* flat
   ! arrays directly and therefore does not need this marshal; the
   ! two input paths converge here.
   ! ------------------------------------------------------------------
   !
   use sfincs_data
   !
   implicit none
   !
   integer :: i, n, ierr_parse
   character(len=256) :: errmsg
   !
   if (.not. allocated(src_structures)) then
      !
      nr_src_structures = 0
      return
      !
   endif
   !
   n = size(src_structures)
   nr_src_structures = n
   !
   if (n <= 0) then
      !
      deallocate(src_structures)
      return
      !
   endif
   !
   call allocate_struc_flat_arrays(n)
   !
   ! Copy scalar / vector per-structure fields and parse rule
   ! expressions into the shared rule_* stream via add_rule.
   !
   do i = 1, n
      !
      ! String fields: truncation warning if longer than struc_name_len.
      !
      if (allocated(src_structures(i)%name)) then
         !
         if (len(src_structures(i)%name) > struc_name_len) then
            !
            write(logstr,'(a,i0,a,i0,a)')' Warning ! src_structure name length > ', struc_name_len, &
                 ' at entry ', i, '; truncating'
            call write_log(logstr, 0)
            !
         endif
         !
         struc_name(i) = src_structures(i)%name
         !
      endif
      !
      struc_type(i)   = int(src_structures(i)%structure_type, 1)
      !
      ! struc_status is runtime-only (not on the TOML type); leave it at
      ! the default of 0 (open) set by allocate_struc_flat_arrays.
      !
      struc_src_1_x(i) = src_structures(i)%src_1_x
      struc_src_1_y(i) = src_structures(i)%src_1_y
      struc_src_2_x(i) = src_structures(i)%src_2_x
      struc_src_2_y(i) = src_structures(i)%src_2_y
      !
      ! obs_1 / obs_2 default to the corresponding src_* when the TOML
      ! reader did not see the key (tracked via has_obs_1 / has_obs_2).
      ! This lets 0.0 remain a legal coordinate value.
      !
      if (src_structures(i)%has_obs_1) then
         !
         struc_obs_1_x(i) = src_structures(i)%obs_1_x
         struc_obs_1_y(i) = src_structures(i)%obs_1_y
         !
      else
         !
         struc_obs_1_x(i) = src_structures(i)%src_1_x
         struc_obs_1_y(i) = src_structures(i)%src_1_y
         !
      endif
      !
      if (src_structures(i)%has_obs_2) then
         !
         struc_obs_2_x(i) = src_structures(i)%obs_2_x
         struc_obs_2_y(i) = src_structures(i)%obs_2_y
         !
      else
         !
         struc_obs_2_x(i) = src_structures(i)%src_2_x
         struc_obs_2_y(i) = src_structures(i)%src_2_y
         !
      endif
      !
      struc_q(i)                 = src_structures(i)%q
      struc_qmax(i)              = src_structures(i)%qmax
      struc_flow_coef(i)         = src_structures(i)%flow_coef
      struc_width(i)             = src_structures(i)%width
      struc_sill_elevation(i)    = src_structures(i)%sill_elevation
      struc_mannings_n(i)        = src_structures(i)%mannings_n
      struc_opening_duration(i)  = src_structures(i)%opening_duration
      struc_closing_duration(i)  = src_structures(i)%closing_duration
      !
      ! Parse rule expressions. Missing / empty strings leave the
      ! rule_id at 0, which the evaluator interprets as "never fires".
      ! Stash the source string for the init-time log summary.
      !
      if (allocated(src_structures(i)%rule_open)) then
         !
         call add_rule(src_structures(i)%rule_open, &
              struc_rule_open(i), ierr_parse, errmsg)
         !
         if (ierr_parse /= 0) then
            !
            write(logstr,'(a,a,a,a)')' Error ! src_structure "', trim(struc_name(i)), &
                 '" rules.open parse failed: ', trim(errmsg)
            call write_log(logstr, 1)
            call stop_sfincs(trim(logstr), -1)
            !
         endif
         !
         struc_rule_open_src(i) = src_structures(i)%rule_open
         !
      endif
      !
      if (allocated(src_structures(i)%rule_close)) then
         !
         call add_rule(src_structures(i)%rule_close, &
              struc_rule_close(i), ierr_parse, errmsg)
         !
         if (ierr_parse /= 0) then
            !
            write(logstr,'(a,a,a,a)')' Error ! src_structure "', trim(struc_name(i)), &
                 '" rules.close parse failed: ', trim(errmsg)
            call write_log(logstr, 1)
            call stop_sfincs(trim(logstr), -1)
            !
         endif
         !
         struc_rule_close_src(i) = src_structures(i)%rule_close
         !
      endif
      !
   enddo
   !
   ! Shrink the shared rule stream to exactly the concatenated length.
   !
   call finalize_rule_storage()
   !
   ! Shared post-processing.
   !
   call finalize_src_structures_state()
   !
   ! Drop the derived-type array; flat arrays carry all runtime state now.
   !
   deallocate(src_structures)
   !
   end subroutine
   !
   !
   subroutine update_src_structures(t, dt, tloop)
   !
   ! Compute discharges through each drainage structure, accumulate them
   ! into qsrc(np) (intake: -qq, outfall: +qq), and store per-structure
   ! signed discharge in qstruc(nr_src_structures) for his output.
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
   integer :: istruc, nmin, nmout, nm_o1, nm_o2
   real*4  :: qq, qqmax, elapsed, z1r, z2r
   real*4  :: frac, wdt, mng, zsill, dist, dzds, hgate, qq0
   logical :: open_fires, close_fires
   !
   if (nr_src_structures <= 0) return
   !
   call system_clock(count0, count_rate, count_max)
   !
   !$acc parallel loop present( z_volume, zs, zb, qsrc, qstruc, &
   !$acc                        struc_nm_in, struc_nm_out, &
   !$acc                        struc_nm_obs_1, struc_nm_obs_2, &
   !$acc                        struc_type, &
   !$acc                        struc_q, struc_qmax, struc_flow_coef, &
   !$acc                        struc_width, struc_sill_elevation, &
   !$acc                        struc_mannings_n, &
   !$acc                        struc_opening_duration, struc_closing_duration, &
   !$acc                        struc_distance, struc_status, struc_fraction_open, &
   !$acc                        struc_t_state, &
   !$acc                        struc_rule_open, struc_rule_close, &
   !$acc                        rule_opcode, rule_atom, rule_cmp, rule_threshold, &
   !$acc                        rule_start, rule_length ) &
   !$acc              private( nmin, nmout, nm_o1, nm_o2, qq, qqmax, elapsed, &
   !$acc                       z1r, z2r, frac, wdt, mng, zsill, dist, dzds, hgate, qq0, &
   !$acc                       open_fires, close_fires )
   !$omp parallel do &
   !$omp   private( nmin, nmout, nm_o1, nm_o2, qq, qqmax, elapsed, &
   !$omp            z1r, z2r, frac, wdt, mng, zsill, dist, dzds, hgate, qq0, &
   !$omp            open_fires, close_fires ) &
   !$omp   schedule ( static )
   do istruc = 1, nr_src_structures
      !
      nmin  = struc_nm_in(istruc)
      nmout = struc_nm_out(istruc)
      qqmax = struc_qmax(istruc)
      !
      if (nmin > 0 .and. nmout > 0) then
         !
         select case(struc_type(istruc))
            !
            case(structure_pump)
               !
               qq = struc_q(istruc)
               !
            case(structure_culvert)
               !
               ! Bidirectional: Q = flow_coef * sign(dh) * sqrt(|dh|)
               !
               if (zs(nmin) > zs(nmout)) then
                  !
                  qq =  struc_flow_coef(istruc) * sqrt(zs(nmin)  - zs(nmout))
                  !
               else
                  !
                  qq = -struc_flow_coef(istruc) * sqrt(zs(nmout) - zs(nmin))
                  !
               endif
               !
               qq = sign(min(abs(qq), qqmax), qq)
               !
            case(structure_check_valve)
               !
               ! One-way: flow only when z(in) > z(out); clipped to [0, qmax].
               !
               qq = struc_flow_coef(istruc) * sqrt(max(0.0, zs(nmin) - zs(nmout)))
               qq = min(qq, qqmax)
               !
            case(structure_gate)
               !
               ! Rule-driven state machine + bidirectional culvert-style
               ! flow, scaled by the momentary open fraction.
               !
               ! Status codes: 0=closed, 1=open, 2=opening, 3=closing.
               ! Opening/closing re-use the rule-evaluation branch only in
               ! the terminal (0) and (1) states; transient states (2, 3)
               ! advance purely on elapsed time so they cannot thrash.
               !
               nm_o1 = struc_nm_obs_1(istruc)
               nm_o2 = struc_nm_obs_2(istruc)
               !
               if (nm_o1 > 0) then
                  !
                  z1r = real(zs(nm_o1), 4)
                  !
               else
                  !
                  z1r = 0.0
                  !
               endif
               !
               if (nm_o2 > 0) then
                  !
                  z2r = real(zs(nm_o2), 4)
                  !
               else
                  !
                  z2r = 0.0
                  !
               endif
               !
               select case (int(struc_status(istruc)))
                  !
                  case (0)
                     !
                     ! closed - look for an open trigger
                     !
                     open_fires = evaluate_rule(struc_rule_open(istruc), z1r, z2r)
                     !
                     if (open_fires) then
                        !
                        struc_status(istruc)  = 2
                        struc_t_state(istruc) = real(t, 4)
                        !
                     endif
                     !
                  case (1)
                     !
                     ! open - look for a close trigger
                     !
                     close_fires = evaluate_rule(struc_rule_close(istruc), z1r, z2r)
                     !
                     if (close_fires) then
                        !
                        struc_status(istruc)  = 3
                        struc_t_state(istruc) = real(t, 4)
                        !
                     endif
                     !
                  case (2)
                     !
                     ! opening - advance on elapsed time; do not re-check rules
                     !
                     elapsed = real(t, 4) - struc_t_state(istruc)
                     !
                     if (struc_opening_duration(istruc) <= 0.0 .or. &
                         elapsed >= struc_opening_duration(istruc)) then
                        !
                        struc_status(istruc)        = 1
                        struc_fraction_open(istruc) = 1.0
                        !
                     else
                        !
                        struc_fraction_open(istruc) = elapsed / struc_opening_duration(istruc)
                        !
                     endif
                     !
                  case (3)
                     !
                     ! closing - advance on elapsed time; do not re-check rules
                     !
                     elapsed = real(t, 4) - struc_t_state(istruc)
                     !
                     if (struc_closing_duration(istruc) <= 0.0 .or. &
                         elapsed >= struc_closing_duration(istruc)) then
                        !
                        struc_status(istruc)        = 0
                        struc_fraction_open(istruc) = 0.0
                        !
                     else
                        !
                        struc_fraction_open(istruc) = 1.0 - elapsed / struc_closing_duration(istruc)
                        !
                     endif
                     !
               end select
               !
               ! Flow uses the src pair (nmin/nmout), not the obs pair.
               ! Bates et al. (2010) inertial formulation, per unit width:
               !    q^{n+1} = (q^n - g*h*(dzs/ds)*dt) /
               !              (1 + g*n^2*dt*|q^n| / h^{7/3})
               ! with h = max(max(zs_in, zs_out) - zsill, 0).
               ! Multiply by width * fraction_open to get the structure
               ! discharge. qstruc(istruc) holds q from the previous step
               ! in full (signed, m^3/s) discharge form, so convert via
               ! width * fraction_open to get qq0 in per-unit-width units.
               ! Sign convention: qq > 0 means flow nmin -> nmout, matching
               ! dzds = (zs_out - zs_in)/dist (positive downstream level
               ! -> negative dzds -> positive qq).
               !
               frac  = struc_fraction_open(istruc)
               wdt   = struc_width(istruc)
               mng   = struc_mannings_n(istruc)
               zsill = struc_sill_elevation(istruc)
               dist  = struc_distance(istruc)
               !
               dzds  = (real(zs(nmout), 4) - real(zs(nmin), 4)) / dist
               hgate = max(max(real(zs(nmin), 4), real(zs(nmout), 4)) - zsill, 0.0)
               !
               if (hgate > 0.0 .and. frac > 0.0) then
                  !
                  qq0 = qstruc(istruc) / (wdt * max(frac, 0.001))
                  qq  = (qq0 - g * hgate * dzds * dt) / &
                        (1.0 + g * mng * mng * dt * abs(qq0) / hgate**(7.0/3.0))
                  qq  = qq * wdt * frac
                  !
               else
                  !
                  qq = 0.0
                  !
               endif
               !
               qq = sign(min(abs(qq), qqmax), qq)
               !
         end select
         !
         ! Relaxation: blend new and previous discharge to damp oscillations.
         ! Gates use the Bates (2010) inertial form which already carries
         ! its own temporal inertia via qq0; additional blending would
         ! double-damp and suppress the dynamic response.
         !
         if (struc_type(istruc) /= structure_gate) then
            !
            qq = 1.0 / (structure_relax / dt) * qq + (1.0 - (1.0 / (structure_relax / dt))) * qstruc(istruc)
            !
         endif
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
         qstruc(istruc) = qq
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
   !    name = "south_tide_gate"     ! required, string (sole identifier)
   !    type = "gate"                ! required, one of pump/check_valve/culvert/gate
   !    src_1_x = ... ; src_1_y = ... ; src_2_x = ... ; src_2_y = ...
   !    obs_1_x = ... ; obs_1_y = ... ; obs_2_x = ... ; obs_2_y = ...
   !    q = ...                      ! pump discharge
   !    qmax = ...                   ! max discharge magnitude (safety clamp)
   !    width = ... ; sill_elevation = ... ; mannings_n = ...
   !    opening_duration = ... ; closing_duration = ...
   !    flow_coef = ...              ! culvert / check_valve flow coefficient
   !    rules.open  = "(z1<0.5 | z2-z1>0.05) & z2<1.5"   ! optional trigger expr
   !    rules.close = "z2>2.0"                           ! optional trigger expr
   !
   ! Per-type required keys (enforced on parse):
   !    pump        : name, src_1_x, src_1_y, src_2_x, src_2_y, q
   !    culvert     : name, src_1_x, src_1_y, src_2_x, src_2_y, flow_coef
   !    check_valve : name, src_1_x, src_1_y, src_2_x, src_2_y
   !    gate        : name, src_1_x, src_1_y, src_2_x, src_2_y, width, sill_elevation
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
   type(toml_array), pointer        :: arr_structs
   type(toml_table), pointer        :: tbl_struct, tbl_rules
   character(len=:), allocatable    :: name_str, type_str, rule_str
   integer                          :: n_struct, i, stat
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
      ! Required name string (presence enforced by check_required below,
      ! so that the missing-key error path flows through a single place).
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
            call check_required(tbl_struct, [ character(len=16) :: &
                 'name', 'src_1_x', 'src_1_y', 'src_2_x', 'src_2_y', 'q' ], i, ierr)
            !
         case (structure_check_valve)
            !
            call check_required(tbl_struct, [ character(len=16) :: &
                 'name', 'src_1_x', 'src_1_y', 'src_2_x', 'src_2_y' ], i, ierr)
            !
         case (structure_culvert)
            !
            call check_required(tbl_struct, [ character(len=16) :: &
                 'name', 'src_1_x', 'src_1_y', 'src_2_x', 'src_2_y', 'flow_coef' ], i, ierr)
            !
         case (structure_gate)
            !
            call check_required(tbl_struct, [ character(len=16) :: &
                 'name', 'src_1_x', 'src_1_y', 'src_2_x', 'src_2_y', &
                 'width', 'sill_elevation' ], i, ierr)
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
      ! Coordinates - src pair is required (enforced above). obs pair
      ! defaults to src in the marshal when the key is absent; track
      ! presence here so the marshal can distinguish "user gave (0,0)"
      ! from "user gave nothing".
      !
      call get_value(tbl_struct, 'src_1_x', structures(i)%src_1_x, 0.0, stat=stat)
      call get_value(tbl_struct, 'src_1_y', structures(i)%src_1_y, 0.0, stat=stat)
      call get_value(tbl_struct, 'src_2_x', structures(i)%src_2_x, 0.0, stat=stat)
      call get_value(tbl_struct, 'src_2_y', structures(i)%src_2_y, 0.0, stat=stat)
      !
      structures(i)%has_obs_1 = tbl_struct%has_key('obs_1_x') .or. tbl_struct%has_key('obs_1_y')
      structures(i)%has_obs_2 = tbl_struct%has_key('obs_2_x') .or. tbl_struct%has_key('obs_2_y')
      !
      call get_value(tbl_struct, 'obs_1_x', structures(i)%obs_1_x, 0.0, stat=stat)
      call get_value(tbl_struct, 'obs_1_y', structures(i)%obs_1_y, 0.0, stat=stat)
      call get_value(tbl_struct, 'obs_2_x', structures(i)%obs_2_x, 0.0, stat=stat)
      call get_value(tbl_struct, 'obs_2_y', structures(i)%obs_2_y, 0.0, stat=stat)
      !
      ! Named physical parameters. Defaults are picked to avoid NaN in
      ! arithmetic and to match the legacy-reader fallbacks.
      !
      call get_value(tbl_struct, 'q',                structures(i)%q,                0.0,    stat=stat)
      call get_value(tbl_struct, 'qmax',             structures(i)%qmax,             1.0e30, stat=stat)
      call get_value(tbl_struct, 'width',            structures(i)%width,            0.0,    stat=stat)
      call get_value(tbl_struct, 'sill_elevation',   structures(i)%sill_elevation,   0.0,    stat=stat)
      call get_value(tbl_struct, 'mannings_n',       structures(i)%mannings_n,       0.024,  stat=stat)
      call get_value(tbl_struct, 'opening_duration', structures(i)%opening_duration, 600.0,  stat=stat)
      call get_value(tbl_struct, 'closing_duration', structures(i)%closing_duration, 600.0,  stat=stat)
      call get_value(tbl_struct, 'flow_coef',        structures(i)%flow_coef,        1.0,    stat=stat)
      !
      ! Optional rules sub-table with string "open" / "close" expressions.
      ! Absent sub-table, or absent keys within it, leaves the rule strings
      ! unallocated on the derived type; marshal treats that as "no trigger".
      !
      nullify(tbl_rules)
      call get_value(tbl_struct, 'rules', tbl_rules, requested=.false., stat=stat)
      !
      if (associated(tbl_rules)) then
         !
         if (allocated(rule_str)) deallocate(rule_str)
         call get_value(tbl_rules, 'open', rule_str, stat=stat)
         if (allocated(rule_str)) structures(i)%rule_open = rule_str
         !
         if (allocated(rule_str)) deallocate(rule_str)
         call get_value(tbl_rules, 'close', rule_str, stat=stat)
         if (allocated(rule_str)) structures(i)%rule_close = rule_str
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
   subroutine check_required(table, keys, seq_index, ierr)
   !
   ! Verify that every key in "keys" is present in the TOML table. Missing
   ! keys are reported to the log (naming the structure by its 1-based
   ! sequence index, since "name" itself may be the missing key) and ierr
   ! is set non-zero. Presence is checked via has_key so that a legal
   ! value of 0.0 is not mistaken for "missing".
   !
   use tomlf
   !
   implicit none
   !
   type(toml_table), pointer,    intent(in)    :: table
   character(len=*),             intent(in)    :: keys(:)
   integer,                      intent(in)    :: seq_index
   integer,                      intent(inout) :: ierr
   !
   integer :: k
   !
   do k = 1, size(keys)
      !
      if (.not. table%has_key(trim(keys(k)))) then
         !
         write(logstr,'(a,i0,a,a,a)')' Error ! Structure #', seq_index, &
              ' is missing required key "', trim(keys(k)), '"'
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
   !
   !
   subroutine initialize_gate_status()
   !
   ! Seed each gate's status and fraction_open from the initial water-level
   ! field. Called right after the marshal, by which point zs(:) has already
   ! been populated by initialize_hydro -> set_initial_conditions. For
   ! non-gate structures the defaults from allocate_struc_flat_arrays
   ! (status=0=closed, fraction_open=1.0) already encode "no-op"; we only
   ! touch rows where struc_type == structure_gate.
   !
   ! Status encoding: 0=closed, 1=open, 2=opening, 3=closing.
   !
   use sfincs_data
   !
   implicit none
   !
   integer :: istruc, nm1, nm2
   real    :: z1, z2
   logical :: open_fires, close_fires
   character(len=16) :: status_str
   !
   if (nr_src_structures <= 0) return
   !
   do istruc = 1, nr_src_structures
      !
      if (struc_type(istruc) /= structure_gate) cycle
      !
      nm1 = struc_nm_obs_1(istruc)
      nm2 = struc_nm_obs_2(istruc)
      !
      if (nm1 > 0) then
         !
         z1 = real(zs(nm1), 4)
         !
      else
         !
         z1 = 0.0
         !
      endif
      !
      if (nm2 > 0) then
         !
         z2 = real(zs(nm2), 4)
         !
      else
         !
         z2 = 0.0
         !
      endif
      !
      open_fires  = evaluate_rule(struc_rule_open(istruc),  z1, z2)
      close_fires = evaluate_rule(struc_rule_close(istruc), z1, z2)
      !
      if (open_fires .and. .not. close_fires) then
         !
         struc_status(istruc)        = 1
         struc_fraction_open(istruc) = 1.0
         status_str                  = 'open'
         !
      elseif (.not. open_fires .and. close_fires) then
         !
         struc_status(istruc)        = 0
         struc_fraction_open(istruc) = 0.0
         status_str                  = 'closed'
         !
      elseif (open_fires .and. close_fires) then
         !
         struc_status(istruc)        = 1
         struc_fraction_open(istruc) = 1.0
         status_str                  = 'open'
         write(logstr,'(a,a,a)')'Warning ! gate ', trim(struc_name(istruc)), &
              ': both open and close rules fire at init; keeping gate open'
         call write_log(logstr, 0)
         !
      else
         !
         struc_status(istruc)        = 0
         struc_fraction_open(istruc) = 0.0
         status_str                  = 'closed'
         !
      endif
      !
      ! Transition timer is only consulted after a transition triggers;
      ! seed with t0 so the first rule-driven transition has a sane baseline.
      !
      struc_t_state(istruc) = t0
      !
      write(logstr,'(a,a,a,a)')'gate ', trim(struc_name(istruc)), &
           ' initialised status=', trim(status_str)
      call write_log(logstr, 0)
      !
   enddo
   !
   end subroutine
   !
   !
   subroutine write_src_structures_log_summary()
   !
   ! Emit a one-block-per-structure description of every parsed src
   ! structure to the log file. Intended for operator review at init
   ! time; not printed to stdout.
   !
   implicit none
   !
   integer :: i
   character(len=32) :: type_str
   !
   if (nr_src_structures <= 0) return
   !
   call write_log('------------------------------------------', 0)
   call write_log('Flow control structures', 0)
   call write_log('------------------------------------------', 0)
   !
   write(logstr,'(a,i0,a)')'Added ', nr_src_structures, ' flow control structures'
   call write_log(logstr, 0)
   call write_log('', 0)
   !
   do i = 1, nr_src_structures
      !
      select case (int(struc_type(i)))
         !
         case (structure_pump)
            !
            type_str = 'pump'
            !
         case (structure_culvert)
            !
            type_str = 'culvert'
            !
         case (structure_check_valve)
            !
            type_str = 'check_valve'
            !
         case (structure_gate)
            !
            type_str = 'gate'
            !
         case default
            !
            type_str = 'unknown'
            !
      end select
      !
      write(logstr,'(a,i0,a)')'Structure ', i, ':'
      call write_log(logstr, 0)
      !
      write(logstr,'(a,a)')'  name:        ', trim(struc_name(i))
      call write_log(logstr, 0)
      !
      write(logstr,'(a,a)')'  type:        ', trim(type_str)
      call write_log(logstr, 0)
      !
      write(logstr,'(a,f0.3,a,f0.3,a)')'  src_1:       (', struc_src_1_x(i), ', ', struc_src_1_y(i), ')'
      call write_log(logstr, 0)
      !
      write(logstr,'(a,f0.3,a,f0.3,a)')'  src_2:       (', struc_src_2_x(i), ', ', struc_src_2_y(i), ')'
      call write_log(logstr, 0)
      !
      ! obs coords are meaningful for culvert / check_valve / gate.
      !
      if (struc_type(i) == structure_culvert   .or. &
          struc_type(i) == structure_check_valve .or. &
          struc_type(i) == structure_gate) then
         !
         write(logstr,'(a,f0.3,a,f0.3,a)')'  obs_1:       (', struc_obs_1_x(i), ', ', struc_obs_1_y(i), ')'
         call write_log(logstr, 0)
         !
         write(logstr,'(a,f0.3,a,f0.3,a)')'  obs_2:       (', struc_obs_2_x(i), ', ', struc_obs_2_y(i), ')'
         call write_log(logstr, 0)
         !
      endif
      !
      if (struc_type(i) == structure_pump) then
         !
         write(logstr,'(a,f0.4,a)')'  discharge:   ', struc_q(i), ' (m3/s)'
         call write_log(logstr, 0)
         !
      endif
      !
      if (struc_type(i) == structure_culvert   .or. &
          struc_type(i) == structure_check_valve .or. &
          struc_type(i) == structure_gate) then
         !
         write(logstr,'(a,es12.4,a)')'  qmax:        ', struc_qmax(i), ' (m3/s)'
         call write_log(logstr, 0)
         !
      endif
      !
      if (struc_type(i) == structure_culvert .or. &
          struc_type(i) == structure_check_valve) then
         !
         write(logstr,'(a,f0.4)')'  flow_coef:   ', struc_flow_coef(i)
         call write_log(logstr, 0)
         !
      endif
      !
      if (struc_type(i) == structure_gate) then
         !
         write(logstr,'(a,f0.4,a)')'  width:       ', struc_width(i), ' (m)'
         call write_log(logstr, 0)
         !
         write(logstr,'(a,f0.4,a)')'  sill_elev.:  ', struc_sill_elevation(i), ' (m)'
         call write_log(logstr, 0)
         !
         write(logstr,'(a,f0.4)')'  mannings_n:  ', struc_mannings_n(i)
         call write_log(logstr, 0)
         !
         write(logstr,'(a,f0.2,a)')'  opening:     ', struc_opening_duration(i), ' (s)'
         call write_log(logstr, 0)
         !
         write(logstr,'(a,f0.2,a)')'  closing:     ', struc_closing_duration(i), ' (s)'
         call write_log(logstr, 0)
         !
      endif
      !
      if (struc_rule_open(i) > 0) then
         !
         if (len_trim(struc_rule_open_src(i)) > 0) then
            !
            write(logstr,'(a,a,a)')'  rules.open:  "', trim(struc_rule_open_src(i)), '"'
            !
         else
            !
            write(logstr,'(a)')'  rules.open:  (set)'
            !
         endif
         !
         call write_log(logstr, 0)
         !
      endif
      !
      if (struc_rule_close(i) > 0) then
         !
         if (len_trim(struc_rule_close_src(i)) > 0) then
            !
            write(logstr,'(a,a,a)')'  rules.close: "', trim(struc_rule_close_src(i)), '"'
            !
         else
            !
            write(logstr,'(a)')'  rules.close: (set)'
            !
         endif
         !
         call write_log(logstr, 0)
         !
      endif
      !
      call write_log('', 0)
      !
   enddo
   !
   end subroutine
   !
   !
   subroutine convert_legacy_to_toml(legacy_path, toml_path, ierr)
   !
   ! Transcribe a legacy fixed-column drn file into a TOML sibling file,
   ! so that downstream code only has to consume the TOML schema. One
   ! [[src_structure]] block is emitted per non-blank, non-comment line
   ! of the legacy file. Water-level-triggered gates (legacy dtype 4) are
   ! converted to TOML gate blocks with synthesised rule expressions.
   ! Schedule-triggered gates (legacy dtype 5) are refused; the new rule
   ! grammar is water-level-only and has no time atom.
   !
   ! The converter is deliberately minimal: no coord sanity checks, no
   ! duplicate-name detection, no preservation of comments. It exists only
   ! to remove the parallel legacy parsing machinery that used to live
   ! in this module.
   !
   implicit none
   !
   character(len=*), intent(in)  :: legacy_path
   character(len=*), intent(in)  :: toml_path
   integer,          intent(out) :: ierr
   !
   integer            :: u_in, u_out, stat, n_struct, dtype
   real*4             :: x2, y2, x1, y1, par
   real*4             :: g_width, g_sill, g_mann, g_zmin, g_zmax, g_tcls
   character(len=512) :: line, trimmed
   character(len=32)  :: name_str
   character(len=16)  :: type_name, par_name
   character(len=13)  :: zmin_str, zmax_str
   character(len=128) :: rule_open_str, rule_close_str
   !
   ierr     = 0
   n_struct = 0
   u_in     = 501
   u_out    = 502
   !
   open(u_in,  file=trim(legacy_path), status='old',     action='read',  iostat=stat)
   !
   if (stat /= 0) then
      !
      write(logstr,'(a,a,a)')' Error ! Could not open legacy drn file "', trim(legacy_path), '" for reading'
      call write_log(logstr, 1)
      ierr = 1
      return
      !
   endif
   !
   open(u_out, file=trim(toml_path),   status='replace', action='write', iostat=stat)
   !
   if (stat /= 0) then
      !
      write(logstr,'(a,a,a)')' Error ! Could not open TOML output file "', trim(toml_path), '" for writing'
      close(u_in)
      call write_log(logstr, 1)
      ierr = 1
      return
      !
   endif
   !
   write(u_out,'(a)') '# Auto-generated from legacy drn file by SFINCS.'
   write(u_out,'(a)') '# Do not edit; edit the legacy source or rewrite as TOML directly.'
   write(u_out,'(a)') ''
   !
   do while (.true.)
      !
      read(u_in,'(a)', iostat=stat) line
      !
      if (stat /= 0) exit
      !
      ! Skip blank / comment lines.
      !
      trimmed = adjustl(line)
      !
      if (len_trim(trimmed) == 0) cycle
      if (trimmed(1:1) == '#' .or. trimmed(1:1) == '!') cycle
      !
      ! Columns: x2, y2, x1, y1, dtype, par.
      ! (legacy snk -> src_2; legacy src -> src_1).
      !
      read(line, *, iostat=stat) x2, y2, x1, y1, dtype, par
      !
      if (stat /= 0) then
         !
         write(logstr,'(a,a,a)')' Error ! Could not parse legacy drn line in "', trim(legacy_path), '"'
         call write_log(logstr, 1)
         write(logstr,'(a,a)')'        line: ', trim(line)
         call write_log(logstr, 1)
         close(u_in)
         close(u_out)
         ierr = 1
         return
         !
      endif
      !
      ! Branch on dtype. Gates (4, 5) and unknown codes set ierr and bail.
      !
      select case (dtype)
         !
         case (1)
            !
            type_name = 'pump'
            par_name  = 'q'
            !
         case (2)
            !
            type_name = 'culvert'
            par_name  = 'flow_coef'
            !
         case (3)
            !
            type_name = 'check_valve'
            par_name  = 'flow_coef'
            !
         case (4)
            !
            ! Water-level-triggered gate. Legacy columns past dtype:
            !   width, sill_elevation, mannings_n, zmin, zmax, t_close.
            ! Re-read the whole line to pull those extra columns.
            !
            read(line, *, iostat=stat) x2, y2, x1, y1, dtype, &
                 g_width, g_sill, g_mann, g_zmin, g_zmax, g_tcls
            !
            if (stat /= 0) then
               !
               write(logstr,'(a,a,a)')' Error ! Could not parse legacy dtype-4 gate line in "', trim(legacy_path), '"'
               call write_log(logstr, 1)
               write(logstr,'(a,a)')'        line: ', trim(line)
               call write_log(logstr, 1)
               close(u_in)
               close(u_out)
               ierr = 1
               return
               !
            endif
            !
            ! Synthesise rule strings with the legacy numeric values baked in.
            ! Grammar accepts '<', '>', '&', '|' only (no '<=' / '>=').
            !
            write(zmin_str,'(es13.6)') g_zmin
            write(zmax_str,'(es13.6)') g_zmax
            write(rule_open_str, '(a,a,a,a)') 'z1>', trim(adjustl(zmin_str)), ' & z1<', trim(adjustl(zmax_str))
            write(rule_close_str,'(a,a,a,a)') 'z1<', trim(adjustl(zmin_str)), ' | z1>', trim(adjustl(zmax_str))
            !
            n_struct = n_struct + 1
            !
            if (g_zmin >= g_zmax) then
               !
               write(logstr,'(a,i0,a)')' Warning ! legacy gate entry ', n_struct, ': zmin >= zmax, open rule will never fire'
               call write_log(logstr, 0)
               !
            endif
            !
            write(name_str,'(a,i0)') 'legacy_', n_struct
            !
            write(u_out,'(a)')         '[[src_structure]]'
            write(u_out,'(a,a,a)')     'name             = "', trim(name_str), '"'
            write(u_out,'(a)')         'type             = "gate"'
            write(u_out,'(a,es14.6)')  'src_1_x          = ', x1
            write(u_out,'(a,es14.6)')  'src_1_y          = ', y1
            write(u_out,'(a,es14.6)')  'src_2_x          = ', x2
            write(u_out,'(a,es14.6)')  'src_2_y          = ', y2
            write(u_out,'(a,es14.6)')  'width            = ', g_width
            write(u_out,'(a,es14.6)')  'sill_elevation   = ', g_sill
            write(u_out,'(a,es14.6)')  'mannings_n       = ', g_mann
            write(u_out,'(a,es14.6)')  'opening_duration = ', g_tcls
            write(u_out,'(a,es14.6)')  'closing_duration = ', g_tcls
            write(u_out,'(a,a,a)')     'rules.open       = "', trim(rule_open_str),  '"'
            write(u_out,'(a,a,a)')     'rules.close      = "', trim(rule_close_str), '"'
            write(u_out,'(a)')         ''
            !
            cycle
            !
         case (5)
            !
            write(logstr,'(a)')' Error ! legacy schedule-triggered gate (dtype 5) not supported - rewrite as TOML with rule-based triggers or file an issue to add a time atom to the rule grammar'
            call write_log(logstr, 1)
            close(u_in)
            close(u_out)
            ierr = 1
            return
            !
         case default
            !
            write(logstr,'(a,i0,a)')' Error ! unknown drainage_type ', dtype, ' in legacy drn file'
            call write_log(logstr, 1)
            close(u_in)
            close(u_out)
            ierr = 1
            return
            !
      end select
      !
      n_struct = n_struct + 1
      write(name_str,'(a,i0)') 'legacy_', n_struct
      !
      write(u_out,'(a)')         '[[src_structure]]'
      write(u_out,'(a,a,a)')     'name    = "', trim(name_str), '"'
      write(u_out,'(a,a,a)')     'type    = "', trim(type_name),'"'
      write(u_out,'(a,es14.6)')  'src_1_x = ', x1
      write(u_out,'(a,es14.6)')  'src_1_y = ', y1
      write(u_out,'(a,es14.6)')  'src_2_x = ', x2
      write(u_out,'(a,es14.6)')  'src_2_y = ', y2
      write(u_out,'(a,a,a,es14.6)') trim(par_name), repeat(' ', max(1, 7 - len_trim(par_name))), '= ', par
      write(u_out,'(a)')         ''
      !
   enddo
   !
   close(u_in)
   close(u_out)
   !
   write(logstr,'(a,a,a,a,a)')' Converted legacy drn file "', trim(legacy_path), &
        '" to TOML "', trim(toml_path), '"'
   call write_log(logstr, 0)
   !
   end subroutine

end module
