module sfincs_src_structures
   !
   ! Point structures that move water between two grid cells by user-specified
   ! rules rather than by momentum conservation:
   !    type 1 - pump           (fixed discharge)
   !    type 2 - culvert_simple (bidirectional, optional direction filter)
   !    type 3 - culvert        (physics-based pipe flow with entrance /
   !                             friction / exit losses, bidirectional,
   !                             optional direction filter)
   !    type 4 - gate           (rule-driven state machine, bidirectional)
   !
   ! Legacy TOML alias accepted by the parser:
   !    "check_valve" -> culvert_simple + direction = "positive"
   !
   ! Orifice behaviour is not a first-class type; use type = "culvert"
   ! with submergence_ratio = 0.0 to reproduce it.
   !
   ! These used to live in sfincs_discharges.f90 alongside the river point
   ! discharges read from src/dis/netsrcdis. They have been split out so that
   ! each module has a single responsibility.
   !
   ! Runtime handoff to the continuity module is via the cell-wise qsrc(np)
   ! array (in sfincs_data): this module accumulates qq on endpoint-1
   ! (src_struc_nm_s1) and endpoint-2 (src_struc_nm_s2) cells. Per-structure
   ! signed discharge is also stored in src_struc_q_now(nr_src_structures)
   ! for his output.
   !
   ! Sign convention: a positive qq means flow from nm_s1 to nm_s2.
   ! No direction is baked into the endpoint names themselves; for pumps,
   ! endpoint 1 (nm_s1) is the intake and endpoint 2 (nm_s2) is the discharge,
   ! and the pump logic enforces qq >= 0. All other structure types are
   ! bidirectional and the sign of qq carries the flow direction.
   !
   ! Concurrency: qsrc updates use atomic because two structures (or a river
   ! source and a structure) can land in the same cell.
   !
   ! Subroutines:
   !
   !   initialize_src_structures()
   !     Main entry point. Detects legacy vs TOML, dispatches through the
   !     TOML reader, flattens into src_struc_* arrays, resolves grid-cell
   !     indices, and seeds rule-driven gate statuses from the initial zs.
   !     Called from sfincs_lib at init time.
   !
   !   update_src_structures(t, dt)
   !     Advances the open/close state machine for rule-driven structures,
   !     evaluates the per-type flux formula, and accumulates signed
   !     discharges into qsrc and src_struc_q_now. Called from update_continuity
   !     (sfincs_continuity) once per time step, after update_discharges.
   !
   !   read_toml_src_structures(filename, structures, ierr)
   !     Parse a TOML drn file into an allocatable t_src_structure(:) array.
   !     Validates required per-type keys; returns ierr /= 0 on failure.
   !     Called from initialize_src_structures (this module).
   !
   !   check_required(table, keys, seq_index, ierr)
   !     Helper for read_toml_src_structures: verifies that every key in a
   !     required-key list is present in a given TOML table. Called from
   !     read_toml_src_structures (this module).
   !
   !   parse_structure_type(str, code, ierr)
   !     Translate a TOML "type" string to one of the structure_* codes.
   !     Called from read_toml_src_structures (this module).
   !
   !   parse_direction(str, code, ierr)
   !     Translate a TOML "direction" string to one of the direction_* codes.
   !     Called from read_toml_src_structures (this module).
   !
   !   to_lower(str) result(lower)
   !     Return a lowercase copy of a string (ASCII). Called from
   !     parse_structure_type, parse_direction, and convert_legacy_to_toml
   !     (all in this module).
   !
   !   write_src_structures_log_summary()
   !     Emit a one-block-per-structure human-readable description to the
   !     log file; called from initialize_src_structures (this module) once
   !     at init time after marshalling.
   !
   !   convert_legacy_to_toml(legacy_path, toml_path, ierr)
   !     Transcribe a legacy fixed-column .drn file into a TOML sibling so
   !     the downstream code only has to consume the TOML schema. Called
   !     from initialize_src_structures (this module) when the drn file
   !     fails TOML probing.
   !
   use sfincs_log
   use sfincs_error
   use sfincs_rule_expression, only: add_rule, evaluate_rule, finalize_rule_storage
   !
   private :: parse_structure_type, parse_direction, to_lower, check_required
   private :: convert_legacy_to_toml
   private :: write_src_structures_log_summary
   !
   ! Structure type codes
   !
   integer, parameter :: structure_pump           = 1
   integer, parameter :: structure_culvert_simple = 2
   integer, parameter :: structure_culvert        = 3
   integer, parameter :: structure_gate           = 4
   !
   ! Direction filter codes (culvert_simple / culvert). Controls which sign
   ! of discharge is allowed through the structure. Default is "both".
   !
   integer, parameter :: direction_both     = 1
   integer, parameter :: direction_positive = 2
   integer, parameter :: direction_negative = 3
   !
   ! Pump reduction curve depth (m). Pump discharge is scaled by
   ! min(1, h_up/pump_reduction_depth) so the pump cannot pump the intake
   ! cell dry. Fixed constant, not user-tunable.
   !
   real*4, parameter :: pump_reduction_depth = 0.1
   !
   ! Derived type for the TOML-based src structure input.
   !
   ! Gate open/close triggers are described by small boolean expressions
   ! in strings (e.g. "(z1<0.5 | z2-z1>0.05) & z2<1.5"). Those strings
   ! live here as raw characters on the derived type; the parser runs
   ! during marshalling and emits bytecode into the shared rule_*
   ! streams owned by the sfincs_rule_expression module.
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
      ! Direction filter (one of direction_both / _positive / _negative).
      ! Honoured only for culvert_simple and culvert; other types ignore it.
      !
      integer :: direction
      !
      ! Geometry - x_s1/y_s1 and x_s2/y_s2 define the two endpoint cell
      ! coordinates; x_o1/y_o1 and x_o2/y_o2 are optional observation-point
      ! coordinates and default to the endpoint coordinates in the marshal
      ! when the TOML reader did not see the keys (tracked via has_o1 / has_o2).
      !
      real :: x_s1, y_s1
      real :: x_s2, y_s2
      real :: x_o1, y_o1
      real :: x_o2, y_o2
      logical :: has_o1
      logical :: has_o2
      !
      ! Parameters
      !
      ! q                 - pump discharge
      ! width             - gate / culvert width
      ! sill_elevation    - gate sill elevation
      ! mannings_n        - gate Manning's n
      ! opening_duration  - time (s) to go from closed to fully open
      ! closing_duration  - time (s) to go from open to fully closed
      ! flow_coef         - culvert_simple / check_valve / culvert flow coefficient
      ! height            - culvert pipe height (m, rectangular cross-section)
      ! invert_1          - culvert bed elevation at endpoint 1 (m)
      ! invert_2          - culvert bed elevation at endpoint 2 (m)
      ! submergence_ratio - culvert submergence threshold h_dn/h_up (-)
      !
      real :: q
      real :: width
      real :: sill_elevation
      real :: mannings_n
      real :: opening_duration
      real :: closing_duration
      real :: flow_coef
      !
      ! Detailed-culvert geometry + submergence threshold
      !
      real :: height
      real :: invert_1
      real :: invert_2
      real :: submergence_ratio
      !
      ! Gate control rule expressions (raw strings; parsed by marshal).
      ! Either or both may be unallocated, meaning "no trigger for this action".
      !
      character(len=:), allocatable :: rule_open
      character(len=:), allocatable :: rule_close
      !
   end type t_src_structure
   !
   ! Module-level storage for structures parsed from a TOML input file.
   ! Populated by the dispatcher and flattened into the flat arrays below
   ! by the marshal.
   !
   type(t_src_structure), allocatable :: src_structures(:)   ! intermediate derived-type array; flattened + deallocated by marshal_src_structures_to_flat_arrays on the toml path (gpu deep-copy avoidance).
   !
   ! Module-level runtime state for src structures (moved from sfincs_data).
   ! Populated by the legacy reader or by marshal_src_structures_to_flat_arrays
   ! from the TOML path; consumed by update_src_structures and the his output.
   ! Public so downstream modules (sfincs_openacc, sfincs_output, sfincs_ncoutput,
   ! sfincs_lib) can reference them.
   !
   ! Meta / name
   !
   integer, parameter :: src_struc_name_len = 128          ! max length of struct name strings
   character(len=src_struc_name_len), dimension(:), allocatable, public :: src_struc_name
   !
   ! Kind / state
   !
   integer*1, dimension(:), allocatable, public :: src_struc_type
   integer,   dimension(:), allocatable, public :: src_struc_direction   ! direction_* code; honoured by culvert_simple and culvert
   integer*1, dimension(:), allocatable, public :: src_struc_status
   real*4,    dimension(:), allocatable, public :: src_struc_distance
   real*4,    dimension(:), allocatable, public :: src_struc_fraction_open
   !
   ! Input file path (sfincs.inp keyword 'drnfile'); 'none' when no drainage
   ! structures file is supplied.
   !
   character(len=256), public :: drnfile
   !
   ! Cell mapping
   !
   integer, public :: nr_src_structures
   integer*4, dimension(:), allocatable, public :: src_struc_nm_s1     ! (nr_src_structures) endpoint-1 cell indices
   integer*4, dimension(:), allocatable, public :: src_struc_nm_s2     ! (nr_src_structures) endpoint-2 cell indices
   integer*4, dimension(:), allocatable, public :: src_struc_nm_o1     ! (nr_src_structures) obs-1 cell indices (gate rule inputs; defaults to endpoint-1 cell)
   integer*4, dimension(:), allocatable, public :: src_struc_nm_o2     ! (nr_src_structures) obs-2 cell indices (gate rule inputs; defaults to endpoint-2 cell)
   !
   ! Gate transition timer (simulation time at which current status was entered).
   ! Only meaningful for structure_gate; ignored for other types.
   !
   real*4,    dimension(:), allocatable, public :: src_struc_t_state
   !
   ! Coordinates
   !
   real*4, dimension(:), allocatable, public :: src_struc_x_s1, src_struc_y_s1
   real*4, dimension(:), allocatable, public :: src_struc_x_s2, src_struc_y_s2
   real*4, dimension(:), allocatable, public :: src_struc_x_o1, src_struc_y_o1
   real*4, dimension(:), allocatable, public :: src_struc_x_o2, src_struc_y_o2
   !
   ! Named parameters
   !
   real*4, dimension(:), allocatable, public :: src_struc_q                  ! pump discharge
   real*4, dimension(:), allocatable, public :: src_struc_flow_coef          ! culvert_simple / check_valve / culvert flow coefficient
   real*4, dimension(:), allocatable, public :: src_struc_width              ! gate / culvert width
   real*4, dimension(:), allocatable, public :: src_struc_sill_elevation     ! gate sill elevation
   real*4, dimension(:), allocatable, public :: src_struc_mannings_n         ! gate / culvert Manning's n
   real*4, dimension(:), allocatable, public :: src_struc_opening_duration   ! gate opening duration (s)
   real*4, dimension(:), allocatable, public :: src_struc_closing_duration   ! gate closing duration (s)
   !
   ! Detailed-culvert geometry
   !
   real*4, dimension(:), allocatable, public :: src_struc_height             ! culvert pipe height (m)
   real*4, dimension(:), allocatable, public :: src_struc_invert_1           ! culvert bed elevation at endpoint 1 (m)
   real*4, dimension(:), allocatable, public :: src_struc_invert_2           ! culvert bed elevation at endpoint 2 (m)
   !
   ! Detailed-culvert submergence threshold
   !
   real*4, dimension(:), allocatable, public :: src_struc_submergence_ratio  ! culvert submergence threshold h_dn/h_up (-)
   !
   ! Runtime state
   !
   real*4, dimension(:), allocatable, public :: src_struc_q_now               ! (nr_src_structures) signed discharge this step per structure, mirrors the qsrc pattern
   !
   ! Per-structure rule ids into the registry owned by sfincs_rule_expression.
   ! A rule_id of 0 means "no rule; never fires".
   !
   ! src_struc_rule_open_src / src_struc_rule_close_src hold the raw source strings
   ! (for log emission only); these do not need to travel to GPU.
   !
   integer, dimension(:), allocatable, public :: src_struc_rule_open       ! (nr_src_structures) rule_id for open action, 0 = no rule
   integer, dimension(:), allocatable, public :: src_struc_rule_close      ! (nr_src_structures) rule_id for close action, 0 = no rule
   !
   integer, parameter :: src_struc_rule_src_len = 256
   character(len=src_struc_rule_src_len), dimension(:), allocatable, public :: src_struc_rule_open_src
   character(len=src_struc_rule_src_len), dimension(:), allocatable, public :: src_struc_rule_close_src
   !
contains
   !
   !-----------------------------------------------------------------------------------------------------!
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
      ! After parsing, the derived-type src_structures(:) array is flattened
      ! into the src_struc_* 1D arrays (the runtime's sole state representation),
      ! grid-cell indices and distances are resolved, a descriptive block is
      ! written to the log, and gate statuses are seeded from the initial
      ! water-level field.
      !
      ! Called from: sfincs_lib (once, at init time).
      !
      use sfincs_data
      use quadtree
      use tomlf, only : toml_table, toml_error, toml_load
      !
      implicit none
      !
      ! Dispatcher locals
      !
      type(toml_table), allocatable :: probe_top
      type(toml_error), allocatable :: probe_err
      integer                       :: ierr_toml, ierr_conv
      logical                       :: ok, is_toml
      character(len=512)            :: toml_path
      !
      ! Marshal locals
      !
      integer                       :: i, ierr_parse
      character(len=256)            :: errmsg
      !
      ! Cell-index / distance locals
      !
      integer                       :: istruc, nmq
      real*4                        :: x_s1_tmp, y_s1_tmp, x_s2_tmp, y_s2_tmp
      !
      ! Gate-status seeding locals
      !
      integer                       :: nm_o1, nm_o2
      real                          :: zs_o1, zs_o2
      logical                       :: open_fires, close_fires
      character(len=16)             :: status_str
      !
      drainage_structures = .false.
      !
      if (drnfile(1:4) == 'none') return
      !
      ! Existence check
      !
      ok = check_file_exists(drnfile, 'Drainage points drn file', .true.)
      !
      ! Probe TOML / convert legacy / re-read TOML
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
         ! to the TOML reader. The converter derives its own output path from
         ! drnfile.
         !
         call convert_legacy_to_toml(drnfile, toml_path, ierr_conv)
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
      ! Marshal src_structures(:) -> src_struc_* flat arrays.
      !
      ! The runtime reads all src-structure state from flat per-struct
      ! arrays (the src_struc_* family: src_struc_type, src_struc_q, src_struc_flow_coef, ...).
      ! The TOML reader, however, naturally produces a derived-type array
      ! src_structures(:) of t_src_structure, which carries allocatable
      ! components: character(len=:), allocatable :: name, plus the rule
      ! expression strings.
      !
      ! nvfortran's openacc implicit deep-copy of derived types that
      ! contain allocatable components has been unreliable in practice:
      ! pushing a type(...), allocatable :: arr(:) with nested allocatables
      ! to the device tends to produce runtime issues. Flat arrays of
      ! primitive types (real, integer, fixed-length character) copy
      ! cleanly across !$acc enter data copyin(...), so we keep the live
      ! runtime state in those.
      !
      ! The marshal is the one-shot bridge: toml -> src_structures(:)
      ! -> src_struc_* flat arrays -> deallocate(src_structures). After it
      ! runs, nothing of the derived-type array survives, so no gpu
      ! region ever sees a problematic allocatable-in-derived-type.
      !
      if (.not. allocated(src_structures)) then
         !
         nr_src_structures = 0
         !
         call write_src_structures_log_summary()
         !
         return
         !
      endif
      !
      nr_src_structures = size(src_structures)
      !
      if (nr_src_structures <= 0) then
         !
         deallocate(src_structures)
         !
         call write_src_structures_log_summary()
         !
         return
         !
      endif
      !
      drainage_structures = .true.
      !
      ! Allocate flat arrays to size nr_src_structures and seed defaults.
      !
      allocate(src_struc_nm_s1(nr_src_structures))
      allocate(src_struc_nm_s2(nr_src_structures))
      allocate(src_struc_nm_o1(nr_src_structures))
      allocate(src_struc_nm_o2(nr_src_structures))
      allocate(src_struc_q_now(nr_src_structures))
      allocate(src_struc_type(nr_src_structures))
      allocate(src_struc_direction(nr_src_structures))
      allocate(src_struc_distance(nr_src_structures))
      allocate(src_struc_status(nr_src_structures))
      allocate(src_struc_fraction_open(nr_src_structures))
      allocate(src_struc_t_state(nr_src_structures))
      allocate(src_struc_name(nr_src_structures))
      allocate(src_struc_x_s1(nr_src_structures))
      allocate(src_struc_y_s1(nr_src_structures))
      allocate(src_struc_x_s2(nr_src_structures))
      allocate(src_struc_y_s2(nr_src_structures))
      allocate(src_struc_x_o1(nr_src_structures))
      allocate(src_struc_y_o1(nr_src_structures))
      allocate(src_struc_x_o2(nr_src_structures))
      allocate(src_struc_y_o2(nr_src_structures))
      allocate(src_struc_q(nr_src_structures))
      allocate(src_struc_flow_coef(nr_src_structures))
      allocate(src_struc_width(nr_src_structures))
      allocate(src_struc_sill_elevation(nr_src_structures))
      allocate(src_struc_mannings_n(nr_src_structures))
      allocate(src_struc_opening_duration(nr_src_structures))
      allocate(src_struc_closing_duration(nr_src_structures))
      allocate(src_struc_height(nr_src_structures))
      allocate(src_struc_invert_1(nr_src_structures))
      allocate(src_struc_invert_2(nr_src_structures))
      allocate(src_struc_submergence_ratio(nr_src_structures))
      allocate(src_struc_rule_open(nr_src_structures))
      allocate(src_struc_rule_close(nr_src_structures))
      allocate(src_struc_rule_open_src(nr_src_structures))
      allocate(src_struc_rule_close_src(nr_src_structures))
      !
      src_struc_rule_open      = 0
      src_struc_rule_close     = 0
      src_struc_rule_open_src  = ' '
      src_struc_rule_close_src = ' '
      !
      src_struc_nm_s1          = 0
      src_struc_nm_s2          = 0
      src_struc_nm_o1          = 0
      src_struc_nm_o2          = 0
      src_struc_q_now          = 0.0
      src_struc_type           = 0
      src_struc_direction      = direction_both
      src_struc_distance       = 0.0
      src_struc_fraction_open  = 1.0   ! default "fully open": structures without rules bypass the state machine and use this as a no-op multiplier in the common-tail scaling
      src_struc_status         = 1     ! 0=closed, 1=open, 2=opening, 3=closing; default open (see above). Rule-driven structures overwrite this in the init-time seeding below.
      src_struc_t_state        = 0.0
      src_struc_name           = ' '
      src_struc_x_s1           = 0.0
      src_struc_y_s1           = 0.0
      src_struc_x_s2           = 0.0
      src_struc_y_s2           = 0.0
      src_struc_x_o1           = 0.0
      src_struc_y_o1           = 0.0
      src_struc_x_o2           = 0.0
      src_struc_y_o2           = 0.0
      src_struc_q              = 0.0
      src_struc_flow_coef      = 1.0
      src_struc_width          = 0.0
      src_struc_sill_elevation = 0.0
      src_struc_mannings_n     = 0.024
      src_struc_opening_duration = 600.0
      src_struc_closing_duration = 600.0
      src_struc_height           = 0.0
      src_struc_invert_1         = 0.0
      src_struc_invert_2         = 0.0
      src_struc_submergence_ratio = 0.667
      !
      ! Copy scalar / coord / string / parameter fields from src_structures(:)
      ! into the flat arrays, and parse rule source strings via add_rule.
      !
      do i = 1, nr_src_structures
         !
         ! String fields: truncation warning if longer than src_struc_name_len.
         !
         if (allocated(src_structures(i)%name)) then
            !
            if (len(src_structures(i)%name) > src_struc_name_len) then
               !
               write(logstr,'(a,i0,a,i0,a)')' Warning ! src_structure name length > ', src_struc_name_len, &
                    ' at entry ', i, '; truncating'
               call write_log(logstr, 0)
               !
            endif
            !
            src_struc_name(i) = src_structures(i)%name
            !
         endif
         !
         src_struc_type(i)      = int(src_structures(i)%structure_type, 1)
         src_struc_direction(i) = src_structures(i)%direction
         !
         ! src_struc_status is runtime-only (not on the TOML type); leave it at
         ! the default of 0 (closed) set above.
         !
         src_struc_x_s1(i) = src_structures(i)%x_s1
         src_struc_y_s1(i) = src_structures(i)%y_s1
         src_struc_x_s2(i) = src_structures(i)%x_s2
         src_struc_y_s2(i) = src_structures(i)%y_s2
         !
         ! obs 1 / obs 2 default to the corresponding endpoint when the TOML
         ! reader did not see the key (tracked via has_o1 / has_o2).
         ! This lets 0.0 remain a legal coordinate value.
         !
         if (src_structures(i)%has_o1) then
            !
            src_struc_x_o1(i) = src_structures(i)%x_o1
            src_struc_y_o1(i) = src_structures(i)%y_o1
            !
         else
            !
            src_struc_x_o1(i) = src_structures(i)%x_s1
            src_struc_y_o1(i) = src_structures(i)%y_s1
            !
         endif
         !
         if (src_structures(i)%has_o2) then
            !
            src_struc_x_o2(i) = src_structures(i)%x_o2
            src_struc_y_o2(i) = src_structures(i)%y_o2
            !
         else
            !
            src_struc_x_o2(i) = src_structures(i)%x_s2
            src_struc_y_o2(i) = src_structures(i)%y_s2
            !
         endif
         !
         src_struc_q(i)                 = src_structures(i)%q
         src_struc_flow_coef(i)         = src_structures(i)%flow_coef
         src_struc_width(i)             = src_structures(i)%width
         src_struc_sill_elevation(i)    = src_structures(i)%sill_elevation
         src_struc_mannings_n(i)        = src_structures(i)%mannings_n
         src_struc_opening_duration(i)  = src_structures(i)%opening_duration
         src_struc_closing_duration(i)  = src_structures(i)%closing_duration
         src_struc_height(i)            = src_structures(i)%height
         src_struc_invert_1(i)          = src_structures(i)%invert_1
         src_struc_invert_2(i)          = src_structures(i)%invert_2
         src_struc_submergence_ratio(i) = src_structures(i)%submergence_ratio
         !
         ! Parse rule expressions. Missing / empty strings leave the
         ! rule_id at 0, which the evaluator interprets as "never fires".
         ! Stash the source string for the init-time log summary.
         !
         if (allocated(src_structures(i)%rule_open)) then
            !
            call add_rule(src_structures(i)%rule_open, &
                 src_struc_rule_open(i), ierr_parse, errmsg)
            !
            if (ierr_parse /= 0) then
               !
               write(logstr,'(a,a,a,a)')' Error ! src_structure "', trim(src_struc_name(i)), &
                    '" rules_open parse failed: ', trim(errmsg)
               call write_log(logstr, 1)
               call stop_sfincs(trim(logstr), -1)
               !
            endif
            !
            src_struc_rule_open_src(i) = src_structures(i)%rule_open
            !
         endif
         !
         if (allocated(src_structures(i)%rule_close)) then
            !
            call add_rule(src_structures(i)%rule_close, &
                 src_struc_rule_close(i), ierr_parse, errmsg)
            !
            if (ierr_parse /= 0) then
               !
               write(logstr,'(a,a,a,a)')' Error ! src_structure "', trim(src_struc_name(i)), &
                    '" rules_close parse failed: ', trim(errmsg)
               call write_log(logstr, 1)
               call stop_sfincs(trim(logstr), -1)
               !
            endif
            !
            src_struc_rule_close_src(i) = src_structures(i)%rule_close
            !
         endif
         !
      enddo
      !
      ! Shrink the shared rule bytecode stream to exactly the concatenated
      ! length (also allocates zero-length arrays when no rules were seen).
      !
      call finalize_rule_storage()
      !
      ! Drop the derived-type array; flat arrays carry all runtime state now.
      !
      deallocate(src_structures)
      !
      ! Resolve cell-index lookups (src_struc_nm_s1 / _s2 / _o1 / _o2)
      ! and centre-to-centre distance from coordinate pairs.
      !
      do istruc = 1, nr_src_structures
         !
         nmq = find_quadtree_cell(src_struc_x_s1(istruc), src_struc_y_s1(istruc))
         if (nmq > 0) src_struc_nm_s1(istruc) = index_sfincs_in_quadtree(nmq)
         !
         nmq = find_quadtree_cell(src_struc_x_s2(istruc), src_struc_y_s2(istruc))
         if (nmq > 0) src_struc_nm_s2(istruc) = index_sfincs_in_quadtree(nmq)
         !
         ! obs cell indices feed the gate rule evaluator. The marshal has
         ! already defaulted x_o1/y_o1 and x_o2/y_o2 to the endpoint
         ! coordinates when the TOML reader did not see the keys, so this
         ! lookup gives us obs-1 == endpoint-1 and obs-2 == endpoint-2 for
         ! those cases without extra branching.
         !
         nmq = find_quadtree_cell(src_struc_x_o1(istruc), src_struc_y_o1(istruc))
         if (nmq > 0) src_struc_nm_o1(istruc) = index_sfincs_in_quadtree(nmq)
         !
         nmq = find_quadtree_cell(src_struc_x_o2(istruc), src_struc_y_o2(istruc))
         if (nmq > 0) src_struc_nm_o2(istruc) = index_sfincs_in_quadtree(nmq)
         !
         if (src_struc_nm_s1(istruc) > 0 .and. src_struc_nm_s2(istruc) > 0) then
            !
            x_s1_tmp = z_xz(src_struc_nm_s1(istruc))
            y_s1_tmp = z_yz(src_struc_nm_s1(istruc))
            x_s2_tmp = z_xz(src_struc_nm_s2(istruc))
            y_s2_tmp = z_yz(src_struc_nm_s2(istruc))
            src_struc_distance(istruc) = sqrt( (x_s2_tmp - x_s1_tmp)**2 + (y_s2_tmp - y_s1_tmp)**2 )
            !
         endif
         !
      enddo
      !
      if (any(src_struc_nm_s1 == 0) .or. any(src_struc_nm_s2 == 0)) then
         !
         write(logstr,'(a)') 'Warning ! For some source-structure endpoints no matching active grid cell was found!'
         call write_log(logstr, 0)
         write(logstr,'(a)') 'Warning ! These points will be skipped, please check your input!'
         call write_log(logstr, 0)
         !
      endif
      !
      ! Write the per-structure descriptive block to the log file.
      ! Emitted before the gate-status seeding so the per-gate init status
      ! lines trail the structure block they annotate.
      !
      call write_src_structures_log_summary()
      !
      ! Initial-status seeding for rule-driven structures.
      !
      ! zs(:) has already been populated by initialize_domain -> initialize_hydro
      ! -> set_initial_conditions by the time we get here, so obs-point lookups
      ! against zs are valid. For structures with no rule expressions the defaults
      ! assigned above (status=1=open, fraction_open=1.0) already encode "no-op":
      ! the state machine is skipped at runtime and the common-tail scaling by
      ! fraction_open is a 1.0 multiply.
      !
      ! Status encoding: 0=closed, 1=open, 2=opening, 3=closing.
      !
      do istruc = 1, nr_src_structures
         !
         ! Skip structures without rules - keep the "always open" defaults.
         !
         if (src_struc_rule_open(istruc) <= 0 .and. src_struc_rule_close(istruc) <= 0) cycle
         !
         nm_o1 = src_struc_nm_o1(istruc)
         nm_o2 = src_struc_nm_o2(istruc)
         !
         if (nm_o1 > 0) then
            !
            zs_o1 = real(zs(nm_o1), 4)
            !
         else
            !
            zs_o1 = 0.0
            !
         endif
         !
         if (nm_o2 > 0) then
            !
            zs_o2 = real(zs(nm_o2), 4)
            !
         else
            !
            zs_o2 = 0.0
            !
         endif
         !
         open_fires  = evaluate_rule(src_struc_rule_open(istruc),  zs_o1, zs_o2)
         close_fires = evaluate_rule(src_struc_rule_close(istruc), zs_o1, zs_o2)
         !
         if (open_fires .and. .not. close_fires) then
            !
            src_struc_status(istruc)        = 1
            src_struc_fraction_open(istruc) = 1.0
            status_str                  = 'open'
            !
         elseif (.not. open_fires .and. close_fires) then
            !
            src_struc_status(istruc)        = 0
            src_struc_fraction_open(istruc) = 0.0
            status_str                  = 'closed'
            !
         elseif (open_fires .and. close_fires) then
            !
            src_struc_status(istruc)        = 1
            src_struc_fraction_open(istruc) = 1.0
            status_str                  = 'open'
            write(logstr,'(a,a,a,a)')'Warning ! structure ', trim(src_struc_name(istruc)), &
                 ': both open and close rules fire at init; keeping structure open'
            call write_log(logstr, 0)
            !
         else
            !
            src_struc_status(istruc)        = 0
            src_struc_fraction_open(istruc) = 0.0
            status_str                  = 'closed'
            !
         endif
         !
         ! Transition timer is only consulted after a transition triggers;
         ! seed with t0 so the first rule-driven transition has a sane baseline.
         !
         src_struc_t_state(istruc) = t0
         !
         write(logstr,'(a,a,a,a)')'structure ', trim(src_struc_name(istruc)), &
              ' initialised status=', trim(status_str)
         call write_log(logstr, 0)
         !
      enddo
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine update_src_structures(t, dt)
      !
      ! Compute discharges through each drainage structure, accumulate them
      ! into qsrc(np) (endpoint 1: -qq, endpoint 2: +qq), and store
      ! per-structure signed discharge in src_struc_q_now(nr_src_structures)
      ! for his output. Sign convention: qq > 0 means flow from endpoint 1
      ! to endpoint 2.
      !
      ! Called AFTER update_discharges, which zeros qsrc first.
      !
      ! Atomic updates on qsrc(nm) guard against two structures (or a river
      ! and a structure) writing to the same cell under parallel execution.
      !
      ! Called from: update_continuity (sfincs_continuity), once per time step.
      !
      use sfincs_data
      use sfincs_timers
      !
      implicit none
      !
      real*8  :: t
      real*4  :: dt
      !
      integer :: istruc, nm_s1, nm_s2, nm_o1, nm_o2
      real*4  :: qq, elapsed, zs_o1, zs_o2
      real*4  :: frac, wdt, mng, zsill, dist, dzds, hgate, qq0, alpha
      real*4  :: dh, a_eff
      real*4  :: h_up, h_dn, qq_sign
      logical :: open_fires, close_fires
      !
      if (nr_src_structures <= 0) return
      !
      call timer_start('drainage structures')
      !
      !$acc parallel loop present( z_volume, zs, zb, qsrc, src_struc_q_now, &
      !$acc                        src_struc_nm_s1, src_struc_nm_s2, &
      !$acc                        src_struc_nm_o1, src_struc_nm_o2, &
      !$acc                        src_struc_type, src_struc_direction, &
      !$acc                        src_struc_q, src_struc_flow_coef, &
      !$acc                        src_struc_width, src_struc_sill_elevation, &
      !$acc                        src_struc_mannings_n, &
      !$acc                        src_struc_opening_duration, src_struc_closing_duration, &
      !$acc                        src_struc_height, &
      !$acc                        src_struc_invert_1, src_struc_invert_2, &
      !$acc                        src_struc_submergence_ratio, &
      !$acc                        src_struc_distance, src_struc_status, src_struc_fraction_open, &
      !$acc                        src_struc_t_state, &
      !$acc                        src_struc_rule_open, src_struc_rule_close, &
      !$acc                        rule_opcode, rule_atom, rule_cmp, rule_threshold, &
      !$acc                        rule_start, rule_length ) &
      !$acc              private( nm_s1, nm_s2, nm_o1, nm_o2, qq, elapsed, &
      !$acc                       zs_o1, zs_o2, frac, wdt, mng, zsill, dist, dzds, hgate, qq0, alpha, &
      !$acc                       dh, a_eff, &
      !$acc                       h_up, h_dn, qq_sign, &
      !$acc                       open_fires, close_fires )
      !$omp parallel do &
      !$omp   private( nm_s1, nm_s2, nm_o1, nm_o2, qq, elapsed, &
      !$omp            zs_o1, zs_o2, frac, wdt, mng, zsill, dist, dzds, hgate, qq0, alpha, &
      !$omp            dh, a_eff, &
      !$omp            h_up, h_dn, qq_sign, &
      !$omp            open_fires, close_fires ) &
      !$omp   schedule ( static )
      do istruc = 1, nr_src_structures
         !
         nm_s1 = src_struc_nm_s1(istruc)
         nm_s2 = src_struc_nm_s2(istruc)
         !
         if (nm_s1 > 0 .and. nm_s2 > 0) then
            !
            ! Open/close rule state machine (any structure type, any status).
            !
            ! Only runs if the user provided at least one of rules_open /
            ! rules_close. Structures without rules stay at the init-time
            ! defaults (status=1=open, fraction_open=1.0), which turns the
            ! common-tail scaling below into a no-op.
            !
            ! Status codes: 0=closed, 1=open, 2=opening, 3=closing.
            ! Transient states 2 and 3 advance purely on elapsed time so the
            ! state machine cannot thrash; rule evaluation happens in the
            ! terminal states 0 and 1 only. Obs points feed the rule inputs
            ! and default to the src pair in the marshal.
            !
            if (src_struc_rule_open(istruc) > 0 .or. src_struc_rule_close(istruc) > 0) then
               !
               nm_o1 = src_struc_nm_o1(istruc)
               nm_o2 = src_struc_nm_o2(istruc)
               !
               if (nm_o1 > 0) then
                  !
                  zs_o1 = real(zs(nm_o1), 4)
                  !
               else
                  !
                  zs_o1 = 0.0
                  !
               endif
               !
               if (nm_o2 > 0) then
                  !
                  zs_o2 = real(zs(nm_o2), 4)
                  !
               else
                  !
                  zs_o2 = 0.0
                  !
               endif
               !
               select case (int(src_struc_status(istruc)))
                  !
                  case (0)
                     !
                     ! closed - look for an open trigger
                     !
                     open_fires = evaluate_rule(src_struc_rule_open(istruc), zs_o1, zs_o2)
                     !
                     if (open_fires) then
                        !
                        src_struc_status(istruc)  = 2
                        src_struc_t_state(istruc) = real(t, 4)
                        !
                     endif
                     !
                  case (1)
                     !
                     ! open - look for a close trigger
                     !
                     close_fires = evaluate_rule(src_struc_rule_close(istruc), zs_o1, zs_o2)
                     !
                     if (close_fires) then
                        !
                        src_struc_status(istruc)  = 3
                        src_struc_t_state(istruc) = real(t, 4)
                        !
                     endif
                     !
                  case (2)
                     !
                     ! opening - advance on elapsed time; do not re-check rules
                     !
                     elapsed = real(t, 4) - src_struc_t_state(istruc)
                     !
                     if (src_struc_opening_duration(istruc) <= 0.0 .or. &
                         elapsed >= src_struc_opening_duration(istruc)) then
                        !
                        src_struc_status(istruc)        = 1
                        src_struc_fraction_open(istruc) = 1.0
                        !
                     else
                        !
                        src_struc_fraction_open(istruc) = elapsed / src_struc_opening_duration(istruc)
                        !
                     endif
                     !
                  case (3)
                     !
                     ! closing - advance on elapsed time; do not re-check rules
                     !
                     elapsed = real(t, 4) - src_struc_t_state(istruc)
                     !
                     if (src_struc_closing_duration(istruc) <= 0.0 .or. &
                         elapsed >= src_struc_closing_duration(istruc)) then
                        !
                        src_struc_status(istruc)        = 0
                        src_struc_fraction_open(istruc) = 0.0
                        !
                     else
                        !
                        src_struc_fraction_open(istruc) = 1.0 - elapsed / src_struc_closing_duration(istruc)
                        !
                     endif
                     !
               end select
               !
            endif
            !
            ! Per-type flux formula. Produces a raw signed discharge qq in
            ! m^3/s, before the common-tail scaling by fraction_open and
            ! direction filter.
            !
            select case(src_struc_type(istruc))
               !
               case(structure_pump)
                  !
                  ! Pump endpoint mapping: endpoint 1 (nm_s1) is the intake,
                  ! endpoint 2 (nm_s2) is the discharge. The pump enforces
                  ! qq >= 0 (no reverse flow); the sign convention in the
                  ! common tail below then sends water from nm_s1 to nm_s2.
                  !
                  qq = src_struc_q(istruc)
                  !
                  ! Reduction curve: scale by upstream depth so the pump cannot
                  ! pump the intake cell dry. pump_reduction_depth is a module-level
                  ! constant (see top of module); not user-tunable.
                  !
                  ! Turn this off for now. Does not work with subgrid.
                  !
                  !h_up = max(real(zs(nm_s1), 4) - zb(nm_s1), 0.0)
                  !qq   = qq * min(1.0, h_up / pump_reduction_depth)
                  !
               case(structure_culvert_simple)
                  !
                  ! Bidirectional: Q = flow_coef * sign(dh) * sqrt(|dh|).
                  ! The legacy "check_valve" alias maps to direction_positive
                  ! in the parser; the direction filter in the common tail
                  ! below restricts the allowed sign when requested.
                  !
                  if (zs(nm_s1) > zs(nm_s2)) then
                     !
                     qq =  src_struc_flow_coef(istruc) * sqrt(zs(nm_s1) - zs(nm_s2))
                     !
                  else
                     !
                     qq = -src_struc_flow_coef(istruc) * sqrt(zs(nm_s2) - zs(nm_s1))
                     !
                  endif
                  !
               case(structure_culvert)
                  !
                  ! Regime-aware culvert. The controlling sill is the higher
                  ! of the two inverts (flow cannot pass until the upstream
                  ! water level reaches it). Upstream / downstream are picked
                  ! by the water-level difference, so the structure is
                  ! bidirectional and the direction filter in the common tail
                  ! below restricts the sign when requested.
                  !
                  ! Two regimes, selected by h_dn/h_up against the user-set
                  ! submergence_ratio threshold (default 2/3 = 0.667, the
                  ! standard broad-crested-weir / Villemonte value):
                  !
                  !    submerged (h_dn/h_up >= threshold):
                  !       qq = flow_coef * a_eff * sqrt(2 g |dh|)
                  !    free / inlet-controlled (h_dn/h_up < threshold):
                  !       qq = flow_coef * a_eff * sqrt(2 g  h_up)
                  !
                  ! The flow area a_eff = width * min(h_up, height) caps at
                  ! the barrel height, so a deeply-submerged inlet can't
                  ! give unbounded discharge.
                  !
                  zsill = max(src_struc_invert_1(istruc), src_struc_invert_2(istruc))
                  !
                  dh = real(zs(nm_s1), 4) - real(zs(nm_s2), 4)
                  !
                  if (dh >= 0.0) then
                     !
                     h_up    = max(real(zs(nm_s1), 4) - zsill, 0.0)
                     h_dn    = max(real(zs(nm_s2), 4) - zsill, 0.0)
                     qq_sign =  1.0
                     !
                  else
                     !
                     h_up    = max(real(zs(nm_s2), 4) - zsill, 0.0)
                     h_dn    = max(real(zs(nm_s1), 4) - zsill, 0.0)
                     qq_sign = -1.0
                     !
                  endif
                  !
                  if (h_up <= 0.0) then
                     !
                     qq = 0.0
                     !
                  else
                     !
                     a_eff = src_struc_width(istruc) * min(h_up, src_struc_height(istruc))
                     !
                     if (h_dn / h_up >= src_struc_submergence_ratio(istruc)) then
                        !
                        qq = qq_sign * src_struc_flow_coef(istruc) * a_eff * sqrt(2.0 * g * abs(dh))
                        !
                     else
                        !
                        qq = qq_sign * src_struc_flow_coef(istruc) * a_eff * sqrt(2.0 * g * h_up)
                        !
                     endif
                     !
                  endif
                  !
               case(structure_gate)
                  !
                  ! Bidirectional culvert-style flow. Flow uses the src pair
                  ! (nm_s1/nm_s2), not the obs pair. Bates et al. (2010)
                  ! inertial formulation, per unit width:
                  !    q^{n+1} = (q^n - g*h*(dzs/ds)*dt) /
                  !              (1 + g*n^2*dt*|q^n| / h^{7/3})
                  ! with h = max(max(zs_s1, zs_s2) - zsill, 0).
                  ! Multiply by width to get the full structure discharge;
                  ! scaling by fraction_open happens in the common tail.
                  ! src_struc_q_now(istruc) holds the previous step's discharge
                  ! after the full common-tail scaling (width*fraction_open),
                  ! so unscale by (width*fraction_open) to recover qq0 in
                  ! per-unit-width form. Sign convention: qq > 0 means flow
                  ! nm_s1 -> nm_s2, matching dzds = (zs_s2 - zs_s1)/dist.
                  !
                  frac  = src_struc_fraction_open(istruc)
                  wdt   = src_struc_width(istruc)
                  mng   = src_struc_mannings_n(istruc)
                  zsill = src_struc_sill_elevation(istruc)
                  dist  = src_struc_distance(istruc)
                  !
                  dzds  = (real(zs(nm_s2), 4) - real(zs(nm_s1), 4)) / dist
                  hgate = max(max(real(zs(nm_s1), 4), real(zs(nm_s2), 4)) - zsill, 0.0)
                  !
                  if (hgate > 0.0 .and. frac > 0.0) then
                     !
                     qq0 = src_struc_q_now(istruc) / (wdt * max(frac, 0.001))
                     qq  = (qq0 - g * hgate * dzds * dt) / &
                           (1.0 + g * mng * mng * dt * abs(qq0) / hgate**(7.0/3.0))
                     qq  = qq * wdt
                     qq  = src_struc_flow_coef(istruc) * qq
                     !
                  else
                     !
                     qq = 0.0
                     !
                  endif
                  !
            end select
            !
            ! Common tail: scale by fraction_open (state-machine output) and
            ! apply the direction filter. Structures without rules sit at
            ! fraction_open=1.0 so the scaling is a no-op; structures with
            ! direction_both (the default) see the filter as a no-op too.
            !
            qq = qq * src_struc_fraction_open(istruc)
            !
            if (src_struc_direction(istruc) == direction_positive .and. qq < 0.0) qq = 0.0
            if (src_struc_direction(istruc) == direction_negative .and. qq > 0.0) qq = 0.0
            !
            ! Relaxation: blend new and previous discharge to damp oscillations.
            ! structure_relax is a dimensionless step count: alpha = 1/N damps
            ! the discharge response over roughly N time steps. Typical 1-10.
            !
            alpha = 1.0 / structure_relax
            qq = alpha * qq + (1.0 - alpha) * src_struc_q_now(istruc)
            !
            ! Limit discharge by available volume in the donor cell (endpoint 1
            ! for qq > 0, endpoint 2 for qq < 0).
            !
            if (subgrid) then
               !
               if (qq > 0.0) then
                  !
                  qq = min(qq,  max(z_volume(nm_s1), 0.0) / dt)
                  !
               else
                  !
                  qq = max(qq, -max(z_volume(nm_s2), 0.0) / dt)
                  !
               endif
               !
            else
               !
               if (qq > 0.0) then
                  !
                  qq = min(qq,  max((zs(nm_s1) - zb(nm_s1)) * cell_area(z_flags_iref(nm_s1)), 0.0) / dt)
                  !
               else
                  !
                  qq = max(qq, -max((zs(nm_s2) - zb(nm_s2)) * cell_area(z_flags_iref(nm_s2)), 0.0) / dt)
                  !
               endif
               !
            endif
            !
            src_struc_q_now(istruc) = qq
            !
            ! Accumulate into cell-wise qsrc. Atomic guards against multiple
            ! structures (or a river and a structure) in the same cell. Sign
            ! convention qq > 0 means flow nm_s1 -> nm_s2, so qq is subtracted
            ! at endpoint 1 and added at endpoint 2.
            !
            !$acc atomic update
            !$omp atomic
            qsrc(nm_s1) = qsrc(nm_s1) - qq
            !$acc atomic update
            !$omp atomic
            qsrc(nm_s2) = qsrc(nm_s2) + qq
            !
         endif
         !
      enddo
      !$omp end parallel do
      !$acc end parallel loop
      !
      call timer_stop('drainage structures')
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
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
      !    type = "gate"                ! required, one of pump/culvert_simple/gate/culvert
      !                                 ! legacy alias: "check_valve" -> culvert_simple + direction="positive"
      !                                 ! note: "culvert" now resolves to the detailed-culvert physics type;
      !                                 !       users wanting the lumped one-coefficient form must say
      !                                 !       "culvert_simple" explicitly. Orifice behaviour is recoverable
      !                                 !       as "culvert" with submergence_ratio = 0.0.
      !    direction = "both"           ! optional, culvert_simple/culvert only
      !                                 ! one of "both" (default), "positive", "negative"
      !                                 ! positive: allow flow src_1 -> src_2 only
      !                                 ! negative: allow flow src_2 -> src_1 only
      !    src_1 = [x, y] ; src_2 = [x, y]
      !    obs_1 = [x, y] ; obs_2 = [x, y]
      !    q = ...                      ! pump discharge
      !    width = ... ; sill_elevation = ... ; mannings_n = ...
      !    opening_duration = ... ; closing_duration = ...
      !    flow_coef = ...              ! culvert_simple / culvert flow coefficient
      !    height = ...                                 ! culvert pipe height (m)
      !    invert_1 = ... ; invert_2 = ...              ! culvert invert elevations at src_1/src_2 ends
      !    submergence_ratio = ...                      ! culvert submergence threshold h_dn/h_up (-)
      !    rules_open  = "(z1<0.5 | z2-z1>0.05) & z2<1.5"   ! optional trigger expr
      !    rules_close = "z2>2.0"                           ! optional trigger expr
      !
      ! Per-type required keys (enforced on parse):
      !    pump           : name, src_1, src_2, q
      !    culvert_simple : name, src_1, src_2, flow_coef
      !    gate           : name, src_1, src_2, width, sill_elevation
      !    culvert        : name, src_1, src_2,
      !                     width, height, invert_1, invert_2
      !                     (optional: flow_coef=0.6, submergence_ratio=0.667)
      !
      ! On success, structures is allocated to the exact number of entries
      ! (can be 0). On any I/O or parse failure, structures is left
      ! unallocated and ierr is non-zero.
      !
      ! This routine does not modify module state; it is the caller's job to
      ! decide what to do with the parsed array.
      !
      ! Called from: initialize_src_structures (this module).
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
      type(toml_table), pointer        :: tbl_struct
      character(len=:), allocatable    :: name_str, type_str, rule_str, dir_str, type_str_lc
      integer                          :: n_struct, i, stat, ierr_parse
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
         call parse_structure_type(type_str, structures(i)%structure_type, ierr_parse)
         !
         if (ierr_parse /= 0) then
            !
            ierr = ierr_parse
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
                    'name', 'q' ], i, ierr)
               call check_required_coord_pair(tbl_struct, 'src_1', i, ierr)
               call check_required_coord_pair(tbl_struct, 'src_2', i, ierr)
               !
            case (structure_culvert_simple)
               !
               call check_required(tbl_struct, [ character(len=16) :: &
                    'name', 'flow_coef' ], i, ierr)
               call check_required_coord_pair(tbl_struct, 'src_1', i, ierr)
               call check_required_coord_pair(tbl_struct, 'src_2', i, ierr)
               !
            case (structure_culvert)
               !
               call check_required(tbl_struct, [ character(len=16) :: &
                    'name', 'width', 'height', 'invert_1', 'invert_2' ], i, ierr)
               call check_required_coord_pair(tbl_struct, 'src_1', i, ierr)
               call check_required_coord_pair(tbl_struct, 'src_2', i, ierr)
               !
            case (structure_gate)
               !
               call check_required(tbl_struct, [ character(len=16) :: &
                    'name', 'width', 'sill_elevation' ], i, ierr)
               call check_required_coord_pair(tbl_struct, 'src_1', i, ierr)
               call check_required_coord_pair(tbl_struct, 'src_2', i, ierr)
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
         call read_coord_pair(tbl_struct, 'src_1', structures(i)%x_s1, structures(i)%y_s1, i, ierr)
         call read_coord_pair(tbl_struct, 'src_2', structures(i)%x_s2, structures(i)%y_s2, i, ierr)
         !
         structures(i)%has_o1 = tbl_struct%has_key('obs_1')
         structures(i)%has_o2 = tbl_struct%has_key('obs_2')
         !
         call read_coord_pair(tbl_struct, 'obs_1', structures(i)%x_o1, structures(i)%y_o1, i, ierr)
         call read_coord_pair(tbl_struct, 'obs_2', structures(i)%x_o2, structures(i)%y_o2, i, ierr)
         !
         ! Named physical parameters. Defaults are picked to avoid NaN in
         ! arithmetic and to match the legacy-reader fallbacks.
         !
         call get_value(tbl_struct, 'q',                structures(i)%q,                0.0,    stat=stat)
         call get_value(tbl_struct, 'width',            structures(i)%width,            0.0,    stat=stat)
         call get_value(tbl_struct, 'sill_elevation',   structures(i)%sill_elevation,   0.0,    stat=stat)
         !
         ! opening_duration / closing_duration default depends on type: gate keeps
         ! its historical 600 s (legacy "dtype 4" gates always had finite ramp
         ! durations), pump / culvert_simple / culvert default to 0 s (instant
         ! open/close when a rule fires; skips the transient states 2 and 3).
         !
         if (structures(i)%structure_type == structure_gate) then
            !
            call get_value(tbl_struct, 'opening_duration', structures(i)%opening_duration, 600.0, stat=stat)
            call get_value(tbl_struct, 'closing_duration', structures(i)%closing_duration, 600.0, stat=stat)
            !
         else
            !
            call get_value(tbl_struct, 'opening_duration', structures(i)%opening_duration, 0.0, stat=stat)
            call get_value(tbl_struct, 'closing_duration', structures(i)%closing_duration, 0.0, stat=stat)
            !
         endif
         !
         ! flow_coef default differs by type: 1.0 for culvert_simple (legacy
         ! lumped one-coefficient form), 0.6 for the detailed culvert
         ! (standard orifice discharge coefficient).
         !
         if (structures(i)%structure_type == structure_culvert) then
            !
            call get_value(tbl_struct, 'flow_coef', structures(i)%flow_coef, 0.6, stat=stat)
            !
         else
            !
            call get_value(tbl_struct, 'flow_coef', structures(i)%flow_coef, 1.0, stat=stat)
            !
         endif
         !
         ! mannings_n (gate only). Default 0.024 for concrete-lined gate sill.
         !
         call get_value(tbl_struct, 'mannings_n', structures(i)%mannings_n, 0.024, stat=stat)
         !
         ! Detailed-culvert geometry + submergence threshold. Geometry keys
         ! are required (enforced above); submergence_ratio defaults to 2/3
         ! (0.667), the standard broad-crested-weir / Villemonte value.
         !
         call get_value(tbl_struct, 'height',            structures(i)%height,            0.0,   stat=stat)
         call get_value(tbl_struct, 'invert_1',          structures(i)%invert_1,          0.0,   stat=stat)
         call get_value(tbl_struct, 'invert_2',          structures(i)%invert_2,          0.0,   stat=stat)
         call get_value(tbl_struct, 'submergence_ratio', structures(i)%submergence_ratio, 0.667, stat=stat)
         !
         ! Optional direction filter (culvert_simple / culvert). Default is
         ! direction_both. Unknown strings are a hard error.
         !
         structures(i)%direction = direction_both
         !
         if (allocated(dir_str)) deallocate(dir_str)
         call get_value(tbl_struct, 'direction', dir_str, stat=stat)
         !
         if (allocated(dir_str)) then
            !
            call parse_direction(dir_str, structures(i)%direction, ierr_parse)
            !
            if (ierr_parse /= 0) then
               !
               ierr = ierr_parse
               write(logstr,'(a,a,a,i0)')' Error ! Unknown direction "', trim(dir_str), &
                    '" in src_structure entry ', i
               call write_log(logstr, 1)
               call cleanup_on_error()
               return
               !
            endif
            !
         endif
         !
         ! Legacy alias side-effect: "check_valve" always pins direction_positive
         ! regardless of any explicit direction key. Detect on the lowered type
         ! string so "Check_Valve" etc. are handled identically.
         !
         type_str_lc = to_lower(type_str)
         !
         if (type_str_lc == 'check_valve') then
            !
            structures(i)%direction = direction_positive
            !
         endif
         !
         ! Optional rules_open / rules_close string expressions. Absent keys
         ! leave the rule strings unallocated on the derived type; marshal
         ! treats that as "no trigger".
         !
         if (allocated(rule_str)) deallocate(rule_str)
         call get_value(tbl_struct, 'rules_open', rule_str, stat=stat)
         if (allocated(rule_str)) structures(i)%rule_open = rule_str
         !
         if (allocated(rule_str)) deallocate(rule_str)
         call get_value(tbl_struct, 'rules_close', rule_str, stat=stat)
         if (allocated(rule_str)) structures(i)%rule_close = rule_str
         !
      enddo
      !
      contains
         !
         subroutine cleanup_on_error()
            !
            ! Internal helper for the parse loop: drop the partially-filled
            ! structures(:) array so the caller always sees it unallocated on
            ! error exit. Trivial deallocator.
            !
            ! Called from: read_toml_src_structures (host routine) on every
            ! error bail-out path.
            !
            if (allocated(structures)) deallocate(structures)
            !
         end subroutine
         !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine check_required(table, keys, seq_index, ierr)
      !
      ! Verify that every key in "keys" is present in the TOML table. Missing
      ! keys are reported to the log (naming the structure by its 1-based
      ! sequence index, since "name" itself may be the missing key) and ierr
      ! is set non-zero. Presence is checked via has_key so that a legal
      ! value of 0.0 is not mistaken for "missing".
      !
      ! Called from: read_toml_src_structures (this module), once per
      ! structure entry in the per-type required-field validation block.
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
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine check_required_coord_pair(table, key_base, seq_index, ierr)
      !
      ! Verify that a coordinate pair "<key_base> = [x, y]" is present in the
      ! TOML table. Emits a single missing-key error when absent.
      !
      ! Called from: read_toml_src_structures (this module), once per required
      ! coordinate pair (src_1, src_2) in the per-type validation block.
      !
      use tomlf
      !
      implicit none
      !
      type(toml_table), pointer, intent(in)    :: table
      character(len=*),          intent(in)    :: key_base
      integer,                   intent(in)    :: seq_index
      integer,                   intent(inout) :: ierr
      !
      if (.not. table%has_key(trim(key_base))) then
         !
         write(logstr,'(a,i0,a,a,a)')' Error ! Structure #', seq_index, &
              ' is missing required coordinate pair "', trim(key_base), ' = [x, y]"'
         call write_log(logstr, 1)
         ierr = 1
         !
      endif
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine read_coord_pair(table, key_base, x, y, seq_index, ierr)
      !
      ! Read a coordinate pair "<key_base> = [x, y]" from a TOML table.
      !
      ! If the key is absent, x and y are left at 0.0 and no error is raised
      ! here — presence of required pairs is enforced separately by
      ! check_required_coord_pair.
      !
      ! Called from: read_toml_src_structures (this module), once per
      ! coordinate pair (src_1, src_2, obs_1, obs_2) per structure entry.
      !
      use tomlf
      !
      implicit none
      !
      type(toml_table), pointer, intent(in)    :: table
      character(len=*),          intent(in)    :: key_base
      real,                      intent(out)   :: x, y
      integer,                   intent(in)    :: seq_index
      integer,                   intent(inout) :: ierr
      !
      type(toml_array), pointer :: arr
      integer                   :: n, stat
      !
      x = 0.0
      y = 0.0
      !
      if (.not. table%has_key(trim(key_base))) return
      !
      nullify(arr)
      call get_value(table, trim(key_base), arr, requested=.false., stat=stat)
      !
      if (.not. associated(arr)) then
         !
         write(logstr,'(a,a,a,i0,a)')' Error ! Key "', trim(key_base), &
              '" in src_structure #', seq_index, ' must be a 2-element array [x, y]'
         call write_log(logstr, 1)
         ierr = 1
         return
         !
      endif
      !
      n = len(arr)
      !
      if (n /= 2) then
         !
         write(logstr,'(a,a,a,i0,a,i0,a)')' Error ! Key "', trim(key_base), &
              '" in src_structure #', seq_index, ' must have exactly 2 elements (got ', n, ')'
         call write_log(logstr, 1)
         ierr = 1
         return
         !
      endif
      !
      call get_value(arr, 1, x, stat=stat)
      call get_value(arr, 2, y, stat=stat)
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine parse_structure_type(str, code, ierr)
      !
      ! Translate a TOML "type" string to one of the structure_* codes.
      !
      ! Legacy alias accepted (quietly, no warning):
      !    "check_valve" -> structure_culvert_simple
      !                     (caller is responsible for pinning direction_positive)
      !
      ! Note: "culvert" now resolves to structure_culvert (the detailed
      ! physics-based pipe-flow type). Users wanting the lumped one-coefficient
      ! form must say "culvert_simple" explicitly.
      !
      ! Called from: read_toml_src_structures (this module), once per entry
      ! to resolve the required "type" key.
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
         case ('culvert_simple', 'check_valve')
            !
            code = structure_culvert_simple
            !
         case ('gate')
            !
            code = structure_gate
            !
         case ('culvert')
            !
            code = structure_culvert
            !
         case default
            !
            ierr = 1
            !
      end select
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine parse_direction(str, code, ierr)
      !
      ! Translate a TOML "direction" string to one of the direction_* codes.
      ! Accepts "both" / "positive" / "negative" case-insensitively.
      !
      ! Called from: read_toml_src_structures (this module) when an optional
      ! "direction" key is present on a structure entry.
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
         case ('both')
            !
            code = direction_both
            !
         case ('positive')
            !
            code = direction_positive
            !
         case ('negative')
            !
            code = direction_negative
            !
         case default
            !
            ierr = 1
            !
      end select
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   function to_lower(str) result(lower)
      !
      ! Return a lowercase copy of str (ASCII only).
      !
      ! Called from: parse_structure_type, parse_direction, and
      ! convert_legacy_to_toml (all in this module).
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
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine write_src_structures_log_summary()
      !
      ! Emit a one-block-per-structure description of every parsed src
      ! structure to the log file. Intended for operator review at init
      ! time; not printed to stdout.
      !
      ! Called from: initialize_src_structures (this module), once after
      ! the marshal runs and cell indices have been resolved.
      !
      implicit none
      !
      integer :: i
      character(len=32) :: type_str, dir_str
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
         select case (int(src_struc_type(i)))
            !
            case (structure_pump)
               !
               type_str = 'pump'
               !
            case (structure_culvert_simple)
               !
               type_str = 'culvert_simple'
               !
            case (structure_culvert)
               !
               type_str = 'culvert'
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
         write(logstr,'(a22,1x,a)')              '  name               :',              trim(src_struc_name(i))
         call write_log(logstr, 0)
         !
         write(logstr,'(a22,1x,a)')              '  type               :',              trim(type_str)
         call write_log(logstr, 0)
         !
         write(logstr,'(a22,1x,a,a,a,a,a)')      '  src_1              :',             '(', trim(fmt_real(src_struc_x_s1(i), 3)), ', ', trim(fmt_real(src_struc_y_s1(i), 3)), ')'
         call write_log(logstr, 0)
         !
         write(logstr,'(a22,1x,a,a,a,a,a)')      '  src_2              :',             '(', trim(fmt_real(src_struc_x_s2(i), 3)), ', ', trim(fmt_real(src_struc_y_s2(i), 3)), ')'
         call write_log(logstr, 0)
         !
         ! obs coords are meaningful for culvert_simple / gate.
         !
         if (src_struc_type(i) == structure_culvert_simple .or. &
             src_struc_type(i) == structure_gate) then
            !
            write(logstr,'(a22,1x,a,a,a,a,a)')   '  obs_1              :',              '(', trim(fmt_real(src_struc_x_o1(i), 3)), ', ', trim(fmt_real(src_struc_y_o1(i), 3)), ')'
            call write_log(logstr, 0)
            !
            write(logstr,'(a22,1x,a,a,a,a,a)')   '  obs_2              :',              '(', trim(fmt_real(src_struc_x_o2(i), 3)), ', ', trim(fmt_real(src_struc_y_o2(i), 3)), ')'
            call write_log(logstr, 0)
            !
         endif
         !
         if (src_struc_type(i) == structure_pump) then
            !
            write(logstr,'(a22,1x,a,a)')            '  discharge          :',       trim(fmt_real(src_struc_q(i), 4)), ' (m3/s)'
            call write_log(logstr, 0)
            !
         endif
         !
         if (src_struc_type(i) == structure_culvert_simple) then
            !
            write(logstr,'(a22,1x,a)')              '  flow_coef          :',       trim(fmt_real(src_struc_flow_coef(i), 4))
            call write_log(logstr, 0)
            !
         endif
         !
         ! Direction filter (culvert_simple / culvert)
         !
         if (src_struc_type(i) == structure_culvert_simple .or. &
             src_struc_type(i) == structure_culvert) then
            !
            select case (src_struc_direction(i))
               !
               case (direction_both)
                  !
                  dir_str = 'both'
                  !
               case (direction_positive)
                  !
                  dir_str = 'positive'
                  !
               case (direction_negative)
                  !
                  dir_str = 'negative'
                  !
               case default
                  !
                  dir_str = 'unknown'
                  !
            end select
            !
            write(logstr,'(a22,1x,a)')              '  direction          :',      trim(dir_str)
            call write_log(logstr, 0)
            !
         endif
         !
         if (src_struc_type(i) == structure_culvert) then
            !
            write(logstr,'(a22,1x,a,a)')            '  width              :',              trim(fmt_real(src_struc_width(i), 4)),             ' (m)'
            call write_log(logstr, 0)
            !
            write(logstr,'(a22,1x,a,a)')            '  height             :',             trim(fmt_real(src_struc_height(i), 4)),            ' (m)'
            call write_log(logstr, 0)
            !
            write(logstr,'(a22,1x,a,a)')            '  invert_1           :',           trim(fmt_real(src_struc_invert_1(i), 4)),          ' (m)'
            call write_log(logstr, 0)
            !
            write(logstr,'(a22,1x,a,a)')            '  invert_2           :',           trim(fmt_real(src_struc_invert_2(i), 4)),          ' (m)'
            call write_log(logstr, 0)
            !
            write(logstr,'(a22,1x,a)')              '  flow_coef          :',          trim(fmt_real(src_struc_flow_coef(i), 4))
            call write_log(logstr, 0)
            !
            write(logstr,'(a22,1x,a)')              '  submergence_ratio  :',  trim(fmt_real(src_struc_submergence_ratio(i), 4))
            call write_log(logstr, 0)
            !
         endif
         !
         if (src_struc_type(i) == structure_gate) then
            !
            write(logstr,'(a22,1x,a,a)')            '  width              :',              trim(fmt_real(src_struc_width(i), 4)),             ' (m)'
            call write_log(logstr, 0)
            !
            write(logstr,'(a22,1x,a,a)')            '  sill_elevation     :',     trim(fmt_real(src_struc_sill_elevation(i), 4)),    ' (m)'
            call write_log(logstr, 0)
            !
            write(logstr,'(a22,1x,a)')              '  mannings_n         :',         trim(fmt_real(src_struc_mannings_n(i), 4))
            call write_log(logstr, 0)
            !
            write(logstr,'(a22,1x,a,a)')            '  opening_duration   :',   trim(fmt_real(src_struc_opening_duration(i), 2)),  ' (s)'
            call write_log(logstr, 0)
            !
            write(logstr,'(a22,1x,a,a)')            '  closing_duration   :',   trim(fmt_real(src_struc_closing_duration(i), 2)),  ' (s)'
            call write_log(logstr, 0)
            !
         endif
         !
         if (src_struc_rule_open(i) > 0) then
            !
            if (len_trim(src_struc_rule_open_src(i)) > 0) then
               !
               write(logstr,'(a22,1x,a,a,a)')       '  rules_open         :',       '"', trim(src_struc_rule_open_src(i)), '"'
               !
            else
               !
               write(logstr,'(a22,1x,a)')           '  rules_open         :',       '(set)'
               !
            endif
            !
            call write_log(logstr, 0)
            !
         endif
         !
         if (src_struc_rule_close(i) > 0) then
            !
            if (len_trim(src_struc_rule_close_src(i)) > 0) then
               !
               write(logstr,'(a22,1x,a,a,a)')       '  rules_close        :',      '"', trim(src_struc_rule_close_src(i)), '"'
               !
            else
               !
               write(logstr,'(a22,1x,a)')           '  rules_close        :',      '(set)'
               !
            endif
            !
            call write_log(logstr, 0)
            !
         endif
         !
         ! Opening/closing durations. For gate structures these are always
         ! printed (above); for other types only print if rules are set and
         ! the duration is non-zero (non-default).
         !
         if (src_struc_type(i) /= structure_gate) then
            !
            if ((src_struc_rule_open(i)  > 0 .or. src_struc_rule_close(i) > 0) .and. &
                (src_struc_opening_duration(i) > 0.0 .or. src_struc_closing_duration(i) > 0.0)) then
               !
               write(logstr,'(a22,1x,a,a)')         '  opening_duration   :',  trim(fmt_real(src_struc_opening_duration(i), 2)),  ' (s)'
               call write_log(logstr, 0)
               !
               write(logstr,'(a22,1x,a,a)')         '  closing_duration   :',  trim(fmt_real(src_struc_closing_duration(i), 2)),  ' (s)'
               call write_log(logstr, 0)
               !
            endif
            !
         endif
         !
         call write_log('', 0)
         !
      enddo
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
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
      ! The output path is derived from legacy_path: if it ends in ".drn"
      ! (case-insensitive) the suffix ".toml" is inserted before the ".drn",
      ! otherwise ".toml" is appended. The resolved path is returned in
      ! toml_path for the caller to feed into the TOML reader.
      !
      ! The converter is deliberately minimal: no coord sanity checks, no
      ! duplicate-name detection, no preservation of comments. It exists only
      ! to remove the parallel legacy parsing machinery that used to live
      ! in this module.
      !
      ! Called from: initialize_src_structures (this module) when toml-f
      ! rejects the drn file on the initial probe.
      !
      implicit none
      !
      character(len=*), intent(in)  :: legacy_path
      character(len=*), intent(out) :: toml_path
      integer,          intent(out) :: ierr
      !
      integer            :: u_in, u_out, stat, n_struct, dtype
      integer            :: len_in, ext_pos
      real*4             :: x2, y2, x1, y1, par
      real*4             :: g_width, g_sill, g_mann, g_zmin, g_zmax, g_tcls
      character(len=512) :: line, trimmed
      character(len=32)  :: name_str
      character(len=16)  :: type_name, par_name, dir_name
      character(len=13)  :: zmin_str, zmax_str
      character(len=128) :: rule_open_str, rule_close_str
      !
      ierr     = 0
      n_struct = 0
      u_in     = 501
      u_out    = 502
      !
      ! Derive the TOML sibling path from legacy_path. If legacy_path ends
      ! in ".drn" (case-insensitive), insert ".toml" before the extension;
      ! else append ".toml".
      !
      len_in  = len_trim(legacy_path)
      ext_pos = 0
      !
      if (len_in >= 4) then
         !
         if (to_lower(legacy_path(len_in-3:len_in)) == '.drn') then
            !
            ext_pos = len_in - 3
            !
         endif
         !
      endif
      !
      if (ext_pos > 0) then
         !
         toml_path = legacy_path(1:ext_pos-1) // '.toml' // legacy_path(ext_pos:len_in)
         !
      else
         !
         toml_path = legacy_path(1:len_in) // '.toml'
         !
      endif
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
         ! Columns: x1, y1, x2, y2, dtype, par.
         ! (legacy xsnk=intake -> src_1; legacy xsrc=outfall -> src_2).
         !
         read(line, *, iostat=stat) x1, y1, x2, y2, dtype, par
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
         ! dir_name is left blank unless dtype pins a direction filter; a blank
         ! dir_name causes the emitter below to skip the direction key entirely,
         ! which reads back as direction_both (the default).
         !
         dir_name = ''
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
               ! legacy culvert -> bidirectional culvert_simple
               !
               type_name = 'culvert_simple'
               par_name  = 'flow_coef'
               !
            case (3)
               !
               ! legacy check_valve -> culvert_simple with direction = "positive"
               !
               type_name = 'culvert_simple'
               par_name  = 'flow_coef'
               dir_name  = 'positive'
               !
            case (4)
               !
               ! Water-level-triggered gate. Legacy columns past dtype:
               !   width, sill_elevation, mannings_n, zmin, zmax, t_close.
               ! Re-read the whole line to pull those extra columns.
               !
               read(line, *, iostat=stat) x1, y1, x2, y2, dtype, &
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
               write(u_out,'(a,es14.6,a,es14.6,a)') 'src_1            = [', x1, ', ', y1, ']'
               write(u_out,'(a,es14.6,a,es14.6,a)') 'src_2            = [', x2, ', ', y2, ']'
               write(u_out,'(a,es14.6)')  'width            = ', g_width
               write(u_out,'(a,es14.6)')  'sill_elevation   = ', g_sill
               write(u_out,'(a,es14.6)')  'mannings_n       = ', g_mann
               write(u_out,'(a,es14.6)')  'opening_duration = ', g_tcls
               write(u_out,'(a,es14.6)')  'closing_duration = ', g_tcls
               write(u_out,'(a,a,a)')     'rules_open       = "', trim(rule_open_str),  '"'
               write(u_out,'(a,a,a)')     'rules_close      = "', trim(rule_close_str), '"'
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
         !
         if (len_trim(dir_name) > 0) then
            !
            write(u_out,'(a,a,a)')  'direction = "', trim(dir_name), '"'
            !
         endif
         !
         write(u_out,'(a,es14.6,a,es14.6,a)') 'src_1   = [', x1, ', ', y1, ']'
         write(u_out,'(a,es14.6,a,es14.6,a)') 'src_2   = [', x2, ', ', y2, ']'
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
   !
end module
