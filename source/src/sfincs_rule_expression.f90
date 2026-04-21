module sfincs_rule_expression
   !
   ! Boolean rule mini-language used by gate-like src_structures to decide
   ! when to open or close. The grammar is:
   !
   !    expr     := or_expr
   !    or_expr  := and_expr ( '|' and_expr )*
   !    and_expr := comp     ( '&' comp     )*
   !    comp     := '(' expr ')' | atom cmp_op number
   !    atom     := 'z1' | 'z2' | 'z2-z1'          (case-insensitive)
   !    cmp_op   := '<' | '>' | '<=' | '>=' | '=' | '=='
   !    number   := real literal
   !
   ! Precedence: paren > comp > '&' > '|'. Left-associative.
   !
   ! Each rule is compiled to a reverse-polish bytecode stream in four
   ! parallel module-level arrays (opcode / atom / cmp / threshold) and
   ! registered in a small parallel (start, length) registry indexed by
   ! integer rule_id. Callers keep only that rule_id per-structure; the
   ! module owns both the op stream (so all rules across all structures
   ! concatenate into one buffer) and the registry.
   !
   ! The evaluator is a small fixed-depth logical stack machine that is
   ! ACC-safe (no allocations, no strings, no I/O).
   !
   implicit none
   !
   private
   !
   ! Public bytecode storage. These four parallel arrays hold the
   ! concatenated op streams for every rule that has been parsed. They
   ! are public so sfincs_openacc can name them in !$acc directives.
   !
   public :: rule_opcode, rule_atom, rule_cmp, rule_threshold, rule_n_ops
   !
   ! Public rule registry. rule_start(id) / rule_length(id) index into
   ! the op streams above. n_rules is the highest rule_id currently
   ! allocated; entries at index 0 are not used (rule_id == 0 is the
   ! "no rule" sentinel).
   !
   public :: rule_start, rule_length, n_rules
   !
   ! Public API.
   !
   public :: add_rule, evaluate_rule, finalize_rule_storage
   !
   ! ---------------------------------------------------------------
   ! Opcodes.
   !
   integer, parameter :: op_cmp               = 1
   integer, parameter :: op_and               = 2
   integer, parameter :: op_or                = 3
   !
   ! Atom codes (the z value being compared).
   !
   integer, parameter :: atom_z1              = 1
   integer, parameter :: atom_z2              = 2
   integer, parameter :: atom_z2_minus_z1     = 3
   !
   ! Comparator codes.
   !
   integer, parameter :: cmp_lt               = 1
   integer, parameter :: cmp_gt               = 2
   integer, parameter :: cmp_le               = 3
   integer, parameter :: cmp_ge               = 4
   integer, parameter :: cmp_eq               = 5
   !
   ! Parser / evaluator capacity limits.
   !
   integer, parameter :: expr_stack_max        = 16   ! evaluator logical-stack depth
   integer, parameter :: expr_ops_max_per_rule = 64   ! max bytecode length per rule (parser rejects longer)
   integer, parameter :: expr_tokens_max       = 128  ! max tokens in a single rule string
   !
   ! ---------------------------------------------------------------
   ! Shared bytecode storage.
   !
   ! add_rule appends into these arrays and auto-grows them. Handles
   ! (rule_start, rule_length) in the registry index into them.
   !
   integer, allocatable, save :: rule_opcode(:)
   integer, allocatable, save :: rule_atom(:)
   integer, allocatable, save :: rule_cmp(:)
   real,    allocatable, save :: rule_threshold(:)
   integer,              save :: rule_n_ops    = 0
   integer,              save :: rule_capacity = 0
   !
   ! Rule registry. Indexed by rule_id in [1 .. n_rules]. rule_id == 0
   ! is the "never fires" sentinel and has no registry entry.
   !
   integer, allocatable, save :: rule_start(:)
   integer, allocatable, save :: rule_length(:)
   integer,              save :: n_rules                = 0
   integer,              save :: rule_registry_capacity = 0
   !
   integer, parameter :: initial_capacity          = 256
   integer, parameter :: initial_registry_capacity = 16
   !
contains
   !
   !
   subroutine add_rule(src, rule_id, ierr, errmsg)
   !
   ! Parse the boolean expression in src, append its bytecode to the
   ! module-level rule_* arrays, and register a new rule entry that
   ! points at it. Returns the new rule_id. An empty or whitespace-only
   ! src returns rule_id = 0 (the "never fires" sentinel; no registry
   ! entry is created). On parse failure ierr /= 0, rule_id = 0, and
   ! errmsg carries the diagnostic.
   !
   implicit none
   !
   character(len=*),           intent(in)  :: src
   integer,                    intent(out) :: rule_id
   integer,                    intent(out) :: ierr
   character(len=*), optional, intent(out) :: errmsg
   !
   ! Per-call scratch buffers for the parsed bytecode; sized to the
   ! parser's per-rule cap. Copied into the module storage on success.
   !
   integer           :: ops_buf   (expr_ops_max_per_rule)
   integer           :: atoms_buf (expr_ops_max_per_rule)
   integer           :: cmps_buf  (expr_ops_max_per_rule)
   real              :: thr_buf   (expr_ops_max_per_rule)
   integer           :: nops, new_start
   character(len=256):: local_errmsg
   !
   character(len=len(src)) :: src_nospace
   integer                 :: ic, ip, jp
   !
   rule_id = 0
   ierr    = 0
   if (present(errmsg)) errmsg = ''
   !
   ! Empty / whitespace-only source: never fires, no registry entry.
   !
   if (len_trim(src) == 0) return
   !
   ! Strip all whitespace (space, tab, LF, CR) so callers can write
   ! 'z1 < 0.5' or 'z2 - z1 > 0.05' as freely as 'z1<0.5' / 'z2-z1>0.05'.
   !
   jp = 0
   do ip = 1, len(src)
      !
      ic = iachar(src(ip:ip))
      !
      if (ic /= iachar(' ') .and. ic /= 9 .and. ic /= 10 .and. ic /= 13) then
         !
         jp = jp + 1
         src_nospace(jp:jp) = src(ip:ip)
         !
      endif
      !
   enddo
   !
   if (jp == 0) return
   !
   call parse_rule_expression(src_nospace(1:jp), ops_buf, atoms_buf, cmps_buf, thr_buf, &
        nops, ierr, local_errmsg)
   !
   if (ierr /= 0) then
      !
      if (present(errmsg)) errmsg = local_errmsg
      return
      !
   endif
   !
   if (nops <= 0) then
      !
      ierr = 1
      if (present(errmsg)) errmsg = 'rule parse produced no ops'
      return
      !
   endif
   !
   ! Ensure the op stream has room for the new ops.
   !
   call grow_rule_storage(rule_n_ops + nops)
   !
   new_start = rule_n_ops + 1
   !
   rule_opcode   (new_start : new_start + nops - 1) = ops_buf  (1:nops)
   rule_atom     (new_start : new_start + nops - 1) = atoms_buf(1:nops)
   rule_cmp      (new_start : new_start + nops - 1) = cmps_buf (1:nops)
   rule_threshold(new_start : new_start + nops - 1) = thr_buf  (1:nops)
   !
   rule_n_ops = rule_n_ops + nops
   !
   ! Register the new rule and return its id.
   !
   call grow_rule_registry(n_rules + 1)
   !
   n_rules             = n_rules + 1
   rule_start(n_rules) = new_start
   rule_length(n_rules)= nops
   rule_id             = n_rules
   !
   end subroutine
   !
   !
   subroutine finalize_rule_storage()
   !
   ! Shrink the rule_* op streams to exactly rule_n_ops and the registry
   ! to exactly n_rules. If nothing was ever allocated (no rules parsed)
   ! everything is allocated to size 0 so downstream openacc directives
   ! can reference the arrays safely.
   !
   implicit none
   !
   integer, allocatable :: tmp_i(:)
   real,    allocatable :: tmp_r(:)
   !
   ! Op streams.
   !
   if (.not. allocated(rule_opcode)) then
      !
      allocate(rule_opcode(0))
      allocate(rule_atom(0))
      allocate(rule_cmp(0))
      allocate(rule_threshold(0))
      rule_capacity = 0
      rule_n_ops    = 0
      !
   else if (rule_capacity /= rule_n_ops) then
      !
      allocate(tmp_i(rule_n_ops))
      if (rule_n_ops > 0) tmp_i = rule_opcode(1:rule_n_ops)
      call move_alloc(tmp_i, rule_opcode)
      !
      allocate(tmp_i(rule_n_ops))
      if (rule_n_ops > 0) tmp_i = rule_atom(1:rule_n_ops)
      call move_alloc(tmp_i, rule_atom)
      !
      allocate(tmp_i(rule_n_ops))
      if (rule_n_ops > 0) tmp_i = rule_cmp(1:rule_n_ops)
      call move_alloc(tmp_i, rule_cmp)
      !
      allocate(tmp_r(rule_n_ops))
      if (rule_n_ops > 0) tmp_r = rule_threshold(1:rule_n_ops)
      call move_alloc(tmp_r, rule_threshold)
      !
      rule_capacity = rule_n_ops
      !
   endif
   !
   ! Registry.
   !
   if (.not. allocated(rule_start)) then
      !
      allocate(rule_start(0))
      allocate(rule_length(0))
      rule_registry_capacity = 0
      n_rules                = 0
      !
   else if (rule_registry_capacity /= n_rules) then
      !
      allocate(tmp_i(n_rules))
      if (n_rules > 0) tmp_i = rule_start(1:n_rules)
      call move_alloc(tmp_i, rule_start)
      !
      allocate(tmp_i(n_rules))
      if (n_rules > 0) tmp_i = rule_length(1:n_rules)
      call move_alloc(tmp_i, rule_length)
      !
      rule_registry_capacity = n_rules
      !
   endif
   !
   end subroutine
   !
   !
   pure function evaluate_rule(rule_id, z1, z2) result(fired)
   !
   ! Fixed-depth stack machine that evaluates a compiled rule against
   ! the two water levels z1 (intake) and z2 (outfall). A rule_id of 0
   ! short-circuits to .false. ("never fires").
   !
   !$acc routine seq
   !
   implicit none
   !
   integer, intent(in) :: rule_id
   real,    intent(in) :: z1, z2
   logical             :: fired
   !
   logical :: stack(expr_stack_max)
   integer :: sp, k, idx, rs, rl
   real    :: zval
   logical :: a, b
   !
   fired = .false.
   !
   if (rule_id <= 0) return
   !
   rs = rule_start(rule_id)
   rl = rule_length(rule_id)
   !
   if (rl <= 0) return
   !
   sp = 0
   !
   do k = 1, rl
      !
      idx = rs + k - 1
      !
      select case (rule_opcode(idx))
         !
         case (op_cmp)
            !
            select case (rule_atom(idx))
               !
               case (atom_z1)
                  !
                  zval = z1
                  !
               case (atom_z2)
                  !
                  zval = z2
                  !
               case (atom_z2_minus_z1)
                  !
                  zval = z2 - z1
                  !
               case default
                  !
                  zval = 0.0
                  !
            end select
            !
            if (sp >= expr_stack_max) return
            sp = sp + 1
            !
            select case (rule_cmp(idx))
               !
               case (cmp_lt)
                  !
                  stack(sp) = zval < rule_threshold(idx)
                  !
               case (cmp_gt)
                  !
                  stack(sp) = zval > rule_threshold(idx)
                  !
               case (cmp_le)
                  !
                  stack(sp) = zval <= rule_threshold(idx)
                  !
               case (cmp_ge)
                  !
                  stack(sp) = zval >= rule_threshold(idx)
                  !
               case (cmp_eq)
                  !
                  stack(sp) = zval == rule_threshold(idx)
                  !
               case default
                  !
                  stack(sp) = .false.
                  !
            end select
            !
         case (op_and)
            !
            b = stack(sp)
            a = stack(sp - 1)
            sp = sp - 1
            stack(sp) = a .and. b
            !
         case (op_or)
            !
            b = stack(sp)
            a = stack(sp - 1)
            sp = sp - 1
            stack(sp) = a .or. b
            !
      end select
      !
   enddo
   !
   if (sp >= 1) fired = stack(1)
   !
   end function
   !
   !
   subroutine grow_rule_storage(min_capacity)
   !
   ! Ensure rule_capacity >= min_capacity. On first growth, allocates to
   ! max(initial_capacity, min_capacity). On subsequent growth, doubles
   ! until the requested capacity fits. Existing contents are preserved.
   !
   implicit none
   !
   integer, intent(in) :: min_capacity
   !
   integer :: new_capacity
   integer, allocatable :: tmp_i(:)
   real,    allocatable :: tmp_r(:)
   !
   if (.not. allocated(rule_opcode)) then
      !
      new_capacity = max(initial_capacity, min_capacity)
      allocate(rule_opcode   (new_capacity))
      allocate(rule_atom     (new_capacity))
      allocate(rule_cmp      (new_capacity))
      allocate(rule_threshold(new_capacity))
      rule_capacity = new_capacity
      return
      !
   endif
   !
   if (min_capacity <= rule_capacity) return
   !
   new_capacity = max(2 * rule_capacity, min_capacity)
   !
   allocate(tmp_i(new_capacity))
   if (rule_n_ops > 0) tmp_i(1:rule_n_ops) = rule_opcode(1:rule_n_ops)
   call move_alloc(tmp_i, rule_opcode)
   !
   allocate(tmp_i(new_capacity))
   if (rule_n_ops > 0) tmp_i(1:rule_n_ops) = rule_atom(1:rule_n_ops)
   call move_alloc(tmp_i, rule_atom)
   !
   allocate(tmp_i(new_capacity))
   if (rule_n_ops > 0) tmp_i(1:rule_n_ops) = rule_cmp(1:rule_n_ops)
   call move_alloc(tmp_i, rule_cmp)
   !
   allocate(tmp_r(new_capacity))
   if (rule_n_ops > 0) tmp_r(1:rule_n_ops) = rule_threshold(1:rule_n_ops)
   call move_alloc(tmp_r, rule_threshold)
   !
   rule_capacity = new_capacity
   !
   end subroutine
   !
   !
   subroutine grow_rule_registry(min_capacity)
   !
   ! Ensure rule_registry_capacity >= min_capacity. On first growth,
   ! allocates to max(initial_registry_capacity, min_capacity). On
   ! subsequent growth, doubles until the requested capacity fits.
   ! Existing contents are preserved.
   !
   implicit none
   !
   integer, intent(in) :: min_capacity
   !
   integer :: new_capacity
   integer, allocatable :: tmp_i(:)
   !
   if (.not. allocated(rule_start)) then
      !
      new_capacity = max(initial_registry_capacity, min_capacity)
      allocate(rule_start (new_capacity))
      allocate(rule_length(new_capacity))
      rule_registry_capacity = new_capacity
      return
      !
   endif
   !
   if (min_capacity <= rule_registry_capacity) return
   !
   new_capacity = max(2 * rule_registry_capacity, min_capacity)
   !
   allocate(tmp_i(new_capacity))
   if (n_rules > 0) tmp_i(1:n_rules) = rule_start(1:n_rules)
   call move_alloc(tmp_i, rule_start)
   !
   allocate(tmp_i(new_capacity))
   if (n_rules > 0) tmp_i(1:n_rules) = rule_length(1:n_rules)
   call move_alloc(tmp_i, rule_length)
   !
   rule_registry_capacity = new_capacity
   !
   end subroutine
   !
   !
   subroutine parse_rule_expression(src, ops, atoms, cmps, thresholds, nops, ierr, errmsg)
   !
   ! Recursive-descent parser that compiles a rule string to reverse-
   ! polish bytecode in four parallel arrays. op_cmp entries use all
   ! three of atoms / cmps / thresholds; op_and / op_or use only opcode.
   !
   implicit none
   !
   character(len=*), intent(in)  :: src
   integer,          intent(out) :: ops(:)
   integer,          intent(out) :: atoms(:)
   integer,          intent(out) :: cmps(:)
   real,             intent(out) :: thresholds(:)
   integer,          intent(out) :: nops
   integer,          intent(out) :: ierr
   character(len=*), intent(out) :: errmsg
   !
   ! Token kinds:
   !   1 = ident (z1/z2/z2-z1)    payload: atom code in tok_atom
   !   2 = number                 payload: real in tok_num
   !   3 = lparen
   !   4 = rparen
   !   5 = and
   !   6 = or
   !   7 = lt
   !   8 = gt
   !   9 = le
   !  10 = ge
   !  11 = eq
   !
   integer :: tok_kind(expr_tokens_max)
   integer :: tok_atom(expr_tokens_max)
   real    :: tok_num (expr_tokens_max)
   integer :: tok_pos (expr_tokens_max)
   integer :: n_tokens, ip
   !
   nops   = 0
   ierr   = 0
   errmsg = ''
   !
   call tokenize_rule(src, tok_kind, tok_atom, tok_num, tok_pos, n_tokens, ierr, errmsg)
   !
   if (ierr /= 0) return
   !
   if (n_tokens == 0) then
      !
      ierr   = 1
      errmsg = 'empty rule expression'
      return
      !
   endif
   !
   ip = 1
   !
   call parse_or_expr(tok_kind, tok_atom, tok_num, tok_pos, n_tokens, ip, &
        ops, atoms, cmps, thresholds, nops, ierr, errmsg)
   !
   if (ierr /= 0) return
   !
   if (ip <= n_tokens) then
      !
      ierr = 1
      write(errmsg,'(a,i0)') 'unexpected trailing token at position ', tok_pos(ip)
      return
      !
   endif
   !
   end subroutine
   !
   !
   subroutine tokenize_rule(src, tok_kind, tok_atom, tok_num, tok_pos, n_tokens, ierr, errmsg)
   !
   ! One-pass tokenizer. Emits fixed-size parallel arrays: kind, atom
   ! code, number value, source position. Token kinds match the
   ! parameters above in parse_rule_expression.
   !
   implicit none
   !
   character(len=*), intent(in)    :: src
   integer,          intent(out)   :: tok_kind(:)
   integer,          intent(out)   :: tok_atom(:)
   real,             intent(out)   :: tok_num(:)
   integer,          intent(out)   :: tok_pos(:)
   integer,          intent(out)   :: n_tokens
   integer,          intent(inout) :: ierr
   character(len=*), intent(inout) :: errmsg
   !
   integer, parameter :: tok_ident  = 1
   integer, parameter :: tok_number = 2
   integer, parameter :: tok_lparen = 3
   integer, parameter :: tok_rparen = 4
   integer, parameter :: tok_and    = 5
   integer, parameter :: tok_or     = 6
   integer, parameter :: tok_lt     = 7
   integer, parameter :: tok_gt     = 8
   integer, parameter :: tok_le     = 9
   integer, parameter :: tok_ge     = 10
   integer, parameter :: tok_eq     = 11
   !
   integer :: pos, slen, start, kstart, ic, atom_code, iostat_read
   character(len=:), allocatable :: lower
   character(len=32) :: num_buf
   logical :: matched
   !
   lower    = to_lower_local(src)
   slen     = len(lower)
   pos      = 1
   n_tokens = 0
   !
   do while (pos <= slen)
      !
      ! Skip whitespace.
      !
      ic = iachar(lower(pos:pos))
      !
      if (ic == iachar(' ') .or. ic == 9 .or. ic == 10 .or. ic == 13) then
         !
         pos = pos + 1
         cycle
         !
      endif
      !
      if (n_tokens >= expr_tokens_max) then
         !
         ierr = 1
         write(errmsg,'(a,i0,a)') 'too many tokens (>', expr_tokens_max, ') in rule expression'
         return
         !
      endif
      !
      start = pos
      !
      ! Single-character tokens.
      !
      matched = .true.
      !
      select case (lower(pos:pos))
         !
         case ('(')
            !
            n_tokens = n_tokens + 1
            tok_kind(n_tokens) = tok_lparen
            tok_pos (n_tokens) = start
            pos = pos + 1
            !
         case (')')
            !
            n_tokens = n_tokens + 1
            tok_kind(n_tokens) = tok_rparen
            tok_pos (n_tokens) = start
            pos = pos + 1
            !
         case ('&')
            !
            n_tokens = n_tokens + 1
            tok_kind(n_tokens) = tok_and
            tok_pos (n_tokens) = start
            pos = pos + 1
            !
         case ('|')
            !
            n_tokens = n_tokens + 1
            tok_kind(n_tokens) = tok_or
            tok_pos (n_tokens) = start
            pos = pos + 1
            !
         case ('<')
            !
            n_tokens = n_tokens + 1
            tok_pos (n_tokens) = start
            if (pos + 1 <= slen .and. lower(min(pos+1,slen):min(pos+1,slen)) == '=') then
               !
               tok_kind(n_tokens) = tok_le
               pos = pos + 2
               !
            else
               !
               tok_kind(n_tokens) = tok_lt
               pos = pos + 1
               !
            endif
            !
         case ('>')
            !
            n_tokens = n_tokens + 1
            tok_pos (n_tokens) = start
            if (pos + 1 <= slen .and. lower(min(pos+1,slen):min(pos+1,slen)) == '=') then
               !
               tok_kind(n_tokens) = tok_ge
               pos = pos + 2
               !
            else
               !
               tok_kind(n_tokens) = tok_gt
               pos = pos + 1
               !
            endif
            !
         case ('=')
            !
            n_tokens = n_tokens + 1
            tok_kind(n_tokens) = tok_eq
            tok_pos (n_tokens) = start
            if (pos + 1 <= slen .and. lower(min(pos+1,slen):min(pos+1,slen)) == '=') then
               !
               pos = pos + 2
               !
            else
               !
               pos = pos + 1
               !
            endif
            !
         case default
            !
            matched = .false.
            !
      end select
      !
      if (matched) cycle
      !
      ! Number: optional sign is not part of the grammar (z2-z1 is a
      ! fixed atom, not arithmetic). A leading '-' or '+' is only
      ! treated as a number's sign when the next char is digit or dot.
      !
      ic = iachar(lower(pos:pos))
      !
      if (ic == iachar('-') .or. ic == iachar('+') .or. &
          ic == iachar('.') .or. (ic >= iachar('0') .and. ic <= iachar('9'))) then
         !
         if (lower(pos:pos) == '-' .or. lower(pos:pos) == '+') then
            !
            if (pos + 1 > slen) then
               !
               ierr = 1
               write(errmsg,'(a,i0)') 'trailing sign without digits at position ', pos
               return
               !
            endif
            !
            ic = iachar(lower(pos+1:pos+1))
            !
            if (.not. (ic == iachar('.') .or. (ic >= iachar('0') .and. ic <= iachar('9')))) then
               !
               ierr = 1
               write(errmsg,'(a,a,a,i0)') 'unexpected character "', lower(pos:pos), &
                    '" at position ', pos
               return
               !
            endif
            !
         endif
         !
         kstart = pos
         !
         ! Leading sign.
         !
         if (lower(pos:pos) == '-' .or. lower(pos:pos) == '+') pos = pos + 1
         !
         ! Integer part.
         !
         do while (pos <= slen)
            !
            ic = iachar(lower(pos:pos))
            if (.not. (ic >= iachar('0') .and. ic <= iachar('9'))) exit
            pos = pos + 1
            !
         enddo
         !
         ! Fractional part.
         !
         if (pos <= slen) then
            !
            if (lower(pos:pos) == '.') then
               !
               pos = pos + 1
               !
               do while (pos <= slen)
                  !
                  ic = iachar(lower(pos:pos))
                  if (.not. (ic >= iachar('0') .and. ic <= iachar('9'))) exit
                  pos = pos + 1
                  !
               enddo
               !
            endif
            !
         endif
         !
         ! Exponent.
         !
         if (pos <= slen) then
            !
            if (lower(pos:pos) == 'e') then
               !
               pos = pos + 1
               !
               if (pos <= slen) then
                  !
                  if (lower(pos:pos) == '+' .or. lower(pos:pos) == '-') pos = pos + 1
                  !
               endif
               !
               do while (pos <= slen)
                  !
                  ic = iachar(lower(pos:pos))
                  if (.not. (ic >= iachar('0') .and. ic <= iachar('9'))) exit
                  pos = pos + 1
                  !
               enddo
               !
            endif
            !
         endif
         !
         num_buf = lower(kstart:pos-1)
         !
         read(num_buf, *, iostat=iostat_read) tok_num(n_tokens + 1)
         !
         if (iostat_read /= 0) then
            !
            ierr = 1
            write(errmsg,'(a,a,a,i0)') 'invalid number "', trim(num_buf), &
                 '" at position ', kstart
            return
            !
         endif
         !
         n_tokens = n_tokens + 1
         tok_kind(n_tokens) = tok_number
         tok_pos (n_tokens) = kstart
         !
         cycle
         !
      endif
      !
      ! Identifiers: z1, z2, z2-z1. The 'z2-z1' atom contains a '-',
      ! which would otherwise be eaten by the number path; we match it
      ! as a longest-match-first prefix here.
      !
      kstart = pos
      !
      select case (lower(pos:pos))
         !
         case ('z')
            !
            if (pos + 4 <= slen) then
               !
               if (lower(pos:pos+4) == 'z2-z1') then
                  !
                  atom_code = atom_z2_minus_z1
                  n_tokens = n_tokens + 1
                  tok_kind(n_tokens) = tok_ident
                  tok_atom(n_tokens) = atom_code
                  tok_pos (n_tokens) = kstart
                  pos = pos + 5
                  cycle
                  !
               endif
               !
            endif
            !
            if (pos + 1 <= slen) then
               !
               if (lower(pos:pos+1) == 'z1') then
                  !
                  atom_code = atom_z1
                  n_tokens = n_tokens + 1
                  tok_kind(n_tokens) = tok_ident
                  tok_atom(n_tokens) = atom_code
                  tok_pos (n_tokens) = kstart
                  pos = pos + 2
                  cycle
                  !
               elseif (lower(pos:pos+1) == 'z2') then
                  !
                  atom_code = atom_z2
                  n_tokens = n_tokens + 1
                  tok_kind(n_tokens) = tok_ident
                  tok_atom(n_tokens) = atom_code
                  tok_pos (n_tokens) = kstart
                  pos = pos + 2
                  cycle
                  !
               endif
               !
            endif
            !
            ierr = 1
            write(errmsg,'(a,i0)') 'unknown z-identifier at position ', pos
            return
            !
         case default
            !
            ierr = 1
            write(errmsg,'(a,a,a,i0)') 'unexpected character "', lower(pos:pos), &
                 '" at position ', pos
            return
            !
      end select
      !
   enddo
   !
   end subroutine
   !
   !
   recursive subroutine parse_or_expr(tok_kind, tok_atom, tok_num, tok_pos, n_tokens, ip, &
                                      ops, atoms, cmps, thresholds, nops, ierr, errmsg)
   !
   ! Parse an or-expression. Emits a parse_and_expr, then while the
   ! next token is 'or', consumes it, parses another and-expression,
   ! and emits op_or.
   !
   implicit none
   !
   integer,          intent(in)    :: tok_kind(:), tok_atom(:), tok_pos(:)
   real,             intent(in)    :: tok_num(:)
   integer,          intent(in)    :: n_tokens
   integer,          intent(inout) :: ip
   integer,          intent(inout) :: ops(:), atoms(:), cmps(:)
   real,             intent(inout) :: thresholds(:)
   integer,          intent(inout) :: nops
   integer,          intent(inout) :: ierr
   character(len=*), intent(inout) :: errmsg
   !
   integer, parameter :: tok_or = 6
   !
   call parse_and_expr(tok_kind, tok_atom, tok_num, tok_pos, n_tokens, ip, &
        ops, atoms, cmps, thresholds, nops, ierr, errmsg)
   if (ierr /= 0) return
   !
   do while (ip <= n_tokens)
      !
      if (tok_kind(ip) /= tok_or) exit
      !
      ip = ip + 1
      !
      call parse_and_expr(tok_kind, tok_atom, tok_num, tok_pos, n_tokens, ip, &
           ops, atoms, cmps, thresholds, nops, ierr, errmsg)
      if (ierr /= 0) return
      !
      if (nops >= expr_ops_max_per_rule) then
         !
         ierr = 1
         errmsg = 'rule expression too long (op buffer full)'
         return
         !
      endif
      !
      nops = nops + 1
      ops(nops)        = op_or
      atoms(nops)      = 0
      cmps(nops)       = 0
      thresholds(nops) = 0.0
      !
   enddo
   !
   end subroutine
   !
   !
   recursive subroutine parse_and_expr(tok_kind, tok_atom, tok_num, tok_pos, n_tokens, ip, &
                                       ops, atoms, cmps, thresholds, nops, ierr, errmsg)
   !
   ! Parse an and-expression. Emits a parse_comp, then while the next
   ! token is 'and', consumes it, parses another comp, and emits op_and.
   !
   implicit none
   !
   integer,          intent(in)    :: tok_kind(:), tok_atom(:), tok_pos(:)
   real,             intent(in)    :: tok_num(:)
   integer,          intent(in)    :: n_tokens
   integer,          intent(inout) :: ip
   integer,          intent(inout) :: ops(:), atoms(:), cmps(:)
   real,             intent(inout) :: thresholds(:)
   integer,          intent(inout) :: nops
   integer,          intent(inout) :: ierr
   character(len=*), intent(inout) :: errmsg
   !
   integer, parameter :: tok_and = 5
   !
   call parse_comp(tok_kind, tok_atom, tok_num, tok_pos, n_tokens, ip, &
        ops, atoms, cmps, thresholds, nops, ierr, errmsg)
   if (ierr /= 0) return
   !
   do while (ip <= n_tokens)
      !
      if (tok_kind(ip) /= tok_and) exit
      !
      ip = ip + 1
      !
      call parse_comp(tok_kind, tok_atom, tok_num, tok_pos, n_tokens, ip, &
           ops, atoms, cmps, thresholds, nops, ierr, errmsg)
      if (ierr /= 0) return
      !
      if (nops >= expr_ops_max_per_rule) then
         !
         ierr = 1
         errmsg = 'rule expression too long (op buffer full)'
         return
         !
      endif
      !
      nops = nops + 1
      ops(nops)        = op_and
      atoms(nops)      = 0
      cmps(nops)       = 0
      thresholds(nops) = 0.0
      !
   enddo
   !
   end subroutine
   !
   !
   recursive subroutine parse_comp(tok_kind, tok_atom, tok_num, tok_pos, n_tokens, ip, &
                                   ops, atoms, cmps, thresholds, nops, ierr, errmsg)
   !
   ! Parse either a parenthesised expression or a leaf comparison
   ! "atom <|> number". Emits op_cmp for the leaf case.
   !
   implicit none
   !
   integer,          intent(in)    :: tok_kind(:), tok_atom(:), tok_pos(:)
   real,             intent(in)    :: tok_num(:)
   integer,          intent(in)    :: n_tokens
   integer,          intent(inout) :: ip
   integer,          intent(inout) :: ops(:), atoms(:), cmps(:)
   real,             intent(inout) :: thresholds(:)
   integer,          intent(inout) :: nops
   integer,          intent(inout) :: ierr
   character(len=*), intent(inout) :: errmsg
   !
   integer, parameter :: tok_ident  = 1
   integer, parameter :: tok_number = 2
   integer, parameter :: tok_lparen = 3
   integer, parameter :: tok_rparen = 4
   integer, parameter :: tok_lt     = 7
   integer, parameter :: tok_gt     = 8
   integer, parameter :: tok_le     = 9
   integer, parameter :: tok_ge     = 10
   integer, parameter :: tok_eq     = 11
   !
   integer :: atom_code, cmp_code
   !
   if (ip > n_tokens) then
      !
      ierr   = 1
      errmsg = 'unexpected end of expression'
      return
      !
   endif
   !
   if (tok_kind(ip) == tok_lparen) then
      !
      ip = ip + 1
      !
      call parse_or_expr(tok_kind, tok_atom, tok_num, tok_pos, n_tokens, ip, &
           ops, atoms, cmps, thresholds, nops, ierr, errmsg)
      if (ierr /= 0) return
      !
      if (ip > n_tokens) then
         !
         ierr   = 1
         errmsg = 'missing closing ")"'
         return
         !
      endif
      !
      if (tok_kind(ip) /= tok_rparen) then
         !
         ierr = 1
         write(errmsg,'(a,i0)') 'expected ")" at position ', tok_pos(ip)
         return
         !
      endif
      !
      ip = ip + 1
      return
      !
   endif
   !
   ! Leaf: atom cmp_op number.
   !
   if (tok_kind(ip) /= tok_ident) then
      !
      ierr = 1
      write(errmsg,'(a,i0)') 'expected atom (z1/z2/z2-z1) at position ', tok_pos(ip)
      return
      !
   endif
   !
   atom_code = tok_atom(ip)
   ip = ip + 1
   !
   if (ip > n_tokens) then
      !
      ierr   = 1
      errmsg = 'expected comparator after atom'
      return
      !
   endif
   !
   select case (tok_kind(ip))
      !
      case (tok_lt)
         !
         cmp_code = cmp_lt
         !
      case (tok_gt)
         !
         cmp_code = cmp_gt
         !
      case (tok_le)
         !
         cmp_code = cmp_le
         !
      case (tok_ge)
         !
         cmp_code = cmp_ge
         !
      case (tok_eq)
         !
         cmp_code = cmp_eq
         !
      case default
         !
         ierr = 1
         write(errmsg,'(a,i0)') 'expected comparator ("<", ">", "<=", ">=", "=") at position ', tok_pos(ip)
         return
         !
   end select
   !
   ip = ip + 1
   !
   if (ip > n_tokens) then
      !
      ierr   = 1
      errmsg = 'expected number after comparator'
      return
      !
   endif
   !
   if (tok_kind(ip) /= tok_number) then
      !
      ierr = 1
      write(errmsg,'(a,i0)') 'expected numeric threshold at position ', tok_pos(ip)
      return
      !
   endif
   !
   if (nops >= expr_ops_max_per_rule) then
      !
      ierr = 1
      errmsg = 'rule expression too long (op buffer full)'
      return
      !
   endif
   !
   nops = nops + 1
   ops(nops)        = op_cmp
   atoms(nops)      = atom_code
   cmps(nops)       = cmp_code
   thresholds(nops) = tok_num(ip)
   !
   ip = ip + 1
   !
   end subroutine
   !
   !
   function to_lower_local(str) result(lower)
   !
   ! Return a lowercase copy of str (ASCII only). Local to this module
   ! so rule-parsing doesn't depend on sfincs_src_structures for a case
   ! fold; the trivial duplication is worth the decoupling.
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
end module
