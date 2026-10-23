!***************************************************************************************************
! xnet_controls.f90
! This file contains modules and subroutines to control the execution of XNet.
!***************************************************************************************************

#include "xnet_macros.fh"

Module xnet_controls
  !-------------------------------------------------------------------------------------------------
  ! This module contains the controls used by XNet and routines to read and validate them.
  !-------------------------------------------------------------------------------------------------
  Use, Intrinsic :: iso_fortran_env, Only: lun_stderr=>error_unit, lun_stdin=>input_unit, &
    & lun_stdout=>output_unit
  Use xnet_types, Only: dp
  Use xnet_util, Only: last_nonblank_index
  Implicit None

  Character(*), Parameter :: standalone_controls_file = 'controls.nml'
  Integer, Parameter :: max_includes_per_file = 16 ! Prevent accidental unbounded include fan-out
  Integer, Parameter :: max_include_depth = 16     ! Bound recursive input and detect runaway nesting
  Integer, Parameter :: controls_path_length = 1024 ! Bounded storage for controls-file paths

  Type :: xnet_controls_t
    ! The initial values below are intentionally invalid.  Use set_xnet_controls_defaults before use.
    ! Problem Description
    Character(80) :: description(3) = ' '

    ! Job Controls
    Integer :: szone = 0      ! Starting zone
    Integer :: nzone = 0      ! Number of zones
    Integer :: iweak0 = 0     ! >0: strong and weak; =0: no weak; <0: weak only
    Integer :: iscrn = 0      ! If =0, screening is ignored
    Integer :: iprocess = 0   ! If >0, process nuclear data at run time

    ! Zone Batching Controls
    Integer :: nzbatchmx = 0  ! Maximum number of zones in a batch

    ! Integration Controls
    Integer :: isolv = 0                ! Integration method (1=BE, 3=BDF)
    Integer :: kstmx = 0                ! Maximum number of timesteps before exit
    Integer :: kitmx = 0                ! Maximum iterations within a timestep
    Integer :: ijac = 0                 ! Jacobian rebuild interval after the first iteration
    Integer :: iconvc = 0               ! Convergence condition (0=mass conservation)
    Real(dp) :: changemx = -1.0_dp      ! Relative abundance change used to choose the timestep
    Real(dp) :: yacc = -1.0_dp          ! Minimum abundance used in timestep determination
    Real(dp) :: tolm = -1.0_dp          ! Maximum network mass error
    Real(dp) :: tolc = -1.0_dp          ! Iterative convergence limit
    Real(dp) :: ymin = -1.0_dp          ! Abundances below this value are set to zero
    Real(dp) :: tdel_maxmult = -1.0_dp  ! Maximum timestep growth factor

    ! Self-heating Controls
    Integer :: iheat = 0                ! If >0, couple the network implicitly to temperature
    Real(dp) :: changemxt = -1.0_dp     ! Relative temperature change used to choose the timestep
    Real(dp) :: tolt9 = -1.0_dp         ! Iterative temperature convergence limit

    ! NSE Initial Conditions Controls
    Real(dp) :: t9nse = -1.0_dp         ! Temperature in GK above which NSE initial conditions are used

    ! Neutrino Controls
    Integer :: ineutrino = 0            ! If >0, include neutrino capture reactions

    ! Output Controls
    Integer :: idiag = 0                ! Diagnostic output level
    Integer :: itsout = 0               ! Per-timestep output level
    Character(80) :: ev_file_base = ' '  ! ASCII output filename base
    Character(80) :: bin_file_base = ' ' ! Binary output filename base
    Integer :: nnucout = -1             ! Number of species in condensed output
    Character(5), Allocatable :: output_nuclei(:)

    ! Input Controls used by the standalone driver
    Character(80) :: data_dir = ' '
    Character(80), Allocatable :: inab_files(:)
    Character(80), Allocatable :: thermo_files(:)
  End Type xnet_controls_t

  ! Problem Description
  Character(80) :: descript(3)

  ! Job Controls
  Integer              :: szone        ! Starting zone
  Integer              :: nzone        ! Number of zones
  Integer, Allocatable :: zone_id(:,:) ! Current zone index map (index, batch)
  Integer, Allocatable :: iweak(:)     ! Per-zone weak-reaction flag
  Integer              :: iweak0       ! Input weak-reaction flag
  Integer              :: iscrn        ! If =0, screening is ignored
  Integer              :: iprocess     ! If >0, process nuclear data at run time

  ! Zone Batching Controls
  Integer :: nzevolve                 ! Number of zones evolved concurrently
  Integer :: nzbatchmx                ! Maximum number of zones in a batch
  Integer :: nzbatch                  ! Active number of zones in a batch
  Integer :: szbatch                  ! Starting zone for a batch
  Integer :: zb_offset                ! Zone-batch offset in nzevolve
  Integer :: zb_lo                    ! Lower zone-batch extent in nzevolve
  Integer :: zb_hi                    ! Upper zone-batch extent in nzevolve
  Logical, Allocatable, Target :: lzactive(:) ! Active-zone mask
  !$omp threadprivate(nzbatch,szbatch,zb_offset,zb_lo,zb_hi)

  ! Integration Controls
  Integer :: isolv                   ! Integration method (1=BE, 3=BDF)
  Integer :: kstmx                   ! Maximum number of timesteps before exit
  Integer :: kitmx                   ! Maximum iterations within a timestep
  Integer :: ijac                    ! Jacobian rebuild interval after the first iteration
  Integer :: iconvc                  ! Convergence condition (0=mass conservation)
  Real(dp) :: changemx               ! Relative abundance change used to choose the timestep
  Real(dp) :: yacc                   ! Minimum abundance used in timestep determination
  Real(dp) :: tolc                   ! Iterative convergence limit
  Real(dp) :: tolm                   ! Maximum network mass error
  Real(dp) :: ymin                   ! Abundances below this value are set to zero
  Real(dp) :: tdel_maxmult           ! Maximum timestep growth factor
  Integer, Allocatable :: kmon(:,:)  ! Solver behavior monitors
  Integer, Allocatable :: ktot(:,:)  ! Cumulative solver behavior monitors
  Real(dp) :: t9min = 0.01_dp        ! Temperature minimum for strong reactions

  ! Self-heating Controls
  Integer :: iheat                   ! If >0, couple the network implicitly to temperature
  Real(dp) :: changemxt              ! Relative temperature change used to choose the timestep
  Real(dp) :: tolt9                  ! Iterative temperature convergence limit

  ! NSE Initial Conditions Controls
  Real(dp) :: t9nse                  ! Temperature in GK above which NSE is used initially

  ! Neutrino Controls
  Integer :: ineutrino               ! If >0, include neutrino capture reactions

  ! Output Controls
  Integer :: nnucout                 ! Number of species in condensed output
  Integer :: idiag                   ! Diagnostic output level
  Integer :: itsout                  ! Per-timestep output level
  Character(4) :: nnucout_string     ! Condensed-output format field
  Integer, Allocatable :: inucout(:) ! Condensed-output species indices
  Integer, Allocatable :: lun_ts(:)  ! Per-zone time-series logical units
  Integer, Allocatable :: lun_ev(:)  ! Per-zone ASCII-output logical units
  Character(5), Allocatable :: output_nuc(:) ! Condensed-output species names
  Character(80) :: ev_file_base, bin_file_base ! Output filename bases
  Integer :: lun_diag                ! Per-thread diagnostic logical unit
  !$omp threadprivate(lun_diag)

  ! Input Controls
  Character(80) :: inab_file_base, thermo_file_base ! Input filename bases
  Character(80), Allocatable :: inab_file(:)        ! Per-zone initial-abundance files
  Character(80), Allocatable :: thermo_file(:)      ! Per-zone thermodynamic histories
  Integer :: lun_th, lun_ab                         ! Input logical units
  Integer, Allocatable :: iaux(:)                   ! Auxiliary-nucleus logical units
  !$omp threadprivate(lun_th,lun_ab)

  ! Execution and runtime bookkeeping
  Integer :: myid, nproc, tid, nthread ! Task and thread identifiers and counts
  Character(Len=1) :: sweep            ! Current hydrodynamic sweep: x, y, or z
  !$omp threadprivate(tid,sweep)

  ! Accelerator kernels use these controls on the device.
  !XDIR XDECLARE_VAR(iheat,iscrn,iconvc,ymin)

Contains

  Subroutine set_xnet_controls_defaults(controls)
    !-----------------------------------------------------------------------------------------------
    ! Set the default XNet controls.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Output variables
    Type(xnet_controls_t), Intent(out) :: controls

    controls = xnet_controls_t()
    Allocate(controls%output_nuclei(0),controls%inab_files(0),controls%thermo_files(0))
    Include 'controls.defaults'

    Return
  End Subroutine set_xnet_controls_defaults

  Subroutine validate_xnet_controls(controls,ierr,message)
    !-----------------------------------------------------------------------------------------------
    ! Check the XNet controls for valid values.
    !-----------------------------------------------------------------------------------------------
    Use, Intrinsic :: ieee_arithmetic, Only: ieee_is_finite
    Implicit None

    ! Input variables
    Type(xnet_controls_t), Intent(in) :: controls

    ! Output variables
    Integer, Intent(out) :: ierr
    Character(*), Intent(out) :: message

    ierr = 1
    message = ' '

    If ( .not. Allocated(controls%inab_files) .or. .not. Allocated(controls%thermo_files) .or. &
      & .not. Allocated(controls%output_nuclei) ) Then
      message = 'XNet control arrays must be allocated'
    ElseIf ( .not. ieee_is_finite(controls%changemx) .or. &
      & .not. ieee_is_finite(controls%yacc) .or. &
      & .not. ieee_is_finite(controls%tolm) .or. &
      & .not. ieee_is_finite(controls%tolc) .or. &
      & .not. ieee_is_finite(controls%ymin) .or. &
      & .not. ieee_is_finite(controls%tdel_maxmult) .or. &
      & .not. ieee_is_finite(controls%changemxt) .or. &
      & .not. ieee_is_finite(controls%tolt9) .or. &
      & .not. ieee_is_finite(controls%t9nse) ) Then
      message = 'XNet real-valued controls must be finite'
    ElseIf ( controls%nzone < 1 ) Then
      message = 'nzone must be positive'
    ElseIf ( controls%nnucout < 0 ) Then
      message = 'nnucout must be nonnegative'
    ElseIf ( controls%nnucout > Size(controls%output_nuclei) ) Then
      message = 'nnucout exceeds the number of supplied output_nuclei'
    ElseIf ( controls%nnucout > 0 .and. &
      & Any(Len_Trim(controls%output_nuclei(:controls%nnucout)) == 0) ) Then
      message = 'one output_nuclei entry is required for each requested output species'
    ElseIf ( controls%szone < 1 .or. controls%szone > controls%nzone ) Then
      message = 'szone must select an existing zone'
    ElseIf ( controls%isolv /= 1 .and. controls%isolv /= 3 ) Then
      message = 'isolv must be 1 (BE) or 3 (BDF)'
    ElseIf ( controls%ijac < 1 ) Then
      message = 'ijac must be positive'
    ElseIf ( controls%nzbatchmx < 1 .or. controls%kitmx < 1 .or. controls%kstmx < 1 ) Then
      message = 'batch, iteration, and step limits must be positive'
    ElseIf ( controls%tdel_maxmult <= 0.0_dp ) Then
      message = 'tdel_maxmult must be positive'
    ElseIf ( controls%changemx <= 0.0_dp .or. controls%yacc < 0.0_dp ) Then
      message = 'changemx must be positive and yacc nonnegative'
    ElseIf ( controls%changemxt <= 0.0_dp .or. controls%tolt9 <= 0.0_dp ) Then
      message = 'self-heating change and convergence limits must be positive'
    ElseIf ( controls%t9nse < 0.0_dp ) Then
      message = 't9nse must be nonnegative'
    ElseIf ( controls%tolc <= 0.0_dp .or. controls%tolm <= 0.0_dp .or. controls%ymin < 0.0_dp ) Then
      message = 'integration tolerances must be positive and ymin nonnegative'
    Else
      ierr = 0
    EndIf

    Return
  End Subroutine validate_xnet_controls

  Subroutine validate_standalone_controls(controls,ierr,message)
    !-----------------------------------------------------------------------------------------------
    ! Check the additional controls required by standalone XNet.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Type(xnet_controls_t), Intent(in) :: controls

    ! Output variables
    Integer, Intent(out) :: ierr
    Character(*), Intent(out) :: message

    ! Local variables
    Integer :: input_count

    Call validate_xnet_controls(controls,ierr,message)
    If ( ierr /= 0 ) Return

    ierr = 1
    If ( Len_Trim(controls%data_dir) == 0 ) Then
      message = 'data_dir is required for the standalone driver'
    ElseIf ( Size(controls%inab_files) == 0 .or. Size(controls%thermo_files) == 0 ) Then
      message = 'inab_files and thermo_files are required for the standalone driver'
    ElseIf ( Max(last_nonblank_index(controls%inab_files), &
      & last_nonblank_index(controls%thermo_files)) > controls%nzone ) Then
      message = 'the number of input file pairs cannot exceed nzone'
    ElseIf ( is_hdf_thermo_file(controls%thermo_files(1)) ) Then
      ! HDF5 input files may contain several zones.
      input_count = Max(last_nonblank_index(controls%inab_files), &
        & last_nonblank_index(controls%thermo_files))
      If ( Any(Len_Trim(controls%inab_files(:input_count)) == 0) .or. &
        & Any(Len_Trim(controls%thermo_files(:input_count)) == 0) ) Then
        message = 'each HDF5 input entry needs both inab_files and thermo_files'
      Else
        ierr = 0
      EndIf
    ElseIf ( Size(controls%inab_files) < controls%nzone .or. &
      & Size(controls%thermo_files) < controls%nzone ) Then
      message = 'one inab_files and thermo_files entry is required per zone'
    ElseIf ( Any(Len_Trim(controls%inab_files(:controls%nzone)) == 0) .or. &
      & Any(Len_Trim(controls%thermo_files(:controls%nzone)) == 0) ) Then
      message = 'one inab_files and thermo_files entry is required per zone'
    Else
      ierr = 0
    EndIf

    Return
  End Subroutine validate_standalone_controls

  Subroutine read_xnet_controls(controls,filename,ierr,message)
    !-----------------------------------------------------------------------------------------------
    ! Read, validate, and broadcast XNet controls.
    !-----------------------------------------------------------------------------------------------
    Use xnet_parallel, Only: parallel_bcast, parallel_IOProcessor
    Use xnet_util, Only: normalize_path
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: filename

    ! Output variables
    Type(xnet_controls_t), Intent(out) :: controls
    Integer, Intent(out) :: ierr
    Character(*), Intent(out) :: message

    ! Local variables
    Character(controls_path_length) :: normalized_filename
    Character(controls_path_length) :: ancestor_paths(max_include_depth)

    Call set_xnet_controls_defaults(controls)
    ierr = 0
    message = ' '
    ancestor_paths = ' '

    If ( parallel_IOProcessor() ) Then
      Call normalize_path(filename,normalized_filename,message)
      If ( Len_Trim(message) == 0 ) Then
        Call read_controls_file(controls,normalized_filename,1,ancestor_paths,message)
      EndIf
      If ( Len_Trim(message) == 0 ) Call normalize_xnet_controls(controls,message)
      If ( Len_Trim(message) == 0 ) Call validate_standalone_controls(controls,ierr,message)
      If ( Len_Trim(message) /= 0 ) ierr = 1
    EndIf

    Call broadcast_xnet_controls(controls)
    Call parallel_bcast(ierr)
    Call parallel_bcast(message)

    Return
  End Subroutine read_xnet_controls

  Recursive Subroutine read_controls_file(controls,filename,depth,ancestor_paths,message)
    !-----------------------------------------------------------------------------------------------
    ! Read one controls file and then its included files in order.
    !-----------------------------------------------------------------------------------------------
    Use xnet_util, Only: normalize_path, string_lc
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: filename
    Integer, Intent(in) :: depth

    ! Input/Output variables
    Type(xnet_controls_t), Intent(inout) :: controls
    Character(controls_path_length), Intent(inout) :: ancestor_paths(:)

    ! Output variables
    Character(*), Intent(out) :: message

    ! Namelist staging variables
    Character(controls_path_length) :: include_files(max_includes_per_file+1)
    Character(80) :: description(3), ev_file_base, bin_file_base, data_dir
    Character(80), Allocatable :: inab_files(:), thermo_files(:)
    Character(5), Allocatable :: output_nuclei(:)
    Integer :: szone, nzone, iweak0, iscrn, iprocess, nzbatchmx
    Integer :: isolv, kstmx, kitmx, ijac, iconvc
    Integer :: iheat, ineutrino, idiag, itsout, nnucout
    Real(dp) :: changemx, yacc, tolm, tolc, ymin, tdel_maxmult
    Real(dp) :: changemxt, tolt9, t9nse
    Namelist /xnet_controls/ description, szone, nzone, iweak0, iscrn, iprocess, nzbatchmx, isolv, &
      & kstmx, kitmx, ijac, iconvc, changemx, yacc, tolm, tolc, ymin, tdel_maxmult, iheat, &
      & changemxt, tolt9, t9nse, ineutrino, idiag, itsout, ev_file_base, bin_file_base, data_dir, &
      & nnucout, output_nuclei, inab_files, thermo_files, include_files

    ! Local variables
    Character(controls_path_length) :: include_path, normalized_filename
    Character(:), Allocatable :: line
    Integer :: lun, ierr, i, staging_size, input_count, output_count, file_size
    Integer :: largest_inab_index, largest_thermo_index, largest_output_index, largest_include_index

    message = ' '
    Call normalize_path(filename,normalized_filename,message)
    If ( Len_Trim(message) /= 0 ) Return

    If ( depth > max_include_depth ) Then
      message = 'Controls include depth limit exceeded: '//Trim(normalized_filename)
      Return
    ElseIf ( depth > 1 ) Then
      If ( Any(ancestor_paths(:depth-1) == normalized_filename) ) Then
        message = 'Controls include cycle: '//Trim(normalized_filename)
        Return
      EndIf
    EndIf
    ancestor_paths(depth) = normalized_filename

    Open(newunit=lun,file=Trim(normalized_filename),status='old',action='read',iostat=ierr)
    If ( ierr /= 0 ) Then
      message = 'Failed to open controls file: '//Trim(normalized_filename)
      Return
    EndIf

    Inquire(unit=lun,size=file_size,iostat=ierr)
    If ( ierr /= 0 ) Then
      message = 'Failed to inspect controls file: '//Trim(normalized_filename)
      Close(lun)
      Return
    EndIf
    ! Use the file size to avoid truncating a long record during the array-index scan.
    Allocate(Character(Max(1,file_size)) :: line,stat=ierr)
    If ( ierr /= 0 ) Then
      message = 'Failed to allocate controls input buffer'
      Close(lun)
      Return
    EndIf

    ! Find the largest indices needed for the namelist arrays.
    largest_inab_index = 0
    largest_thermo_index = 0
    largest_output_index = 0
    largest_include_index = 0
    Do
      Read(lun,'(a)',iostat=ierr) line
      If ( ierr < 0 ) Exit
      If ( ierr > 0 ) Then
        message = 'Failed to inspect controls file: '//Trim(normalized_filename)
        Close(lun)
        Return
      EndIf
      Call string_lc(line)
      Call inspect_array_control(line,'inab_files',largest_inab_index,ierr)
      If ( ierr == 0 ) Call inspect_array_control(line,'thermo_files',largest_thermo_index,ierr)
      If ( ierr == 0 ) Call inspect_array_control(line,'output_nuclei',largest_output_index,ierr)
      If ( ierr == 0 ) Call inspect_array_control(line,'include_files',largest_include_index,ierr)
      If ( ierr /= 0 ) Then
        message = 'Dynamically sized controls require explicit positive indices in '// &
          & Trim(normalized_filename)
        Close(lun)
        Return
      ElseIf ( largest_include_index > max_includes_per_file+1 ) Then
        Write(message,'(a,i0,a)') 'Controls file has more than ',max_includes_per_file, &
          & ' direct includes: '//Trim(normalized_filename)
        Close(lun)
        Return
      EndIf
    EndDo
    Rewind(lun)

    staging_size = Max(1,Size(controls%inab_files),Size(controls%thermo_files), &
      & Size(controls%output_nuclei),largest_inab_index,largest_thermo_index,largest_output_index)
    Allocate(inab_files(staging_size),thermo_files(staging_size), &
      & output_nuclei(staging_size),stat=ierr)
    If ( ierr /= 0 ) Then
      message = 'Failed to allocate controls arrays'
      Close(lun)
      Return
    EndIf

    ! Begin with values read from earlier files.  Namelist input changes only named fields.
    include_files = ' '
    description = controls%description
    szone = controls%szone
    nzone = controls%nzone
    iweak0 = controls%iweak0
    iscrn = controls%iscrn
    iprocess = controls%iprocess
    nzbatchmx = controls%nzbatchmx
    isolv = controls%isolv
    kstmx = controls%kstmx
    kitmx = controls%kitmx
    ijac = controls%ijac
    iconvc = controls%iconvc
    changemx = controls%changemx
    yacc = controls%yacc
    tolm = controls%tolm
    tolc = controls%tolc
    ymin = controls%ymin
    tdel_maxmult = controls%tdel_maxmult
    iheat = controls%iheat
    changemxt = controls%changemxt
    tolt9 = controls%tolt9
    t9nse = controls%t9nse
    ineutrino = controls%ineutrino
    idiag = controls%idiag
    itsout = controls%itsout
    ev_file_base = controls%ev_file_base
    bin_file_base = controls%bin_file_base
    data_dir = controls%data_dir
    nnucout = controls%nnucout
    inab_files = ' '
    thermo_files = ' '
    output_nuclei = ' '
    If ( Size(controls%inab_files) > 0 ) Then
      inab_files(:Size(controls%inab_files)) = controls%inab_files
    EndIf
    If ( Size(controls%thermo_files) > 0 ) Then
      thermo_files(:Size(controls%thermo_files)) = controls%thermo_files
    EndIf
    If ( Size(controls%output_nuclei) > 0 ) Then
      output_nuclei(:Size(controls%output_nuclei)) = controls%output_nuclei
    EndIf

    Read(lun,nml=xnet_controls,iostat=ierr)
    Close(lun)
    If ( ierr /= 0 ) Then
      Write(message,'(a,i0,a)') 'Malformed or unknown xnet_controls namelist, iostat=',ierr, &
        & ' in '//Trim(normalized_filename)
      Return
    ElseIf ( Len_Trim(include_files(max_includes_per_file+1)) /= 0 ) Then
      Write(message,'(a,i0,a)') 'Controls file has more than ',max_includes_per_file, &
        & ' direct includes: '//Trim(normalized_filename)
      Return
    EndIf

    ! Copy the namelist values to the XNet controls.
    controls%description = description
    controls%szone = szone
    controls%nzone = nzone
    controls%iweak0 = iweak0
    controls%iscrn = iscrn
    controls%iprocess = iprocess
    controls%nzbatchmx = nzbatchmx
    controls%isolv = isolv
    controls%kstmx = kstmx
    controls%kitmx = kitmx
    controls%ijac = ijac
    controls%iconvc = iconvc
    controls%changemx = changemx
    controls%yacc = yacc
    controls%tolm = tolm
    controls%tolc = tolc
    controls%ymin = ymin
    controls%tdel_maxmult = tdel_maxmult
    controls%iheat = iheat
    controls%changemxt = changemxt
    controls%tolt9 = tolt9
    controls%t9nse = t9nse
    controls%ineutrino = ineutrino
    controls%idiag = idiag
    controls%itsout = itsout
    controls%ev_file_base = ev_file_base
    controls%bin_file_base = bin_file_base
    controls%data_dir = data_dir
    controls%nnucout = nnucout

    input_count = Max(last_nonblank_index(inab_files),last_nonblank_index(thermo_files))
    output_count = last_nonblank_index(output_nuclei)
    Call resize_input_controls(controls,input_count,message)
    If ( Len_Trim(message) /= 0 ) Return
    Call resize_output_controls(controls,output_count,message)
    If ( Len_Trim(message) /= 0 ) Return
    If ( input_count > 0 ) Then
      controls%inab_files = inab_files(:input_count)
      controls%thermo_files = thermo_files(:input_count)
    EndIf
    If ( output_count > 0 ) controls%output_nuclei = output_nuclei(:output_count)

    ! Later include files override values read earlier.
    Do i = 1, max_includes_per_file
      If ( Len_Trim(include_files(i)) == 0 ) Cycle
      Call resolve_include(normalized_filename,include_files(i),include_path,message)
      If ( Len_Trim(message) /= 0 ) Return
      Call read_controls_file(controls,include_path,depth+1,ancestor_paths,message)
      If ( Len_Trim(message) /= 0 ) Return
    EndDo
    ancestor_paths(depth) = ' '

    Return
  End Subroutine read_controls_file

  Subroutine inspect_array_control(line,name,largest_index,ierr)
    !-----------------------------------------------------------------------------------------------
    ! Find the largest explicit index used for an array control.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: line, name

    ! Input/Output variables
    Integer, Intent(inout) :: largest_index

    ! Output variables
    Integer, Intent(out) :: ierr

    ! Local variables
    Character :: quote
    Integer :: close_parenthesis, name_end, name_position, offset, search_from, index_value
    Integer :: first_nonblank, line_length, line_position
    Logical :: skip_character

    ierr = 0
    search_from = 1
    line_length = Len_Trim(line)
    Do
      name_position = 0
      quote = ' '
      skip_character = .False.

      ! Find the next occurrence outside quoted strings and comments.
      Do line_position = search_from, line_length
        If ( skip_character ) Then
          skip_character = .False.
          Cycle
        ElseIf ( quote /= ' ' ) Then
          If ( line(line_position:line_position) == quote ) Then
            If ( line_position < line_length ) Then
              If ( line(line_position+1:line_position+1) == quote ) Then
                skip_character = .True.
              Else
                quote = ' '
              EndIf
            Else
              quote = ' '
            EndIf
          EndIf
          Cycle
        ElseIf ( line(line_position:line_position) == "'" .or. &
          & line(line_position:line_position) == '"' ) Then
          quote = line(line_position:line_position)
          Cycle
        ElseIf ( line(line_position:line_position) == '!' ) Then
          Return
        EndIf

        name_end = line_position+Len_Trim(name)-1
        If ( name_end > line_length ) Cycle
        If ( line(line_position:name_end) /= Trim(name) ) Cycle
        If ( line_position > 1 ) Then
          If ( Index('abcdefghijklmnopqrstuvwxyz0123456789_', &
            & line(line_position-1:line_position-1)) /= 0 ) Cycle
        EndIf
        If ( name_end < line_length ) Then
          If ( Index('abcdefghijklmnopqrstuvwxyz0123456789_',line(name_end+1:name_end+1)) /= 0 ) Cycle
        EndIf
        name_position = line_position
        Exit
      EndDo
      If ( name_position == 0 ) Return

      first_nonblank = name_position+Len_Trim(name)
      Do offset = first_nonblank, line_length
        first_nonblank = offset
        If ( line(offset:offset) /= ' ' ) Exit
      EndDo
      ! Array controls must use explicit indices so the array can be allocated before the read.
      If ( first_nonblank > line_length ) Then
        ierr = 1
        Return
      ElseIf ( line(first_nonblank:first_nonblank) /= '(' ) Then
        ierr = 1
        Return
      EndIf
      If ( first_nonblank == line_length ) Then
        ierr = 1
        Return
      EndIf
      close_parenthesis = Index(line(first_nonblank+1:line_length),')')
      If ( close_parenthesis == 0 ) Then
        ierr = 1
        Return
      EndIf
      close_parenthesis = first_nonblank+close_parenthesis
      Read(line(first_nonblank+1:close_parenthesis-1),*,iostat=ierr) index_value
      If ( ierr /= 0 .or. index_value < 1 ) Then
        ierr = 1
        Return
      EndIf
      largest_index = Max(largest_index,index_value)
      search_from = close_parenthesis+1
      If ( search_from > line_length ) Return
    EndDo

    Return
  End Subroutine inspect_array_control

  Subroutine normalize_xnet_controls(controls,message)
    !-----------------------------------------------------------------------------------------------
    ! Expand standalone input filenames as needed for the selected zones.
    !-----------------------------------------------------------------------------------------------
    Use xnet_util, Only: name_ordered
    Implicit None

    ! Input/Output variables
    Type(xnet_controls_t), Intent(inout) :: controls

    ! Output variables
    Character(*), Intent(out) :: message

    ! Local variables
    Character(80) :: inab_file_base, thermo_file_base
    Integer :: izone, last_inab_index, last_input_index, last_thermo_index
    Logical :: hdf_input

    message = ' '
    If ( controls%nzone < 1 ) Return

    ! Reject a file list longer than the number of zones before resizing it.
    last_inab_index = last_nonblank_index(controls%inab_files)
    last_thermo_index = last_nonblank_index(controls%thermo_files)
    If ( Max(last_inab_index,last_thermo_index) > controls%nzone ) Then
      message = 'the number of input file pairs cannot exceed nzone'
      Return
    EndIf

    Call resize_input_controls(controls,controls%nzone,message)
    If ( Len_Trim(message) /= 0 ) Return

    ! Check that each entry contains a complete file pair.
    last_input_index = 0
    Do izone = 1, controls%nzone
      If ( Len_Trim(controls%inab_files(izone)) == 0 .neqv. &
        & Len_Trim(controls%thermo_files(izone)) == 0 ) Then
        message = 'each input file entry needs both inab_files and thermo_files'
        Return
      EndIf
      If ( Len_Trim(controls%inab_files(izone)) /= 0 ) last_input_index = izone
    EndDo

    hdf_input = is_hdf_thermo_file(controls%thermo_files(1))
    ! HDF5 files may contain several zones, so do not append zone numbers.
    If ( hdf_input ) Return

    If ( last_input_index == 1 .and. controls%nzone > 1 ) Then
      ! Expand a single ASCII input-file pair to one pair per zone.
      inab_file_base = controls%inab_files(1)
      thermo_file_base = controls%thermo_files(1)
      Do izone = 1, controls%nzone
        controls%inab_files(izone) = Trim(inab_file_base)
        controls%thermo_files(izone) = Trim(thermo_file_base)
        Call name_ordered(controls%inab_files(izone),izone,controls%nzone)
        Call name_ordered(controls%thermo_files(izone),izone,controls%nzone)
      EndDo
    ElseIf ( last_input_index /= controls%nzone ) Then
      message = 'one input file pair or one pair per zone is required'
    EndIf

    Return
  End Subroutine normalize_xnet_controls

  Logical Function is_hdf_thermo_file(filename)
    !-----------------------------------------------------------------------------------------------
    ! Check whether a thermodynamic input filename has an HDF5 extension.
    !-----------------------------------------------------------------------------------------------
    Use xnet_util, Only: string_lc
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: filename

    ! Local variables
    Character(Len(filename)) :: lowercase_filename
    Integer :: filename_length

    lowercase_filename = filename
    Call string_lc(lowercase_filename)
    filename_length = Len_Trim(lowercase_filename)
    is_hdf_thermo_file = .False.
    If ( filename_length >= 3 ) Then
      is_hdf_thermo_file = lowercase_filename(filename_length-2:filename_length) == '.h5'
    EndIf
    If ( .not. is_hdf_thermo_file .and. filename_length >= 4 ) Then
      is_hdf_thermo_file = lowercase_filename(filename_length-3:filename_length) == '.hdf'
    EndIf
    If ( .not. is_hdf_thermo_file .and. filename_length >= 5 ) Then
      is_hdf_thermo_file = lowercase_filename(filename_length-4:filename_length) == '.hdf5'
    EndIf

    Return
  End Function is_hdf_thermo_file

  Subroutine resolve_include(parent,child,path,message)
    !-----------------------------------------------------------------------------------------------
    ! Resolve an included controls filename relative to its parent file.
    !-----------------------------------------------------------------------------------------------
    Use xnet_util, Only: normalize_path
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: parent, child

    ! Output variables
    Character(controls_path_length), Intent(out) :: path
    Character(*), Intent(out) :: message

    ! Local variables
    Character(2*controls_path_length+1) :: combined_path
    Integer :: slash

    combined_path = ' '
    If ( Len_Trim(child) == 0 ) Then
      message = 'Controls include path is empty'
      Return
    ElseIf ( child(1:1) == '/' ) Then
      combined_path = Trim(child)
    Else
      slash = Scan(Trim(parent),'/',back=.True.)
      If ( slash > 0 ) Then
        combined_path = parent(:slash)//Trim(child)
      Else
        combined_path = Trim(child)
      EndIf
    EndIf

    Call normalize_path(Trim(combined_path),path,message)

    Return
  End Subroutine resolve_include

  Subroutine resize_input_controls(controls,new_size,message)
    !-----------------------------------------------------------------------------------------------
    ! Resize the abundance and thermodynamic filename arrays.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Integer, Intent(in) :: new_size

    ! Input/Output variables
    Type(xnet_controls_t), Intent(inout) :: controls

    ! Output variables
    Character(*), Intent(out) :: message

    ! Local variables
    Character(80), Allocatable :: new_inab_files(:), new_thermo_files(:)
    Integer :: copy_size, ierr

    message = ' '
    If ( new_size < 0 ) Then
      message = 'Controls input-array size must be nonnegative'
      Return
    ElseIf ( Size(controls%inab_files) == new_size .and. &
      & Size(controls%thermo_files) == new_size ) Then
      Return
    EndIf

    Allocate(new_inab_files(new_size),new_thermo_files(new_size),stat=ierr)
    If ( ierr /= 0 ) Then
      message = 'Failed to allocate standalone input filenames'
      Return
    EndIf
    new_inab_files = ' '
    new_thermo_files = ' '
    copy_size = Min(new_size,Size(controls%inab_files),Size(controls%thermo_files))
    If ( copy_size > 0 ) Then
      new_inab_files(:copy_size) = controls%inab_files(:copy_size)
      new_thermo_files(:copy_size) = controls%thermo_files(:copy_size)
    EndIf
    Call Move_Alloc(new_inab_files,controls%inab_files)
    Call Move_Alloc(new_thermo_files,controls%thermo_files)

    Return
  End Subroutine resize_input_controls

  Subroutine resize_output_controls(controls,new_size,message)
    !-----------------------------------------------------------------------------------------------
    ! Resize the output-nuclei array.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Integer, Intent(in) :: new_size

    ! Input/Output variables
    Type(xnet_controls_t), Intent(inout) :: controls

    ! Output variables
    Character(*), Intent(out) :: message

    ! Local variables
    Character(5), Allocatable :: new_output_nuclei(:)
    Integer :: copy_size, ierr

    message = ' '
    If ( new_size < 0 ) Then
      message = 'Controls output-array size must be nonnegative'
      Return
    ElseIf ( Size(controls%output_nuclei) == new_size ) Then
      Return
    EndIf

    Allocate(new_output_nuclei(new_size),stat=ierr)
    If ( ierr /= 0 ) Then
      message = 'Failed to allocate output_nuclei'
      Return
    EndIf
    new_output_nuclei = ' '
    copy_size = Min(new_size,Size(controls%output_nuclei))
    If ( copy_size > 0 ) new_output_nuclei(:copy_size) = controls%output_nuclei(:copy_size)
    Call Move_Alloc(new_output_nuclei,controls%output_nuclei)

    Return
  End Subroutine resize_output_controls

  Subroutine apply_xnet_controls(controls,ierr,message)
    ! Apply validated controls and allocate reusable execution controls, without standalone filenames.
    Implicit None
    Type(xnet_controls_t), Intent(in) :: controls
    Integer, Intent(out) :: ierr
    Character(*), Intent(out) :: message

    Call validate_xnet_controls(controls,ierr,message)
    If ( ierr /= 0 ) Return
    If ( Allocated(zone_id) ) Then
      ierr = 1
      message = 'XNet execution controls are already allocated'
      Return
    EndIf

    descript = controls%description
    szone = controls%szone
    nzone = controls%nzone
    iweak0 = controls%iweak0
    iscrn = controls%iscrn
    iprocess = controls%iprocess
    nzbatchmx = controls%nzbatchmx
    isolv = controls%isolv
    kstmx = controls%kstmx
    kitmx = controls%kitmx
    ijac = controls%ijac
    iconvc = controls%iconvc
    changemx = controls%changemx
    yacc = controls%yacc
    tolm = controls%tolm
    tolc = controls%tolc
    ymin = controls%ymin
    tdel_maxmult = controls%tdel_maxmult
    iheat = controls%iheat
    changemxt = controls%changemxt
    tolt9 = controls%tolt9
    t9nse = controls%t9nse
    ineutrino = controls%ineutrino
    idiag = controls%idiag
    itsout = controls%itsout
    ev_file_base = controls%ev_file_base
    bin_file_base = controls%bin_file_base
    nnucout = controls%nnucout

    ! The BDF integrator does not need the additional Backward-Euler change limits.
    If ( isolv == 3 ) Then
      changemx = 1.0e10_dp
      changemxt = 1.0e10_dp
    EndIf


    nzevolve = nzbatchmx*nthread
    Allocate(zone_id(3,nzevolve),lzactive(nzevolve),iweak(nzevolve),lun_ev(nzevolve),lun_ts(nzevolve),iaux(nzevolve))
    Allocate(kmon(5,nzevolve),ktot(5,nzevolve),output_nuc(nnucout),inucout(nnucout))
    zone_id = 0
    lzactive = .false.
    iweak = iweak0
    iaux = 0
    lun_ev = 0
    lun_ts = 0
    kmon = 0
    ktot = 0
    inucout = 0
    If ( nnucout > 0 ) output_nuc = controls%output_nuclei(:nnucout)
    Write(nnucout_string,'(i4)') nnucout
    nnucout_string = Adjustl(nnucout_string)

    !$omp parallel default(shared)
    zb_offset = (tid-1)*nzbatchmx
    zb_lo = zb_offset+1
    zb_hi = zb_offset+nzbatchmx
    !$omp end parallel

    !XDIR XUPDATE XASYNC(tid) &
    !XDIR XDEVICE(iheat,iscrn,iconvc,ymin)
    !XDIR XENTER_DATA XASYNC(tid) &
    !XDIR XCREATE(lzactive,iweak,kmon,ktot)
    Return
  End Subroutine apply_xnet_controls

  Subroutine apply_standalone_controls(controls,data_dir)
    ! Add standalone filenames to the common execution controls.
    Use xnet_util, Only: xnet_terminate
    Implicit None
    Type(xnet_controls_t), Intent(in) :: controls
    Character(*), Intent(out) :: data_dir
    Character(256) :: message
    Integer :: ierr

    Call validate_standalone_controls(controls,ierr,message)
    If ( ierr /= 0 ) Call xnet_terminate(Trim(message))
    Call apply_xnet_controls(controls,ierr,message)
    If ( ierr /= 0 ) Call xnet_terminate(Trim(message))
    data_dir = controls%data_dir
    Allocate(inab_file(nzone),thermo_file(nzone))
    inab_file = controls%inab_files(:nzone)
    thermo_file = controls%thermo_files(:nzone)
    Return
  End Subroutine apply_standalone_controls
  Subroutine read_controls(data_dir)
    !-----------------------------------------------------------------------------------------------
    ! Read and apply the standalone XNet controls.
    !-----------------------------------------------------------------------------------------------
    Use xnet_util, Only: xnet_terminate
    Implicit None

    ! Output variables
    Character(80), Intent(out) :: data_dir

    ! Local variables
    Type(xnet_controls_t) :: controls
    Character(256) :: message
    Integer :: ierr

    Call read_xnet_controls(controls,standalone_controls_file,ierr,message)
    If ( ierr /= 0 ) Call xnet_terminate(Trim(message))
    Call apply_standalone_controls(controls,data_dir)

    Return
  End Subroutine read_controls

  Character(160) Function escape_namelist_string(value)
    !-----------------------------------------------------------------------------------------------
    ! Escape apostrophes in a namelist character value.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: value

    ! Local variables
    Integer :: i, output_position

    escape_namelist_string = ' '
    output_position = 1
    Do i = 1, Len_Trim(value)
      If ( value(i:i) == "'" ) Then
        escape_namelist_string(output_position:output_position+1) = "''"
        output_position = output_position + 2
      Else
        escape_namelist_string(output_position:output_position) = value(i:i)
        output_position = output_position + 1
      EndIf
    EndDo

    Return
  End Function escape_namelist_string

  Subroutine write_controls(lun_out,data_dir)
    !-----------------------------------------------------------------------------------------------
    ! Write the current XNet controls as a namelist.
    !-----------------------------------------------------------------------------------------------
    Use xnet_parallel, Only: parallel_IOProcessor
    Implicit None

    ! Input variables
    Integer, Intent(in) :: lun_out
    Character(*), Intent(in) :: data_dir

    ! Local variables
    Integer :: i, izone

    If ( idiag >= 0 .or. ( idiag >= -1 .and. parallel_IOProcessor() ) ) Then
      Write(lun_out,'(a)') '! Resolved/effective XNet controls for this run'
      Write(lun_out,'(a)') '&xnet_controls'
      Do i = 1, 3
        Write(lun_out,'(a,i0,a,a,a)') '  description(',i,") = '", &
          & Trim(escape_namelist_string(descript(i))),"',"
      EndDo
      Write(lun_out,'(a,i0,a)') '  szone = ',szone,','
      Write(lun_out,'(a,i0,a)') '  nzone = ',nzone,','
      Write(lun_out,'(a,i0,a)') '  iweak0 = ',iweak0,','
      Write(lun_out,'(a,i0,a)') '  iscrn = ',iscrn,','
      Write(lun_out,'(a,i0,a)') '  iprocess = ',iprocess,','
      Write(lun_out,'(a,i0,a)') '  nzbatchmx = ',nzbatchmx,','
      Write(lun_out,'(a,i0,a)') '  isolv = ',isolv,','
      Write(lun_out,'(a,i0,a)') '  kstmx = ',kstmx,','
      Write(lun_out,'(a,i0,a)') '  kitmx = ',kitmx,','
      Write(lun_out,'(a,i0,a)') '  ijac = ',ijac,','
      Write(lun_out,'(a,i0,a)') '  iconvc = ',iconvc,','
      Write(lun_out,'(a,es24.16,a)') '  changemx = ',changemx,','
      Write(lun_out,'(a,es24.16,a)') '  yacc = ',yacc,','
      Write(lun_out,'(a,es24.16,a)') '  tolm = ',tolm,','
      Write(lun_out,'(a,es24.16,a)') '  tolc = ',tolc,','
      Write(lun_out,'(a,es24.16,a)') '  ymin = ',ymin,','
      Write(lun_out,'(a,es24.16,a)') '  tdel_maxmult = ',tdel_maxmult,','
      Write(lun_out,'(a,i0,a)') '  iheat = ',iheat,','
      Write(lun_out,'(a,es24.16,a)') '  changemxt = ',changemxt,','
      Write(lun_out,'(a,es24.16,a)') '  tolt9 = ',tolt9,','
      Write(lun_out,'(a,es24.16,a)') '  t9nse = ',t9nse,','
      Write(lun_out,'(a,i0,a)') '  ineutrino = ',ineutrino,','
      Write(lun_out,'(a,i0,a)') '  idiag = ',idiag,','
      Write(lun_out,'(a,i0,a)') '  itsout = ',itsout,','
      Write(lun_out,'(3a)') "  ev_file_base = '",Trim(escape_namelist_string(ev_file_base)),"',"
      Write(lun_out,'(3a)') "  bin_file_base = '",Trim(escape_namelist_string(bin_file_base)),"',"
      Write(lun_out,'(a,i0,a)') '  nnucout = ',nnucout,','
      Do i = 1, nnucout
        Write(lun_out,'(a,i0,a,a,a)') '  output_nuclei(',i,") = '", &
          & Trim(escape_namelist_string(output_nuc(i))),"',"
      EndDo
      Write(lun_out,'(3a)') "  data_dir = '",Trim(escape_namelist_string(data_dir)),"',"
      Do izone = 1, nzone
        Write(lun_out,'(a,i0,a,a,a)') '  inab_files(',izone,") = '", &
          & Trim(escape_namelist_string(inab_file(izone))),"',"
        Write(lun_out,'(a,i0,a,a,a)') '  thermo_files(',izone,") = '", &
          & Trim(escape_namelist_string(thermo_file(izone))),"',"
      EndDo
      Write(lun_out,'(a)') '/'
    EndIf

    Return
  End Subroutine write_controls

  Subroutine broadcast_xnet_controls(controls)
    !-----------------------------------------------------------------------------------------------
    ! Broadcast the XNet controls to all MPI ranks.
    !-----------------------------------------------------------------------------------------------
    Use xnet_parallel, Only: parallel_bcast, parallel_IOProcessor
    Implicit None

    ! Input/Output variables
    Type(xnet_controls_t), Intent(inout) :: controls

    ! Local variables
    Integer :: input_count, output_count

    If ( parallel_IOProcessor() ) Then
      input_count = Size(controls%inab_files)
      output_count = Size(controls%output_nuclei)
    Else
      input_count = 0
      output_count = 0
    EndIf

    Call parallel_bcast(input_count)
    Call parallel_bcast(output_count)

    If ( .not. parallel_IOProcessor() ) Then
      If ( Allocated(controls%inab_files) ) Deallocate(controls%inab_files)
      If ( Allocated(controls%thermo_files) ) Deallocate(controls%thermo_files)
      If ( Allocated(controls%output_nuclei) ) Deallocate(controls%output_nuclei)
      Allocate(controls%inab_files(input_count),controls%thermo_files(input_count))
      Allocate(controls%output_nuclei(output_count))
    EndIf

    Call parallel_bcast(controls%description)
    Call parallel_bcast(controls%szone)
    Call parallel_bcast(controls%nzone)
    Call parallel_bcast(controls%iweak0)
    Call parallel_bcast(controls%iscrn)
    Call parallel_bcast(controls%iprocess)
    Call parallel_bcast(controls%nzbatchmx)
    Call parallel_bcast(controls%isolv)
    Call parallel_bcast(controls%kstmx)
    Call parallel_bcast(controls%kitmx)
    Call parallel_bcast(controls%ijac)
    Call parallel_bcast(controls%iconvc)
    Call parallel_bcast(controls%changemx)
    Call parallel_bcast(controls%yacc)
    Call parallel_bcast(controls%tolm)
    Call parallel_bcast(controls%tolc)
    Call parallel_bcast(controls%ymin)
    Call parallel_bcast(controls%tdel_maxmult)
    Call parallel_bcast(controls%iheat)
    Call parallel_bcast(controls%changemxt)
    Call parallel_bcast(controls%tolt9)
    Call parallel_bcast(controls%t9nse)
    Call parallel_bcast(controls%ineutrino)
    Call parallel_bcast(controls%idiag)
    Call parallel_bcast(controls%itsout)
    Call parallel_bcast(controls%ev_file_base)
    Call parallel_bcast(controls%bin_file_base)
    Call parallel_bcast(controls%data_dir)
    Call parallel_bcast(controls%nnucout)
    If ( output_count > 0 ) Call parallel_bcast(controls%output_nuclei)
    If ( input_count > 0 ) Then
      Call parallel_bcast(controls%inab_files)
      Call parallel_bcast(controls%thermo_files)
    EndIf

    Return
  End Subroutine broadcast_xnet_controls
  Subroutine controls_finalize
    ! Release XNet-owned CPU allocations.
    Implicit None
    If ( Allocated(zone_id) ) Deallocate (zone_id)
    If ( Allocated(iweak) ) Deallocate (iweak)
    If ( Allocated(lzactive) ) Deallocate (lzactive)
    If ( Allocated(kmon) ) Deallocate (kmon)
    If ( Allocated(ktot) ) Deallocate (ktot)
    If ( Allocated(inucout) ) Deallocate (inucout)
    If ( Allocated(lun_ts) ) Deallocate (lun_ts)
    If ( Allocated(lun_ev) ) Deallocate (lun_ev)
    If ( Allocated(output_nuc) ) Deallocate (output_nuc)
    If ( Allocated(inab_file) ) Deallocate (inab_file)
    If ( Allocated(thermo_file) ) Deallocate (thermo_file)
    If ( Allocated(iaux) ) Deallocate (iaux)
    Return
  End Subroutine controls_finalize

End Module xnet_controls
