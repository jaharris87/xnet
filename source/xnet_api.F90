! Supported single-instance constant-condition CPU embedding procedures.
#include "xnet_macros.fh"

Module xnet_api
  Use xnet_types, Only: dp
  Use xnet_controls, Only: xnet_controls_t, set_xnet_controls_defaults
  Implicit None
  Private

  Integer, Parameter, Public :: xnet_success = 0, xnet_invalid = 1
  Integer, Parameter, Public :: xnet_lifecycle_error = 2, xnet_unsupported = 3
  Public :: dp, xnet_controls_t, set_xnet_controls_defaults
  Public :: xnet_initialize, xnet_burn_batch, xnet_finalize

  Type, Public :: xnet_network_t
    Character(5), Allocatable :: names(:)
    Integer, Allocatable :: mass_number(:), charge_number(:)
    Real(dp), Allocatable :: binding_energy_mev(:), mass_excess_mev(:)
  End Type xnet_network_t

  Type, Public :: xnet_timing_t
    Real(dp) :: seconds = 0.0_dp, solve = 0.0_dp, jacobian = 0.0_dp
    Real(dp) :: derivative = 0.0_dp, rates = 0.0_dp, screening = 0.0_dp
  End Type xnet_timing_t

  ! Distinguish this lifecycle from applications using the private module state directly.
  Logical :: api_initialized = .false.

Contains

  Subroutine xnet_initialize(controls,network,ierr,message)
    ! Load a preprocessed network once. Reader, EOS and solver failures may still terminate.
#ifdef XNET_EMBEDDING_CPU
    Use nuclear_data, Only: ny, nname, ia, iz, be, mex, read_nuclear_data
    Use reaction_data, Only: read_reaction_data
    Use xnet_abundances, Only: abundances_initialize, y
    Use xnet_conditions, Only: conditions_initialize
    Use xnet_controls, Only: validate_xnet_controls, apply_xnet_controls, zone_id, nzevolve, &
      & myid, nproc, tid, nthread, lun_diag, lun_stdout, isolv, sweep
    Use xnet_eos, Only: eos_initialize
    Use xnet_flux, Only: flux_init
    Use xnet_integrate_be, Only: be_init
    Use xnet_integrate_bdf, Only: bdf_init
    Use xnet_jacobian, Only: read_jacobian_data
    Use xnet_match, Only: read_match_data
    Use xnet_nnu, Only: tmevnu, fluxcms, nnuspec
    Use xnet_nse, Only: nse_initialize
    Use xnet_screening, Only: screening_init
#endif
    Implicit None
    Type(xnet_controls_t), Intent(in) :: controls
    Type(xnet_network_t), Intent(out) :: network
    Integer, Intent(out) :: ierr
    Character(*), Intent(out) :: message
#ifdef XNET_EMBEDDING_CPU
    Character(80) :: data_description
    Character(10), Parameter :: required_files(5) = &
      & [Character(10) :: 'sunet', 'netwinv', 'nets3', 'nets4', 'match_data']
    Integer :: i, input_unit, input_status
    Logical :: exists

    ierr = xnet_invalid
    message = ' '
    If ( api_initialized .or. Allocated(y) ) Then
      ierr = xnet_lifecycle_error
      message = 'XNet is already initialized'
      Return
    EndIf
    Call validate_xnet_controls(controls,ierr,message)
    If ( ierr /= 0 ) Return
    ierr = xnet_invalid
    If ( Len_Trim(controls%data_dir) == 0 ) Then
      message = 'data_dir is required'
      Return
    EndIf
    If ( controls%iprocess /= 0 .or. controls%ineutrino /= 0 .or. controls%itsout > 0 ) Then
      ierr = xnet_unsupported
      message = 'Embedded CPU calls require preprocessed data, no neutrino capture histories, and no time-series files'
      Return
    EndIf
    ! Check file access before installing state. File contents still use the existing fatal readers.
    Do i = 1, Size(required_files)
      Inquire(file=Trim(controls%data_dir)//'/'//Trim(required_files(i)),exist=exists)
      If ( .not. exists ) Then
        message = 'Missing network file: '//Trim(required_files(i))
        Return
      EndIf
      Open(newunit=input_unit,file=Trim(controls%data_dir)//'/'//Trim(required_files(i)), &
        & status='old',action='read',access='stream',form='unformatted',iostat=input_status)
      If ( input_status /= 0 ) Then
        message = 'Cannot read network file: '//Trim(required_files(i))
        Return
      EndIf
      Close(input_unit)
    EndDo
    myid = 0
    nproc = 1
    tid = 1
    nthread = 1
    sweep = 'x'
    lun_diag = lun_stdout
    Call apply_xnet_controls(controls,ierr,message)
    If ( ierr /= 0 ) Return
    ! Both real consumers carry diagnostic integer tuples, with three or four entries.
    Deallocate (zone_id)
    Allocate (zone_id(4,nzevolve))
    zone_id = 0
    Call read_nuclear_data(controls%data_dir,data_description)
    Call read_reaction_data(controls%data_dir)
    Call read_jacobian_data(controls%data_dir)
    Call read_match_data(controls%data_dir)
    Call screening_init
    Call flux_init
    Call eos_initialize
    If ( isolv == 1 ) Call be_init
    If ( isolv == 3 ) Call bdf_init
    Call nse_initialize
    Call abundances_initialize
    Call conditions_initialize(2)
    Allocate (tmevnu(2,nnuspec,nzevolve),fluxcms(2,nnuspec,nzevolve))
    tmevnu = 0.0_dp
    fluxcms = 0.0_dp

    Allocate (network%names(ny),network%mass_number(ny),network%charge_number(ny), &
      & network%binding_energy_mev(ny),network%mass_excess_mev(ny))
    network%names = nname(1:ny)
    network%mass_number = ia
    network%charge_number = iz
    network%binding_energy_mev = be
    network%mass_excess_mev = mex
    api_initialized = .true.
    ierr = xnet_success
#else
    ierr = xnet_unsupported
    message = 'Embedding currently supports dense serial CPU builds only'
#endif
    Return
  End Subroutine xnet_initialize

  Subroutine xnet_burn_batch(dt_s,active,zone_id,density,temperature,ye,mass_fraction, &
    & restart_dt_s,neutrino_loss,counters,timing,ierr,message,aux_x,aux_a,aux_z)
    ! Constant density, with constant temperature unless self-heating is selected.
#ifdef XNET_EMBEDDING_CPU
    Use, Intrinsic :: ieee_arithmetic, Only: ieee_is_finite
    Use nuclear_data, Only: ny, aa
    Use xnet_abundances, Only: y, ystart, xext, aext, zext
    Use xnet_conditions, Only: t, t9, rho, ye_work=>ye, tdel, tdel_old, nt, nstart, &
      & tstart, tstop, tdelstart, t9start, rhostart, yestart, nh, th, t9h, rhoh, yeh
    Use xnet_controls, Only: nzbatchmx, nzevolve, szbatch, nzbatch, lzactive, iheat, &
      & zone_work=>zone_id, iaux, ktot
    Use xnet_evolve, Only: full_net
#endif
    Implicit None
    Real(dp), Intent(in) :: dt_s, density(:)
    Logical, Intent(in) :: active(:)
    Integer, Intent(in) :: zone_id(:,:)
    Real(dp), Intent(inout) :: temperature(:), ye(:), mass_fraction(:,:), restart_dt_s(:)
    Real(dp), Intent(out) :: neutrino_loss(:)
    Integer, Intent(out) :: counters(:,:)
    Type(xnet_timing_t), Intent(out) :: timing
    Integer, Intent(out) :: ierr
    Character(*), Intent(out) :: message
    Real(dp), Optional, Intent(in) :: aux_x(:), aux_a(:), aux_z(:)
#ifdef XNET_EMBEDDING_CPU
    Integer :: n, lane, kstep
    Real(dp), Allocatable :: loss_work(:)
    Real(dp) :: total_mass
    Type(xnet_timing_t) :: before

    timing = xnet_timing_t()
    neutrino_loss = 0.0_dp
    counters = 0
    message = ' '
    ierr = xnet_lifecycle_error
    If ( .not. api_initialized ) Then
      message = 'Initialize XNet before burning'
      Return
    EndIf
    ierr = xnet_invalid
    n = Size(active)
    If ( n > nzbatchmx .or. Size(density) /= n .or. Size(temperature) /= n .or. &
      & Size(ye) /= n .or. Size(restart_dt_s) /= n .or. Size(neutrino_loss) /= n .or. &
      & Size(mass_fraction,1) /= ny .or. Size(mass_fraction,2) /= n .or. &
      & Size(counters,1) /= 5 .or. Size(counters,2) /= n .or. Size(zone_id,2) /= n .or. &
      & (Size(zone_id,1) /= 3 .and. Size(zone_id,1) /= 4) ) Then
      message = 'Inconsistent batch dimensions or batch larger than nzbatchmx'
      Return
    EndIf
    If ( .not. ieee_is_finite(dt_s) ) Then
      message = 'Burn duration must be finite and positive'
      Return
    EndIf
    If ( dt_s <= 0.0_dp ) Then
      message = 'Burn duration must be finite and positive'
      Return
    EndIf
    If ( (Present(aux_x) .neqv. Present(aux_a)) .or. (Present(aux_x) .neqv. Present(aux_z)) ) Then
      message = 'Supply auxiliary mass fraction, A and Z together'
      Return
    EndIf
    If ( Present(aux_x) ) Then
      If ( Size(aux_x) /= n .or. Size(aux_a) /= n .or. Size(aux_z) /= n ) Then
        message = 'Inconsistent auxiliary dimensions'
        Return
      EndIf
    EndIf
    Do lane = 1, n
      If ( .not. active(lane) ) Cycle
      If ( .not. ieee_is_finite(density(lane)) .or. .not. ieee_is_finite(temperature(lane)) .or. &
        & .not. ieee_is_finite(ye(lane)) .or. .not. ieee_is_finite(restart_dt_s(lane)) .or. &
        & .not. All(ieee_is_finite(mass_fraction(:,lane))) ) Then
        message = 'Active-zone inputs must be finite'
        Return
      EndIf
      If ( density(lane) <= 0.0_dp .or. temperature(lane) <= 0.0_dp .or. &
        & ye(lane) < 0.0_dp .or. ye(lane) > 1.0_dp .or. restart_dt_s(lane) < 0.0_dp .or. &
        & Any(mass_fraction(:,lane) < 0.0_dp) ) Then
        message = 'Invalid active-zone density, temperature, Ye, composition or restart timestep'
        Return
      EndIf
      total_mass = Sum(mass_fraction(:,lane))
      If ( Present(aux_x) ) Then
        If ( .not. ieee_is_finite(aux_x(lane)) .or. .not. ieee_is_finite(aux_a(lane)) .or. &
          & .not. ieee_is_finite(aux_z(lane)) ) Then
          message = 'Active-zone auxiliary inputs must be finite'
          Return
        EndIf
        If ( aux_x(lane) < 0.0_dp .or. aux_x(lane) > 1.0_dp .or. aux_a(lane) < 0.0_dp .or. &
          & aux_z(lane) < 0.0_dp .or. aux_z(lane) > aux_a(lane) ) Then
          message = 'Invalid auxiliary mass fraction, A or Z'
          Return
        EndIf
        If ( aux_x(lane) > 0.0_dp .and. aux_a(lane) <= 0.0_dp ) Then
          message = 'Positive auxiliary mass requires positive A'
          Return
        EndIf
        total_mass = total_mass + aux_x(lane)
      EndIf
      If ( total_mass <= 0.0_dp ) Then
        message = 'Active-zone composition must contain positive mass'
        Return
      EndIf
    EndDo
    ierr = xnet_success
    If ( .not. Any(active) ) Return

    before = current_timing()
    szbatch = 1
    nzbatch = Count(active)
    lzactive = .false.
    lzactive(1:n) = active
    zone_work = 0
    zone_work(1:Size(zone_id,1),1:n) = zone_id
    ystart = 0.0_dp
    xext = 0.0_dp
    aext = 1.0_dp
    zext = 0.0_dp
    iaux = 0
    ! The first endpoint is index one; zero can reach history element zero in timestep.
    nstart = 1
    tstart = 0.0_dp
    tstop = dt_s
    tdelstart = 0.0_dp
    t9start = 1.0_dp
    rhostart = 1.0_dp
    yestart = 0.5_dp
    Do lane = 1, n
      If ( .not. active(lane) ) Cycle
      ystart(:,lane) = mass_fraction(:,lane) / aa
      t9start(lane) = temperature(lane)*1.0e-9_dp
      rhostart(lane) = density(lane)
      yestart(lane) = ye(lane)
      ! Preserve the CHIMERA startup hint convention without exposing its sign encoding.
      tdelstart(lane) = -restart_dt_s(lane)
      If ( Present(aux_x) ) Then
        If ( aux_x(lane) > 0.0_dp ) Then
          xext(lane) = aux_x(lane)
          aext(lane) = aux_a(lane)
          zext(lane) = aux_z(lane)
          iaux(lane) = 1
        EndIf
      EndIf
    EndDo
    nh = 2
    th(1,:) = tstart
    th(2,:) = tstop
    t9h(1,:) = t9start
    t9h(2,:) = t9start
    rhoh(1,:) = rhostart
    rhoh(2,:) = rhostart
    yeh(1,:) = yestart
    yeh(2,:) = yestart
    tdel = 0.0_dp
    nt = nstart
    t = tstart
    t9 = t9start
    rho = rhostart
    ye_work = yestart
    y = ystart
    Allocate (loss_work(nzevolve))
    Call full_net(kstep,loss_work)
    Do lane = 1, n
      If ( .not. active(lane) ) Cycle
      mass_fraction(:,lane) = y(:,lane)*aa
      If ( iheat > 0 ) temperature(lane) = t9(lane)*1.0e9_dp
      ye(lane) = ye_work(lane)
      restart_dt_s(lane) = tdel_old(lane)
      neutrino_loss(lane) = loss_work(lane)
      counters(:,lane) = ktot(:,lane)
    EndDo
    Deallocate (loss_work)
    timing = current_timing()
    timing%seconds = timing%seconds - before%seconds
    timing%solve = timing%solve - before%solve
    timing%jacobian = timing%jacobian - before%jacobian
    timing%derivative = timing%derivative - before%derivative
    timing%rates = timing%rates - before%rates
    timing%screening = timing%screening - before%screening
#else
    timing = xnet_timing_t()
    neutrino_loss = 0.0_dp
    counters = 0
    ierr = xnet_unsupported
    message = 'Embedding currently supports dense serial CPU builds only'
#endif
    Return
  End Subroutine xnet_burn_batch

  Subroutine xnet_finalize(ierr,message)
    ! CPU cleanup never finalizes caller MPI or accelerator state.
#ifdef XNET_EMBEDDING_CPU
    Use nuclear_data, Only: nuclear_data_finalize
    Use reaction_data, Only: reaction_data_finalize
    Use xnet_abundances, Only: abundances_finalize
    Use xnet_conditions, Only: conditions_finalize
    Use xnet_controls, Only: controls_finalize
    Use xnet_eos, Only: eos_finalize
    Use xnet_ffn, Only: ffn_finalize
    Use xnet_flux, Only: flux_finalize
    Use xnet_integrate_be, Only: be_finalize
    Use xnet_integrate_bdf, Only: bdf_finalize
    Use xnet_jacobian, Only: jacobian_finalize
    Use xnet_match, Only: match_finalize
    Use xnet_nnu, Only: nnu_finalize
    Use xnet_nse, Only: nse_finalize
    Use xnet_screening, Only: screening_finalize
#endif
    Implicit None
    Integer, Intent(out) :: ierr
    Character(*), Intent(out) :: message
    ierr = xnet_success
    message = ' '
#ifdef XNET_EMBEDDING_CPU
    If ( .not. api_initialized ) Return
    Call be_finalize
    Call bdf_finalize
    Call nse_finalize
    Call eos_finalize
    Call jacobian_finalize
    Call screening_finalize
    Call flux_finalize
    Call match_finalize
    Call ffn_finalize
    Call nnu_finalize
    Call reaction_data_finalize
    Call nuclear_data_finalize
    Call abundances_finalize
    Call conditions_finalize
    Call controls_finalize
    api_initialized = .false.
#else
    ierr = xnet_unsupported
    message = 'Embedding currently supports dense serial CPU builds only'
#endif
    Return
  End Subroutine xnet_finalize

#ifdef XNET_EMBEDDING_CPU
  Function current_timing() Result(timing)
    Use xnet_timers, Only: timer_xnet, timer_solve, timer_jacob, timer_deriv, timer_csect, timer_scrn
    Implicit None
    Type(xnet_timing_t) :: timing
    timing%seconds = timer_xnet
    timing%solve = timer_solve
    timing%jacobian = timer_jacob
    timing%derivative = timer_deriv
    timing%rates = timer_csect
    timing%screening = timer_scrn
    Return
  End Function current_timing
#endif

End Module xnet_api
