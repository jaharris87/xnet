Program api_consumer
  ! This external caller imports no XNet implementation module.
  Use, Intrinsic :: ieee_arithmetic, Only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
  Use xnet_api, Only: dp, xnet_controls_t, xnet_network_t, xnet_timing_t, &
    & set_xnet_controls_defaults, xnet_initialize, xnet_burn_batch, xnet_finalize, &
    & xnet_success, xnet_invalid, xnet_lifecycle_error, xnet_unsupported
  Use consumer_fixture, Only: fixture_state, write_state, check
  Implicit None

  Integer, Parameter :: lanes = 4
  Type(xnet_controls_t) :: controls, bad_controls
  Type(xnet_network_t) :: network, rejected_network
  Type(xnet_timing_t) :: timing
  Integer :: ierr, isolv, iheat, i
  Integer :: ids(3,lanes), ids4(4,lanes), counters(5,lanes)
  Character(80) :: data_dir, arg
  Character(512) :: message
  Logical :: active(lanes)
  Real(dp) :: density(lanes),temperature(lanes),ye(lanes),restart(lanes),loss(lanes)
  Real(dp) :: aux_x(lanes),aux_a(lanes),aux_z(lanes)
  Real(dp), Allocatable :: x(:,:), saved_x(:,:)
  Real(dp) :: saved_temperature(lanes),saved_ye(lanes),saved_restart(lanes)
  Real(dp) :: nan_value, inf_value

  Call get_command_argument(1,data_dir)
  Call get_command_argument(2,arg)
  Read(arg,*) isolv
  Call get_command_argument(3,arg)
  Read(arg,*) iheat
  Call check(Len_Trim(data_dir) > 0,'data directory argument')
  Call probe_uninitialized()
  Call xnet_finalize(ierr,message)
  Call expect(xnet_success,'finalize before initialization')

  Call set_xnet_controls_defaults(controls)
  controls%data_dir = data_dir
  controls%nzone = lanes
  controls%nzbatchmx = lanes
  controls%isolv = isolv
  controls%iheat = iheat
  controls%t9nse = 11.0_dp
  controls%idiag = -2
  If ( isolv == 3 ) Then
    controls%kitmx = 10
    controls%iconvc = 3
    controls%ymin = 1.0e-99_dp
  EndIf

  bad_controls = controls
  bad_controls%nzbatchmx = 0
  Call xnet_initialize(bad_controls,rejected_network,ierr,message)
  Call expect(xnet_invalid,'invalid initialization controls')
  bad_controls = controls
  bad_controls%itsout = 1
  Call xnet_initialize(bad_controls,rejected_network,ierr,message)
  Call expect(xnet_unsupported,'unsupported time-series output')
  bad_controls = controls
  bad_controls%iprocess = 1
  Call xnet_initialize(bad_controls,rejected_network,ierr,message)
  Call expect(xnet_unsupported,'unsupported preprocessing')
  bad_controls = controls
  bad_controls%data_dir = 'missing-embedding-test-network'
  Call xnet_initialize(bad_controls,rejected_network,ierr,message)
  Call expect(xnet_invalid,'missing network data')

  Call xnet_initialize(controls,network,ierr,message)
  Call expect(xnet_success,'initialize')
  Call check(Allocated(network%names),'metadata names allocated')
  i = Size(network%names)
  Call check(i > 0,'nonempty network')
  Call check(Size(network%mass_number) == i .and. Size(network%charge_number) == i,'metadata A/Z extent')
  Call check(Size(network%binding_energy_mev) == i .and. Size(network%mass_excess_mev) == i,'metadata energies extent')
  Call check(All(network%mass_number > 0),'metadata mass number')
  Call check(All(network%charge_number >= 0) .and. &
    & All(network%charge_number <= network%mass_number),'metadata charge number')
  Call check(All(ieee_is_finite(network%binding_energy_mev)) .and. &
    & All(ieee_is_finite(network%mass_excess_mev)),'finite metadata energies')
  Do i = 1, Size(network%names)
    Call check(Len_Trim(network%names(i)) > 0,'metadata species name')
    Write(*,'(a,1x,a,4(1x,es26.17e3))') 'META',Trim(Adjustl(network%names(i))), &
      & Real(network%mass_number(i),dp),Real(network%charge_number(i),dp), &
      & network%binding_energy_mev(i),network%mass_excess_mev(i)
  EndDo
  Allocate(x(Size(network%names),lanes),saved_x(Size(network%names),lanes))
  Call fixture_state(network%names,density,temperature,ye,x,restart,aux_x,aux_a,aux_z)
  Do i = 1, lanes
    ids(:,i) = (/ i,7,11 /)
    ids4(:,i) = (/ i,7,11,13 /)
  EndDo
  active = (/ .True.,.False.,.True.,.True. /)

  Call xnet_initialize(controls,rejected_network,ierr,message)
  Call expect(xnet_lifecycle_error,'duplicate initialization')
  Call invalid_requests()
  Call burn_and_check('first',lanes)
  Call burn_and_check('repeat',lanes)
  active = .True.
  Call burn_and_check('short',lanes-1)
  active = .False.
  Call burn_and_check('inactive',lanes)
  Call fixture_state(network%names,density,temperature,ye,x,restart,aux_x,aux_a,aux_z)
  x = x / 0.95_dp
  ye = 0.5_dp
  active = (/ .True.,.False.,.True.,.True. /)
  Call burn_and_check('no_aux',lanes)

  Call xnet_finalize(ierr,message)
  Call expect(xnet_success,'finalize')
  Call xnet_finalize(ierr,message)
  Call expect(xnet_success,'repeated finalize')
  Call probe_uninitialized()

  ! Reinitialize after real burns, then compare a fresh state with the first result.
  Call xnet_initialize(controls,network,ierr,message)
  Call expect(xnet_success,'reinitialize')
  Call fixture_state(network%names,density,temperature,ye,x,restart,aux_x,aux_a,aux_z)
  active = (/ .True.,.False.,.True.,.True. /)
  Call burn_and_check('reinitialize',lanes)
  Call xnet_finalize(ierr,message)
  Call expect(xnet_success,'final finalize')
  Deallocate(x,saved_x)
  Write(*,'(a)') 'API_CHECKS_PASSED'

Contains

  Subroutine expect(expected,label)
    Implicit None
    Integer, Intent(in) :: expected
    Character(*), Intent(in) :: label

    If ( ierr /= expected ) Write(*,'(a,i0,2a)') 'Unexpected status ',ierr,': ',Trim(message)
    Call check(ierr == expected,label)
    If ( expected /= xnet_success ) Call check(Len_Trim(message) > 0,label//' error message')
    Return
  End Subroutine expect

  Subroutine probe_uninitialized()
    Implicit None
    Real(dp) :: d(1), temp(1), ey(1), frac(1,1), hint(1), nu(1)
    Integer :: zone(3,1), count(5,1)
    Logical :: mask(1)

    d = 1.0e7_dp
    temp = 1.7e9_dp
    ey = 0.5_dp
    frac = 1.0_dp
    hint = 0.0_dp
    mask = .True.
    zone = 1
    Call xnet_burn_batch(0.1_dp,mask,zone,d,temp,ey,frac,hint,nu,count,timing,ierr,message)
    Call expect(xnet_lifecycle_error,'burn outside initialized lifetime')
    Return
  End Subroutine probe_uninitialized

  Subroutine remember_state()
    Implicit None
    saved_x = x
    saved_temperature = temperature
    saved_ye = ye
    saved_restart = restart
    Return
  End Subroutine remember_state

  Subroutine unchanged_state(label)
    Implicit None
    Character(*), Intent(in) :: label
    Call check(All(x == saved_x) .and. All(temperature == saved_temperature) .and. &
      & All(ye == saved_ye) .and. All(restart == saved_restart),label//' changed caller state')
    Return
  End Subroutine unchanged_state

  Subroutine invalid_requests()
    Implicit None
    Integer :: bad_ids(2,lanes), bad_counters(4,lanes)
    Real(dp) :: value

    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    inf_value = ieee_value(0.0_dp,ieee_positive_inf)
    Call remember_state()
    Call xnet_burn_batch(-0.1_dp,active,ids,density,temperature,ye,x,restart,loss,counters,timing,ierr,message)
    Call expect(xnet_invalid,'negative burn interval')
    Call unchanged_state('negative burn interval')
    Call xnet_burn_batch(nan_value,active,ids,density,temperature,ye,x,restart,loss,counters,timing,ierr,message)
    Call expect(xnet_invalid,'nonfinite burn interval')
    bad_ids = 1
    Call xnet_burn_batch(0.1_dp,active,bad_ids,density,temperature,ye,x,restart,loss,counters,timing,ierr,message)
    Call expect(xnet_invalid,'zone identifier shape')
    Call xnet_burn_batch(0.1_dp,active,ids,density,temperature,ye,x,restart,loss,bad_counters,timing,ierr,message)
    Call expect(xnet_invalid,'counter shape')
    Call xnet_burn_batch(0.1_dp,active,ids,density,temperature,ye,x(:Size(x,1)-1,:),restart, &
      & loss,counters,timing,ierr,message)
    Call expect(xnet_invalid,'composition shape')
    Call xnet_burn_batch(0.1_dp,active,ids,density,temperature,ye,x,restart,loss,counters,timing,ierr,message,aux_x=aux_x)
    Call expect(xnet_invalid,'incomplete auxiliary inputs')
    value = density(1)
    density(1) = inf_value
    Call xnet_burn_batch(0.1_dp,active,ids,density,temperature,ye,x,restart,loss,counters,timing,ierr,message)
    Call expect(xnet_invalid,'infinite density')
    density(1) = value
    temperature(1) = nan_value
    Call xnet_burn_batch(0.1_dp,active,ids,density,temperature,ye,x,restart,loss,counters,timing,ierr,message)
    Call expect(xnet_invalid,'NaN temperature')
    temperature = saved_temperature
    ye(1) = 1.1_dp
    Call xnet_burn_batch(0.1_dp,active,ids,density,temperature,ye,x,restart,loss,counters,timing,ierr,message)
    Call expect(xnet_invalid,'electron fraction range')
    ye = saved_ye
    x(1,1) = -1.0_dp
    Call xnet_burn_batch(0.1_dp,active,ids,density,temperature,ye,x,restart,loss,counters,timing,ierr,message)
    Call expect(xnet_invalid,'negative mass fraction')
    x = saved_x
    x(1,1) = nan_value
    Call xnet_burn_batch(0.1_dp,active,ids,density,temperature,ye,x,restart,loss,counters,timing,ierr,message)
    Call expect(xnet_invalid,'nonfinite mass fraction')
    x = saved_x
    x(:,1) = 0.0_dp
    Call xnet_burn_batch(0.1_dp,active,ids,density,temperature,ye,x,restart,loss,counters,timing,ierr,message)
    Call expect(xnet_invalid,'zero total composition')
    x = saved_x
    restart(1) = -1.0_dp
    Call xnet_burn_batch(0.1_dp,active,ids,density,temperature,ye,x,restart,loss,counters,timing,ierr,message)
    Call expect(xnet_invalid,'negative restart interval')
    restart = saved_restart
    aux_a(1) = 0.0_dp
    Call xnet_burn_batch(0.1_dp,active,ids,density,temperature,ye,x,restart,loss,counters,timing,ierr,message, &
      & aux_x,aux_a,aux_z)
    Call expect(xnet_invalid,'auxiliary mass number')
    aux_a(1) = 56.0_dp
    Call unchanged_state('invalid calls')
    Return
  End Subroutine invalid_requests

  Subroutine burn_and_check(label,n)
    Implicit None
    Character(*), Intent(in) :: label
    Integer, Intent(in) :: n
    Integer :: j
    Real(dp) :: times(6)

    Call remember_state()
    loss = -1.0_dp
    counters = -1
    If ( label == 'no_aux' ) Then
      Call xnet_burn_batch(0.1_dp,active(:n),ids4(:,:n),density(:n),temperature(:n),ye(:n),x(:,:n), &
        & restart(:n),loss(:n),counters(:,:n),timing,ierr,message)
    Else
      Call xnet_burn_batch(0.1_dp,active(:n),ids(:,:n),density(:n),temperature(:n),ye(:n),x(:,:n), &
        & restart(:n),loss(:n),counters(:,:n),timing,ierr,message,aux_x(:n),aux_a(:n),aux_z(:n))
    EndIf
    Call expect(xnet_success,label//' burn')
    times = (/ timing%seconds,timing%solve,timing%jacobian,timing%derivative,timing%rates,timing%screening /)
    Call check(All(ieee_is_finite(times)) .and. All(times >= 0.0_dp),label//' timing')
    Write(*,'(a,1x,a,6(1x,es26.17e3))') 'TIMING',Trim(label),times
    Do j = 1, n
      If ( active(j) ) Then
        Call check(All(ieee_is_finite(x(:,j))) .and. All(x(:,j) >= 0.0_dp),label//' composition')
        Call check(ieee_is_finite(temperature(j)) .and. temperature(j) > 0.0_dp,label//' temperature')
        Call check(ieee_is_finite(ye(j)) .and. ye(j) >= 0.0_dp .and. ye(j) <= 1.0_dp,label//' Ye')
        Call check(ieee_is_finite(restart(j)) .and. restart(j) >= 0.0_dp,label//' restart')
        Call check(ieee_is_finite(loss(j)) .and. loss(j) >= 0.0_dp,label//' neutrino loss')
        Call check(counters(1,j) > 0 .and. All(counters(:,j) >= 0),label//' counters')
        If ( iheat == 0 ) Call check(Abs(temperature(j)-saved_temperature(j)) <= &
          & 2.0_dp*Epsilon(1.0_dp)*saved_temperature(j),label//' fixed temperature')
      Else
        Call check(All(x(:,j) == saved_x(:,j)) .and. temperature(j) == saved_temperature(j) .and. &
          & ye(j) == saved_ye(j) .and. restart(j) == saved_restart(j),label//' inactive state')
        Call check(loss(j) == 0.0_dp .and. All(counters(:,j) == 0),label//' inactive outputs')
      EndIf
    EndDo
    If ( n < lanes ) Then
      Call check(All(x(:,n+1:) == saved_x(:,n+1:)) .and. &
        & All(temperature(n+1:) == saved_temperature(n+1:)),label//' beyond short batch')
    EndIf
    If ( .not. Any(active(:n)) ) Call check(All(times == 0.0_dp),label//' empty timing')
    Call write_state(label,temperature(:n),ye(:n),x(:,:n),restart(:n),loss(:n),counters(:,:n))
    Return
  End Subroutine burn_and_check

End Program api_consumer
