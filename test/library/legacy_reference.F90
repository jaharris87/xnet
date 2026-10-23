Program legacy_reference
  ! Independent caller of the legacy setup and evolution path.  Keep the controls and workspace
  ! allocation here so comparisons also exercise initialization by the public library interface.
  Use consumer_fixture, Only: fixture_state, write_state, check
  Use nuclear_data, Only: ny, aa, nname, ia, iz, be, mex, read_nuclear_data
  Use reaction_data, Only: read_reaction_data
  Use xnet_abundances, Only: ystart, yo, y, yt, ydot, xext, aext, zext
  Use xnet_conditions, Only: t, tt, to, tdel, tdel_next, tdel_old, t9, t9t, t9o, t9dot, &
    & rho, rhot, rhoo, ye, yet, yeo, cv, etae, detaedt9, nt, ntt, nto, ints, intso, &
    & nstart, tstart, tstop, tdelstart, t9start, rhostart, yestart, nh, th, t9h, rhoh, yeh
  Use xnet_controls, Only: descript, szone, nzone, zone_id, iweak, iweak0, iscrn, iprocess, &
    & nzevolve, nzbatchmx, nzbatch, szbatch, zb_offset, zb_lo, zb_hi, lzactive, isolv, &
    & kstmx, kitmx, ijac, iconvc, changemx, yacc, tolc, tolm, ymin, tdel_maxmult, kmon, ktot, &
    & iheat, changemxt, tolt9, t9nse, ineutrino, nnucout, idiag, itsout, nnucout_string, &
    & inucout, lun_ts, lun_ev, output_nuc, ev_file_base, bin_file_base, lun_diag, lun_stdout, &
    & inab_file, thermo_file, iaux, myid, nproc, tid, nthread, sweep
  Use xnet_eos, Only: eos_initialize
  Use xnet_evolve, Only: full_net
  Use xnet_flux, Only: flux_init
  Use xnet_gpu, Only: gpu_init
  Use xnet_integrate_be, Only: be_init
  Use xnet_integrate_bdf, Only: bdf_init
  Use xnet_jacobian, Only: read_jacobian_data
  Use xnet_match, Only: read_match_data
  Use xnet_nnu, Only: nnuspec, tmevnu, fluxcms
  Use xnet_nse, Only: nse_initialize
  Use xnet_parallel, Only: parallel_initialize, parallel_myproc, parallel_nprocs
  Use xnet_screening, Only: screening_init
  Use xnet_types, Only: dp
  Implicit None

  Integer, Parameter :: lanes = 4, history_points = 2
  Real(dp), Parameter :: duration = 0.1_dp
  Character(1024) :: data_dir
  Character(80) :: data_desc, argument
  Integer :: ierr, i
  Real(dp) :: density(lanes), temperature(lanes), electron_fraction(lanes), restart(lanes)
  Real(dp) :: auxiliary_x(lanes), auxiliary_a(lanes), auxiliary_z(lanes), loss(lanes), denu(lanes)
  Real(dp), Allocatable :: mass_fraction(:,:)
  Integer :: counters(5,lanes)

  Call check(Command_Argument_Count() == 3,'expected data_dir, isolv, and iheat arguments')
  Call Get_Command_Argument(1,data_dir)
  Call Get_Command_Argument(2,argument)
  Read(argument,*,iostat=ierr) isolv
  Call check(ierr == 0,'invalid integrator argument')
  Call check(isolv == 1 .or. isolv == 3,'integrator must be BE or BDF')
  Call Get_Command_Argument(3,argument)
  Read(argument,*,iostat=ierr) iheat
  Call check(ierr == 0,'invalid self-heating argument')
  Call check(iheat == 0 .or. iheat == 1,'self-heating must be off or on')

  Call parallel_initialize()
  myid = parallel_myproc()
  nproc = parallel_nprocs()
  Call check(nproc == 1,'legacy comparison requires serial execution')
  tid = 1
  nthread = 1
  sweep = 'x'
  Call gpu_init()

  ! These are the compiled defaults, followed by the fixed test's execution and output choices.
  ! Assign them directly rather than sharing apply_xnet_controls with the interface under test.
  descript = 'Legacy constant-condition library reference'
  szone = 1
  nzone = lanes
  iweak0 = 1
  iscrn = 1
  iprocess = 0
  nzbatchmx = lanes
  nzevolve = lanes
  nzbatch = lanes
  szbatch = 1
  zb_offset = 0
  zb_lo = 1
  zb_hi = lanes
  kstmx = 9999
  kitmx = 5
  ijac = 1
  iconvc = 0
  changemx = 1.0e-1_dp
  yacc = 1.0e-7_dp
  tolm = 1.0e-6_dp
  tolc = 1.0e-4_dp
  ymin = 1.0e-30_dp
  tdel_maxmult = 2.0_dp
  changemxt = 1.0e-2_dp
  tolt9 = 1.0e-4_dp
  t9nse = 11.0_dp
  ineutrino = 0
  idiag = -2
  itsout = 0
  nnucout = 0
  nnucout_string = '0'
  ev_file_base = 'ev_'
  bin_file_base = 'ts_'
  lun_diag = lun_stdout
  If ( isolv == 3 ) Then
    kitmx = 10
    iconvc = 3
    ymin = 1.0e-99_dp
    changemx = 1.0e10_dp
    changemxt = 1.0e10_dp
  EndIf

  Allocate(zone_id(3,nzevolve),lzactive(nzevolve),iweak(nzevolve),lun_ev(nzevolve),lun_ts(nzevolve))
  Allocate(kmon(5,nzevolve),ktot(5,nzevolve),inab_file(nzone),thermo_file(nzone))
  Allocate(output_nuc(nnucout),inucout(nnucout))
  zone_id = 1
  lzactive = .false.
  iweak = iweak0
  lun_ev = 0
  lun_ts = 0
  kmon = 0
  ktot = 0
  inab_file = ' '
  thermo_file = ' '

  ! The runner prepares the copied network with the production net_setup executable.
  Call read_nuclear_data(Trim(data_dir),data_desc)
  Do i = 1, ny
    Write(*,'(a,1x,a,4(1x,es26.17e3))') 'META',Trim(Adjustl(nname(i))), &
      & Real(ia(i),dp),Real(iz(i),dp),be(i),mex(i)
  EndDo
  Call read_reaction_data(Trim(data_dir))
  Call read_jacobian_data(Trim(data_dir))
  Call read_match_data(Trim(data_dir))
  Call screening_init()
  Call flux_init()
  Call eos_initialize()
  If ( isolv == 1 ) Call be_init()
  If ( isolv == 3 ) Call bdf_init()
  Call nse_initialize()

  ! Allocate the legacy state as the standalone driver does, with just two history endpoints.
  Allocate(y(ny,nzevolve),yo(ny,nzevolve),yt(ny,nzevolve),ydot(ny,nzevolve),ystart(ny,nzevolve))
  Allocate(xext(nzevolve),aext(nzevolve),zext(nzevolve),iaux(nzevolve))
  y = 0.0_dp ; yo = 0.0_dp ; yt = 0.0_dp ; ydot = 0.0_dp ; ystart = 0.0_dp
  xext = 0.0_dp ; aext = 1.0_dp ; zext = 0.0_dp ; iaux = 0

  Allocate(tdel(nzevolve),tdel_next(nzevolve),tdel_old(nzevolve),t(nzevolve),tt(nzevolve),to(nzevolve), &
    & t9(nzevolve),t9t(nzevolve),t9o(nzevolve),t9dot(nzevolve),rho(nzevolve),rhot(nzevolve),rhoo(nzevolve), &
    & ye(nzevolve),yet(nzevolve),yeo(nzevolve),cv(nzevolve),etae(nzevolve),detaedt9(nzevolve), &
    & nt(nzevolve),ntt(nzevolve),nto(nzevolve),ints(nzevolve),intso(nzevolve))
  tdel = 0.0_dp ; tdel_next = 0.0_dp ; tdel_old = 0.0_dp
  t = 0.0_dp ; tt = 0.0_dp ; to = 0.0_dp
  t9 = 0.0_dp ; t9t = 0.0_dp ; t9o = 0.0_dp ; t9dot = 0.0_dp
  rho = 0.0_dp ; rhot = 0.0_dp ; rhoo = 0.0_dp
  ye = 0.0_dp ; yet = 0.0_dp ; yeo = 0.0_dp
  cv = 0.0_dp ; etae = 0.0_dp ; detaedt9 = 0.0_dp
  nt = 0 ; ntt = 0 ; nto = 0 ; ints = 0 ; intso = 0

  Allocate(nstart(nzevolve),tstart(nzevolve),tstop(nzevolve),tdelstart(nzevolve), &
    & t9start(nzevolve),rhostart(nzevolve),yestart(nzevolve),nh(nzevolve), &
    & th(history_points,nzevolve),t9h(history_points,nzevolve), &
    & rhoh(history_points,nzevolve),yeh(history_points,nzevolve), &
    & tmevnu(history_points,nnuspec,nzevolve),fluxcms(history_points,nnuspec,nzevolve))
  nstart = 0 ; tstart = 0.0_dp ; tstop = 0.0_dp ; tdelstart = 0.0_dp
  t9start = 0.0_dp ; rhostart = 0.0_dp ; yestart = 0.0_dp ; nh = 0
  th = 0.0_dp ; t9h = 0.0_dp ; rhoh = 0.0_dp ; yeh = 0.0_dp
  tmevnu = 0.0_dp ; fluxcms = 0.0_dp

  Allocate(mass_fraction(ny,lanes))
  Call fixture_state(nname(1:ny),density,temperature,electron_fraction,mass_fraction,restart, &
    & auxiliary_x,auxiliary_a,auxiliary_z)
  Call reference_burn('first',lanes,[.true.,.false.,.true.,.true.])
  Call reference_burn('repeat',lanes,[.true.,.false.,.true.,.true.])
  Call reference_burn('short',lanes-1,[.true.,.true.,.true.])
  Call reference_burn('inactive',lanes,[.false.,.false.,.false.,.false.])
  Call fixture_state(nname(1:ny),density,temperature,electron_fraction,mass_fraction,restart, &
    & auxiliary_x,auxiliary_a,auxiliary_z)
  mass_fraction = mass_fraction/0.95_dp
  electron_fraction = 0.5_dp
  auxiliary_x = 0.0_dp
  auxiliary_a = 1.0_dp
  auxiliary_z = 0.0_dp
  Call reference_burn('no_aux',lanes,[.true.,.false.,.true.,.true.])
  Deallocate(mass_fraction)

Contains

  Subroutine reference_burn(label,ncall,active)
    ! Reproduce the application-to-legacy-state conversion independently of the public burn routine.
    Implicit None
    Character(*), Intent(in) :: label
    Integer, Intent(in) :: ncall
    Logical, Intent(in) :: active(ncall)
    Integer :: j, kstep

    nzbatch = ncall
    lzactive = .false.
    lzactive(:ncall) = active
    loss = 0.0_dp
    counters = 0
    If ( Any(active) ) Then
      ! Initialize every storage lane: several legacy setup operations also visit inactive lanes.
      Do j = 1, lanes
        ystart(:,j) = mass_fraction(:,j)/aa
        y(:,j) = ystart(:,j)
        xext(j) = auxiliary_x(j)
        aext(j) = auxiliary_a(j)
        zext(j) = auxiliary_z(j)
        iaux(j) = 1
        zone_id(:,j) = [j,1,1]
        tstart(j) = 0.0_dp
        tstop(j) = duration
        tdelstart(j) = -restart(j)
        t9start(j) = temperature(j)*1.0e-9_dp
        rhostart(j) = density(j)
        yestart(j) = electron_fraction(j)
        nh(j) = history_points
        nstart(j) = 1
        th(:,j) = [0.0_dp,duration]
        t9h(:,j) = t9start(j)
        rhoh(:,j) = rhostart(j)
        yeh(:,j) = yestart(j)
        tdel(j) = 0.0_dp
        nt(j) = nstart(j)
        t(j) = tstart(j)
        t9(j) = t9start(j)
        rho(j) = rhostart(j)
        ye(j) = yestart(j)
      EndDo
      Call full_net(kstep,denu(zb_lo:zb_hi))
      Call check(kstep > 0,'legacy evolution took no steps')
      Do j = 1, ncall
        If ( active(j) ) Then
          Call check(t(j) >= tstop(j),'legacy evolution did not reach stop time')
          If ( iheat > 0 ) temperature(j) = t9(j)*1.0e9_dp
          electron_fraction(j) = ye(j)
          mass_fraction(:,j) = aa*y(:,j)
          restart(j) = tdel_old(j)
          loss(j) = denu(j)
          counters(:,j) = ktot(:,j)
        EndIf
      EndDo
    EndIf
    Call write_state(label,temperature(:ncall),electron_fraction(:ncall),mass_fraction(:,:ncall), &
      & restart(:ncall),loss(:ncall),counters(:,:ncall))
    Return
  End Subroutine reference_burn

End Program legacy_reference
