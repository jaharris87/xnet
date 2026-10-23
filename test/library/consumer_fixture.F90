Module consumer_fixture
  ! Shared problem data and text output; this module does not depend on XNet.
  Use, Intrinsic :: iso_fortran_env, Only: real64
  Implicit None
  Private
  Public :: fixture_state, write_state, check

Contains

  Subroutine fixture_state(names,density,temperature,ye,x,restart,aux_x,aux_a,aux_z)
    Implicit None
    Character(*), Intent(in) :: names(:)
    Real(real64), Intent(out) :: density(:),temperature(:),ye(:),x(:,:),restart(:)
    Real(real64), Intent(out) :: aux_x(:),aux_a(:),aux_z(:)
    Integer :: i, j

    x = 0.0_real64
    Do i = 1, Size(names)
      If ( Trim(Adjustl(names(i))) == 'c12' .or. Trim(Adjustl(names(i))) == 'o16' ) x(i,:) = 0.475_real64
    EndDo
    Do j = 1, Size(density)
      density(j) = 1.0e7_real64
      temperature(j) = (1.7_real64 + 0.05_real64 * Real(j-1,real64)) * 1.0e9_real64
    EndDo
    aux_x = 0.05_real64
    aux_a = 56.0_real64
    aux_z = 26.0_real64
    ye = 0.475_real64 + 0.05_real64 * 26.0_real64 / 56.0_real64
    restart = 0.0_real64
    ! Exercise the caller's positive seconds hint as well as automatic initial steps.
    restart(1) = 1.0e-12_real64
    Return
  End Subroutine fixture_state

  Subroutine write_state(label,temperature,ye,x,restart,loss,counters)
    Implicit None
    Character(*), Intent(in) :: label
    Real(real64), Intent(in) :: temperature(:),ye(:),x(:,:),restart(:),loss(:)
    Integer, Intent(in) :: counters(:,:)
    Integer :: j, i

    Write(*,'(a,1x,a,2(1x,i0))') 'STATE',Trim(label),Size(x,1),Size(x,2)
    Do j = 1, Size(x,2)
      Write(*,'(a,1x,i0,4(1x,es26.17e3),5(1x,i0))') &
        & 'ZONE',j,temperature(j),ye(j),restart(j),loss(j),counters(:,j)
      Write(*,'(a,1x,i0,*(1x,es26.17e3))') 'X',j,(x(i,j),i=1,Size(x,1))
    EndDo
    Return
  End Subroutine write_state

  Subroutine check(condition,message)
    Implicit None
    Logical, Intent(in) :: condition
    Character(*), Intent(in) :: message

    If ( .not. condition ) Then
      Write(*,'(a)') 'FAILED: '//Trim(message)
      Stop 1
    EndIf
    Return
  End Subroutine check

End Module consumer_fixture
