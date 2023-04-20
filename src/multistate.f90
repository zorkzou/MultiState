!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! MultiState@GWEV (MS@GWEV): an interface to compute mixed-spin state via the "External" keyword of Gaussian.
!
! Wenli Zou,  Email: qcband@gmail.com
! Institute of Modern Physics, Northwest University, Xi'an, China
!
! Dec. 01, 2022.
!
! Mar. 21, 2023. Bug fix. The default chi was in au by mistake.
! Apr. 20, 2023. Calculation with -nst 1 may be performed; Bug fix for ONIOM by G16.c.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
program MultiState
  implicit real(kind=8) (a-h,o-z)
  integer,parameter :: iGIN=41, iGOU=42, iCTP=50, MaxStat=9
  character*200 :: ctmp, tag
  logical :: ifchs
  allocatable   :: iINP(:), iFCH(:), SO_chs(:), SO_dlt(:), IZAG(:), XYZG(:), Llist(:)

  allocate(iINP(MaxStat), iFCH(MaxStat), SO_chs(MaxStat*(MaxStat-1)/2), SO_dlt(MaxStat))

!---------------------------------------------------------------------------------------------------------------------------------
! 1. set port numbers:
!    iINP(i)   Gaussian input file for states i
!    iFCH(i)   Gaussian fchk file for states i
!---------------------------------------------------------------------------------------------------------------------------------
  do i = 1, MaxStat
    iINP(i) = 60 + i
    iFCH(i) = 80 + i
  end do

!---------------------------------------------------------------------------------------------------------------------------------
! 2. read arguments
!---------------------------------------------------------------------------------------------------------------------------------
  call RdArg(Imode,iGIN,iGOU,iCTP,iINP,iFCH,MaxStat,NStat,ifchs,SO_chi,SO_chs,SO_dlt,ctmp)

!---------------------------------------------------------------------------------------------------------------------------------
! 3. read *.EIn
!---------------------------------------------------------------------------------------------------------------------------------
  call RdEIn_Line1(iGIN,natom,natall,nder)
  allocate(IZAG(natom), XYZG(3*natom), Llist(natall))
  call RdEIn_XYZ(iGIN,natom,natall,IZAG,XYZG,Llist)

!---------------------------------------------------------------------------------------------------------------------------------
! 4.
! Imode = 0: Generate Gaussian's *.gjf files for scalar states
! Imode = 1: Compute energy and energy derivatives of the mixed-spin ground state, and print Gaussian's *.EOu file
!---------------------------------------------------------------------------------------------------------------------------------
  if(Imode == 0) then
    call GenGjf(iCTP,iINP,NStat,natom,nder,IZAG,XYZG,ctmp,tag)
  else
    call GenEou(iGOU,iFCH,NStat,natom,natall,nder,XYZG,ifchs,SO_chi,SO_chs,SO_dlt,Llist,ctmp,tag)
  end if

  deallocate(iINP, iFCH, SO_chs, SO_dlt, IZAG, XYZG, Llist)

end program MultiState

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Generate Gaussian's *.gjf files for scalar states
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine GenGjf(iCTP,iINP,NStat,natom,nder,IZAG,XYZG,ctmp,tag)
  Implicit Real*8(A-H,O-Z)
  dimension :: iINP(NStat), IZAG(natom), XYZG(3*natom)
  character*200 :: ctmp,tag

  do i = 1, NStat
    call WrtGjf(iCTP,iINP(i),i,natom,nder,IZAG,XYZG,ctmp,tag)
  end do

  Return
End Subroutine GenGjf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Compute energy and energy derivatives of the mixed-spin ground state, and print Gaussian's *.EOu file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine GenEou(iGOU,iFCH,NStat,natom,natall,nder,XYZG,ifchs,SO_chi,SO_chs,SO_dlt,Llist,ctmp,tag)
  Implicit Real*8(A-H,O-Z)
  real(kind=8),parameter :: con_cm2au = 219474.63137d0
  dimension :: iFCH(NStat), XYZG(3*natom), SO_chs(*), SO_dlt(NStat), Llist(natall)
  character*200 :: ctmp, tag
  logical :: ifchs, lrot
  allocatable   :: lrot(:), XYZF(:), ENEi(:), GRDi(:,:), FCMi(:,:), HMAT(:,:), DMAT(:,:), EMIX(:), GRDm(:), FCMm(:), Rmat(:),   &
    Scr1(:), Scr2(:)

  if(NStat == 1) then
    write(*,"(1x,29('='),/,' Single Spin State Calculation',/,1x,29('='))")
  else
    write(*,"(1x,32('='),/,' Multiple-State Spin-Mixing Model',/,1x,32('='))")
  end if

!---------------------------------------------------------------------------------------------------------------------------------
! set array length
!---------------------------------------------------------------------------------------------------------------------------------
  ! a special case: nder = 0 if natom = 1
  if(natom == 1) nder = 0
  na3 = 3*natom
  ntt = na3*(na3+1)/2
  nss = na3*na3

  if(nder == 0) then
    lengrd = 1
    lenfcm = 1
    lenftt = 1
    lenrot = 1
  else if(nder == 1) then
    lengrd = na3
    lenfcm = 16
    lenftt = 1
    lenrot = 9
  else if(nder == 2) then
    lengrd = na3
    lenfcm = nss
    lenftt = ntt
    lenrot = 9
  end if
  lscr1= max(lenfcm,NStat*NStat*2)

!---------------------------------------------------------------------------------------------------------------------------------
! read data from the fchk files
!---------------------------------------------------------------------------------------------------------------------------------
  allocate(XYZF(na3), lrot(NStat), ENEi(NStat), GRDi(lengrd,NStat), FCMi(lenftt,NStat), HMAT(NStat,NStat), DMAT(NStat,NStat),  &
    EMIX(NStat), GRDm(lengrd), FCMm(lenftt), Rmat(lenrot), Scr1(lscr1), Scr2(NStat))
  XYZF = 0.0d0
  GRDi = 0.0d0
  GRDm = 0.0d0
  FCMi = 0.0d0
  FCMm = 0.0d0
  lrot =.false.

  ! energies
  do i = 1, NStat
    call RdFchkEne(iFCH(i),ENEi(i),ctmp,tag)
  end do

  ! gradients
  if(nder > 0) then
    do i = 1, NStat
      call RdFchkGRD(iFCH(i),natom,XYZF,GRDi(1,i),ctmp,tag)
      call CheckXYZ(natom,lrot(i),XYZF,XYZG,Scr1)
      if(lrot(1) .neqv. lrot(i)) call XError("Standard orientations do not agree.")
    end do
  end if

  ! hessians
  if(nder > 1) then
    do i = 1, NStat
      call RdFchkFCM(iFCH(i),ntt,FCMi(1,i),ctmp,tag)
    end do
  end if

!---------------------------------------------------------------------------------------------------------------------------------
! Construct SO model Hamiltonian, and do diagonalization.
!---------------------------------------------------------------------------------------------------------------------------------
  write(*,"(' Number of spin states: ',i2)") NStat
  if(NStat > 1) then
    if(ifchs) then
      ii = 0
      do i = 1, NStat-1
        do j = i+1, NStat
          ii = ii + 1
          HMAT(i,j) = SO_chs(ii)
          HMAT(j,i) = SO_chs(ii)
        end do
      end do
      write(*,"(' chi values (in cm^-1):')")
      write(*,"(5(4x,f8.1))") con_cm2au*SO_chs(1:ii)
    else
      HMAT = -SO_chi
      write(*,"(' chi (in cm^-1): ',/,4x,f8.1)") con_cm2au*SO_chi
    end if
  end if

  write(*,"(' Energies of spin states:',/,1x,47('-'),/,' No.',14x,'E(calc.)',12x,'E(shifted)',/,1x,47('-'))")
  do i = 1, NStat
    HMAT(i,i) = ENEi(i) + SO_dlt(i)
    write(*,"(i4,2f22.9)") i, ENEi(i), HMAT(i,i)
  end do
  write(*,"(1x,47('-'))")

  if(NStat == 1) then
    EMIX(1) = ENEi(1)
    DMAT(1,1) = 1.0d0
  else
    call DCopy(NStat*NStat,HMAT,1,DMAT,1)
    call DSYEV('V','L', NStat, DMAT, NStat, EMIX, Scr1, lscr1, ist)
      if(ist /= 0) call XError("Error in sub. main: diagnolization failed.")
    
    write(*,"(' Energies and weights of mixed-spin states:',/,1x,66('-'),/,' No.',16x,'E(mix)',11x,'Weights',/,1x,66('-'))")
    do i = 1, NStat
      write(*,"(i4,f22.9,5x,4(2x,f6.1,'%'))") i, EMIX(i), (DMAT(j,i)*DMAT(j,i)*1.0d2, j=1,min(NStat,4))
      if(NStat > 4) write (*,"(30x,4(2x,f6.1,'%'))")      (DMAT(j,i)*DMAT(j,i)*1.0d2, j=5,NStat)
    end do
    write(*,"(1x,66('-'))")
  end if

!---------------------------------------------------------------------------------------------------------------------------------
! Compute rotation matrix, if necessary
!---------------------------------------------------------------------------------------------------------------------------------
  if(lrot(1)) call qrotmol(0,6,natom,XYZG,XYZF,Rmat,lscr1,Scr1)

!---------------------------------------------------------------------------------------------------------------------------------
! Calculate gradients and hessians of the mixed-spin ground state
!---------------------------------------------------------------------------------------------------------------------------------

  if(nder > 0) then
    if(NStat == 1) then
      GRDm(:) = GRDi(:,1)
    else
      do i = 1, na3
        do j = 1, NStat
          GRDm(i) = GRDm(i) + DMAT(j,1)*DMAT(j,1)*GRDi(i,j)
        end do
      end do
    end if
    if(lrot(1)) call rotvec(natom,-1,Rmat,GRDm,Scr1)
  end if

  if(nder > 1) then
    if(NStat == 1) then
      FCMm(:) = FCMi(:,1)
    else
      ! d^T * H_m,n * d
      do i = 1, ntt
        do j = 1, NStat
          FCMm(i) = FCMm(i) + DMAT(j,1)*DMAT(j,1)*FCMi(i,j)
        end do
      end do

      ! + 2 d^T * H_m * D * q_n
      ii = 0
      Do nu = 1, na3
        ! q --> Scr1
        call qcalc(na3,NStat,nu,GRDi,DMAT,EMIX,Scr1,Scr2)
        ! q' = D * q --> Scr2
        call dgemm('N','N',NStat,1,NStat,1.0d0,DMAT,NStat,Scr1,NStat,0.0d0,Scr2,NStat)
        ! (2 d) .* q' -- > Scr1
        do i = 1, NStat
          Scr1(i) = 2.0d0 * DMAT(i,1) * Scr2(i)
        end do
        Do mu = 1, nu
          ii = ii + 1
          do i = 1, NStat
            FCMm(ii) = FCMm(ii) + Scr1(i) * GRDi(mu,i)
          end do
        end Do
      end Do
    end if

    if(lrot(1)) then
      call LT2Sqr(na3,FCMm,Scr1)
      call rothess(natom,-1,Rmat,Scr1,FCMm)
      call Sqr2LT(na3,Scr1,FCMm)
    end if
  end if

!---------------------------------------------------------------------------------------------------------------------------------
! Print Gaussian's *.EOu file
!---------------------------------------------------------------------------------------------------------------------------------
  if(natom == natall) then
    call WrtEOU(iGOU,nder,na3,ntt,EMIX(1),GRDm,FCMm)
  else   ! ONIOM by g16.c
    call WrtEOU_ONIOM(iGOU,natom,natall,Llist,nder,EMIX(1),GRDm,FCMm)
  end if

  deallocate(XYZF, lrot, ENEi, GRDi, FCMi, HMAT, DMAT, EMIX, GRDm, FCMm, Rmat, Scr1, Scr2)

  Return
End Subroutine GenEou

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! q vector calculation
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine qcalc(na3,nstat,ixyz,GRDi,DMAT,EMIX,Qvec,Scr1)
  Implicit Real*8(A-H,O-Z)
  dimension :: GRDi(na3,nstat), DMAT(nstat,nstat), EMIX(nstat), Qvec(nstat), Scr1(nstat)

  ! H_mu * d --> Scr1
  do i = 1, nstat
    Scr1(i) = GRDi(ixyz,i) * DMAT(i,1)
  end do

  ! (H_mu * d)^T * D --> Qvec
  call dgemm('N','N',1,nstat,nstat,1.0d0,Scr1,1,DMAT,nstat,0.0d0,Qvec,1)

  do i = 1, nstat
    x = EMIX(1) - EMIX(i)
    if(x > 1.0d-5) then
      Qvec(i) = Qvec(i) / x
    else
      Qvec(i) = 0.0d0
    end if
  end do

  Return
End Subroutine qcalc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Write Gaussian's *.EOu file for ONIOM job by G16.c
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine WrtEOU_ONIOM(iGOU,natom,natall,Llist,nder,Energy,GRD,FCM)
  Implicit Real*8(A-H,O-Z)
  dimension :: Llist(natall), GRD(*), FCM(*)
  allocatable   :: GRDall(:), FCMall(:)

  na3 = natall*3
  ntt = na3*(na3+1)/2

  rewind(iGOU)
  if(nder == 0) then
    ! energy, dipole moment
    write(iGOU,"(4d20.12)") Energy, (0.0d0, i=1,3)
    ! gradients
    write(iGOU,"(3d20.12)") (0.0d0, i=1,na3)
    ! polarizability (l.t. part)
    write(iGOU,"(3d20.12)") (0.0d0, i=1,6)
    ! apt
    write(iGOU,"(3d20.12)") (0.0d0, i=1,3*na3)
    ! hessian matrix (l.t. part)
    write(iGOU,"(3d20.12)") (0.0d0, i=1,ntt)

  else if(nder == 1) then
    allocate(GRDall(na3))
    ! energy, dipole moment
    write(iGOU,"(4d20.12)") Energy, (0.0d0, i=1,3)
    ! gradients
    call FillGRD(natom,natall,Llist,GRD,GRDall)
    write(iGOU,"(3d20.12)") (GRDall(i), i=1,na3)
    ! polarizability (l.t. part)
    write(iGOU,"(3d20.12)") (0.0d0, i=1,6)
    ! apt
    write(iGOU,"(3d20.12)") (0.0d0, i=1,3*na3)
    ! hessian matrix (l.t. part)
    write(iGOU,"(3d20.12)") (0.0d0, i=1,ntt)
    deallocate(GRDall)

  else if(nder == 2) then
    allocate(GRDall(na3),FCMall(ntt))
    ! energy, dipole moment
    write(iGOU,"(4d20.12)") Energy, (0.0d0, i=1,3)
    ! gradients
    call FillGRD(natom,natall,Llist,GRD,GRDall)
    write(iGOU,"(3d20.12)") (GRDall(i), i=1,na3)
    ! polarizability (l.t. part)
    write(iGOU,"(3d20.12)") (0.0d0, i=1,6)
    ! apt
    write(iGOU,"(3d20.12)") (0.0d0, i=1,3*na3)
    ! hessian matrix (l.t. part)
    call FillFCM(natom,natall,Llist,FCM,FCMall)
    write(iGOU,"(3d20.12)") (FCMall(i), i=1,ntt)
    deallocate(GRDall,FCMall)

  end if

  Return
End Subroutine WrtEOU_ONIOM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Write Gaussian's *.EOu file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine WrtEOU(iGOU,nder,na3,ntt,Energy,GRD,FCM)
  Implicit Real*8(A-H,O-Z)
  dimension :: GRD(*), FCM(*)

  rewind(iGOU)
  if(nder == 0) then
    ! energy, dipole moment
    write(iGOU,"(4d20.12)") Energy, (0.0d0, i=1,3)
    ! gradients
    write(iGOU,"(3d20.12)") (0.0d0, i=1,na3)
    ! polarizability (l.t. part)
    write(iGOU,"(3d20.12)") (0.0d0, i=1,6)
    ! apt
    write(iGOU,"(3d20.12)") (0.0d0, i=1,3*na3)
    ! hessian matrix (l.t. part)
    write(iGOU,"(3d20.12)") (0.0d0, i=1,ntt)

  else if(nder == 1) then
    ! energy, dipole moment
    write(iGOU,"(4d20.12)") Energy, (0.0d0, i=1,3)
    ! gradients
    write(iGOU,"(3d20.12)") (GRD(i), i=1,na3)
    ! polarizability (l.t. part)
    write(iGOU,"(3d20.12)") (0.0d0, i=1,6)
    ! apt
    write(iGOU,"(3d20.12)") (0.0d0, i=1,3*na3)
    ! hessian matrix (l.t. part)
    write(iGOU,"(3d20.12)") (0.0d0, i=1,ntt)

  else if(nder == 2) then
    ! energy, dipole moment
    write(iGOU,"(4d20.12)") Energy, (0.0d0, i=1,3)
    ! gradients
    write(iGOU,"(3d20.12)") (GRD(i), i=1,na3)
    ! polarizability (l.t. part)
    write(iGOU,"(3d20.12)") (0.0d0, i=1,6)
    ! apt
    write(iGOU,"(3d20.12)") (0.0d0, i=1,3*na3)
    ! hessian matrix (l.t. part)
    write(iGOU,"(3d20.12)") (FCM(i), i=1,ntt)

  end if

  Return
End Subroutine WrtEOU

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! fill the force constant matrix
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine FillFCM(n1,n2,Llist,f1,f2)
  Implicit Real*8(A-H,O-Z)
  dimension :: Llist(n2), f1(*), f2(*)
  allocatable :: fm1(:,:,:,:), fm2(:,:,:,:)

  allocate(fm1(3,n1,3,n1), fm2(3,n2,3,n2))

  call LT2Sqr(n1*3,f1,fm1)

  fm2 = 0.0d0

  i1 = 0
  do i2 = 1, n2
    if(Llist(i2) == 0) cycle
    i1 = i1 + 1

    j1 = 0
    do j2 = 1, n2
      if(Llist(j2) == 0) cycle
      j1 = j1 + 1
      fm2(:,j2,:,i2) = fm1(:,j1,:,i1)
    end do
  end do

  call Sqr2LT(n2*3,fm2,f2)

  deallocate(fm1, fm2)

  Return
End Subroutine FillFCM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! fill the gradient array
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine FillGRD(n1,n2,Llist,g1,g2)
  Implicit Real*8(A-H,O-Z)
  dimension :: Llist(n2), g1(3,n1), g2(3,n2)

  g2 = 0.0d0

  i1 = 0
  do i2 = 1, n2
    if(Llist(i2) == 0) cycle
    i1 = i1 + 1
    g2(:,i2) = g1(:,i1)
  end do

  Return
End Subroutine FillGRD

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! This subroutine determines the rotation matrix for the best fit superimposition of two molecular Coordinates.
!
! The basic method was described in S. K. Kearsley, Acta Cryst. A45, 208 (1989), and coded in the PDBSUP program by B.Rupp and
! S.Parkin at LLNL (1996).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine qrotmol(iprint,iout,natm,cartarget,carprobe,rmat,lscr,scr)
  implicit real(kind=8) (a-h,o-z)
  dimension cartarget(3,natm), carprobe(3,natm), rmat(3,3), scr(lscr)
  allocatable :: dxp(:,:), dxm(:,:), qmat(:,:), eval(:), evec(:)

  allocate(dxp(3,natm), dxm(3,natm), qmat(4,4), eval(4), evec(16))

  ! translate the molecules to their geometric centers
  call shift(natm,scr,cartarget)
  call shift(natm,scr,carprobe)
  if(iprint > 0) then
    write(iout,"(' Cartesian coordinates of target molecule:')")
    do i=1,natm
      write(iout,"(3f20.12)")cartarget(:,i)
    end do
    write(iout,"(' Cartesian coordinates of probe molecule:')")
    do i=1,natm
      write(iout,"(3f20.12)")carprobe(:,i)
    end do
  end if

  ! coordinate differences: plus and minus
  dxp = carprobe + cartarget
  dxm = carprobe - cartarget

  ! construct Kearsley's Q-matrix
  call conqmt(iprint,iout,natm,dxp,dxm,qmat)

  ! solve Q * L = L * A
  call DCopy(16,qmat,1,evec,1)
  call DSYEV('V','L', 4, evec, 4, eval, scr, lscr, ist)
    if(ist /= 0) call XError("Error in sub. qrotmol: diagnolization failed.")

  ! construct the best fit rotation matrix using the eigenvectors
  call conrot(iprint,iout,evec,rmat)

  deallocate(dxp, dxm, qmat, eval, evec)

  return

  contains

  ! shift the coordinates to the geometric center
  Subroutine Shift(n,o,xyz)
    implicit real(kind=8) (a-h,o-z)
    dimension :: xyz(3,*), o(3)

    o = 0.0d0
    do i=1,n
      o = o + xyz(:,i)
    end do
    o = o /dble(n)
    do i=1,n
      xyz(:,i) = xyz(:,i) - o
    end do

    return
  end Subroutine Shift

  ! Construct Kearsley's Q-matrix. See the last equation in
  ! S.K.Kearsley, On the orthogonal transformation used for structural comparisons, Acta Cryst. A45, 208 (1989)
  subroutine conqmt(iprint,iout,n,xp,xm,q)
    implicit real(kind=8) (a-h,o-z)
    dimension         :: xp(3,n), xm(3,n), q(4,4)

    q = 0.0d0

    do i = 1, n
      q(1,1) = q(1,1) + xm(1,i)*xm(1,i) + xm(2,i)*xm(2,i) + xm(3,i)*xm(3,i)
      q(2,2) = q(2,2) + xp(2,i)*xp(2,i) + xp(3,i)*xp(3,i) + xm(1,i)*xm(1,i)
      q(3,3) = q(3,3) + xp(1,i)*xp(1,i) + xp(3,i)*xp(3,i) + xm(2,i)*xm(2,i)
      q(4,4) = q(4,4) + xp(1,i)*xp(1,i) + xp(2,i)*xp(2,i) + xm(3,i)*xm(3,i)

      q(2,1) = q(2,1) + xp(2,i)*xm(3,i) - xm(2,i)*xp(3,i)
      q(3,1) = q(3,1) + xm(1,i)*xp(3,i) - xp(1,i)*xm(3,i)
      q(4,1) = q(4,1) + xp(1,i)*xm(2,i) - xm(1,i)*xp(2,i)
      q(3,2) = q(3,2) + xm(1,i)*xm(2,i) - xp(1,i)*xp(2,i)
      q(4,2) = q(4,2) + xm(1,i)*xm(3,i) - xp(1,i)*xp(3,i)
      q(4,3) = q(4,3) + xm(2,i)*xm(3,i) - xp(2,i)*xp(3,i)

      q(1,2) = q(2,1)
      q(1,3) = q(3,1)
      q(1,4) = q(4,1)
      q(2,3) = q(3,2)
      q(2,4) = q(4,2)
      q(3,4) = q(4,3)
    end do

    if(iprint > 0) then
      write(iout,"(/,' Q-matrix:')")
      write(iout,"(4x,4f16.6)") q
    end if

    return
  end subroutine conqmt

  ! Constructing the best fit rotation matrix. See the 2nd equation in
  ! S.K.Kearsley, On the orthogonal transformation used for structural comparisons, Acta Cryst. A45, 208 (1989)
  subroutine conrot(iprint,iout,q,r)
    implicit real(kind=8) (a-h,o-z)
    dimension         :: q(4), r(3,3)

    r(1,1) = q(1)*q(1) + q(2)*q(2) - q(3)*q(3) - q(4)*q(4)
    r(2,1) = 2 * (q(2)*q(3) + q(1)*q(4))
    r(3,1) = 2 * (q(2)*q(4) - q(1)*q(3))
    r(1,2) = 2 * (q(2)*q(3) - q(1)*q(4))
    r(2,2) = q(1)*q(1) - q(2)*q(2) + q(3)*q(3) - q(4)*q(4)
    r(3,2) = 2 * (q(3)*q(4) + q(1)*q(2))
    r(1,3) = 2 * (q(2)*q(4) + q(1)*q(3))
    r(2,3) = 2 * (q(3)*q(4) - q(1)*q(2))
    r(3,3) = q(1)*q(1) - q(2)*q(2) - q(3)*q(3) + q(4)*q(4)

    if(iprint > 0) then
      write(iout,"(/,' Rotation matrix for the best fit:')")
      do i = 1, 3
        write(iout,"(4x,3f16.6)") r(i,:)
      end do
    end if

    return
  end subroutine conrot

end subroutine qrotmol

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read Hessian matrix from fchk
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine RdFchkFCM(ifchk,ntt,FCM,ctmp,tag)
  Implicit Real*8(A-H,O-Z)
  dimension :: FCM(ntt)
  character*80 :: ctmp,tag

  !rewind(ifchk)
  tag='Cartesian Force Constants                  R   N= '
  do while(.true.)
    read(ifchk,"(a80)",iostat=ist) ctmp
      if(ist /= 0) call XError("No Hessian matrix found in the fchk file.")
    if(index(ctmp,tag(1:50)) > 0) then
      read(ifchk,"(5e16.8)")(FCM(i),i=1,NTT)
      exit
    end if
  end do

  Return
End Subroutine RdFchkFCM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! lower triangular matrix --> symmetric square matrix
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine LT2Sqr(N,T,S)
  implicit real(kind=8) (a-h,o-z)
  real(kind=8) :: T(*),S(N,N)

  k=0
  do i=1,N
    do j=1,i-1
      k=k+1
      S(j,i)=T(k)
      S(i,j)=T(k)
    end do
    k=k+1
    S(i,i)=T(k)
  end do

  return
end subroutine LT2Sqr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! symmetric square matrix --> lower triangular matrix
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Sqr2LT(N,S,T)
 implicit real(kind=8) (a-h,o-z)
 real(kind=8) :: T(*),S(N,N)

 ii=0
 Do i=1,N
   Do j=1,i
     ii=ii+1
     T(ii)=(S(j,i)+S(i,j))*0.5d0
   end Do
 end Do

 return
end subroutine Sqr2LT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! It performs a 3x3 rotation operation of Cartesian vectors:
! mode >=0: XYZ' = XYZ * ROT
!      < 0: XYZ' = XYZ * ROT^-1 = XYZ * ROT^T
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rotvec(NV,mode,ROT,XYZ,SCR)
  implicit double precision (a-h,o-z)
  dimension XYZ(3,NV),ROT(3,3),SCR(3)

  if(mode >= 0) then
    do i=1,NV
      SCR = XYZ(:,i)
      XYZ(:,i) = 0.0d0
      do j=1,3
        do k=1,3
          XYZ(j,i)=XYZ(j,i)+ROT(k,j)*SCR(k)
        end do
      end do
    end do
  else
    do i=1,NV
      SCR = XYZ(:,i)
      XYZ(:,i) = 0.0d0
      do j=1,3
        do k=1,3
          XYZ(j,i)=XYZ(j,i)+ROT(j,k)*SCR(k)
        end do
      end do
    end do
  end if

  return
end subroutine rotvec

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Rotate the Hessian matrix
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rothess(natom,mode,ROT,FCM,SCR)
  implicit double precision (a-h,o-z)
  dimension FCM(3,natom,3,natom),ROT(9),SCR(3,3,2)

  do i=1,natom
    do j=1,natom
      SCR(:,:,1) = FCM(:,j,:,i)
      call rotmat(1,mode,ROT,SCR(1,1,1),SCR(1,1,2))
      FCM(:,j,:,i) = SCR(:,:,1)
    end do
  end do

  return
end subroutine rothess

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! It performs a 3x3 rotation operation of Cartesian matrices:
! mode >=0: M' = ROT^T * M * ROT
!      < 0: M' = ROT * M * ROT^T
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rotmat(NM,mode,ROT,AM,SCR)
  implicit double precision (a-h,o-z)
  dimension AM(9,NM),ROT(9),SCR(9)

  if(mode >= 0) then
    do i=1,NM
      call dgemm('n','n',3,3,3,1.0d0,AM(1,i),3,ROT,3,0.0d0,SCR,3)
      call dgemm('t','n',3,3,3,1.0d0,ROT,3,SCR,3,0.0d0,AM(1,i),3)
    end do
  else
    do i=1,NM
      call dgemm('n','t',3,3,3,1.0d0,AM(1,i),3,ROT,3,0.0d0,SCR,3)
      call dgemm('n','n',3,3,3,1.0d0,ROT,3,SCR,3,0.0d0,AM(1,i),3)
    end do
  end if

  return
end subroutine rotmat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! check (translated) Cartesian coordinates
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine CheckXYZ(natom,lrot,xyz1,xyz2,Scr)
  Implicit Real*8(A-H,O-Z)
  dimension :: xyz1(3,natom), xyz2(3,natom), Scr(3,3)
  logical :: lrot

  lrot = .false.

  ! O1: center of geometry-1
  Scr(:,1) = 0.0d0
  do i = 1, natom
    Scr(:,1) = Scr(:,1) + xyz1(:,i)
  end do
  Scr(:,1) = Scr(:,1) / dble(natom)

  ! O2: center of geometry-2
  Scr(:,2) = 0.0d0
  do i = 1, natom
    Scr(:,2) = Scr(:,2) + xyz2(:,i)
  end do
  Scr(:,2) = Scr(:,2) / dble(natom)

  ! check distances:
  ! A_i' = A_i - O1, B_i' = B_i - O2, C_i = B_i' - A_i'
  do i = 1, natom
    Scr(:,3) = (xyz2(:,i) - Scr(:,2)) - (xyz1(:,i) - Scr(:,1))
    x = sqrt( DDOT(3,Scr(1,3),1,Scr(1,3),1) )
    if(x > 1.0d-4) then
      ! call XError("The geometry has been rotated.")
      lrot = .true.
      exit
    end if
  end do

  Return
End Subroutine CheckXYZ

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read Cartesian coordinates and gradients from CFour's GRD file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine RdFchkGRD(ifchk,natom,XYZ,GRD,ctmp,tag)
  Implicit Real*8(A-H,O-Z)
  dimension :: XYZ(3*natom), GRD(3*natom)
  character*80 :: ctmp,tag

  rewind(ifchk)
  tag="Current cartesian coordinates              R   N= "
  do while(.true.)
    read(ifchk,"(a80)",iostat=ist) ctmp
      if(ist /= 0) call XError("No Cartesian coordinates found in the fchk file.")
    if(index(ctmp,tag(1:50)) > 0) then
      read(ifchk,"(5e16.8)") XYZ
      exit
    end if
  end do

  tag="Cartesian Gradient                         R   N= "
  do while(.true.)
    read(ifchk,"(a80)",iostat=ist) ctmp
      if(ist /= 0) call XError("No Cartesian gradients found in the fchk file.")
    if(index(ctmp,tag(1:50)) > 0) then
      read(ifchk,"(5e16.8)") GRD
      exit
    end if
  end do

  Return
End Subroutine RdFchkGRD

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read total energy from the fchk file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine RdFchkEne(ifchk,ENE,ctmp,tag)
  Implicit Real*8(A-H,O-Z)
  character*80 :: ctmp,tag

  ENE = 0.0d0

  ! read total energy
  rewind(ifchk)
  tag="Total Energy                               R"
  do while(.true.)
    read(ifchk,"(a80)",iostat=ist) ctmp
      if(ist /= 0) call XError("No total energy found in the fchk file.")
    if(index(ctmp,tag(1:44)) > 0) then
      read(ctmp(45:),*) ENE
      exit
    end if
  end do

  Return
End Subroutine RdFchkEne

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Write *.gjf file for istate
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine WrtGjf(iCTP,iINP,istate,natom,nder,IZA,XYZ,ctmp,tag)
  implicit real(kind=8) (a-h,o-z)
  parameter(au2ang=0.52917720859d0)
  dimension :: IZA(natom), XYZ(3,natom)
  character*200 :: ctmp,tag

  tag = "$SPi$GRADi$FREQi"
  write(tag( 4: 4),"(i1)") istate
  write(tag(10:10),"(i1)") istate
  write(tag(16:16),"(i1)") istate

  rewind(iINP)
  rewind(iCTP)
  do while(.true.)
    read(iCTP,"(a200)",iostat=ist) ctmp
    if(ist /= 0) then
      if(nder == 0) then
        ctmp = "No $SPi found in the templet."
        write(ctmp(7:7),"(i1)") istate
        call XError(ctmp(1:29))
      else if(nder == 1) then
        ctmp = "No $GRADi found in the templet."
        write(ctmp(9:9),"(i1)") istate
        call XError(ctmp(1:31))
      else
        ctmp = "No $FREQi found in the templet."
        write(ctmp(9:9),"(i1)") istate
        call XError(ctmp(1:31))
      end if
    end if

    if(ctmp(1:1) == '!') cycle

    istr=nonspace(ctmp)
    if(ctmp(istr:istr) == '$') call charl2u(ctmp)
    if(nder == 0 .and. index(ctmp,tag( 1: 4)) > 0) exit
    if(nder == 1 .and. index(ctmp,tag( 5:10)) > 0) exit
    if(nder == 2 .and. index(ctmp,tag(11:16)) > 0) exit
  end do

  loop1: do while(.true.)
    read(iCTP,"(a200)",iostat=ist) ctmp
      if(ist /= 0) call XError("No *END_OF_INPUT found in the templet.")

    if(ctmp(1:1) == '!') cycle

    istr=nonspace(ctmp)
    if(ctmp(istr:istr) == '*') call charl2u(ctmp)
    if(index(ctmp,"*BEFORE_GEOM") > 0) then
      loop2: do while(.true.)
        read(iCTP,"(a200)",iostat=ist) ctmp
          if(ist /= 0) call XError("Please check *BEFORE_GEOM in the templet.")
        istr=nonspace(ctmp)
        if(ctmp(istr:istr) == '*') call charl2u(ctmp)
        if(index(ctmp,"*END_OF_INPUT") > 0) exit loop2
        write(iINP,"(a)") trim(ctmp)
      end do loop2

      ! write Cartesian geometry
      do i=1,natom
        call ElemZA(ctmp,IZA(i))
        write(iINP,"(a3,3f20.12)") ctmp(1:3), XYZ(:,i)*au2ang
      end do
      write(iINP,*)
    end if

    istr=nonspace(ctmp)
    if(ctmp(istr:istr) == '*') call charl2u(ctmp)
    if(index(ctmp,"*AFTER_GEOM") > 0) then
      loop3: do while(.true.)
        read(iCTP,"(a200)",iostat=ist) ctmp
          if(ist /= 0) call XError("Please check *AFTER_GEOM in the templet.")
        istr=nonspace(ctmp)
        if(ctmp(istr:istr) == '*') call charl2u(ctmp)
        if(index(ctmp,"*END_OF_INPUT") > 0) exit loop1
        write(iINP,"(a)") trim(ctmp)
      end do loop3
    end if

  end do loop1

  write(iINP,"(//)")

  return
end subroutine WrtGjf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read Cartesian coordinates (in a.u.) from Gaussian's *.EIn. Fixed for G16.c ONIOM calculation.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine RdEIn_XYZ(iGIN,natom,natall,IZA,XYZ,Llist)
  Implicit Real*8(A-H,O-Z)
  dimension :: IZA(natom), XYZ(3,natom), Llist(natall)
  allocatable :: Scr(:)

  IZA = 0
  XYZ = 0.0d0
  Llist = 0

  allocate(Scr(3))

  rewind(iGIN)
  read(iGIN,*)
  i1 = 0
  do i = 1, natall
    read(iGIN,*,iostat=ist) iza1, Scr
    if(ist /= 0) call XError("Please check the Cartesian coordinates in *.EIn.")
    if(iza1 > 0) then
      i1 = i1 + 1
      IZA(i1) = iza1
      XYZ(:,i1) = Scr
      Llist(i) = 1
    end if
  end do

  deallocate(Scr)

  if(i1 /= natom) call XError("Wrong natom in *.EIn.")

  Return
End Subroutine RdEIn_XYZ

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read arguments:
!
! Job type:
! -GEN        Imode = 0 (default). Generate input files for states. Required options and files:
!      -GIN   Gaussian's *.EIn file name
!      -CTP   templet file name
!      -NST   Number of states (NStat) to be mixed, which can be 1, 2 (default), 3, ..., MaxStat
!      -INi   Gaussian input file for states i, i=1,2,..
!
! -MIX        Imode = 1. Compute energy, gradients, and/or hessians of the mixed-spin ground state.
!      -GIN   Gaussian's *.EIn file name
!      -GOU   Gaussian's *.EOu file name
!      -NST   Number of states (NStat) to be mixed, which can be 1, 2 (default), 3, ..., MaxStat
!      -FCi   Gaussian fchk file for states i, i=1,2,..
!      -CHI   Empirical SO constant (in cm^-1). Default: 400 cm^-1
!      -CHS   Specify each SO constant separately (in cm^-1)
!      -SHi   Energy shifts for states i, i=1,2,.. (in cm^-1)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdArg(imode,iGIN,iGOU,iCTP,iINP,iFCH,MaxStat,NStat,ifchs,SO_chi,SO_chs,SO_dlt,ctmp)
  implicit real(kind=8) (a-h,o-z)
  real(kind=8),parameter :: con_cm2au = 219474.63137d0
  dimension :: iINP(MaxStat), iFCH(MaxStat), SO_chs(*), SO_dlt(MaxStat)
  character :: ctmp*200
  logical :: ifchs, ifopen

  imode = 0
  NStat = 2
  SO_chi= 4.0d2 / con_cm2au
  SO_dlt= 0.0d0
  ifchs = .false.

  i = 0
  do while(.true.)
    i = i + 1
    call Get1Cmd(i,.True.,ctmp)
    istr=nonspace(ctmp)
    iend=len_trim(ctmp)

    if(iend == 0) then
      exit
    else if(ctmp(istr:iend) == '-MIX') then
      imode = 1
    else if(ctmp(istr:iend) == '-GEN') then
      imode = 0
    else if(ctmp(istr:iend) == '-NST') then
      i = i + 1
      call Get1Cmd(i,.False.,ctmp)
      istr=nonspace(ctmp)
      iend=len_trim(ctmp)
      if(iend > 0) then
        read(ctmp(istr:iend),*,iostat=ist) NStat
        if(ist /= 0) call XError("Cannot read NStat!")
        if(NStat < 1 .or. NStat > MaxStat) call XError("NStat is out of range.")
      else
        call XError("NStat is not provided!")
      end if
    else if(ctmp(istr:iend) == '-GIN') then
      i = i + 1
      call Get1Cmd(i,.False.,ctmp)
      istr=nonspace(ctmp)
      iend=len_trim(ctmp)
      if(iend == 0) then
        call XError("No file name provided for -GIN.")
      else
        open(iGIN,file=ctmp(istr:iend),status='old',iostat=ist)
        if(ist /= 0) call XError("Cannot open old file for -GIN.")
      end if
    else if(ctmp(istr:iend) == '-GOU') then
      i = i + 1
      call Get1Cmd(i,.False.,ctmp)
      istr=nonspace(ctmp)
      iend=len_trim(ctmp)
      if(iend == 0) then
        call XError("No file name provided for -GOU.")
      else
        open(iGOU,file=ctmp(istr:iend),status='replace',iostat=ist)
        if(ist /= 0) call XError("Cannot open file for -GOU.")
      end if
    else if(ctmp(istr:iend) == '-CTP') then
      i = i + 1
      call Get1Cmd(i,.False.,ctmp)
      istr=nonspace(ctmp)
      iend=len_trim(ctmp)
      if(iend == 0) then
        call XError("No file name provided for -CTP.")
      else
        open(iCTP,file=ctmp(istr:iend),status='old',iostat=ist)
        if(ist /= 0) call XError("Cannot open file for -CTP.")
      end if
    else if(ctmp(istr:iend-1) == '-IN') then
      read(ctmp(iend:iend),*,iostat=ist) iport
      if(ist /= 0) call XError("Cannot read port number from -INi!")
      if(iport < 1 .or. iport > MaxStat) call XError("The port number in -INi is out of range.")
      i = i + 1
      call Get1Cmd(i,.False.,ctmp)
      istr=nonspace(ctmp)
      iend=len_trim(ctmp)
      if(iend == 0) then
        call XError("No file name provided for -INi.")
      else
        open(iINP(iport),file=ctmp(istr:iend),status='replace',iostat=ist)
        if(ist /= 0) call XError("Cannot open input file for -INi.")
      end if
    else if(ctmp(istr:iend-1) == '-FC') then
      read(ctmp(iend:iend),*,iostat=ist) iport
      if(ist /= 0) call XError("Cannot read port number from -FCi!")
      if(iport < 1 .or. iport > MaxStat) call XError("The port number in -FCi is out of range.")
      i = i + 1
      call Get1Cmd(i,.False.,ctmp)
      istr=nonspace(ctmp)
      iend=len_trim(ctmp)
      if(iend == 0) then
        call XError("No file name provided for -FCi.")
      else
        open(iFCH(iport),file=ctmp(istr:iend),status='old',iostat=ist)
        if(ist /= 0) call XError("Cannot open fchk file for -FCi.")
      end if
    else if(ctmp(istr:iend) == '-CHI') then
      i = i + 1
      call Get1Cmd(i,.False.,ctmp)
      istr=nonspace(ctmp)
      iend=len_trim(ctmp)
      if(iend > 0) then
        read(ctmp(istr:iend),*,iostat=ist) SO_chi
        if(ist /= 0) call XError("Cannot read SO_chi!")
      else
        call XError("The SO_chi value is not provided!")
      end if
      ! cm-1 --> au
      SO_chi = max(1.0d1,abs(SO_chi)) / con_cm2au
    else if(ctmp(istr:iend) == '-CHS') then
      ifchs = .true.
      do j = 1, NStat*(NStat-1)/2
        i = i + 1
        call Get1Cmd(i,.False.,ctmp)
        istr=nonspace(ctmp)
        iend=len_trim(ctmp)
        if(iend > 0) then
          read(ctmp(istr:iend),*,iostat=ist) SO_chs(j)
          if(ist /= 0) call XError("Cannot read an element of the SO_chs array!")
        else
          call XError("The array elements of SO_chs are insufficient!")
        end if
        ! cm-1 --> au
        SO_chs(j) = SO_chs(j) / con_cm2au
      end do
    else if(ctmp(istr:iend-1) == '-SH') then
      read(ctmp(iend:iend),*,iostat=ist) j
      if(ist /= 0) call XError("Cannot read the index number from -SHi!")
      if(j < 1 .or. j > MaxStat) call XError("The index number in -SHi is out of range.")
      i = i + 1
      call Get1Cmd(i,.False.,ctmp)
      istr=nonspace(ctmp)
      iend=len_trim(ctmp)
      if(iend > 0) then
        read(ctmp(istr:iend),*,iostat=ist) SO_dlt(j)
        if(ist /= 0) call XError("Cannot read delta!")
      else
        call XError("The delta value is not provided for -SHi!")
      end if
      ! cm-1 --> au
      SO_dlt(j) = SO_dlt(j) / con_cm2au
    else
      call XError("Unknown argument "//ctmp(istr:iend))
    end if
  end do

  if(Imode == 1) then
    inquire(unit=iGIN,opened=ifopen)
      if(.NOT. ifopen) call XError("-GIN is not defined for -MIX!")
    inquire(unit=iGOU,opened=ifopen)
      if(.NOT. ifopen) call XError("-GOU is not defined for -MIX!")
    do i = 1, NStat
      inquire(unit=iFCH(i),opened=ifopen)
      if(.NOT. ifopen) then
        write(*,"(/,' The fchk file name of state ',i1,' is required for -MIX.',/)") i
        call XError("A fchk file has not been provided by -FCi!")
      end if
    end do
    if(NStat == 1) then
      SO_chi= 0.0d0
      SO_dlt= 0.0d0
      ifchs = .false.
    end if
  else
    inquire(unit=iGIN,opened=ifopen)
      if(.NOT.ifopen) call XError("-GIN is not defined for -GEN!")
    inquire(unit=iCTP,opened=ifopen)
      if(.NOT.ifopen) call XError("-CTP is not defined for -GEN!")
    do i = 1, NStat
      inquire(unit=iINP(i),opened=ifopen)
      if(.NOT. ifopen) then
        write(*,"(/,' The input file name of state ',i1,' is required for -GEN.',/)") i
        call XError("A input file has not been provided by -INi!")
      end if
    end do
  end if

  return
end subroutine RdArg

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read an argument, in upper case if ifUP
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine Get1Cmd(i,ifUP,ctmp)
  Implicit Real*8(A-H,O-Z)
  logical ifUP
  character*200 ctmp

  call get_command_argument(i,ctmp)
  ! call getarg(i,ctmp)
  if(ifUP) call charl2u(ctmp)

  Return
End Subroutine Get1Cmd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! position of the first non-space character in a string.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nonspace(string)
  implicit double precision (a-h,o-z)
  character*(*) string

  length=LEN_TRIM(string)
  if(length <= 1) then
    i=length
  else
    do i=1,length
      if(string(i:i) /= ' ') exit
    end do
  endif

  nonspace=i

  return
end function nonspace

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! cha --> CHA
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine charl2u(cha)
  implicit real(kind=8) (a-h,o-z)
  character*(*) :: cha
  character*1  :: L2U

  do i=1,len_trim(cha)
    cha(i:i)=L2U(cha(i:i))
  end do

  return
end subroutine charl2u

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! CHA --> cha
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine charu2l(cha)
  implicit real(kind=8) (a-h,o-z)
  character*(*) :: cha
  character*1  :: U2L

  do i=1,len_trim(cha)
    cha(i:i)=U2L(cha(i:i))
  end do

  return
end subroutine charu2l

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! L --> l
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U2L(letter)
  implicit real(kind=8) (a-h,o-z)
  character*1 :: letter,U2L

  if((ichar(letter) >= 65).and.(ichar(letter) <= 90))then
    U2L=char(ichar(letter)+32)
  else
    U2L=letter
  endif

  return
end function U2L

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! l --> L
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L2U(letter)
  implicit real(kind=8) (a-h,o-z)
  character*1 :: letter,L2U

  if( ichar(letter) >= 97 .and. ichar(letter) <= 122 )then
    L2U=char(ichar(letter)-32)
  else
    L2U=letter
  endif

  return
end function L2U

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read an error message and stop
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine XError(inf)
  implicit real(kind=8) (a-h,o-z)
  character*(*) :: inf

  write(*,"(/,' *** Error! ',a)")trim(inf)

  stop

  return
end subroutine XError

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read the first line of Gaussian's *.EIn.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine RdEIn_Line1(iGIN,natom,natall,nder)
  Implicit Real*8(A-H,O-Z)

  rewind(iGIN)

  natall=-1
  nder=-1

  read(iGIN,*,iostat=ist) natall,nder
  if(ist /= 0) call XError("Please check the first line in *.EIn.")
  if(natall < 1) call XError("Natall < 1.")
  if(nder < 0 .or. nder > 2) call XError("Nder is out of range.")

  ! G16.c's ONIOM calculation also prints dummy atoms, which has to be fixed.
  natom = 0
  do i=1, natall
    read(iGIN,*) i1
    if(i1 > 0) natom = natom + 1
  end do

  Return
End Subroutine RdEIn_Line1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Returns element symbol "el" for nuclear charge iza.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ElemZA(el,iza)
  implicit real(kind=8) (a-h,o-z)
  parameter (maxza=120)
  character*3 :: el,atomlib(maxza)
  data (atomlib(i),i=1,maxza) / &
   'H  ','HE ','LI ','BE ','B  ','C  ','N  ','O  ','F  ','NE ',   'NA ','MG ','AL ','SI ','P  ','S  ','CL ','AR ','K  ','CA ', &
   'SC ','TI ','V  ','CR ','MN ','FE ','CO ','NI ','CU ','ZN ',   'GA ','GE ','AS ','SE ','BR ','KR ','RB ','SR ','Y  ','ZR ', &
   'NB ','MO ','TC ','RU ','RH ','PD ','AG ','CD ','IN ','SN ',   'SB ','TE ','I  ','XE ','CS ','BA ','LA ','CE ','PR ','ND ', &
   'PM ','SM ','EU ','GD ','TB ','DY ','HO ','ER ','TM ','YB ',   'LU ','HF ','TA ','W  ','RE ','OS ','IR ','PT ','AU ','HG ', &
   'TL ','PB ','BI ','PO ','AT ','RN ','FR ','RA ','AC ','TH ',   'PA ','U  ','NP ','PU ','AM ','CM ','BK ','CF ','ES ','FM ', &
   'MD ','NO ','LR ','RF ','DB ','SG ','BH ','HS ','MT ','DS ',   'RG ','CN ','NH ','FL ','MC ','LV ','TS ','OG ','119','120'/
  save atomlib

  el = "???"
  if(iza > 0 .and. iza <= maxza) el = adjustl(atomlib(iza))
  if(iza > 0 .and. iza <= 118) call charu2l(el(2:2))

  return
end subroutine ElemZA

!--- END

