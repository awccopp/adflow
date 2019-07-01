subroutine test(speed1, speed2, block_size)

  use constants
  use blockPointers, only : bnx=>nx, bny=>ny, bnz=>nz, &
       bil=>il, bjl=>jl, bkl=>kl, &
       bie=>ie, bje=>je, bke=>ke, &
       bib=>ib, bjb=>jb, bkb=>kb, &
       bw=>w, bp=>p, bgamma=>gamma, &
       bx=>x, brlv=>rlv, brev=>rev, bvol=>vol, bVolRef=>volRef, bd2wall=>d2wall, &
       biblank=>iblank, bPorI=>porI, bPorJ=>porJ, bPorK=>porK, bdw=>dw, bfw=>fw
  use utils
  use block
  use surfaceFamilies, only : fullFamList
  use flowVarRefState
  use communication, only : myid, adflow_comm_world
  use blockette, only : blocketteRes, BS
  use inputDiscretization, only : useBlockettes
  implicit none

  real(kind=realType), intent(out) :: speed1, speed2
  integer(kind=intType), intent(out) :: block_size

  ! Misc
  integer(kind=intType) :: i, j, k, l, iSize, nn, nIter, ii,jj,kk, iiter, loopCount, ierr
  real(kind=realType) ::  timeA, timeB, tmp, tmp2, norm
  integer(kind=intType) :: omp_get_thread_num
  logical :: useBlockettesSave
  logical :: updateVars

  ! First call the residual routine to update all intermed. variables
  call blocketteRes(useUpdateDt=.True.)

  ! Disable blockette code for the first run
  useBlockettesSave = useBlockettes
  useBlockettes = .False.

  niter = 5
  block_size = BS
  call mpi_barrier(adflow_comm_world, ierr)
  timeA = mpi_wtime()

  do i=1, nIter
     ! Call blockette code without vectorized routines
     call blocketteRes(useUpdateDt=.True., useUpdateVars=updateVars)
  end do
  call mpi_barrier(adflow_comm_world, ierr)
  timeB = mpi_wtime()

  iSize = 0
  do nn=1, nDom
     call setPointers(nn, 1, 1)
     iSize = iSize + bnx*bny*bnz
     !print *,'dw:', dw(2,2,2,1), dw(2,2,2,6)
     !print *, 'sum:', myid, sum(bdw(2:bil, 2:bjl, 2:bkl, 1:nw))
  end do

  speed1 = iSize*nIter/(timeB-timeA)/1e6

  call calcResNorm(norm)

  if (myid == 0) then
    print *,'speed1:', speed1, timeB-timeA, 0, isize, norm ,bnx, bny, bnz
  end if

  ! Enable the vectorized code for the second run
  useBlockettes = .True.

  call mpi_barrier(adflow_comm_world, ierr)
  ! Test the blockette code
  timeA = mpi_wtime()
  do i=1, nIter
     ! Call blockette code with vectorized routines
     call blocketteRes(useUpdateDt=.True., useUpdateVars=updateVars)
  end do
  call mpi_barrier(adflow_comm_world, ierr)
  timeB = mpi_wtime()

  !print *,'dw:', dw(2,2,2,1), dw(2,2,2,6)
  !print *, 'sum:', myid, sum(bdw(2:bil, 2:bjl, 2:bkl, 1:nw))

  speed2 = iSize*nIter/(timeB-timeA)/1e6
  call calcResNorm(norm)

  if (myid == 0) then
    print *,'speed2:',speed2 , timeB-timeA, loopCount, isize, norm, BS
  end if
  !print *, 'loop count:', loopCount, il, jl, kl

  ! Restore the useBlockettes option
  useBlockettes = useBlockettesSave

end subroutine test

subroutine test_b(speed1, speed2, block_size)

  use constants
  use blockPointers, only : bnx=>nx, bny=>ny, bnz=>nz, &
       bil=>il, bjl=>jl, bkl=>kl, &
       bie=>ie, bje=>je, bke=>ke, &
       bib=>ib, bjb=>jb, bkb=>kb, &
       bw=>w, bp=>p, bgamma=>gamma, &
       bx=>x, brlv=>rlv, brev=>rev, bvol=>vol, bVolRef=>volRef, bd2wall=>d2wall, &
       biblank=>iblank, bPorI=>porI, bPorJ=>porJ, bPorK=>porK, bdw=>dw, bfw=>fw
  use utils
  use block
  use surfaceFamilies, only : fullFamList
  use flowVarRefState
  use communication, only : myid, adflow_comm_world
  use blockette, only : blocketteRes, BS
  use inputDiscretization, only : useBlockettes
  use adjointAPI, only : computeMatrixFreeProductBwdFast
  use ADjointVars , only: nCellsLocal
#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
#else
  implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#endif

  Vec vecX, vecY


  real(kind=realType), intent(out) :: speed1, speed2
  integer(kind=intType), intent(out) :: block_size

  ! Misc
  integer(kind=intType) :: i, j, k, l, iSize, nn, nIter, ii,jj,kk, iiter, loopCount, sps
  real(kind=realType) ::  timeA, timeB, tmp, tmp2, norm
  integer(kind=intType) :: omp_get_thread_num
  logical :: useBlockettesSave
  logical :: updateVars

  integer(kind=intType) ::ierr, nDimW

  real(kind=realType), pointer :: dwb_pointer(:)
  real(kind=realType), pointer :: wb_pointer(:)

  ! First call the residual routine to update all intermed. variables
  call blocketteRes(useUpdateDt=.True.)

  ! We need to create a dummy AD seed, and an output vector to read the output
  nDimW = nw * nCellsLocal(1_intTYpe) * 1_intType

  ! print *,nDimW

  ! create dummy arrays
  call VecCreate(ADFLOW_COMM_WORLD, vecX, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecSetSizes(vecX, nDimW, PETSC_DECIDE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecSetBlockSize(vecX, nw, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecSetType(vecX, VECMPI, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  !  Create duplicates
  call VecDuplicate(vecX, vecY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! copy the residual to dwbar
  call VecGetArrayF90(vecX, dwb_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! print *,'setting the dummy seeds'

  ii = 0
  do nn=1,nDom
    call setPointers_d(nn, 1, 1)
    do k=2, bkl
       do j=2, bjl
          do i=2, bil
             do l=1, 6
                ii = ii + 1
                ! print *, i,j,k,l,ii, bdw(i,j,k,l)
                dwb_pointer(ii) = bdw(i, j, k, l)
             end do
          end do
       end do
    end do
  end do

  ! print *,'flag1'

  call VecRestoreArrayF90(Vecx, dwb_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! now we can run tests like the usual code
  call VecGetArrayReadF90(vecX, dwb_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecGetArrayF90(VecY, wb_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Disable blockette code for the first run
  useBlockettesSave = useBlockettes
  useBlockettes = .False.

  ! print *,'flag2'

  niter = 5
  block_size = BS
  call mpi_barrier(adflow_comm_world, ierr)
  timeA = mpi_wtime()

  do i=1, nIter
     ! Call blockette code without vectorized routines
     call computeMatrixFreeProductBwdFast(dwb_pointer, wb_pointer, size(dwb_pointer))
  end do
  call mpi_barrier(adflow_comm_world, ierr)
  timeB = mpi_wtime()

  iSize = 0
  do nn=1, nDom
     call setPointers(nn, 1, 1)
     iSize = iSize + bnx*bny*bnz
     !print *,'dw:', dw(2,2,2,1), dw(2,2,2,6)
     !print *, 'sum:', myid, sum(bdw(2:bil, 2:bjl, 2:bkl, 1:nw))
  end do

  speed1 = iSize*nIter/(timeB-timeA)/1e6
  ! print *,'flag3'

  ! call calcDWDNorm(norm)
  call calcADNorm(norm)

  if (myid == 0) then
    print *,'speed1:', speed1, timeB-timeA, 0, isize, norm ,bnx, bny, bnz
  end if

  ! Enable the vectorized code for the second run
  useBlockettes = .True.

  ! ! we need to re-set the AD seeds
  ! ii = 0
  ! do nn=1,nDom
  !   call setPointers_d(nn, 1, 1)
  !   do k=2, bkl
  !      do j=2, bjl
  !         do i=2, bil
  !            do l=1, 6
  !               ii = ii + 1
  !               ! print *, i,j,k,l,ii, bdw(i,j,k,l)
  !               dwb_pointer(ii) = bdw(i, j, k, l)
  !            end do
  !         end do
  !      end do
  !   end do
  ! end do

  call mpi_barrier(adflow_comm_world, ierr)
  ! Test the blockette code
  timeA = mpi_wtime()
  do i=1, nIter
     ! Call blockette code with vectorized routines
     call computeMatrixFreeProductBwdFast(dwb_pointer, wb_pointer, size(dwb_pointer))
  end do
  call mpi_barrier(adflow_comm_world, ierr)
  timeB = mpi_wtime()

  call VecRestoreArrayF90(vecY, wb_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecRestoreArrayF90(Vecx, dwb_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  !print *,'dw:', dw(2,2,2,1), dw(2,2,2,6)
  !print *, 'sum:', myid, sum(bdw(2:bil, 2:bjl, 2:bkl, 1:nw))

  speed2 = iSize*nIter/(timeB-timeA)/1e6
  ! call calcDWDNorm(norm)
  call calcADNorm(norm)

  ! We need to also compare the norms of pd, rlvd, revd, radid, radjd, radkd

  if (myid == 0) then
    print *,'speed2:',speed2 , timeB-timeA, 0, isize, norm, BS
  end if
  !print *, 'loop count:', loopCount, il, jl, kl

  ! Restore the useBlockettes option
  useBlockettes = useBlockettesSave

end subroutine test_b

subroutine calcResNorm(norm)

  use constants
    use blockPointers, only : bnx=>nx, bny=>ny, bnz=>nz, &
         bil=>il, bjl=>jl, bkl=>kl, &
         bie=>ie, bje=>je, bke=>ke, &
         bib=>ib, bjb=>jb, bkb=>kb, &
         bw=>w, bp=>p, bgamma=>gamma, &
         bx=>x, brlv=>rlv, brev=>rev, bvol=>vol, bVolRef=>volRef, bd2wall=>d2wall, &
         biblank=>iblank, bPorI=>porI, bPorJ=>porJ, bPorK=>porK, bdw=>dw, bfw=>fw
  implicit none

  integer(kind=intType) :: i, j, k, l
  real(kind=realType), intent(out) :: norm

  ! Calculate the norm
  norm = 0.0_realType
  ! do k = 1, bnz
  !   do j = 1, bny
  !     do i = 1, bnx
  do k=2, bkl
     do j=2, bjl
        do i=2, bil
           do l = 1, 6
             norm = norm + bdw(i,j,k,l)*bdw(i,j,k,l)
        end do
      end do
    end do
  end do

  norm = sqrt(norm)/(bnx*bny*bnz)

end subroutine calcResNorm

subroutine calcADNorm(norm)

    use constants
    use block, only : nDom

    use blockPointers, only : bnx=>nx, bny=>ny, bnz=>nz, &
         bil=>il, bjl=>jl, bkl=>kl, &
         bie=>ie, bje=>je, bke=>ke, &
         bib=>ib, bjb=>jb, bkb=>kb, &
         bw=>w, bp=>p, bgamma=>gamma, &
         bx=>x, brlv=>rlv, brev=>rev, bvol=>vol, bVolRef=>volRef, bd2wall=>d2wall, &
         biblank=>iblank, bPorI=>porI, bPorJ=>porJ, bPorK=>porK, bdw=>dw, bfw=>fw, bwd=>wd
    use utils, only : setPointers_d, setPointers

    implicit none

    integer(kind=intType) :: i, j, k, l, sps, nn
    real(kind=realType), intent(out) :: norm

    norm = 0.0_realType
    do nn=1,nDom
       do sps=1,1
         call setPointers_d(nn, 1, sps)
         do k=2, bkl
            do j=2, bjl
               do i=2, bil
                   do l=1, 6
                      norm = norm + bwd(i, j, k, l)*bwd(i, j, k, l)
                   end do
                end do
             end do
          end do
       end do
    end do

    norm = sqrt(norm)/(bnx*bny*bnz)
end subroutine calcADNorm

subroutine calcDWDNorm(norm)

    use constants
    use block, only : nDom

    use blockPointers, only : bnx=>nx, bny=>ny, bnz=>nz, &
         bil=>il, bjl=>jl, bkl=>kl, &
         bie=>ie, bje=>je, bke=>ke, &
         bib=>ib, bjb=>jb, bkb=>kb, &
         bw=>w, bp=>p, bgamma=>gamma, &
         bx=>x, brlv=>rlv, brev=>rev, bvol=>vol, bVolRef=>volRef, bd2wall=>d2wall, &
         biblank=>iblank, bPorI=>porI, bPorJ=>porJ, bPorK=>porK, bdw=>dw, bfw=>fw, bdwd=>dwd
    use utils, only : setPointers_d, setPointers

    implicit none

    integer(kind=intType) :: i, j, k, l, sps, nn
    real(kind=realType), intent(out) :: norm

    norm = 0.0_realType
    do nn=1,nDom
       do sps=1,1
         call setPointers_d(nn, 1, sps)
         do k=2, bkl
            do j=2, bjl
               do i=2, bil
                   do l=1, 6
                      norm = norm + bdwd(i, j, k, l)*bdwd(i, j, k, l)
                   end do
                end do
             end do
          end do
       end do
    end do

    norm = sqrt(norm)/(bnx*bny*bnz)
end subroutine calcDWDNorm

subroutine orig

    use constants
    use communication, only : adflow_comm_world
    use BCRoutines, only : applyallBC_block
    use turbbcRoutines, only : applyallTurbBCthisblock, bcTurbTreatment
    use iteration, only : currentLevel, rfil
    use inputAdjoint,  only : viscPC
    use flowVarRefState, only : nwf, nw
    use blockPointers, only : nDom, il, jl, kl
    use flowVarRefState, only : viscous
    use inputPhysics , only : turbProd, equationMode, equations, turbModel
    use inputDiscretization, only : lowSpeedPreconditioner, lumpedDiss, spaceDiscr, useAPproxWallDistance
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use initializeFlow, only : referenceState
    use section, only: sections, nSections
    use monitor, only : timeUnsteadyRestart
    use sa, only : saSource, saViscous, saResScale, qq
    use haloExchange, only : exchangeCoor, whalo2
    use wallDistance, only : updateWallDistancesQuickly
    use solverUtils, only : timeStep_block
    use flowUtils, only : allNodalGradients, computeLamViscosity, computePressureSimple, &
         computeSpeedOfSoundSquared, adjustInflowAngle
    use fluxes, only : inviscidDissFluxScalarApprox, inviscidDissFluxMatrixApprox, &
         inviscidUpwindFlux, inviscidDissFluxScalar, inviscidDissFluxMatrix, &
         viscousFlux, viscousFluxApprox, inviscidCentralFlux
    use utils, only : setPointers, EChk
    use turbUtils, only : turbAdvection, computeEddyViscosity
    use residuals, only : initRes_block
    use surfaceIntegrations, only : integrateSurfaces
    use adjointExtra, only : volume_block, metric_block, boundaryNormals,&
         xhalo_block, sumdwandfw, resScale
    use oversetData, only : oversetPresent
    use inputOverset, only : oversetUpdateMode
    use oversetCommUtilities, only : updateOversetConnectivity
    implicit none

    ! Working Variables
    integer(kind=intType) :: ierr, nn, sps, fSize
    real(kind=realType), dimension(nSections) :: t
    real(kind=realType), dimension(nLocalValues, nTimeIntervalsSpectral) :: localVal, globalVal

    ! Zero out the accumulation array for forces and other integrated quantities
    localVal = zero

    call initRes_block(1, nw, nn, sps)
    allocate(qq(2:il, 2:jl, 2:kl))
    call saSource
    call turbAdvection(1_intType, 1_intType, itu1-1, qq)
    !call unsteadyTurbTerm(1_intType, 1_intType, itu1-1, qq)
    call saViscous
    call saResScale
    deallocate(qq)

    rFil = one
    ! Compute the mean flow residuals
    call timeStep_block(.false.)
    call inviscidCentralFlux
    call inviscidDissFluxScalar
    call computeSpeedOfSoundSquared
    call allNodalGradients
    call viscousFlux
    call sumDwAndFw
    call resScale

  end subroutine orig

! subroutine superCode(w, x, vol, P, rlv, rev, gamma, porI, porJ, porK, d2wall, &
!      volRef, dw, fw, iblank)

!   use constants
!   use blockette, only : BS
!   use sa, only : cv13, kar2Inv, cw36, cb3Inv
!   use flowvarRefState, only : timeRef, nw, nwf, nt1, nt2, eddyModel, gammaInf, pInfCorr, &
!        viscous, rhoInf
!   use inputDiscretization, only : adis, vis4, vis2
!   use section, only : sections
!   use blockPointers, only : sectionID
!   use paramTurb
!   use inputPhysics, only :useft2SA, useRotationSA, turbProd, equations, useQCR,&
!        prandtl, prandtlturb
!   use turbMod, only : secondOrd
!   use utils, only : getCorrectForK
!   use iteration, only : rFil
!   use inputIteration, only : turbResScale
!   implicit none

!   ! Dummy Block dimensions
!   integer(kind=intType), parameter :: nx=BS, ny=BS, nz=BS
!   integer(kind=intType), parameter :: il=BS+1, jl=BS+1, kl=BS+1
!   integer(kind=intType), parameter :: ie=BS+2, je=BS+2, ke=BS+2

!   !Inputs
!   real(kind=realType), dimension(0:ie, 0:je, 0:ke, 1:nw), intent(in) :: w
!   real(kind=realType), dimension(0:ie, 0:je, 0:ke, 3), intent(in) :: x
!   real(kind=realType), dimension(1:ie, 1:je, 1:ke), intent(in) :: &
!        rlv, rev, vol
!   real(kind=realType), dimension(2:il, 2:jl, 2:kl) :: volRef, d2wall
!   real(kind=realType), dimension(0:ie, 0:je, 0:ke), intent(in) :: P, gamma
!   integer(kind=porType), dimension(1:il, 2:jl, 2:kl), intent(in) :: porI
!   integer(kind=porType), dimension(2:il, 1:jl, 2:kl), intent(in) :: porJ
!   integer(kind=porType), dimension(2:il, 2:jl, 1:kl), intent(in) :: porK
!   integer(kind=intType), dimension(2:il, 2:jl, 2:kl), intent(in) :: iblank

!   ! Input/ouput
!   real(kind=realType), dimension(1:ie, 1:je, 1:ke, 1:nw), intent(inout) :: fw

!   ! Ouputs
!   real(kind=realType), dimension(1:ie, 1:je, 1:ke, 1:nw), intent(out) :: dw

!   ! Intermediate working arrays
!   real(kind=realType), dimension(0:il, 1:jl, 1:jl, 3) :: sI
!   real(kind=realType), dimension(1:il, 0:jl, 1:kl, 3) :: sJ
!   real(kind=realType), dimension(1:il, 1:jl, 0:kl, 3) :: sK
!   real(kind=realType), dimension(1:ie, 1:je, 1:ke) :: radI, radJ, radK
!   real(kind=realType),dimension(1:il,1:jl,1:kl, 3) :: dss
!   real(kind=realType), dimension(0:ie, 0:je, 0:ke) :: ss
!   real(kind=realType), dimension(1:il, 1:jl, 1:kl) :: ux, uy, uz
!   real(kind=realType), dimension(1:il, 1:jl, 1:kl) :: vx, vy, vz
!   real(kind=realType), dimension(1:il, 1:jl, 1:kl) :: wx, wy, wz
!   real(kind=realType), dimension(1:il, 1:jl, 1:kl) :: qx, qy, qz
!   real(kind=realType), dimension(1:il, 1:jl, 1:kl) :: aa
!   ! Integers
!   integer(kind=intType) :: i, j, k, n, m, l, ii, jj

!   ! Reals
!   real(kind=realType) :: v1(3), v2(3)

!   ! Variables for sa Souce
!   real(kind=realType) :: fv1, fv2, ft2
!   real(kind=realType) :: sst, nu, dist2Inv, chi, chi2, chi3
!   real(kind=realType) :: rr, gg, gg6, termFw, fwSa, term1, term2
!   real(kind=realType) :: dfv1, dfv2, dft2, drr, dgg, dfw, sqrtProd
!   real(kind=realType) :: uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz
!   real(kind=realType) :: div2, fact, sxx, syy, szz, sxy, sxz, syz
!   real(kind=realType) :: vortx, vorty, vortz
!   real(kind=realType) :: omegax, omegay, omegaz
!   real(kind=realType) :: strainMag2, strainProd, vortProd
!   real(kind=realType), parameter :: xminn = 1.e-10_realType
!   real(kind=realType), parameter :: f23 = two*third

!   ! Variables for sa Viscous
!   real(kind=realType) :: voli, volmi, volpi, xm, ym, zm, xp, yp, zp
!   real(kind=realType) :: xa, ya, za, ttm, ttp, cnud, cam, cap
!   real(kind=realType) :: nutm, nutp, num, nup, cdm, cdp
!   real(kind=realType) :: c1m, c1p, c10, b1, c1, d1, qs

!   ! Variables for sa Advection
!   real(kind=realType) :: uu, dwt, dwtm1, dwtp1, dwti, dwtj, dwtk
!   integer(kind=intType), parameter :: nAdv=1
!   integer(kind=intType) :: offset

!   ! Variables for spectral Radius
!   real(kind=realType) :: plim, rlim, clim2
!   real(kind=realType) :: cc2, qsi, qsj, qsk, sx, sy, sz, rmu
!   real(kind=realType) :: ri, rj, rk, rij, rjk, rki
!   real(kind=realType) :: vsi, vsj, vsk, rfl, dpi, dpj, dpk
!   real(kind=realType) :: sFace, tmp
!   logical :: doScaling

!   ! Variables for inviscid central flux
!   real(kind=realType) :: qsp, qsm, rqsp, rqsm, porVel, porFlux
!   real(kind=realType) :: pa, vnp, vnm, fs

!   ! Variables for inviscid diss flux scalar
!   real(kind=realType), parameter :: dssMax = 0.25_realType
!   real(kind=realType) :: sslim, rhoi
!   real(kind=realType) :: sfil, fis2, fis4
!   real(kind=realType) :: ppor, rrad, dis2, dis4
!   real(kind=realType) :: ddw1,ddw2,ddw3,ddw4,ddw5

!   ! Variables for speed of sound
!   logical :: correctForK
!   real(kind=realType) :: pp
!   real(kind=realType), parameter :: twoThird = two*third

!   ! Variables for nodal gradients
!   real(kind=realType) :: a2, oVol, uBar, vBar, wBar

!   ! Variables for viscous flux
!   real(kind=realType) :: rFilv, por, mul, mue, mut, heatCoef
!   real(kind=realType) :: gm1, factLamHeat, factTurbHeat
!   real(kind=realType) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z
!   real(kind=realType) :: q_x, q_y, q_z
!   real(kind=realType) :: corr, ssx, ssy, ssz, fracDiv, snrm
!   real(kind=realType) :: tauxx, tauyy, tauzz
!   real(kind=realType) :: tauxy, tauxz, tauyz
!   real(kind=realType) :: exx, eyy, ezz
!   real(kind=realType) :: exy, exz, eyz
!   real(kind=realType) :: Wxx, Wyy, Wzz
!   real(kind=realType) :: Wxy, Wxz, Wyz, Wyx, Wzx, Wzy
!   real(kind=realType) :: den, Ccr1
!   real(kind=realType) :: fmx, fmy, fmz, frhoE

!   ! Variables for final summing
!   integer(kind=intType) :: nTurb

!   ! Note that for dw/fw only the 2:BS+1 values are significant. The
!   ! lower "1" padding is so that we can safely scatter the fluxes
!   ! somwhere without issue.

!   ! ---------------------------------------------
!   !              Metric computation
!   ! ---------------------------------------------

!   do k=1,kl
!      n = k -1
!      do j=1,jl
!         m = j -1
!         do i=0,il

!            ! Determine the two diagonal vectors of the face.
!            v2(1) = x(i,j,k,1) - x(i,m,n,1)
!            v2(2) = x(i,j,k,2) - x(i,m,n,2)
!            v2(3) = x(i,j,k,3) - x(i,m,n,3)


!            v1(1) = x(i,j,n,1) - x(i,m,k,1)
!            v1(2) = x(i,j,n,2) - x(i,m,k,2)
!            v1(3) = x(i,j,n,3) - x(i,m,k,3)

!            ! The face normal, which is the cross product of the two
!            ! diagonal vectors times fact; remember that fact is
!            ! either -0.5 or 0.5.

!            si(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
!            si(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
!            si(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

!         enddo
!      enddo
!   enddo


!   ! Projected areas of cell faces in the j direction.

!   do k=1,kl
!      n = k -1
!      do j=0,jl
!         do i=1,il
!            l = i -1

!            ! Determine the two diagonal vectors of the face.

!            v1(1) = x(i,j,n,1) - x(l,j,k,1)
!            v1(2) = x(i,j,n,2) - x(l,j,k,2)
!            v1(3) = x(i,j,n,3) - x(l,j,k,3)

!            v2(1) = x(l,j,n,1) - x(i,j,k,1)
!            v2(2) = x(l,j,n,2) - x(i,j,k,2)
!            v2(3) = x(l,j,n,3) - x(i,j,k,3)

!            ! The face normal, which is the cross product of the two
!            ! diagonal vectors times fact; remember that fact is
!            ! either -0.5 or 0.5.

!            sj(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
!            sj(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
!            sj(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

!         enddo
!      enddo
!   enddo

!   ! Projected areas of cell faces in the k direction.

!   do k=0,kl
!      do j=1,jl
!         m = j -1
!         do i=1,il
!            l = i -1

!            ! Determine the two diagonal vectors of the face.

!            v1(1) = x(i,j,k,1) - x(l,m,k,1)
!            v1(2) = x(i,j,k,2) - x(l,m,k,2)
!            v1(3) = x(i,j,k,3) - x(l,m,k,3)

!            v2(1) = x(l,j,k,1) - x(i,m,k,1)
!            v2(2) = x(l,j,k,2) - x(i,m,k,2)
!            v2(3) = x(l,j,k,3) - x(i,m,k,3)

!            ! The face normal, which is the cross product of the two
!            ! diagonal vectors times fact; remember that fact is
!            ! either -0.5 or 0.5.

!            sk(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
!            sk(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
!            sk(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

!         enddo
!      enddo
!   enddo

!   ! ---------------------------------------------
!   !                     Init Res
!   ! ---------------------------------------------

!   ! Obviously this needs to be more complex for the actual code.
!   dw = zero

!   ! ---------------------------------------------
!   !                    SA Source Term
!   ! ---------------------------------------------

!   ! Set model constants
!   cv13    = rsaCv1**3
!   kar2Inv = one/(rsaK**2)
!   cw36    = rsaCw3**6
!   cb3Inv  = one/rsaCb3

!   ! Determine the non-dimensional wheel speed of this block.

!   omegax = timeRef*sections(sectionID)%rotRate(1)
!   omegay = timeRef*sections(sectionID)%rotRate(2)
!   omegaz = timeRef*sections(sectionID)%rotRate(3)

!   do k=2, kl
!      do j=2, jl
!         do i=2, il

!            ! Compute the gradient of u in the cell center. Use is made
!            ! of the fact that the surrounding normals sum up to zero,
!            ! such that the cell i,j,k does not give a contribution.
!            ! The gradient is scaled by the factor 2*vol.

!            uux = w(i+1,j,k,ivx)*si(i,j,k,1) - w(i-1,j,k,ivx)*si(i-1,j,k,1) &
!                 + w(i,j+1,k,ivx)*sj(i,j,k,1) - w(i,j-1,k,ivx)*sj(i,j-1,k,1) &
!                 + w(i,j,k+1,ivx)*sk(i,j,k,1) - w(i,j,k-1,ivx)*sk(i,j,k-1,1)
!            uuy = w(i+1,j,k,ivx)*si(i,j,k,2) - w(i-1,j,k,ivx)*si(i-1,j,k,2) &
!                 + w(i,j+1,k,ivx)*sj(i,j,k,2) - w(i,j-1,k,ivx)*sj(i,j-1,k,2) &
!                 + w(i,j,k+1,ivx)*sk(i,j,k,2) - w(i,j,k-1,ivx)*sk(i,j,k-1,2)
!            uuz = w(i+1,j,k,ivx)*si(i,j,k,3) - w(i-1,j,k,ivx)*si(i-1,j,k,3) &
!                 + w(i,j+1,k,ivx)*sj(i,j,k,3) - w(i,j-1,k,ivx)*sj(i,j-1,k,3) &
!                 + w(i,j,k+1,ivx)*sk(i,j,k,3) - w(i,j,k-1,ivx)*sk(i,j,k-1,3)

!            ! Idem for the gradient of v.

!            vvx = w(i+1,j,k,ivy)*si(i,j,k,1) - w(i-1,j,k,ivy)*si(i-1,j,k,1) &
!                 + w(i,j+1,k,ivy)*sj(i,j,k,1) - w(i,j-1,k,ivy)*sj(i,j-1,k,1) &
!                 + w(i,j,k+1,ivy)*sk(i,j,k,1) - w(i,j,k-1,ivy)*sk(i,j,k-1,1)
!            vvy = w(i+1,j,k,ivy)*si(i,j,k,2) - w(i-1,j,k,ivy)*si(i-1,j,k,2) &
!                 + w(i,j+1,k,ivy)*sj(i,j,k,2) - w(i,j-1,k,ivy)*sj(i,j-1,k,2) &
!                 + w(i,j,k+1,ivy)*sk(i,j,k,2) - w(i,j,k-1,ivy)*sk(i,j,k-1,2)
!            vvz = w(i+1,j,k,ivy)*si(i,j,k,3) - w(i-1,j,k,ivy)*si(i-1,j,k,3) &
!                 + w(i,j+1,k,ivy)*sj(i,j,k,3) - w(i,j-1,k,ivy)*sj(i,j-1,k,3) &
!                 + w(i,j,k+1,ivy)*sk(i,j,k,3) - w(i,j,k-1,ivy)*sk(i,j,k-1,3)

!            ! And for the gradient of w.

!            wwx = w(i+1,j,k,ivz)*si(i,j,k,1) - w(i-1,j,k,ivz)*si(i-1,j,k,1) &
!                 + w(i,j+1,k,ivz)*sj(i,j,k,1) - w(i,j-1,k,ivz)*sj(i,j-1,k,1) &
!                 + w(i,j,k+1,ivz)*sk(i,j,k,1) - w(i,j,k-1,ivz)*sk(i,j,k-1,1)
!            wwy = w(i+1,j,k,ivz)*si(i,j,k,2) - w(i-1,j,k,ivz)*si(i-1,j,k,2) &
!                 + w(i,j+1,k,ivz)*sj(i,j,k,2) - w(i,j-1,k,ivz)*sj(i,j-1,k,2) &
!                 + w(i,j,k+1,ivz)*sk(i,j,k,2) - w(i,j,k-1,ivz)*sk(i,j,k-1,2)
!            wwz = w(i+1,j,k,ivz)*si(i,j,k,3) - w(i-1,j,k,ivz)*si(i-1,j,k,3) &
!                 + w(i,j+1,k,ivz)*sj(i,j,k,3) - w(i,j-1,k,ivz)*sj(i,j-1,k,3) &
!                 + w(i,j,k+1,ivz)*sk(i,j,k,3) - w(i,j,k-1,ivz)*sk(i,j,k-1,3)

!            ! Compute the components of the stress tensor.
!            ! The combination of the current scaling of the velocity
!            ! gradients (2*vol) and the definition of the stress tensor,
!            ! leads to the factor 1/(4*vol).

!            fact = fourth/vol(i,j,k)

!            if (turbProd .eq. strain) then

!               sxx = two*fact*uux
!               syy = two*fact*vvy
!               szz = two*fact*wwz

!               sxy = fact*(uuy + vvx)
!               sxz = fact*(uuz + wwx)
!               syz = fact*(vvz + wwy)

!               ! Compute 2/3 * divergence of velocity squared

!               div2 = f23*(sxx+syy+szz)**2

!               ! Compute strain production term

!               strainMag2 = two*(sxy**2 + sxz**2 + syz**2) &
!                    +           sxx**2 + syy**2 + szz**2

!               strainProd = two*strainMag2 - div2

!               sqrtProd = sqrt(strainProd)

!            else if (turbProd .eq. vorticity) then

!               ! Compute the three components of the vorticity vector.
!               ! Substract the part coming from the rotating frame.

!               vortx = two*fact*(wwy - vvz) - two*omegax
!               vorty = two*fact*(uuz - wwx) - two*omegay
!               vortz = two*fact*(vvx - uuy) - two*omegaz

!               ! Compute the vorticity production term

!               vortProd = vortx**2 + vorty**2 + vortz**2

!               ! First take the square root of the production term to
!               ! obtain the correct production term for spalart-allmaras.
!               ! We do this to avoid if statements.

!               sqrtProd = sqrt(vortProd)

!            end if

!            ! Compute the laminar kinematic viscosity, the inverse of
!            ! wall distance squared, the ratio chi (ratio of nuTilde
!            ! and nu) and the functions fv1 and fv2. The latter corrects
!            ! the production term near a viscous wall.

!            nu       = rlv(i,j,k)/w(i,j,k,irho)
!            dist2Inv = one/(d2Wall(i,j,k)**2)
!            chi      = w(i,j,k,itu1)/nu
!            chi2     = chi*chi
!            chi3     = chi*chi2
!            fv1      = chi3/(chi3+cv13)
!            fv2      = one - chi/(one + chi*fv1)

!            ! The function ft2, which is designed to keep a laminar
!            ! solution laminar. When running in fully turbulent mode
!            ! this function should be set to 0.0.

!            if (useft2SA) then
!               ft2 = rsaCt3*exp(-rsaCt4*chi2)
!            else
!               ft2 = zero
!            end if

!            ! Correct the production term to account for the influence
!            ! of the wall.

!            sst = sqrtProd + w(i,j,k,itu1)*fv2*kar2Inv*dist2Inv

!            ! Add rotation term (useRotationSA defined in inputParams.F90)

!            if (useRotationSA) then
!               sst = sst + rsaCrot*min(zero,sqrt(two*strainMag2))
!            end if

!            ! Make sure that this term remains positive
!            ! (the function fv2 is negative between chi = 1 and 18.4,
!            ! which can cause sst to go negative, which is undesirable).

!            sst = max(sst,xminn)

!            ! Compute the function fw. The argument rr is cut off at 10
!            ! to avoid numerical problems. This is ok, because the
!            ! asymptotical value of fw is then already reached.

!            rr     = w(i,j,k,itu1)*kar2Inv*dist2Inv/sst
!            rr     = min(rr,10.0_realType)
!            gg     = rr + rsaCw2*(rr**6 - rr)
!            gg6    = gg**6
!            termFw = ((one + cw36)/(gg6 + cw36))**sixth
!            fwSa   = gg*termFw

!            ! Compute the source term; some terms are saved for the
!            ! linearization. The source term is stored in dvt.

!            term1 = rsaCb1*(one-ft2)*sqrtProd
!            term2 = dist2Inv*(kar2Inv*rsaCb1*((one-ft2)*fv2 + ft2) &
!                 -           rsaCw1*fwSa)

!            dw(i, j, k, idvt) = (term1 + term2*w(i,j,k,itu1))*w(i,j,k,itu1)

!         enddo
!      enddo
!   enddo

!   ! ---------------------------------------------
!   !                    SA Viscous Term
!   ! ---------------------------------------------
!   !
!   !       Viscous terms in k-direction.
!   !
!   do k=2, kl
!      do j=2, jl
!         do i=2, il

!            ! Compute the metrics in zeta-direction, i.e. along the
!            ! line k = constant.

!            voli  = one/vol(i,j,k)
!            volmi = two/(vol(i,j,k) + vol(i,j,k-1))
!            volpi = two/(vol(i,j,k) + vol(i,j,k+1))

!            xm = sk(i,j,k-1,1)*volmi
!            ym = sk(i,j,k-1,2)*volmi
!            zm = sk(i,j,k-1,3)*volmi
!            xp = sk(i,j,k,  1)*volpi
!            yp = sk(i,j,k,  2)*volpi
!            zp = sk(i,j,k,  3)*volpi

!            xa  = half*(sk(i,j,k,1) + sk(i,j,k-1,1))*voli
!            ya  = half*(sk(i,j,k,2) + sk(i,j,k-1,2))*voli
!            za  = half*(sk(i,j,k,3) + sk(i,j,k-1,3))*voli
!            ttm = xm*xa + ym*ya + zm*za
!            ttp = xp*xa + yp*ya + zp*za

!            ! Computation of the viscous terms in zeta-direction; note
!            ! that cross-derivatives are neglected, i.e. the mesh is
!            ! assumed to be orthogonal.
!            ! Furthermore, the grad(nu)**2 has been rewritten as
!            ! div(nu grad(nu)) - nu div(grad nu) to enhance stability.
!            ! The second derivative in zeta-direction is constructed as
!            ! the central difference of the first order derivatives, i.e.
!            ! d^2/dzeta^2 = d/dzeta (d/dzeta k+1/2 - d/dzeta k-1/2).
!            ! In this way the metric can be taken into account.

!            ! Compute the diffusion coefficients multiplying the nodes
!            ! k+1, k and k-1 in the second derivative. Make sure that
!            ! these coefficients are nonnegative.

!            cnud = -rsaCb2*w(i,j,k,itu1)*cb3Inv
!            cam  =  ttm*cnud
!            cap  =  ttp*cnud

!            nutm = half*(w(i,j,k-1,itu1) + w(i,j,k,itu1))
!            nutp = half*(w(i,j,k+1,itu1) + w(i,j,k,itu1))
!            nu   = rlv(i,j,k)/w(i,j,k,irho)
!            num  = half*(rlv(i,j,k-1)/w(i,j,k-1,irho) + nu)
!            nup  = half*(rlv(i,j,k+1)/w(i,j,k+1,irho) + nu)
!            cdm  = (num + (one + rsaCb2)*nutm)*ttm*cb3Inv
!            cdp  = (nup + (one + rsaCb2)*nutp)*ttp*cb3Inv

!            c1m = max(cdm+cam, zero)
!            c1p = max(cdp+cap, zero)
!            c10 = c1m + c1p

!            ! Update the residual for this cell and store the possible
!            ! coefficients for the matrix in b1, c1 and d1.

!            dw(i,j,k,idvt) = dw(i,j,k,idvt)      + c1m*w(i,j,k-1,itu1) &
!                 - c10*w(i,j,k,itu1) + c1p*w(i,j,k+1,itu1)
!         end do
!      enddo
!   enddo
!   !
!   !       Viscous terms in j-direction.
!   !
!   do k=2, kl
!      do j=2, jl
!         do i=2, il

!            ! Compute the metrics in eta-direction, i.e. along the
!            ! line j = constant.

!            voli  = one/vol(i,j,k)
!            volmi = two/(vol(i,j,k) + vol(i,j-1,k))
!            volpi = two/(vol(i,j,k) + vol(i,j+1,k))

!            xm = sj(i,j-1,k,1)*volmi
!            ym = sj(i,j-1,k,2)*volmi
!            zm = sj(i,j-1,k,3)*volmi
!            xp = sj(i,j,  k,1)*volpi
!            yp = sj(i,j,  k,2)*volpi
!            zp = sj(i,j,  k,3)*volpi

!            xa  = half*(sj(i,j,k,1) + sj(i,j-1,k,1))*voli
!            ya  = half*(sj(i,j,k,2) + sj(i,j-1,k,2))*voli
!            za  = half*(sj(i,j,k,3) + sj(i,j-1,k,3))*voli
!            ttm = xm*xa + ym*ya + zm*za
!            ttp = xp*xa + yp*ya + zp*za

!            ! Computation of the viscous terms in eta-direction; note
!            ! that cross-derivatives are neglected, i.e. the mesh is
!            ! assumed to be orthogonal.
!            ! Furthermore, the grad(nu)**2 has been rewritten as
!            ! div(nu grad(nu)) - nu div(grad nu) to enhance stability.
!            ! The second derivative in eta-direction is constructed as
!            ! the central difference of the first order derivatives, i.e.
!            ! d^2/deta^2 = d/deta (d/deta j+1/2 - d/deta j-1/2).
!            ! In this way the metric can be taken into account.

!            ! Compute the diffusion coefficients multiplying the nodes
!            ! j+1, j and j-1 in the second derivative. Make sure that
!            ! these coefficients are nonnegative.

!            cnud = -rsaCb2*w(i,j,k,itu1)*cb3Inv
!            cam  =  ttm*cnud
!            cap  =  ttp*cnud

!            nutm = half*(w(i,j-1,k,itu1) + w(i,j,k,itu1))
!            nutp = half*(w(i,j+1,k,itu1) + w(i,j,k,itu1))
!            nu   = rlv(i,j,k)/w(i,j,k,irho)
!            num  = half*(rlv(i,j-1,k)/w(i,j-1,k,irho) + nu)
!            nup  = half*(rlv(i,j+1,k)/w(i,j+1,k,irho) + nu)
!            cdm  = (num + (one + rsaCb2)*nutm)*ttm*cb3Inv
!            cdp  = (nup + (one + rsaCb2)*nutp)*ttp*cb3Inv

!            c1m = max(cdm+cam, zero)
!            c1p = max(cdp+cap, zero)
!            c10 = c1m + c1p

!            ! Update the residual for this cell and store the possible
!            ! coefficients for the matrix in b1, c1 and d1.

!            dw(i,j,k,idvt) = dw(i,j,k,idvt)      + c1m*w(i,j-1,k,itu1) &
!                 - c10*w(i,j,k,itu1) + c1p*w(i,j+1,k,itu1)

!         enddo
!      enddo
!   enddo
!   !
!   !       Viscous terms in i-direction.
!   !
!   do k=2, kl
!      do j=2, jl
!         do i=2, il

!            ! Compute the metrics in xi-direction, i.e. along the
!            ! line i = constant.

!            voli  = one/vol(i,j,k)
!            volmi = two/(vol(i,j,k) + vol(i-1,j,k))
!            volpi = two/(vol(i,j,k) + vol(i+1,j,k))

!            xm = si(i-1,j,k,1)*volmi
!            ym = si(i-1,j,k,2)*volmi
!            zm = si(i-1,j,k,3)*volmi
!            xp = si(i,  j,k,1)*volpi
!            yp = si(i,  j,k,2)*volpi
!            zp = si(i,  j,k,3)*volpi

!            xa  = half*(si(i,j,k,1) + si(i-1,j,k,1))*voli
!            ya  = half*(si(i,j,k,2) + si(i-1,j,k,2))*voli
!            za  = half*(si(i,j,k,3) + si(i-1,j,k,3))*voli
!            ttm = xm*xa + ym*ya + zm*za
!            ttp = xp*xa + yp*ya + zp*za

!            ! Computation of the viscous terms in xi-direction; note
!            ! that cross-derivatives are neglected, i.e. the mesh is
!            ! assumed to be orthogonal.
!            ! Furthermore, the grad(nu)**2 has been rewritten as
!            ! div(nu grad(nu)) - nu div(grad nu) to enhance stability.
!            ! The second derivative in xi-direction is constructed as
!            ! the central difference of the first order derivatives, i.e.
!            ! d^2/dxi^2 = d/dxi (d/dxi i+1/2 - d/dxi i-1/2).
!            ! In this way the metric can be taken into account.

!            ! Compute the diffusion coefficients multiplying the nodes
!            ! i+1, i and i-1 in the second derivative. Make sure that
!            ! these coefficients are nonnegative.

!            cnud = -rsaCb2*w(i,j,k,itu1)*cb3Inv
!            cam  =  ttm*cnud
!            cap  =  ttp*cnud

!            nutm = half*(w(i-1,j,k,itu1) + w(i,j,k,itu1))
!            nutp = half*(w(i+1,j,k,itu1) + w(i,j,k,itu1))
!            nu   = rlv(i,j,k)/w(i,j,k,irho)
!            num  = half*(rlv(i-1,j,k)/w(i-1,j,k,irho) + nu)
!            nup  = half*(rlv(i+1,j,k)/w(i+1,j,k,irho) + nu)
!            cdm  = (num + (one + rsaCb2)*nutm)*ttm*cb3Inv
!            cdp  = (nup + (one + rsaCb2)*nutp)*ttp*cb3Inv

!            c1m = max(cdm+cam, zero)
!            c1p = max(cdp+cap, zero)
!            c10 = c1m + c1p

!            ! Update the residual for this cell and store the possible
!            ! coefficients for the matrix in b1, c1 and d1.

!            dw(i,j,k,idvt) = dw(i,j,k,idvt)      + c1m*w(i-1,j,k,itu1) &
!                 - c10*w(i,j,k,itu1) + c1p*w(i+1,j,k,itu1)
!         enddo
!      enddo
!   enddo

!   ! ---------------------------------------------
!   !                SA Advection
!   ! ---------------------------------------------
!   offset=itu1-1
!   do k=2, kl
!      do j=2, jl
!         do i=2, il

!            ! Compute the grid velocity if present.
!            ! It is taken as the average of k and k-1,

!            voli = half/vol(i,j,k)

!            ! if( addGridVelocities ) &
!            !      qs = (sFaceK(i,j,k) + sFaceK(i,j,k-1))*voli

!            ! Compute the normal velocity, where the normal direction
!            ! is taken as the average of faces k and k-1.

!            xa = (sk(i,j,k,1) + sk(i,j,k-1,1))*voli
!            ya = (sk(i,j,k,2) + sk(i,j,k-1,2))*voli
!            za = (sk(i,j,k,3) + sk(i,j,k-1,3))*voli

!            uu = xa*w(i,j,k,ivx) + ya*w(i,j,k,ivy) + za*w(i,j,k,ivz) - qs

!            ! Determine the situation we are having here, i.e. positive
!            ! or negative normal velocity.

!            velKdir: if(uu > zero) then

!               ! Velocity has a component in positive k-direction.
!               ! Loop over the number of advection equations.

!               do ii=1,nAdv

!                  ! Set the value of jj such that it corresponds to the
!                  ! turbulent entry in w.

!                  jj = ii + offset

!                  ! Check whether a first or a second order discretization
!                  ! must be used.

!                  if( secondOrd ) then

!                     ! Second order; store the three differences for the
!                     ! discretization of the derivative in k-direction.

!                     dwtm1 = w(i,j,k-1,jj) - w(i,j,k-2,jj)
!                     dwt   = w(i,j,k,  jj) - w(i,j,k-1,jj)
!                     dwtp1 = w(i,j,k+1,jj) - w(i,j,k,  jj)

!                     ! Construct the derivative in this cell center. This
!                     ! is the first order upwind derivative with two
!                     ! nonlinear corrections.

!                     dwtk = dwt

!                     if(dwt*dwtp1 > zero) then
!                        if(abs(dwt) < abs(dwtp1)) then
!                           dwtk = dwtk + half*dwt
!                        else
!                           dwtk = dwtk + half*dwtp1
!                        endif
!                     endif

!                     if(dwt*dwtm1 > zero) then
!                        if(abs(dwt) < abs(dwtm1)) then
!                           dwtk = dwtk - half*dwt
!                        else
!                           dwtk = dwtk - half*dwtm1
!                        endif
!                     endif

!                  else

!                     ! 1st order upwind scheme.

!                     dwtk = w(i,j,k,jj) - w(i,j,k-1,jj)

!                  endif

!                  ! Update the residual. The convective term must be
!                  ! substracted, because it appears on the other side of
!                  ! the equation as the source and viscous terms.

!                  dw(i,j,k,idvt+ii-1) = dw(i,j,k,idvt+ii-1) - uu*dwtk
!               enddo

!            else velKdir

!               ! Velocity has a component in negative k-direction.
!               ! Loop over the number of advection equations
!               do ii=1,nAdv

!                  ! Set the value of jj such that it corresponds to the
!                  ! turbulent entry in w.

!                  jj = ii + offset

!                  ! Check whether a first or a second order discretization
!                  ! must be used.

!                  if( secondOrd ) then

!                     ! Store the three differences for the discretization of
!                     ! the derivative in k-direction.

!                     dwtm1 = w(i,j,k,  jj) - w(i,j,k-1,jj)
!                     dwt   = w(i,j,k+1,jj) - w(i,j,k,  jj)
!                     dwtp1 = w(i,j,k+2,jj) - w(i,j,k+1,jj)

!                     ! Construct the derivative in this cell center. This is
!                     ! the first order upwind derivative with two nonlinear
!                     ! corrections.

!                     dwtk = dwt

!                     if(dwt*dwtp1 > zero) then
!                        if(abs(dwt) < abs(dwtp1)) then
!                           dwtk = dwtk - half*dwt
!                        else
!                           dwtk = dwtk - half*dwtp1
!                        endif
!                     endif

!                     if(dwt*dwtm1 > zero) then
!                        if(abs(dwt) < abs(dwtm1)) then
!                           dwtk = dwtk + half*dwt
!                        else
!                           dwtk = dwtk + half*dwtm1
!                        endif
!                     endif

!                  else

!                     ! 1st order upwind scheme.

!                     dwtk = w(i,j,k+1,jj) - w(i,j,k,jj)

!                  endif

!                  ! Update the residual. The convective term must be
!                  ! substracted, because it appears on the other side
!                  ! of the equation as the source and viscous terms.

!                  dw(i,j,k,idvt+ii-1) = dw(i,j,k,idvt+ii-1) - uu*dwtk
!               end do
!            endif velKdir
!         enddo
!      enddo
!   enddo

!   !
!   !       Upwind discretization of the convective term in j (eta)
!   !       direction. Either the 1st order upwind or the second order
!   !       fully upwind interpolation scheme, kappa = -1, is used in
!   !       combination with the minmod limiter.
!   !       The possible grid velocity must be taken into account.
!   !
!   do k=2, kl
!      do j=2, jl
!         do i=2, il


!            ! Compute the grid velocity if present.
!            ! It is taken as the average of j and j-1,

!            voli = half/vol(i,j,k)
!            ! if( addGridVelocities ) &
!            !      qs = (sFaceJ(i,j,k) + sFaceJ(i,j-1,k))*voli

!            ! Compute the normal velocity, where the normal direction
!            ! is taken as the average of faces j and j-1.

!            xa = (sj(i,j,k,1) + sj(i,j-1,k,1))*voli
!            ya = (sj(i,j,k,2) + sj(i,j-1,k,2))*voli
!            za = (sj(i,j,k,3) + sj(i,j-1,k,3))*voli

!            uu = xa*w(i,j,k,ivx) + ya*w(i,j,k,ivy) + za*w(i,j,k,ivz) - qs

!            ! Determine the situation we are having here, i.e. positive
!            ! or negative normal velocity.

!            velJdir: if(uu > zero) then

!               ! Velocity has a component in positive j-direction.
!               ! Loop over the number of advection equations.
!               do ii=1,nAdv

!                  ! Set the value of jj such that it corresponds to the
!                  ! turbulent entry in w.

!                  jj = ii + offset

!                  ! Check whether a first or a second order discretization
!                  ! must be used.

!                  if( secondOrd ) then

!                     ! Second order; store the three differences for the
!                     ! discretization of the derivative in j-direction.

!                     dwtm1 = w(i,j-1,k,jj) - w(i,j-2,k,jj)
!                     dwt   = w(i,j,  k,jj) - w(i,j-1,k,jj)
!                     dwtp1 = w(i,j+1,k,jj) - w(i,j,  k,jj)

!                     ! Construct the derivative in this cell center. This is
!                     ! the first order upwind derivative with two nonlinear
!                     ! corrections.

!                     dwtj = dwt

!                     if(dwt*dwtp1 > zero) then
!                        if(abs(dwt) < abs(dwtp1)) then
!                           dwtj = dwtj + half*dwt
!                        else
!                           dwtj = dwtj + half*dwtp1
!                        endif
!                     endif

!                     if(dwt*dwtm1 > zero) then
!                        if(abs(dwt) < abs(dwtm1)) then
!                           dwtj = dwtj - half*dwt
!                        else
!                           dwtj = dwtj - half*dwtm1
!                        endif
!                     endif

!                  else

!                     ! 1st order upwind scheme.

!                     dwtj = w(i,j,k,jj) - w(i,j-1,k,jj)

!                  endif

!                  ! Update the residual. The convective term must be
!                  ! substracted, because it appears on the other side of
!                  ! the equation as the source and viscous terms.

!                  dw(i,j,k,idvt+ii-1) = dw(i,j,k,idvt+ii-1) - uu*dwtj
!               enddo

!            else velJdir

!               ! Velocity has a component in negative j-direction.
!               ! Loop over the number of advection equations.
!               do ii=1,nAdv

!                  ! Set the value of jj such that it corresponds to the
!                  ! turbulent entry in w.

!                  jj = ii + offset

!                  ! Check whether a first or a second order discretization
!                  ! must be used.

!                  if( secondOrd ) then

!                     ! Store the three differences for the discretization of
!                     ! the derivative in j-direction.

!                     dwtm1 = w(i,j,  k,jj) - w(i,j-1,k,jj)
!                     dwt   = w(i,j+1,k,jj) - w(i,j,  k,jj)
!                     dwtp1 = w(i,j+2,k,jj) - w(i,j+1,k,jj)

!                     ! Construct the derivative in this cell center. This is
!                     ! the first order upwind derivative with two nonlinear
!                     ! corrections.

!                     dwtj = dwt

!                     if(dwt*dwtp1 > zero) then
!                        if(abs(dwt) < abs(dwtp1)) then
!                           dwtj = dwtj - half*dwt
!                        else
!                           dwtj = dwtj - half*dwtp1
!                        endif
!                     endif

!                     if(dwt*dwtm1 > zero) then
!                        if(abs(dwt) < abs(dwtm1)) then
!                           dwtj = dwtj + half*dwt
!                        else
!                           dwtj = dwtj + half*dwtm1
!                        endif
!                     endif

!                  else

!                     ! 1st order upwind scheme.

!                     dwtj = w(i,j+1,k,jj) - w(i,j,k,jj)

!                  endif

!                  ! Update the residual. The convective term must be
!                  ! substracted, because it appears on the other side
!                  ! of the equation as the source and viscous terms.

!                  dw(i,j,k,idvt+ii-1) = dw(i,j,k,idvt+ii-1) - uu*dwtj
!               enddo
!            endif velJdir
!         enddo
!      enddo
!   enddo
!   !
!   !       Upwind discretization of the convective term in i (xi)
!   !       direction. Either the 1st order upwind or the second order
!   !       fully upwind interpolation scheme, kappa = -1, is used in
!   !       combination with the minmod limiter.
!   !       The possible grid velocity must be taken into account.
!   !
!   qs = zero
!   do k=2, kl
!      do j=2, jl
!         do i=2, il
!            ! Compute the grid velocity if present.
!            ! It is taken as the average of i and i-1,

!            voli = half/vol(i,j,k)
!            ! if( addGridVelocities ) &
!            !         qs = (sFaceI(i,j,k) + sFaceI(i-1,j,k))*voli

!            ! Compute the normal velocity, where the normal direction
!            ! is taken as the average of faces i and i-1.

!            xa = (si(i,j,k,1) + si(i-1,j,k,1))*voli
!            ya = (si(i,j,k,2) + si(i-1,j,k,2))*voli
!            za = (si(i,j,k,3) + si(i-1,j,k,3))*voli

!            uu = xa*w(i,j,k,ivx) + ya*w(i,j,k,ivy) + za*w(i,j,k,ivz) - qs

!            ! Determine the situation we are having here, i.e. positive
!            ! or negative normal velocity.

!            velIdir: if(uu > zero) then

!               ! Velocity has a component in positive i-direction.
!               ! Loop over the number of advection equations.
!               do ii=1,nAdv

!                  ! Set the value of jj such that it corresponds to the
!                  ! turbulent entry in w.

!                  jj = ii + offset

!                  ! Check whether a first or a second order discretization
!                  ! must be used.

!                  if( secondOrd ) then

!                     ! Second order; store the three differences for the
!                     ! discretization of the derivative in i-direction.

!                     dwtm1 = w(i-1,j,k,jj) - w(i-2,j,k,jj)
!                     dwt   = w(i,  j,k,jj) - w(i-1,j,k,jj)
!                     dwtp1 = w(i+1,j,k,jj) - w(i,  j,k,jj)

!                     ! Construct the derivative in this cell center. This is
!                     ! the first order upwind derivative with two nonlinear
!                     ! corrections.

!                     dwti = dwt

!                     if(dwt*dwtp1 > zero) then
!                        if(abs(dwt) < abs(dwtp1)) then
!                           dwti = dwti + half*dwt
!                        else
!                           dwti = dwti + half*dwtp1
!                        endif
!                     endif

!                     if(dwt*dwtm1 > zero) then
!                        if(abs(dwt) < abs(dwtm1)) then
!                           dwti = dwti - half*dwt
!                        else
!                           dwti = dwti - half*dwtm1
!                        endif
!                     endif

!                  else

!                     ! 1st order upwind scheme.

!                     dwti = w(i,j,k,jj) - w(i-1,j,k,jj)

!                  endif

!                  ! Update the residual. The convective term must be
!                  ! substracted, because it appears on the other side of
!                  ! the equation as the source and viscous terms.

!                  dw(i,j,k,idvt+ii-1) = dw(i,j,k,idvt+ii-1) - uu*dwti
!               enddo

!            else velIdir

!               ! Velocity has a component in negative i-direction.
!               ! Loop over the number of advection equations.
!               do ii=1,nAdv

!                  ! Set the value of jj such that it corresponds to the
!                  ! turbulent entry in w.

!                  jj = ii + offset

!                  ! Check whether a first or a second order discretization
!                  ! must be used.

!                  if( secondOrd ) then

!                     ! Second order; store the three differences for the
!                     ! discretization of the derivative in i-direction.

!                     dwtm1 = w(i,  j,k,jj) - w(i-1,j,k,jj)
!                     dwt   = w(i+1,j,k,jj) - w(i,  j,k,jj)
!                     dwtp1 = w(i+2,j,k,jj) - w(i+1,j,k,jj)

!                     ! Construct the derivative in this cell center. This is
!                     ! the first order upwind derivative with two nonlinear
!                     ! corrections.

!                     dwti = dwt

!                     if(dwt*dwtp1 > zero) then
!                        if(abs(dwt) < abs(dwtp1)) then
!                           dwti = dwti - half*dwt
!                        else
!                           dwti = dwti - half*dwtp1
!                        endif
!                     endif

!                     if(dwt*dwtm1 > zero) then
!                        if(abs(dwt) < abs(dwtm1)) then
!                           dwti = dwti + half*dwt
!                        else
!                           dwti = dwti + half*dwtm1
!                        endif
!                     endif

!                  else

!                     ! 1st order upwind scheme.

!                     dwti = w(i+1,j,k,jj) - w(i,j,k,jj)

!                  endif

!                  ! Update the residual. The convective term must be
!                  ! substracted, because it appears on the other side
!                  ! of the equation as the source and viscous terms.

!                  dw(i,j,k,idvt+ii-1) = dw(i,j,k,idvt+ii-1) - uu*dwti

!                  ! Update the central jacobian. First the term which is
!                  ! always present, i.e. -uu.
!               enddo

!            endif velIdir
!         enddo
!      enddo
!   enddo

!   ! ---------------------------------------------
!   !               Spectral Radius
!   ! ---------------------------------------------

!   ! Set the value of plim. To be fully consistent this must have
!   ! the dimension of a pressure. Therefore a fraction of pInfCorr
!   ! is used. Idem for rlim; compute clim2 as well.

!   plim  = 0.001_realType*pInfCorr
!   rlim  = 0.001_realType*rhoInf
!   clim2 = 0.000001_realType*gammaInf*pInfCorr/rhoInf
!   doScaling = .True.

!   ! Initialize sFace to zero. This value will be used if the
!   ! block is not moving.

!   sFace = zero
!   !
!   !           Inviscid contribution, depending on the preconditioner.
!   !           Compute the cell centered values of the spectral radii.
!   !
!   do k=1,kl
!      do j=1,jl
!         do i=1,il

!            ! Compute the velocities and speed of sound squared.

!            uux  = w(i,j,k,ivx)
!            uuy  = w(i,j,k,ivy)
!            uuz  = w(i,j,k,ivz)
!            cc2 = gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho)
!            cc2 = max(cc2,clim2)

!            ! Set the dot product of the grid velocity and the
!            ! normal in i-direction for a moving face. To avoid
!            ! a number of multiplications by 0.5 simply the sum
!            ! is taken.

!            ! if( addGridVelocities ) &
!            !      sFace = sFaceI(i-1,j,k) + sFaceI(i,j,k)

!            ! Spectral radius in i-direction.

!            sx = si(i-1,j,k,1) + si(i,j,k,1)
!            sy = si(i-1,j,k,2) + si(i,j,k,2)
!            sz = si(i-1,j,k,3) + si(i,j,k,3)

!            qsi = uux*sx + uuy*sy + uuz*sz - sFace

!            ri = half*(abs(qsi) &
!                 +       sqrt(cc2*(sx**2 + sy**2 + sz**2)))

!            ! The grid velocity in j-direction.

!            ! if( addGridVelocities ) &
!            !      sFace = sFaceJ(i,j-1,k) + sFaceJ(i,j,k)

!            ! Spectral radius in j-direction.

!            sx = sj(i,j-1,k,1) + sj(i,j,k,1)
!            sy = sj(i,j-1,k,2) + sj(i,j,k,2)
!            sz = sj(i,j-1,k,3) + sj(i,j,k,3)

!            qsj = uux*sx + uuy*sy + uuz*sz - sFace

!            rj = half*(abs(qsj) &
!                 +       sqrt(cc2*(sx**2 + sy**2 + sz**2)))

!            ! The grid velocity in k-direction.

!            ! if( addGridVelocities ) &
!            !      sFace = sFaceK(i,j,k-1) + sFaceK(i,j,k)

!            ! Spectral radius in k-direction.

!            sx = sk(i,j,k-1,1) + sk(i,j,k,1)
!            sy = sk(i,j,k-1,2) + sk(i,j,k,2)
!            sz = sk(i,j,k-1,3) + sk(i,j,k,3)

!            qsk = uux*sx + uuy*sy + uuz*sz - sFace

!            rk = half*(abs(qsk) &
!                 +       sqrt(cc2*(sx**2 + sy**2 + sz**2)))

!            !
!            !           Adapt the spectral radii if directional scaling must be
!            !           applied.
!            !

!            ! Avoid division by zero by clipping radi, radJ and
!            ! radK.

!            ri = max(ri, eps)
!            rj = max(rj, eps)
!            rk = max(rk, eps)

!            ! Compute the scaling in the three coordinate
!            ! directions.

!            rij = (ri/rj)**adis
!            rjk = (rj/rk)**adis
!            rki = (rk/ri)**adis

!            ! Create the scaled versions of the aspect ratios.
!            ! Note that the multiplication is done with radi, radJ
!            ! and radK, such that the influence of the clipping
!            ! is negligible.

!            radi(i,j,k) = ri*(one + one/rij + rki)
!            radJ(i,j,k) = rj*(one + one/rjk + rij)
!            radK(i,j,k) = rk*(one + one/rki + rjk)
!         end do
!      enddo
!   enddo

!   ! ---------------------------------------------
!   !               Inviscid central flux
!   ! ---------------------------------------------

!   ! Initialize sFace to zero. This value will be used if the
!   ! block is not moving.
!   sFace = zero
!   do k=2, kl
!      do j=2, jl
!         do i=1, nx

!            ! Set the dot product of the grid velocity and the
!            ! normal in i-direction for a moving face.

!            !if( addGridVelocities ) sFace = sFaceI(i,j,k)

!            ! Compute the normal velocities of the left and right state.

!            vnp = w(i+1,j,k,ivx)*sI(i,j,k,1) &
!                 + w(i+1,j,k,ivy)*sI(i,j,k,2) &
!                 + w(i+1,j,k,ivz)*sI(i,j,k,3)
!            vnm = w(i,  j,k,ivx)*sI(i,j,k,1) &
!                 + w(i,  j,k,ivy)*sI(i,j,k,2) &
!                 + w(i,  j,k,ivz)*sI(i,j,k,3)
!            ! Set the values of the porosities for this face.
!            ! porVel defines the porosity w.r.t. velocity;
!            ! porFlux defines the porosity w.r.t. the entire flux.
!            ! The latter is only zero for a discontinuous block
!            ! boundary that must be treated conservatively.
!            ! The default value of porFlux is 0.5, such that the
!            ! correct central flux is scattered to both cells.
!            ! In case of a boundFlux the normal velocity is set
!            ! to sFace.

!            porVel  = one
!            porFlux = half
!            if(porI(i,j,k) == noFlux)    porFlux = zero
!            if(porI(i,j,k) == boundFlux) then
!               porVel = zero
!               vnp    = sFace
!               vnm    = sFace
!            endif

!            ! Incorporate porFlux in porVel.

!            porVel = porVel*porFlux

!            ! Compute the normal velocities relative to the grid for
!            ! the face as well as the mass fluxes.

!            qsp = (vnp -sFace)*porVel
!            qsm = (vnm -sFace)*porVel

!            rqsp = qsp*w(i+1,j,k,irho)
!            rqsm = qsm*w(i,  j,k,irho)

!            ! Compute the sum of the pressure multiplied by porFlux.
!            ! For the default value of porFlux, 0.5, this leads to
!            ! the average pressure.

!            pa = porFlux*(p(i+1,j,k) + p(i,j,k))

!            ! Compute the fluxes and scatter them to the cells
!            ! i,j,k and i+1,j,k. Store the density flux in the
!            ! mass flow of the appropriate sliding mesh interface.

!            fs = rqsp + rqsm
!            dw(i+1,j,k,irho) = dw(i+1,j,k,irho) - fs
!            dw(i,  j,k,irho) = dw(i,  j,k,irho) + fs

!            fs = rqsp*w(i+1,j,k,ivx) + rqsm*w(i,j,k,ivx) &
!                 + pa*sI(i,j,k,1)
!            dw(i+1,j,k,imx) = dw(i+1,j,k,imx) - fs
!            dw(i,  j,k,imx) = dw(i,  j,k,imx) + fs

!            fs = rqsp*w(i+1,j,k,ivy) + rqsm*w(i,j,k,ivy) &
!                 + pa*sI(i,j,k,2)
!            dw(i+1,j,k,imy) = dw(i+1,j,k,imy) - fs
!            dw(i,  j,k,imy) = dw(i,  j,k,imy) + fs

!            fs = rqsp*w(i+1,j,k,ivz) + rqsm*w(i,j,k,ivz) &
!                 + pa*sI(i,j,k,3)
!            dw(i+1,j,k,imz) = dw(i+1,j,k,imz) - fs
!            dw(i,  j,k,imz) = dw(i,  j,k,imz) + fs

!            fs = qsp*w(i+1,j,k,irhoE) + qsm*w(i,j,k,irhoE) &
!                 + porFlux*(vnp*p(i+1,j,k) + vnm*p(i,j,k))
!            dw(i+1,j,k,irhoE) = dw(i+1,j,k,irhoE) - fs
!            dw(i,  j,k,irhoE) = dw(i,  j,k,irhoE) + fs
!         enddo
!      enddo
!   enddo

!   sface = zero
!   do k=2,kl
!      do j=1,ny
!         do i=2,il

!            ! Set the dot product of the grid velocity and the
!            ! normal in j-direction for a moving face.

!            !if( addGridVelocities ) sFace = sFaceJ(i,j,k)

!            ! Compute the normal velocities of the left and right state.

!            vnp = w(i,j+1,k,ivx)*sJ(i,j,k,1) &
!                 + w(i,j+1,k,ivy)*sJ(i,j,k,2) &
!                 + w(i,j+1,k,ivz)*sJ(i,j,k,3)
!            vnm = w(i,j,  k,ivx)*sJ(i,j,k,1) &
!                 + w(i,j,  k,ivy)*sJ(i,j,k,2) &
!                 + w(i,j,  k,ivz)*sJ(i,j,k,3)

!            ! Set the values of the porosities for this face.
!            ! porVel defines the porosity w.r.t. velocity;
!            ! porFlux defines the porosity w.r.t. the entire flux.
!            ! The latter is only zero for a discontinuous block
!            ! boundary that must be treated conservatively.
!            ! The default value of porFlux is 0.5, such that the
!            ! correct central flux is scattered to both cells.
!            ! In case of a boundFlux the normal velocity is set
!            ! to sFace.

!            porVel  = one
!            porFlux = half
!            if(porJ(i,j,k) == noFlux)    porFlux = zero
!            if(porJ(i,j,k) == boundFlux) then
!               porVel = zero
!               vnp    = sFace
!               vnm    = sFace
!            endif

!            ! Incorporate porFlux in porVel.

!            porVel = porVel*porFlux

!            ! Compute the normal velocities for the face as well as the
!            ! mass fluxes.

!            qsp = (vnp - sFace)*porVel
!            qsm = (vnm - sFace)*porVel

!            rqsp = qsp*w(i,j+1,k,irho)
!            rqsm = qsm*w(i,j,  k,irho)

!            ! Compute the sum of the pressure multiplied by porFlux.
!            ! For the default value of porFlux, 0.5, this leads to
!            ! the average pressure.

!            pa = porFlux*(p(i,j+1,k) + p(i,j,k))

!            ! Compute the fluxes and scatter them to the cells
!            ! i,j,k and i,j+1,k. Store the density flux in the
!            ! mass flow of the appropriate sliding mesh interface.

!            fs = rqsp + rqsm
!            dw(i,j+1,k,irho) = dw(i,j+1,k,irho) - fs
!            dw(i,j,  k,irho) = dw(i,j,  k,irho) + fs

!            fs = rqsp*w(i,j+1,k,ivx) + rqsm*w(i,j,k,ivx) &
!                 + pa*sJ(i,j,k,1)
!            dw(i,j+1,k,imx) = dw(i,j+1,k,imx) - fs
!            dw(i,j,  k,imx) = dw(i,j,  k,imx) + fs

!            fs = rqsp*w(i,j+1,k,ivy) + rqsm*w(i,j,k,ivy) &
!                 + pa*sJ(i,j,k,2)
!            dw(i,j+1,k,imy) = dw(i,j+1,k,imy) - fs
!            dw(i,j,  k,imy) = dw(i,j,  k,imy) + fs

!            fs = rqsp*w(i,j+1,k,ivz) + rqsm*w(i,j,k,ivz) &
!                 + pa*sJ(i,j,k,3)
!            dw(i,j+1,k,imz) = dw(i,j+1,k,imz) - fs
!            dw(i,j,  k,imz) = dw(i,j,  k,imz) + fs

!            fs = qsp*w(i,j+1,k,irhoE) + qsm*w(i,j,k,irhoE) &
!                 + porFlux*(vnp*p(i,j+1,k) + vnm*p(i,j,k))
!            dw(i,j+1,k,irhoE) = dw(i,j+1,k,irhoE) - fs
!            dw(i,j,  k,irhoE) = dw(i,j,  k,irhoE) + fs
!         enddo
!      enddo
!   enddo

!   sface = zero
!   do k=1,nz
!      do j=2,jl
!         do i=2,il

!            ! Set the dot product of the grid velocity and the
!            ! normal in k-direction for a moving face.

!            !if( addGridVelocities ) sFace = sFaceK(i,j,k)

!            ! Compute the normal velocities of the left and right state.

!            vnp = w(i,j,k+1,ivx)*sK(i,j,k,1) &
!                 + w(i,j,k+1,ivy)*sK(i,j,k,2) &
!                 + w(i,j,k+1,ivz)*sK(i,j,k,3)
!            vnm = w(i,j,k,  ivx)*sK(i,j,k,1) &
!                 + w(i,j,k,  ivy)*sK(i,j,k,2) &
!                 + w(i,j,k,  ivz)*sK(i,j,k,3)

!            ! Set the values of the porosities for this face.
!            ! porVel defines the porosity w.r.t. velocity;
!            ! porFlux defines the porosity w.r.t. the entire flux.
!            ! The latter is only zero for a discontinuous block
!            ! block boundary that must be treated conservatively.
!            ! The default value of porFlux is 0.5, such that the
!            ! correct central flux is scattered to both cells.
!            ! In case of a boundFlux the normal velocity is set
!            ! to sFace.

!            porVel  = one
!            porFlux = half

!            if(porK(i,j,k) == noFlux)    porFlux = zero
!            if(porK(i,j,k) == boundFlux) then
!               porVel = zero
!               vnp    = sFace
!               vnm    = sFace
!            endif

!            ! Incorporate porFlux in porVel.

!            porVel = porVel*porFlux

!            ! Compute the normal velocities for the face as well as the
!            ! mass fluxes.

!            qsp = (vnp - sFace)*porVel
!            qsm = (vnm - sFace)*porVel

!            rqsp = qsp*w(i,j,k+1,irho)
!            rqsm = qsm*w(i,j,k,  irho)

!            ! Compute the sum of the pressure multiplied by porFlux.
!            ! For the default value of porFlux, 0.5, this leads to
!            ! the average pressure.

!            pa = porFlux*(p(i,j,k+1) + p(i,j,k))

!            ! Compute the fluxes and scatter them to the cells
!            ! i,j,k and i,j,k+1. Store the density flux in the
!            ! mass flow of the appropriate sliding mesh interface.

!            fs = rqsp + rqsm
!            dw(i,j,k+1,irho) = dw(i,j,k+1,irho) - fs
!            dw(i,j,k,  irho) = dw(i,j,k,  irho) + fs

!            fs = rqsp*w(i,j,k+1,ivx) + rqsm*w(i,j,k,ivx) &
!                 + pa*sK(i,j,k,1)
!            dw(i,j,k+1,imx) = dw(i,j,k+1,imx) - fs
!            dw(i,j,k,  imx) = dw(i,j,k,  imx) + fs

!            fs = rqsp*w(i,j,k+1,ivy) + rqsm*w(i,j,k,ivy) &
!                 + pa*sK(i,j,k,2)
!            dw(i,j,k+1,imy) = dw(i,j,k+1,imy) - fs
!            dw(i,j,k,  imy) = dw(i,j,k,  imy) + fs

!            fs = rqsp*w(i,j,k+1,ivz) + rqsm*w(i,j,k,ivz) &
!                 + pa*sK(i,j,k,3)
!            dw(i,j,k+1,imz) = dw(i,j,k+1,imz) - fs
!            dw(i,j,k,  imz) = dw(i,j,k,  imz) + fs

!            fs = qsp*w(i,j,k+1,irhoE) + qsm*w(i,j,k,irhoE) &
!                 + porFlux*(vnp*p(i,j,k+1) + vnm*p(i,j,k))
!            dw(i,j,k+1,irhoE) = dw(i,j,k+1,irhoE) - fs
!            dw(i,j,k,  irhoE) = dw(i,j,k,  irhoE) + fs
!         enddo
!      enddo
!   enddo

!   ! ---------------------------------------------
!   !             Inviscid Diss Flux Scalar
!   ! ---------------------------------------------

!   ! Determine the variables used to compute the switch.
!   ! For the inviscid case this is the pressure; for the viscous
!   ! case it is the entropy.

!   select case (equations)
!   case (EulerEquations)

!      ! Inviscid case. Pressure switch is based on the pressure.
!      ! Also set the value of sslim. To be fully consistent this
!      ! must have the dimension of pressure and it is therefore
!      ! set to a fraction of the free stream value.

!      sslim = 0.001_realType*pInfCorr

!      ! Copy the pressure in ss. Only need the entries used in the
!      ! discretization, i.e. not including the corner halo's, but we'll
!      ! just copy all anyway.

!      ss = P
!      !===============================================================

!   case (NSEquations, RANSEquations)

!      ! Viscous case. Pressure switch is based on the entropy.
!      ! Also set the value of sslim. To be fully consistent this
!      ! must have the dimension of entropy and it is therefore
!      ! set to a fraction of the free stream value.

!      sslim = 0.001_realType*pInfCorr/(rhoInf**gammaInf)

!      ! Store the entropy in ss. See above.
!      do k=0, ke
!         do j=0, je
!            do i=0, ie
!               ss(i,j,k) = p(i,j,k)/(w(i,j,k,irho)**gamma(i,j,k))
!            end do
!         end do
!      end do
!   end select

!   ! Compute the pressure sensor for each cell, in each direction:
!   do k=1,kl
!      do j=1,jl
!         do i=1,il
!            dss(i,j,k,1) = abs((ss(i+1,j,k) - two*ss(i,j,k) + ss(i-1,j,k)) &
!                 /     (ss(i+1,j,k) + two*ss(i,j,k) + ss(i-1,j,k) + sslim))

!            dss(i,j,k,2) = abs((ss(i,j+1,k) - two*ss(i,j,k) + ss(i,j-1,k)) &
!                 /     (ss(i,j+1,k) + two*ss(i,j,k) + ss(i,j-1,k) + sslim))

!            dss(i,j,k,3) = abs((ss(i,j,k+1) - two*ss(i,j,k) + ss(i,j,k-1)) &
!                 /     (ss(i,j,k+1) + two*ss(i,j,k) + ss(i,j,k-1) + sslim))
!         end do
!      end do
!   end do

!   ! Set a couple of constants for the scheme.

!   fis2 = rFil*vis2
!   fis4 = rFil*vis4
!   sfil = one - rFil

!   ! Initialize the dissipative residual to a certain times,
!   ! possibly zero, the previously stored value. Owned cells
!   ! only, because the halo values do not matter.

!   fw = sfil*fw
!   !
!   !       Dissipative fluxes in the i-direction.
!   !
!   do k=2,kl
!      do j=2,jl
!         do i=1,nx

!            ! Compute the dissipation coefficients for this face.

!            ppor = zero
!            if(porI(i,j,k) == normalFlux) ppor = half
!            rrad = ppor*(radI(i,j,k) + radI(i+1,j,k))

!            dis2 = fis2*rrad*min(dssMax, max(dss(i,j,k,1), dss(i+1,j,k,1)))
!            dis4 = dim(fis4*rrad, dis2)

!            ! Compute and scatter the dissipative flux.
!            ! Density. Store it in the mass flow of the
!            ! appropriate sliding mesh interface.

!            ddw1 = w(i+1,j,k,irho) - w(i,j,k,irho)
!            fs  = dis2*ddw1 &
!                 - dis4*(w(i+2,j,k,irho) - w(i-1,j,k,irho) - three*ddw1)

!            fw(i+1,j,k,irho) = fw(i+1,j,k,irho) + fs
!            fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

!            ! X-momentum.

!            ddw2 = w(i+1,j,k,ivx)*w(i+1,j,k,irho) - w(i,j,k,ivx)*w(i,j,k,irho)
!            fs  = dis2*ddw2 &
!                 - dis4*(w(i+2,j,k,ivx)*w(i+2,j,k,irho) - w(i-1,j,k,ivx)*w(i-1,j,k,irho) - three*ddw2)

!            fw(i+1,j,k,imx) = fw(i+1,j,k,imx) + fs
!            fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

!            ! Y-momentum.

!            ddw3 = w(i+1,j,k,ivy)*w(i+1,j,k,irho) - w(i,j,k,ivy)*w(i,j,k,irho)
!            fs  = dis2*ddw3 &
!                 - dis4*(w(i+2,j,k,ivy)*w(i+2,j,k,irho) - w(i-1,j,k,ivy)*w(i-1,j,k,irho) - three*ddw3)

!            fw(i+1,j,k,imy) = fw(i+1,j,k,imy) + fs
!            fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

!            ! Z-momentum.

!            ddw4 = w(i+1,j,k,ivz)*w(i+1,j,k,irho) - w(i,j,k,ivz)*w(i,j,k,irho)
!            fs  = dis2*ddw4 &
!                 - dis4*(w(i+2,j,k,ivz)*w(i+2,j,k,irho) - w(i-1,j,k,ivz)*w(i-1,j,k,irho) - three*ddw4)

!            fw(i+1,j,k,imz) = fw(i+1,j,k,imz) + fs
!            fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

!            ! Energy.

!            ddw5 = (w(i+1,j,k,irhoE) + P(i+1,j,K))- (w(i,j,k,irhoE) + P(i,j,k))
!            fs  = dis2*ddw5 &
!                 - dis4*((w(i+2,j,k,irhoE) + P(i+2,j,k)) - (w(i-1,j,k,irhoE) + P(i-1,j,k)) - three*ddw5)

!            fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) + fs
!            fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs
!         end do
!      end do
!   end do
!   !
!   !       Dissipative fluxes in the j-direction.
!   !
!   do k=2,kl
!      do j=1,ny
!         do i=2,il

!            ! Compute the dissipation coefficients for this face.

!            ppor = zero
!            if(porJ(i,j,k) == normalFlux) ppor = half
!            rrad = ppor*(radJ(i,j,k) + radJ(i,j+1,k))

!            dis2 = fis2*rrad*min(dssMax, max(dss(i,j,k,2),dss(i,j+1,k,2)))
!            dis4 = dim(fis4*rrad, dis2)

!            ! Compute and scatter the dissipative flux.
!            ! Density. Store it in the mass flow of the
!            ! appropriate sliding mesh interface.

!            ddw1 = w(i,j+1,k,irho) - w(i,j,k,irho)
!            fs  = dis2*ddw1 &
!                 - dis4*(w(i,j+2,k,irho) - w(i,j-1,k,irho) - three*ddw1)

!            fw(i,j+1,k,irho) = fw(i,j+1,k,irho) + fs
!            fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

!            ! X-momentum.

!            ddw2 = w(i,j+1,k,ivx)*w(i,j+1,k,irho) - w(i,j,k,ivx)*w(i,j,k,irho)
!            fs  = dis2*ddw2 &
!                 - dis4*(w(i,j+2,k,ivx)*w(i,j+2,k,irho) - w(i,j-1,k,ivx)*w(i,j-1,k,irho) - three*ddw2)

!            fw(i,j+1,k,imx) = fw(i,j+1,k,imx) + fs
!            fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

!            ! Y-momentum.

!            ddw3 = w(i,j+1,k,ivy)*w(i,j+1,k,irho) - w(i,j,k,ivy)*w(i,j,k,irho)
!            fs  = dis2*ddw3 &
!                 - dis4*(w(i,j+2,k,ivy)*w(i,j+2,k,irho) - w(i,j-1,k,ivy)*w(i,j-1,k,irho) - three*ddw3)

!            fw(i,j+1,k,imy) = fw(i,j+1,k,imy) + fs
!            fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

!            ! Z-momentum.

!            ddw4 = w(i,j+1,k,ivz)*w(i,j+1,k,irho) - w(i,j,k,ivz)*w(i,j,k,irho)
!            fs  = dis2*ddw4 &
!                 - dis4*(w(i,j+2,k,ivz)*w(i,j+2,k,irho) - w(i,j-1,k,ivz)*w(i,j-1,k,irho) - three*ddw4)

!            fw(i,j+1,k,imz) = fw(i,j+1,k,imz) + fs
!            fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

!            ! Energy.

!            ddw5 = (w(i,j+1,k,irhoE) + P(i,j+1,k)) - (w(i,j,k,irhoE) + P(i,j,k))
!            fs  = dis2*ddw5 &
!                 - dis4*((w(i,j+2,k,irhoE) + P(i,j+2,k)) - (w(i,j-1,k,irhoE) + P(i,j-1,k)) - three*ddw5)

!            fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) + fs
!            fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs
!         end do
!      end do
!   end do
!   !
!   !       Dissipative fluxes in the k-direction.
!   !
!   do k=1,nz
!      do j=2,jl
!         do i=2,il

!            ! Compute the dissipation coefficients for this face.

!            ppor = zero
!            if(porK(i,j,k) == normalFlux) ppor = half
!            rrad = ppor*(radK(i,j,k) + radK(i,j,k+1))

!            dis2 = fis2*rrad*min(dssMax, max(dss(i,j,k,3), dss(i,j,k+1,3)))
!            dis4 = dim(fis4*rrad, dis2)

!            ! Compute and scatter the dissipative flux.
!            ! Density. Store it in the mass flow of the
!            ! appropriate sliding mesh interface.

!            ddw1 = w(i,j,k+1,irho) - w(i,j,k,irho)
!            fs  = dis2*ddw1 &
!                 - dis4*(w(i,j,k+2,irho) - w(i,j,k-1,irho) - three*ddw1)

!            fw(i,j,k+1,irho) = fw(i,j,k+1,irho) + fs
!            fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

!            ! X-momentum.

!            ddw2 = w(i,j,k+1,ivx)*w(i,j,k+1,irho) - w(i,j,k,ivx)*w(i,j,k,irho)
!            fs  = dis2*ddw2 &
!                 - dis4*(w(i,j,k+2,ivx)*w(i,j,k+2,irho) - w(i,j,k-1,ivx)*w(i,j,k-1,irho) - three*ddw2)

!            fw(i,j,k+1,imx) = fw(i,j,k+1,imx) + fs
!            fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

!            ! Y-momentum.

!            ddw3 = w(i,j,k+1,ivy)*w(i,j,k+1,irho) - w(i,j,k,ivy)*w(i,j,k,irho)
!            fs  = dis2*ddw3 &
!                 - dis4*(w(i,j,k+2,ivy)*w(i,j,k+2,irho) - w(i,j,k-1,ivy)*w(i,j,k-1,irho) - three*ddw3)

!            fw(i,j,k+1,imy) = fw(i,j,k+1,imy) + fs
!            fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

!            ! Z-momentum.

!            ddw4 = w(i,j,k+1,ivz)*w(i,j,k+1,irho) - w(i,j,k,ivz)*w(i,j,k,irho)
!            fs  = dis2*ddw4 &
!                 - dis4*(w(i,j,k+2,ivz)*w(i,j,k+2,irho) - w(i,j,k-1,ivz)*w(i,j,k-1,irho) - three*ddw4)

!            fw(i,j,k+1,imz) = fw(i,j,k+1,imz) + fs
!            fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

!            ! Energy.

!            ddw5 = (w(i,j,k+1,irhoE) + P(i,j,k+1)) - (w(i,j,k,irhoE) + P(i,j,k))
!            fs  = dis2*ddw5 &
!                 - dis4*((w(i,j,k+2,irhoE) + P(i,j,k+2)) - (w(i,j,k-1,irhoE) + P(i,j,k-1)) - three*ddw5)

!            fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + fs
!            fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs
!         end do
!      end do
!   end do

!   ! ---------------------------------------------
!   !           Compute the speed of sound squared
!   ! ---------------------------------------------

!   ! Determine if we need to correct for K
!   correctForK = getCorrectForK()

!   if (correctForK) then
!      do k=1,kl
!         do j=1,jl
!            do i=1,il
!               pp = p(i,j,k) - twoThird*w(i,j,k,irho)*w(i,j,k,itu1)
!               aa(i,j,k) = gamma(i,j,k)*pp/w(i,j,k,irho)
!            enddo
!         enddo
!      enddo
!   else
!      do k=1,kl
!         do j=1,jl
!            do i=1,il
!               aa(i,j,k) = gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho)
!            enddo
!         enddo
!      enddo
!   end if

!   ! ---------------------------------------------
!   !           Compute nodal gradients
!   ! ---------------------------------------------

!   ! Zero all nodeal gradients:
!   ux = zero; uy = zero; uz = zero;
!   vx = zero; vy = zero; vz = zero;
!   wx = zero; wy = zero; wz = zero;
!   qx = zero; qy = zero; qz = zero;

!   ! First part. Contribution in the k-direction.
!   ! The contribution is scattered to both the left and right node
!   ! in k-direction.

!   do k=1, kl
!      do j=1, ny
!         do i=1, nx

!            ! Compute 8 times the average normal for this part of
!            ! the control volume. The factor 8 is taken care of later
!            ! on when the division by the volume takes place.

!            sx = sk(i,j,k-1,  1) + sk(i+1,j,k-1,  1) &
!                 + sk(i,j+1,k-1,1) + sk(i+1,j+1,k-1,1) &
!                 + sk(i,j,  k,  1) + sk(i+1,j,  k,  1) &
!                 + sk(i,j+1,k  ,1) + sk(i+1,j+1,k  ,1)
!            sy = sk(i,j,k-1,  2) + sk(i+1,j,k-1,  2) &
!                 + sk(i,j+1,k-1,2) + sk(i+1,j+1,k-1,2) &
!                 + sk(i,j,  k,  2) + sk(i+1,j,  k,  2) &
!                 + sk(i,j+1,k  ,2) + sk(i+1,j+1,k  ,2)
!            sz = sk(i,j,k-1,  3) + sk(i+1,j,k-1,  3) &
!                 + sk(i,j+1,k-1,3) + sk(i+1,j+1,k-1,3) &
!                 + sk(i,j,  k,  3) + sk(i+1,j,  k,  3) &
!                 + sk(i,j+1,k  ,3) + sk(i+1,j+1,k  ,3)

!            ! Compute the average velocities and speed of sound squared
!            ! for this integration point. Node that these variables are
!            ! stored in w(ivx), w(ivy), w(ivz) and p.

!            ubar = fourth*(w(i,j,  k,ivx) + w(i+1,j,  k,ivx) &
!                 +         w(i,j+1,k,ivx) + w(i+1,j+1,k,ivx))
!            vbar = fourth*(w(i,j,  k,ivy) + w(i+1,j,  k,ivy) &
!                 +         w(i,j+1,k,ivy) + w(i+1,j+1,k,ivy))
!            wbar = fourth*(w(i,j,  k,ivz) + w(i+1,j,  k,ivz) &
!                 +         w(i,j+1,k,ivz) + w(i+1,j+1,k,ivz))

!            a2 = fourth*(aa(i,j,k) + aa(i+1,j,k) + aa(i,j+1,k) + aa(i+1,j+1,k))

!            ! Add the contributions to the surface integral to the node
!            ! j-1 and substract it from the node j. For the heat flux it
!            ! is reversed, because the negative of the gradient of the
!            ! speed of sound must be computed.

!            if(k > 1) then
!               ux(i,j,k-1) = ux(i,j,k-1) + ubar*sx
!               uy(i,j,k-1) = uy(i,j,k-1) + ubar*sy
!               uz(i,j,k-1) = uz(i,j,k-1) + ubar*sz

!               vx(i,j,k-1) = vx(i,j,k-1) + vbar*sx
!               vy(i,j,k-1) = vy(i,j,k-1) + vbar*sy
!               vz(i,j,k-1) = vz(i,j,k-1) + vbar*sz

!               wx(i,j,k-1) = wx(i,j,k-1) + wbar*sx
!               wy(i,j,k-1) = wy(i,j,k-1) + wbar*sy
!               wz(i,j,k-1) = wz(i,j,k-1) + wbar*sz

!               qx(i,j,k-1) = qx(i,j,k-1) - a2*sx
!               qy(i,j,k-1) = qy(i,j,k-1) - a2*sy
!               qz(i,j,k-1) = qz(i,j,k-1) - a2*sz
!            endif

!            if(k < kl) then
!               ux(i,j,k) = ux(i,j,k) - ubar*sx
!               uy(i,j,k) = uy(i,j,k) - ubar*sy
!               uz(i,j,k) = uz(i,j,k) - ubar*sz

!               vx(i,j,k) = vx(i,j,k) - vbar*sx
!               vy(i,j,k) = vy(i,j,k) - vbar*sy
!               vz(i,j,k) = vz(i,j,k) - vbar*sz

!               wx(i,j,k) = wx(i,j,k) - wbar*sx
!               wy(i,j,k) = wy(i,j,k) - wbar*sy
!               wz(i,j,k) = wz(i,j,k) - wbar*sz

!               qx(i,j,k) = qx(i,j,k) + a2*sx
!               qy(i,j,k) = qy(i,j,k) + a2*sy
!               qz(i,j,k) = qz(i,j,k) + a2*sz
!            endif
!         end do
!      enddo
!   enddo


!   ! Second part. Contribution in the j-direction.
!   ! The contribution is scattered to both the left and right node
!   ! in j-direction.

!   do k=1, nz
!      do j=1, jl
!         do i=1, nx

!            ! Compute 8 times the average normal for this part of
!            ! the control volume. The factor 8 is taken care of later
!            ! on when the division by the volume takes place.

!            sx = sj(i,j-1,k,  1) + sj(i+1,j-1,k,  1) &
!                 + sj(i,j-1,k+1,1) + sj(i+1,j-1,k+1,1) &
!                 + sj(i,j,  k,  1) + sj(i+1,j,  k,  1) &
!                 + sj(i,j,  k+1,1) + sj(i+1,j,  k+1,1)
!            sy = sj(i,j-1,k,  2) + sj(i+1,j-1,k,  2) &
!                 + sj(i,j-1,k+1,2) + sj(i+1,j-1,k+1,2) &
!                 + sj(i,j,  k,  2) + sj(i+1,j,  k,  2) &
!                 + sj(i,j,  k+1,2) + sj(i+1,j,  k+1,2)
!            sz = sj(i,j-1,k,  3) + sj(i+1,j-1,k,  3) &
!                 + sj(i,j-1,k+1,3) + sj(i+1,j-1,k+1,3) &
!                 + sj(i,j,  k,  3) + sj(i+1,j,  k,  3) &
!                 + sj(i,j,  k+1,3) + sj(i+1,j,  k+1,3)

!            ! Compute the average velocities and speed of sound squared
!            ! for this integration point. Node that these variables are
!            ! stored in w(ivx), w(ivy), w(ivz) and p.

!            ubar = fourth*(w(i,j,k,  ivx) + w(i+1,j,k,  ivx) &
!                 +         w(i,j,k+1,ivx) + w(i+1,j,k+1,ivx))
!            vbar = fourth*(w(i,j,k,  ivy) + w(i+1,j,k,  ivy) &
!                 +         w(i,j,k+1,ivy) + w(i+1,j,k+1,ivy))
!            wbar = fourth*(w(i,j,k,  ivz) + w(i+1,j,k,  ivz) &
!                 +         w(i,j,k+1,ivz) + w(i+1,j,k+1,ivz))

!            a2 = fourth*(aa(i,j,k) + aa(i+1,j,k) + aa(i,j,k+1) + aa(i+1,j,k+1))

!            ! Add the contributions to the surface integral to the node
!            ! j-1 and substract it from the node j. For the heat flux it
!            ! is reversed, because the negative of the gradient of the
!            ! speed of sound must be computed.

!            if(j > 1) then
!               ux(i,j-1,k) = ux(i,j-1,k) + ubar*sx
!               uy(i,j-1,k) = uy(i,j-1,k) + ubar*sy
!               uz(i,j-1,k) = uz(i,j-1,k) + ubar*sz

!               vx(i,j-1,k) = vx(i,j-1,k) + vbar*sx
!               vy(i,j-1,k) = vy(i,j-1,k) + vbar*sy
!               vz(i,j-1,k) = vz(i,j-1,k) + vbar*sz

!               wx(i,j-1,k) = wx(i,j-1,k) + wbar*sx
!               wy(i,j-1,k) = wy(i,j-1,k) + wbar*sy
!               wz(i,j-1,k) = wz(i,j-1,k) + wbar*sz

!               qx(i,j-1,k) = qx(i,j-1,k) - a2*sx
!               qy(i,j-1,k) = qy(i,j-1,k) - a2*sy
!               qz(i,j-1,k) = qz(i,j-1,k) - a2*sz
!            endif

!            if(j < jl) then
!               ux(i,j,k) = ux(i,j,k) - ubar*sx
!               uy(i,j,k) = uy(i,j,k) - ubar*sy
!               uz(i,j,k) = uz(i,j,k) - ubar*sz

!               vx(i,j,k) = vx(i,j,k) - vbar*sx
!               vy(i,j,k) = vy(i,j,k) - vbar*sy
!               vz(i,j,k) = vz(i,j,k) - vbar*sz

!               wx(i,j,k) = wx(i,j,k) - wbar*sx
!               wy(i,j,k) = wy(i,j,k) - wbar*sy
!               wz(i,j,k) = wz(i,j,k) - wbar*sz

!               qx(i,j,k) = qx(i,j,k) + a2*sx
!               qy(i,j,k) = qy(i,j,k) + a2*sy
!               qz(i,j,k) = qz(i,j,k) + a2*sz
!            endif
!         end do
!      enddo
!   enddo
!   !
!   ! Third part. Contribution in the i-direction.
!   ! The contribution is scattered to both the left and right node
!   ! in i-direction.
!   !
!   do k=1,nz
!      do j=1,ny
!         do i=1,il

!            ! Compute 8 times the average normal for this part of
!            ! the control volume. The factor 8 is taken care of later
!            ! on when the division by the volume takes place.

!            sx = si(i-1,j,k,  1) + si(i-1,j+1,k,  1) &
!                 + si(i-1,j,k+1,1) + si(i-1,j+1,k+1,1) &
!                 + si(i,  j,k,  1) + si(i,  j+1,k,  1) &
!                 + si(i,  j,k+1,1) + si(i,  j+1,k+1,1)
!            sy = si(i-1,j,k,  2) + si(i-1,j+1,k,  2) &
!                 + si(i-1,j,k+1,2) + si(i-1,j+1,k+1,2) &
!                 + si(i,  j,k,  2) + si(i,  j+1,k,  2) &
!                 + si(i,  j,k+1,2) + si(i,  j+1,k+1,2)
!            sz = si(i-1,j,k,  3) + si(i-1,j+1,k,  3) &
!                 + si(i-1,j,k+1,3) + si(i-1,j+1,k+1,3) &
!                 + si(i,  j,k,  3) + si(i,  j+1,k,  3) &
!                 + si(i,  j,k+1,3) + si(i,  j+1,k+1,3)

!            ! Compute the average velocities and speed of sound squared
!            ! for this integration point. Node that these variables are
!            ! stored in w(ivx), w(ivy), w(ivz) and p.

!            ubar = fourth*(w(i,j,k,  ivx) + w(i,j+1,k,  ivx) &
!                 +         w(i,j,k+1,ivx) + w(i,j+1,k+1,ivx))
!            vbar = fourth*(w(i,j,k,  ivy) + w(i,j+1,k,  ivy) &
!                 +         w(i,j,k+1,ivy) + w(i,j+1,k+1,ivy))
!            wbar = fourth*(w(i,j,k,  ivz) + w(i,j+1,k,  ivz) &
!                 +         w(i,j,k+1,ivz) + w(i,j+1,k+1,ivz))

!            a2 = fourth*(aa(i,j,k) + aa(i,j+1,k) + aa(i,j,k+1) + aa(i,j+1,k+1))

!            ! Add the contributions to the surface integral to the node
!            ! j-1 and substract it from the node j. For the heat flux it
!            ! is reversed, because the negative of the gradient of the
!            ! speed of sound must be computed.

!            if(i > 1) then
!               ux(i-1,j,k) = ux(i-1,j,k) + ubar*sx
!               uy(i-1,j,k) = uy(i-1,j,k) + ubar*sy
!               uz(i-1,j,k) = uz(i-1,j,k) + ubar*sz

!               vx(i-1,j,k) = vx(i-1,j,k) + vbar*sx
!               vy(i-1,j,k) = vy(i-1,j,k) + vbar*sy
!               vz(i-1,j,k) = vz(i-1,j,k) + vbar*sz

!               wx(i-1,j,k) = wx(i-1,j,k) + wbar*sx
!               wy(i-1,j,k) = wy(i-1,j,k) + wbar*sy
!               wz(i-1,j,k) = wz(i-1,j,k) + wbar*sz

!               qx(i-1,j,k) = qx(i-1,j,k) - a2*sx
!               qy(i-1,j,k) = qy(i-1,j,k) - a2*sy
!               qz(i-1,j,k) = qz(i-1,j,k) - a2*sz
!            endif

!            if(i < il) then
!               ux(i,j,k) = ux(i,j,k) - ubar*sx
!               uy(i,j,k) = uy(i,j,k) - ubar*sy
!               uz(i,j,k) = uz(i,j,k) - ubar*sz

!               vx(i,j,k) = vx(i,j,k) - vbar*sx
!               vy(i,j,k) = vy(i,j,k) - vbar*sy
!               vz(i,j,k) = vz(i,j,k) - vbar*sz

!               wx(i,j,k) = wx(i,j,k) - wbar*sx
!               wy(i,j,k) = wy(i,j,k) - wbar*sy
!               wz(i,j,k) = wz(i,j,k) - wbar*sz

!               qx(i,j,k) = qx(i,j,k) + a2*sx
!               qy(i,j,k) = qy(i,j,k) + a2*sy
!               qz(i,j,k) = qz(i,j,k) + a2*sz
!            endif
!         enddo
!      enddo
!   enddo

!   ! Divide by 8 times the volume to obtain the correct gradients.

!   do k=1,nz
!      do j=1,ny
!         do i=1,nx

!            ! Compute the inverse of 8 times the volume for this node.

!            oVol = one/(vol(i,  j,  k) + vol(i,  j,  k+1) &
!                 +      vol(i+1,j,  k) + vol(i+1,j,  k+1) &
!                 +      vol(i,  j+1,k) + vol(i,  j+1,k+1) &
!                 +      vol(i+1,j+1,k) + vol(i+1,j+1,k+1))

!            ! Compute the correct velocity gradients and "unit" heat
!            ! fluxes. The velocity gradients are stored in ux, etc.

!            ux(i,j,k) = ux(i,j,k)*oVol
!            uy(i,j,k) = uy(i,j,k)*oVol
!            uz(i,j,k) = uz(i,j,k)*oVol

!            vx(i,j,k) = vx(i,j,k)*oVol
!            vy(i,j,k) = vy(i,j,k)*oVol
!            vz(i,j,k) = vz(i,j,k)*oVol

!            wx(i,j,k) = wx(i,j,k)*oVol
!            wy(i,j,k) = wy(i,j,k)*oVol
!            wz(i,j,k) = wz(i,j,k)*oVol

!            qx(i,j,k) = qx(i,j,k)*oVol
!            qy(i,j,k) = qy(i,j,k)*oVol
!            qz(i,j,k) = qz(i,j,k)*oVol
!         end do
!      enddo
!   enddo

!   ! ---------------------------------------------
!   !                   Viscous Flux
!   ! ---------------------------------------------

!   ! Set QCR parameters
!   Ccr1 = 0.3_realType

!   ! The diagonals of the vorticity tensor components are always zero
!   Wxx = zero
!   Wyy = zero
!   Wzz = zero
!   !
!   !         viscous fluxes in the k-direction.
!   !
!   mue = zero
!   do k=1,nz
!      do j=2,jl
!         do i=2,il

!            ! Set the value of the porosity. If not zero, it is set
!            ! to average the eddy-viscosity and to take the factor
!            ! rFilv into account.

!            por = half*rFilv
!            if(porK(i,j,k) == noFlux) por = zero

!            ! Compute the laminar and (if present) the eddy viscosities
!            ! multiplied by the porosity. Compute the factor in front of
!            ! the gradients of the speed of sound squared for the heat
!            ! flux.

!            mul = por*(rlv(i,j,k) + rlv(i,j,k+1))
!            if( eddyModel ) mue = por*(rev(i,j,k) + rev(i,j,k+1))
!            mut = mul + mue

!            gm1          = half*(gamma(i,j,k) + gamma(i,j,k+1)) - one
!            factLamHeat  = one/(prandtl*gm1)
!            factTurbHeat = one/(prandtlTurb*gm1)

!            heatCoef = mul*factLamHeat + mue*factTurbHeat

!            ! Compute the gradients at the face by averaging the four
!            ! nodal values.

!            u_x = fourth*(ux(i-1,j-1,k) + ux(i,j-1,k) &
!                 +         ux(i-1,j,  k) + ux(i,j,  k))
!            u_y = fourth*(uy(i-1,j-1,k) + uy(i,j-1,k) &
!                 +         uy(i-1,j,  k) + uy(i,j,  k))
!            u_z = fourth*(uz(i-1,j-1,k) + uz(i,j-1,k) &
!                 +         uz(i-1,j,  k) + uz(i,j,  k))

!            v_x = fourth*(vx(i-1,j-1,k) + vx(i,j-1,k) &
!                 +         vx(i-1,j,  k) + vx(i,j,  k))
!            v_y = fourth*(vy(i-1,j-1,k) + vy(i,j-1,k) &
!                 +         vy(i-1,j,  k) + vy(i,j,  k))
!            v_z = fourth*(vz(i-1,j-1,k) + vz(i,j-1,k) &
!                 +         vz(i-1,j,  k) + vz(i,j,  k))

!            w_x = fourth*(wx(i-1,j-1,k) + wx(i,j-1,k) &
!                 +         wx(i-1,j,  k) + wx(i,j,  k))
!            w_y = fourth*(wy(i-1,j-1,k) + wy(i,j-1,k) &
!                 +         wy(i-1,j,  k) + wy(i,j,  k))
!            w_z = fourth*(wz(i-1,j-1,k) + wz(i,j-1,k) &
!                 +         wz(i-1,j,  k) + wz(i,j,  k))

!            q_x = fourth*(qx(i-1,j-1,k) + qx(i,j-1,k) &
!                 +         qx(i-1,j,  k) + qx(i,j,  k))
!            q_y = fourth*(qy(i-1,j-1,k) + qy(i,j-1,k) &
!                 +         qy(i-1,j,  k) + qy(i,j,  k))
!            q_z = fourth*(qz(i-1,j-1,k) + qz(i,j-1,k) &
!                 +         qz(i-1,j,  k) + qz(i,j,  k))


!            ! The gradients in the normal direction are corrected, such
!            ! that no averaging takes places here.
!            ! First determine the vector in the direction from the
!            ! cell center k to cell center k+1.

!            ssx = eighth*(x(i-1,j-1,k+1,1) - x(i-1,j-1,k-1,1) &
!                 +         x(i-1,j,  k+1,1) - x(i-1,j,  k-1,1) &
!                 +         x(i,  j-1,k+1,1) - x(i,  j-1,k-1,1) &
!                 +         x(i,  j,  k+1,1) - x(i,  j,  k-1,1))
!            ssy = eighth*(x(i-1,j-1,k+1,2) - x(i-1,j-1,k-1,2) &
!                 +         x(i-1,j,  k+1,2) - x(i-1,j,  k-1,2) &
!                 +         x(i,  j-1,k+1,2) - x(i,  j-1,k-1,2) &
!                 +         x(i,  j,  k+1,2) - x(i,  j,  k-1,2))
!            ssz = eighth*(x(i-1,j-1,k+1,3) - x(i-1,j-1,k-1,3) &
!                 +         x(i-1,j,  k+1,3) - x(i-1,j,  k-1,3) &
!                 +         x(i,  j-1,k+1,3) - x(i,  j-1,k-1,3) &
!                 +         x(i,  j,  k+1,3) - x(i,  j,  k-1,3))

!            ! Determine the length of this vector and create the
!            ! unit normal.

!            snrm  = one/sqrt(ssx*ssx + ssy*ssy + ssz*ssz)
!            ssx = snrm*ssx
!            ssy = snrm*ssy
!            ssz = snrm*ssz

!            ! Correct the gradients.

!            corr = u_x*ssx + u_y*ssy + u_z*ssz        &
!                 - (w(i,j,k+1,ivx) - w(i,j,k,ivx))*snrm
!            u_x  = u_x - corr*ssx
!            u_y  = u_y - corr*ssy
!            u_z  = u_z - corr*ssz

!            corr = v_x*ssx + v_y*ssy + v_z*ssz        &
!                 - (w(i,j,k+1,ivy) - w(i,j,k,ivy))*snrm
!            v_x  = v_x - corr*ssx
!            v_y  = v_y - corr*ssy
!            v_z  = v_z - corr*ssz

!            corr = w_x*ssx + w_y*ssy + w_z*ssz        &
!                 - (w(i,j,k+1,ivz) - w(i,j,k,ivz))*snrm
!            w_x  = w_x - corr*ssx
!            w_y  = w_y - corr*ssy
!            w_z  = w_z - corr*ssz

!            corr = q_x*ssx + q_y*ssy + q_z*ssz &
!                 + (aa(i,j,k+1) - aa(i,j,k))*snrm
!            q_x  = q_x - corr*ssx
!            q_y  = q_y - corr*ssy
!            q_z  = q_z - corr*ssz

!            ! Compute the stress tensor and the heat flux vector.

!            fracDiv = twoThird*(u_x + v_y + w_z)

!            tauxx = mut*(two*u_x - fracDiv)
!            tauyy = mut*(two*v_y - fracDiv)
!            tauzz = mut*(two*w_z - fracDiv)

!            tauxy = mut*(u_y + v_x)
!            tauxz = mut*(u_z + w_x)
!            tauyz = mut*(v_z + w_y)

!            q_x = heatCoef*q_x
!            q_y = heatCoef*q_y
!            q_z = heatCoef*q_z

!            ! Add QCR corrections if necessary
!            if (useQCR) then

!               ! In the QCR formulation, we add an extra term to the shear tensor:
!               !
!               ! tau_ij,QCR = tau_ij - e_ij
!               !
!               ! where, according to TMR website (http://turbmodels.larc.nasa.gov/spalart.html):
!               !
!               ! e_ij = Ccr1*(O_ik*tau_jk + O_jk*tau_ik)
!               !
!               ! We are computing O_ik as follows:
!               !
!               ! O_ik = 2*W_ik/den

!               ! Compute denominator
!               den = sqrt(u_x*u_x + u_y*u_y + u_z*u_z + &
!                    v_x*v_x + v_y*v_y + v_z*v_z + &
!                    w_x*w_x + w_y*w_y + w_z*w_z)

!               ! Denominator should be limited to avoid division by zero in regions with
!               ! no gradients
!               den = max(den, xminn)

!               ! Compute factor that will multiply all tensor components
!               fact = Ccr1/den

!               ! Compute off-diagonal terms of vorticity tensor (we will ommit the 1/2)
!               Wxy = u_y - v_x
!               Wxz = u_z - w_x
!               Wyz = u_y - v_x
!               Wyx = -Wxy
!               Wzx = -Wxz
!               Wzy = -Wyz

!               ! Compute the extra terms of the Boussinesq relation
!               exx = fact*(Wxx*tauxx + Wxy*tauxy + Wxz*tauxz)*two
!               eyy = fact*(Wyx*tauxy + Wyy*tauyy + Wyz*tauyz)*two
!               ezz = fact*(Wzx*tauxz + Wzy*tauyz + Wzz*tauzz)*two

!               exy = fact*(Wxx*tauxy + Wxy*tauyy + Wxz*tauyz + &
!                    Wyx*tauxx + Wyy*tauxy + Wyz*tauxz)
!               exz = fact*(Wxx*tauxz + Wxy*tauyz + Wxz*tauzz + &
!                    Wzx*tauxx + Wzy*tauxy + Wzz*tauxz)
!               eyz = fact*(Wyx*tauxz + Wyy*tauyz + Wyz*tauzz + &
!                    Wzx*tauxy + Wzy*tauyy + Wzz*tauyz)

!               ! Add extra terms
!               tauxx = tauxx - exx
!               tauyy = tauyy - eyy
!               tauzz = tauzz - ezz
!               tauxy = tauxy - exy
!               tauxz = tauxz - exz
!               tauyz = tauyz - eyz

!            end if

!            ! Compute the average velocities for the face. Remember that
!            ! the velocities are stored and not the momentum.

!            ubar = half*(w(i,j,k,ivx) + w(i,j,k+1,ivx))
!            vbar = half*(w(i,j,k,ivy) + w(i,j,k+1,ivy))
!            wbar = half*(w(i,j,k,ivz) + w(i,j,k+1,ivz))

!            ! Compute the viscous fluxes for this k-face.

!            fmx   = tauxx*sk(i,j,k,1) + tauxy*sk(i,j,k,2) &
!                 + tauxz*sk(i,j,k,3)
!            fmy   = tauxy*sk(i,j,k,1) + tauyy*sk(i,j,k,2) &
!                 + tauyz*sk(i,j,k,3)
!            fmz   = tauxz*sk(i,j,k,1) + tauyz*sk(i,j,k,2) &
!                 + tauzz*sk(i,j,k,3)
!            frhoE =         (ubar*tauxx + vbar*tauxy + wbar*tauxz)*sk(i,j,k,1)
!            frhoE = frhoE + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*sk(i,j,k,2)
!            frhoE = frhoE + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*sk(i,j,k,3)
!            frhoE = frhoE -  q_x*sk(i,j,k,1) - q_y*sk(i,j,k,2) - q_z*sk(i,j,k,3)

!            ! Update the residuals of cell k and k+1.

!            fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
!            fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
!            fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
!            fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE

!            fw(i,j,k+1,imx)   = fw(i,j,k+1,imx)   + fmx
!            fw(i,j,k+1,imy)   = fw(i,j,k+1,imy)   + fmy
!            fw(i,j,k+1,imz)   = fw(i,j,k+1,imz)   + fmz
!            fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + frhoE

!            ! Store the stress tensor and the heat flux vector if this
!            ! face is part of a viscous subface. Both the cases k == 1
!            ! and k == kl must be tested.
!         end do
!      enddo
!   end do
!   !
!   !         Viscous fluxes in the j-direction.
!   !
!   do k=2,kl
!      do j=1,ny
!         do i=2,il

!            ! Set the value of the porosity. If not zero, it is set
!            ! to average the eddy-viscosity and to take the factor
!            ! rFilv into account.

!            por = half*rFilv
!            if(porJ(i,j,k) == noFlux) por = zero

!            ! Compute the laminar and (if present) the eddy viscosities
!            ! multiplied by the porosity. Compute the factor in front of
!            ! the gradients of the speed of sound squared for the heat
!            ! flux.

!            mul = por*(rlv(i,j,k) + rlv(i,j+1,k))
!            if( eddyModel ) mue = por*(rev(i,j,k) + rev(i,j+1,k))
!            mut = mul + mue

!            gm1          = half*(gamma(i,j,k) + gamma(i,j+1,k)) - one
!            factLamHeat  = one/(prandtl*gm1)
!            factTurbHeat = one/(prandtlTurb*gm1)

!            heatCoef = mul*factLamHeat + mue*factTurbHeat

!            ! Compute the gradients at the face by averaging the four
!            ! nodal values.

!            u_x = fourth*(ux(i-1,j,k-1) + ux(i,j,k-1) &
!                 +        ux(i-1,j,k  ) + ux(i,j,k  ))
!            u_y = fourth*(uy(i-1,j,k-1) + uy(i,j,k-1) &
!                 +        uy(i-1,j,k  ) + uy(i,j,k  ))
!            u_z = fourth*(uz(i-1,j,k-1) + uz(i,j,k-1) &
!                 +        uz(i-1,j,k  ) + uz(i,j,k  ))

!            v_x = fourth*(vx(i-1,j,k-1) + vx(i,j,k-1) &
!                 +        vx(i-1,j,k  ) + vx(i,j,k  ))
!            v_y = fourth*(vy(i-1,j,k-1) + vy(i,j,k-1) &
!                 +        vy(i-1,j,k  ) + vy(i,j,k  ))
!            v_z = fourth*(vz(i-1,j,k-1) + vz(i,j,k-1) &
!                 +        vz(i-1,j,k  ) + vz(i,j,k  ))

!            w_x = fourth*(wx(i-1,j,k-1) + wx(i,j,k-1) &
!                 +        wx(i-1,j,k  ) + wx(i,j,k  ))
!            w_y = fourth*(wy(i-1,j,k-1) + wy(i,j,k-1) &
!                 +        wy(i-1,j,k  ) + wy(i,j,k  ))
!            w_z = fourth*(wz(i-1,j,k-1) + wz(i,j,k-1) &
!                 +        wz(i-1,j,k  ) + wz(i,j,k  ))

!            q_x = fourth*(qx(i-1,j,k-1) + qx(i,j,k-1) &
!                 +        qx(i-1,j,k  ) + qx(i,j,k  ))
!            q_y = fourth*(qy(i-1,j,k-1) + qy(i,j,k-1) &
!                 +        qy(i-1,j,k  ) + qy(i,j,k  ))
!            q_z = fourth*(qz(i-1,j,k-1) + qz(i,j,k-1) &
!                 +        qz(i-1,j,k  ) + qz(i,j,k  ))

!            ! The gradients in the normal direction are corrected, such
!            ! that no averaging takes places here.
!            ! First determine the vector in the direction from the
!            ! cell center j to cell center j+1.

!            ssx = eighth*(x(i-1,j+1,k-1,1) - x(i-1,j-1,k-1,1) &
!                 +         x(i-1,j+1,k,  1) - x(i-1,j-1,k,  1) &
!                 +         x(i,  j+1,k-1,1) - x(i,  j-1,k-1,1) &
!                 +         x(i,  j+1,k,  1) - x(i,  j-1,k,  1))
!            ssy = eighth*(x(i-1,j+1,k-1,2) - x(i-1,j-1,k-1,2) &
!                 +         x(i-1,j+1,k,  2) - x(i-1,j-1,k,  2) &
!                 +         x(i,  j+1,k-1,2) - x(i,  j-1,k-1,2) &
!                 +         x(i,  j+1,k,  2) - x(i,  j-1,k,  2))
!            ssz = eighth*(x(i-1,j+1,k-1,3) - x(i-1,j-1,k-1,3) &
!                 +         x(i-1,j+1,k,  3) - x(i-1,j-1,k,  3) &
!                 +         x(i,  j+1,k-1,3) - x(i,  j-1,k-1,3) &
!                 +         x(i,  j+1,k,  3) - x(i,  j-1,k,  3))

!            ! Determine the length of this vector and create the
!            ! unit normal.

!            snrm  = one/sqrt(ssx*ssx + ssy*ssy + ssz*ssz)
!            ssx = snrm*ssx
!            ssy = snrm*ssy
!            ssz = snrm*ssz

!            ! Correct the gradients.

!            corr = u_x*ssx + u_y*ssy + u_z*ssz        &
!                 - (w(i,j+1,k,ivx) - w(i,j,k,ivx))*snrm
!            u_x  = u_x - corr*ssx
!            u_y  = u_y - corr*ssy
!            u_z  = u_z - corr*ssz

!            corr = v_x*ssx + v_y*ssy + v_z*ssz        &
!                 - (w(i,j+1,k,ivy) - w(i,j,k,ivy))*snrm
!            v_x  = v_x - corr*ssx
!            v_y  = v_y - corr*ssy
!            v_z  = v_z - corr*ssz

!            corr = w_x*ssx + w_y*ssy + w_z*ssz        &
!                 - (w(i,j+1,k,ivz) - w(i,j,k,ivz))*snrm
!            w_x  = w_x - corr*ssx
!            w_y  = w_y - corr*ssy
!            w_z  = w_z - corr*ssz

!            corr = q_x*ssx + q_y*ssy + q_z*ssz &
!                 + (aa(i,j+1,k) - aa(i,j,k))*snrm
!            q_x  = q_x - corr*ssx
!            q_y  = q_y - corr*ssy
!            q_z  = q_z - corr*ssz

!            ! Compute the stress tensor and the heat flux vector.

!            fracDiv = twoThird*(u_x + v_y + w_z)

!            tauxx = mut*(two*u_x - fracDiv)
!            tauyy = mut*(two*v_y - fracDiv)
!            tauzz = mut*(two*w_z - fracDiv)

!            tauxy = mut*(u_y + v_x)
!            tauxz = mut*(u_z + w_x)
!            tauyz = mut*(v_z + w_y)

!            q_x = heatCoef*q_x
!            q_y = heatCoef*q_y
!            q_z = heatCoef*q_z

!            ! Add QCR corrections if necessary
!            if (useQCR) then

!               ! In the QCR formulation, we add an extra term to the shear tensor:
!               !
!               ! tau_ij,QCR = tau_ij - e_ij
!               !
!               ! where, according to TMR website (http://turbmodels.larc.nasa.gov/spalart.html):
!               !
!               ! e_ij = Ccr1*(O_ik*tau_jk + O_jk*tau_ik)
!               !
!               ! We are computing O_ik as follows:
!               !
!               ! O_ik = 2*W_ik/den

!               ! Compute denominator
!               den = sqrt(u_x*u_x + u_y*u_y + u_z*u_z + &
!                    v_x*v_x + v_y*v_y + v_z*v_z + &
!                    w_x*w_x + w_y*w_y + w_z*w_z)

!               ! Denominator should be limited to avoid division by zero in regions with
!               ! no gradients
!               den = max(den, xminn)

!               ! Compute factor that will multiply all tensor components
!               fact = Ccr1/den

!               ! Compute off-diagonal terms of vorticity tensor (we will ommit the 1/2)
!               Wxy = u_y - v_x
!               Wxz = u_z - w_x
!               Wyz = u_y - v_x
!               Wyx = -Wxy
!               Wzx = -Wxz
!               Wzy = -Wyz

!               ! Compute the extra terms of the Boussinesq relation
!               exx = fact*(Wxx*tauxx + Wxy*tauxy + Wxz*tauxz)*two
!               eyy = fact*(Wyx*tauxy + Wyy*tauyy + Wyz*tauyz)*two
!               ezz = fact*(Wzx*tauxz + Wzy*tauyz + Wzz*tauzz)*two

!               exy = fact*(Wxx*tauxy + Wxy*tauyy + Wxz*tauyz + &
!                    Wyx*tauxx + Wyy*tauxy + Wyz*tauxz)
!               exz = fact*(Wxx*tauxz + Wxy*tauyz + Wxz*tauzz + &
!                    Wzx*tauxx + Wzy*tauxy + Wzz*tauxz)
!               eyz = fact*(Wyx*tauxz + Wyy*tauyz + Wyz*tauzz + &
!                    Wzx*tauxy + Wzy*tauyy + Wzz*tauyz)

!               ! Add extra terms
!               tauxx = tauxx - exx
!               tauyy = tauyy - eyy
!               tauzz = tauzz - ezz
!               tauxy = tauxy - exy
!               tauxz = tauxz - exz
!               tauyz = tauyz - eyz

!            end if

!            ! Compute the average velocities for the face. Remember that
!            ! the velocities are stored and not the momentum.

!            ubar = half*(w(i,j,k,ivx) + w(i,j+1,k,ivx))
!            vbar = half*(w(i,j,k,ivy) + w(i,j+1,k,ivy))
!            wbar = half*(w(i,j,k,ivz) + w(i,j+1,k,ivz))

!            ! Compute the viscous fluxes for this j-face.

!            fmx   = tauxx*sj(i,j,k,1) + tauxy*sj(i,j,k,2) &
!                 + tauxz*sj(i,j,k,3)
!            fmy   = tauxy*sj(i,j,k,1) + tauyy*sj(i,j,k,2) &
!                 + tauyz*sj(i,j,k,3)
!            fmz   = tauxz*sj(i,j,k,1) + tauyz*sj(i,j,k,2) &
!                 + tauzz*sj(i,j,k,3)
!            frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*sj(i,j,k,1) &
!                 + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*sj(i,j,k,2) &
!                 + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*sj(i,j,k,3) &
!                 - q_x*sj(i,j,k,1) - q_y*sj(i,j,k,2) - q_z*sj(i,j,k,3)

!            ! Update the residuals of cell j and j+1.

!            fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
!            fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
!            fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
!            fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE

!            fw(i,j+1,k,imx)   = fw(i,j+1,k,imx)   + fmx
!            fw(i,j+1,k,imy)   = fw(i,j+1,k,imy)   + fmy
!            fw(i,j+1,k,imz)   = fw(i,j+1,k,imz)   + fmz
!            fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) + frhoE
!         enddo
!      enddo
!   enddo
!   !
!   !         Viscous fluxes in the i-direction.
!   !
!   do k=2, kl
!      do j=2, jl
!         do i=1, nx

!            ! Set the value of the porosity. If not zero, it is set
!            ! to average the eddy-viscosity and to take the factor
!            ! rFilv into account.

!            por = half*rFilv
!            if(porI(i,j,k) == noFlux) por = zero

!            ! Compute the laminar and (if present) the eddy viscosities
!            ! multiplied the porosity. Compute the factor in front of
!            ! the gradients of the speed of sound squared for the heat
!            ! flux.

!            mul = por*(rlv(i,j,k) + rlv(i+1,j,k))
!            if( eddyModel ) mue = por*(rev(i,j,k) + rev(i+1,j,k))
!            mut = mul + mue

!            gm1          = half*(gamma(i,j,k) + gamma(i+1,j,k)) - one
!            factLamHeat  = one/(prandtl*gm1)
!            factTurbHeat = one/(prandtlTurb*gm1)

!            heatCoef = mul*factLamHeat + mue*factTurbHeat

!            ! Compute the gradients at the face by averaging the four
!            ! nodal values.

!            u_x = fourth*(ux(i,j-1,k-1) + ux(i,j,k-1) &
!                 +        ux(i,j-1,k  ) + ux(i,j,k  ))
!            u_y = fourth*(uy(i,j-1,k-1) + uy(i,j,k-1) &
!                 +        uy(i,j-1,k  ) + uy(i,j,k  ))
!            u_z = fourth*(uz(i,j-1,k-1) + uz(i,j,k-1) &
!                 +        uz(i,j-1,k  ) + uz(i,j,k  ))

!            v_x = fourth*(vx(i,j-1,k-1) + vx(i,j,k-1) &
!                 +        vx(i,j-1,k  ) + vx(i,j,k  ))
!            v_y = fourth*(vy(i,j-1,k-1) + vy(i,j,k-1) &
!                 +        vy(i,j-1,k  ) + vy(i,j,k  ))
!            v_z = fourth*(vz(i,j-1,k-1) + vz(i,j,k-1) &
!                 +        vz(i,j-1,k  ) + vz(i,j,k  ))

!            w_x = fourth*(wx(i,j-1,k-1) + wx(i,j,k-1) &
!                 +        wx(i,j-1,k  ) + wx(i,j,k  ))
!            w_y = fourth*(wy(i,j-1,k-1) + wy(i,j,k-1) &
!                 +        wy(i,j-1,k  ) + wy(i,j,k  ))
!            w_z = fourth*(wz(i,j-1,k-1) + wz(i,j,k-1) &
!                 +        wz(i,j-1,k  ) + wz(i,j,k  ))

!            q_x = fourth*(qx(i,j-1,k-1) + qx(i,j,k-1) &
!                 +        qx(i,j-1,k  ) + qx(i,j,k  ))
!            q_y = fourth*(qy(i,j-1,k-1) + qy(i,j,k-1) &
!                 +        qy(i,j-1,k  ) + qy(i,j,k  ))
!            q_z = fourth*(qz(i,j-1,k-1) + qz(i,j,k-1) &
!                 +        qz(i,j-1,k  ) + qz(i,j,k  ))

!            ! The gradients in the normal direction are corrected, such
!            ! that no averaging takes places here.
!            ! First determine the vector in the direction from the
!            ! cell center i to cell center i+1.

!            ssx = eighth*(x(i+1,j-1,k-1,1) - x(i-1,j-1,k-1,1) &
!                 +         x(i+1,j-1,k,  1) - x(i-1,j-1,k,  1) &
!                 +         x(i+1,j,  k-1,1) - x(i-1,j,  k-1,1) &
!                 +         x(i+1,j,  k,  1) - x(i-1,j,  k,  1))
!            ssy = eighth*(x(i+1,j-1,k-1,2) - x(i-1,j-1,k-1,2) &
!                 +         x(i+1,j-1,k,  2) - x(i-1,j-1,k,  2) &
!                 +         x(i+1,j,  k-1,2) - x(i-1,j,  k-1,2) &
!                 +         x(i+1,j,  k,  2) - x(i-1,j,  k,  2))
!            ssz = eighth*(x(i+1,j-1,k-1,3) - x(i-1,j-1,k-1,3) &
!                 +         x(i+1,j-1,k,  3) - x(i-1,j-1,k,  3) &
!                 +         x(i+1,j,  k-1,3) - x(i-1,j,  k-1,3) &
!                 +         x(i+1,j,  k,  3) - x(i-1,j,  k,  3))

!            ! Determine the length of this vector and create the
!            ! unit normal.

!            snrm  = one/sqrt(ssx*ssx + ssy*ssy + ssz*ssz)
!            ssx = snrm*ssx
!            ssy = snrm*ssy
!            ssz = snrm*ssz

!            ! Correct the gradients.

!            corr = u_x*ssx + u_y*ssy + u_z*ssz        &
!                 - (w(i+1,j,k,ivx) - w(i,j,k,ivx))*snrm
!            u_x  = u_x - corr*ssx
!            u_y  = u_y - corr*ssy
!            u_z  = u_z - corr*ssz

!            corr = v_x*ssx + v_y*ssy + v_z*ssz        &
!                 - (w(i+1,j,k,ivy) - w(i,j,k,ivy))*snrm
!            v_x  = v_x - corr*ssx
!            v_y  = v_y - corr*ssy
!            v_z  = v_z - corr*ssz

!            corr = w_x*ssx + w_y*ssy + w_z*ssz        &
!                 - (w(i+1,j,k,ivz) - w(i,j,k,ivz))*snrm
!            w_x  = w_x - corr*ssx
!            w_y  = w_y - corr*ssy
!            w_z  = w_z - corr*ssz

!            corr = q_x*ssx + q_y*ssy + q_z*ssz &
!                 + (aa(i+1,j,k) - aa(i,j,k))*snrm
!            q_x  = q_x - corr*ssx
!            q_y  = q_y - corr*ssy
!            q_z  = q_z - corr*ssz

!            ! Compute the stress tensor and the heat flux vector.

!            fracDiv = twoThird*(u_x + v_y + w_z)

!            tauxx = mut*(two*u_x - fracDiv)
!            tauyy = mut*(two*v_y - fracDiv)
!            tauzz = mut*(two*w_z - fracDiv)

!            tauxy = mut*(u_y + v_x)
!            tauxz = mut*(u_z + w_x)
!            tauyz = mut*(v_z + w_y)

!            q_x = heatCoef*q_x
!            q_y = heatCoef*q_y
!            q_z = heatCoef*q_z

!            ! Add QCR corrections if necessary
!            if (useQCR) then

!               ! In the QCR formulation, we add an extra term to the shear tensor:
!               !
!               ! tau_ij,QCR = tau_ij - e_ij
!               !
!               ! where, according to TMR website (http://turbmodels.larc.nasa.gov/spalart.html):
!               !
!               ! e_ij = Ccr1*(O_ik*tau_jk + O_jk*tau_ik)
!               !
!               ! We are computing O_ik as follows:
!               !
!               ! O_ik = 2*W_ik/den

!               ! Compute denominator
!               den = sqrt(u_x*u_x + u_y*u_y + u_z*u_z + &
!                    v_x*v_x + v_y*v_y + v_z*v_z + &
!                    w_x*w_x + w_y*w_y + w_z*w_z)

!               ! Denominator should be limited to avoid division by zero in regions with
!               ! no gradients
!               den = max(den, xminn)

!               ! Compute factor that will multiply all tensor components
!               fact = Ccr1/den

!               ! Compute off-diagonal terms of vorticity tensor (we will ommit the 1/2)
!               Wxy = u_y - v_x
!               Wxz = u_z - w_x
!               Wyz = u_y - v_x
!               Wyx = -Wxy
!               Wzx = -Wxz
!               Wzy = -Wyz

!               ! Compute the extra terms of the Boussinesq relation
!               exx = fact*(Wxx*tauxx + Wxy*tauxy + Wxz*tauxz)*two
!               eyy = fact*(Wyx*tauxy + Wyy*tauyy + Wyz*tauyz)*two
!               ezz = fact*(Wzx*tauxz + Wzy*tauyz + Wzz*tauzz)*two

!               exy = fact*(Wxx*tauxy + Wxy*tauyy + Wxz*tauyz + &
!                    Wyx*tauxx + Wyy*tauxy + Wyz*tauxz)
!               exz = fact*(Wxx*tauxz + Wxy*tauyz + Wxz*tauzz + &
!                    Wzx*tauxx + Wzy*tauxy + Wzz*tauxz)
!               eyz = fact*(Wyx*tauxz + Wyy*tauyz + Wyz*tauzz + &
!                    Wzx*tauxy + Wzy*tauyy + Wzz*tauyz)

!               ! Add extra terms
!               tauxx = tauxx - exx
!               tauyy = tauyy - eyy
!               tauzz = tauzz - ezz
!               tauxy = tauxy - exy
!               tauxz = tauxz - exz
!               tauyz = tauyz - eyz

!            end if

!            ! Compute the average velocities for the face. Remember that
!            ! the velocities are stored and not the momentum.

!            ubar = half*(w(i,j,k,ivx) + w(i+1,j,k,ivx))
!            vbar = half*(w(i,j,k,ivy) + w(i+1,j,k,ivy))
!            wbar = half*(w(i,j,k,ivz) + w(i+1,j,k,ivz))

!            ! Compute the viscous fluxes for this i-face.

!            fmx   = tauxx*si(i,j,k,1) + tauxy*si(i,j,k,2) &
!                 + tauxz*si(i,j,k,3)
!            fmy   = tauxy*si(i,j,k,1) + tauyy*si(i,j,k,2) &
!                 + tauyz*si(i,j,k,3)
!            fmz   = tauxz*si(i,j,k,1) + tauyz*si(i,j,k,2) &
!                 + tauzz*si(i,j,k,3)
!            frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*si(i,j,k,1) &
!                 + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*si(i,j,k,2) &
!                 + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*si(i,j,k,3) &
!                 - q_x*si(i,j,k,1) - q_y*si(i,j,k,2) - q_z*si(i,j,k,3)

!            ! Update the residuals of cell i and i+1.

!            fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
!            fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
!            fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
!            fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE

!            fw(i+1,j,k,imx)   = fw(i+1,j,k,imx)   + fmx
!            fw(i+1,j,k,imy)   = fw(i+1,j,k,imy)   + fmy
!            fw(i+1,j,k,imz)   = fw(i+1,j,k,imz)   + fmz
!            fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) + frhoE

!         enddo
!      enddo
!   enddo

!   ! ---------------------------------------------
!   !                 Sum dw and fw/res scale
!   ! ---------------------------------------------
!   nTurb = nt2-nt1+1
!   do k=2, kl
!      do j=2, jl
!         do i=2, il
!            oVol = one/volRef(i,j,k)*max(real(iblank(i,j,k), realType), zero)
!            dw(i, j, k, 1:nwf) = (dw(i, j, k, 1:nwf) + fw(i, j, k, 1:nwf))*oVol
!            dw(i,j,k, nt1:nt2) = dw(i,j,k,nt1:nt2)*oVol*turbResScale(1:nTurb)

!         end do
!      end do
!   end do


! end subroutine superCode




! ! Try to run superCode
!  timeA = mpi_wtime()

!  do iIter=1, nIter
!     ! Copy out stuff:
!     call setPointers(1, 1, 1)
!     loopCOunt = 0
!     do kk=2, kl-BS, BS
!        do jj=2, jl-BS, BS
!           do ii=2, il-BS, BS
!              loopCount = loopCount + 1
!              ! Double halos
!              do k=0, BS+2
!                 do j=0, BS+2
!                    do i=0, BS+2
!                       ww(i,j,k,:) = w(i+ii,j+jj,k+kk,:)
!                       xx(i,j,k,:) = x(i+ii,j+jj,k+kk,:)
!                       pp(i,j,k) = p(i+ii,j+jj,k+kk)
!                       ggamma(i,j,k) = gamma(i+ii,j+jj,k+kk)
!                    end do
!                 end do
!              end do

!              ! Single Halos
!              do k=1, BS+2
!                 do j=1, BS+2
!                    do i=1, BS+2
!                       rrlv(i,j,k) = rlv(i+ii,j+jj,k+kk)
!                       rrev(i,j,k) = rev(i+ii,j+jj,k+kk)
!                       vvol(i,j,k) = vol(i+ii,j+jj,k+kk)
!                    end do
!                 end do
!              end do

!              ! Node equivalent
!              do k=2, BS
!                 do j=2, BS
!                    do i=1, BS
!                       pporI(i,j,k) = porI(i+ii,j+jj,k+kk)
!                    end do
!                 end do
!              end do

!              do k=2, BS
!                 do j=1, BS
!                    do i=2, BS
!                       pporJ(i,j,k) = porJ(i+ii,j+jj,k+kk)
!                    end do
!                 end do
!              end do

!              do k=1, BS
!                 do j=2, BS
!                    do i=2, BS
!                       pporK(i,j,k) = porK(i+ii,j+jj,k+kk)
!                    end do
!                 end do
!              end do

!              ! No Halos
!              do k=2, BS+1
!                 do j=2, BS+1
!                    do i=2, BS+1
!                       iiblank(i,j,k) = iblank(i+ii,j+jj,k+kk)
!                       dd2wall(i,j,k) = d2wall(i+ii,j+jj,k+kk)
!                       vvolRef(i,j,k) = volRef(i+ii,j+jj,k+kk)
!                    end do
!                 end do
!              end do

!              ! fw would normally have a useful value.
!              ffw = zero

!              call superCode(ww, xx, vvol, PP, rrlv, rrev, ggamma, pporI, pporJ, pporK, dd2wall, &
!                   vvolRef, ddw, ffw, iiblank)
!           end do
!        end do
!     end do
!  end do
!  timeB = mpi_wtime()

!  speed = loopCount*BS**3*nIter/(timeB-timeA)/1e6

!  print *,'speed2:', speed,  timeB-timeA, loopCount

! ! Super block inputs
! real(kind=realType), dimension(0:BS+2, 0:BS+2, 0:BS+2, 1:nw) :: ww
! real(kind=realType), dimension(0:BS+2, 0:BS+2, 0:BS+2, 3) :: xx
! real(kind=realType), dimension(1:BS+2, 1:BS+2, 1:BS+2) :: &
!      rrlv, rrev, vvol
! real(kind=realType), dimension(2:BS+1, 2:BS+1, 2:BS+1) :: vvolRef, dd2wall
! real(kind=realType), dimension(0:BS+2, 0:BS+2, 0:BS+2) :: PP, ggamma
! integer(kind=porType), dimension(1:BS+1, 2:BS+1, 2:BS+1) :: pporI
! integer(kind=porType), dimension(2:BS+1, 1:BS+1, 2:BS+1) :: pporJ
! integer(kind=porType), dimension(2:BS+1, 2:BS+1, 1:BS+1) :: pporK
! integer(kind=intType), dimension(2:BS+1, 2:BS+1, 2:BS+1) :: iiblank

! ! Super block Ouputs
! real(kind=realType), dimension(1:BS+2, 1:BS+2, 1:BS+2, 1:nw) :: ffw, ddw



!     !$OMP& private(bnx, bny, bnz, bil, bjl, bkl, bie, bje, bke, bib, bjb, bkb) &
!     !$OMP& private(singleHaloStart, doubleHaloStart, nodeStart) &
!     !$OMP& private(bw, bp, bgamma, bss, bx, brlv, brev, bvol, baa, bradI, bradJ, bradK) &
!     !$OMP& private(bdss, bvolRef, bd2wall, biblank, bporI, bporJ, bporK, bfw, bdw) &
!     !$OMP& private(bux, buy, buz, bvx, bvy, bvz, bwx, bwy, bwz, bqx, bqy, bqz) & !


     !!$OMP& private(nx, ny, nz, il, jl, kl, ie, je, ke, ib, jb, kb) &
     !!$OMP& private(singleHaloStart, doubleHaloStart, nodeStart) &
     !!$OMP& private(w, p, gamma, ss, x, rlv, rev, vol, aa, radI, radJ, radK) &
     !!$OMP& private(dss, volRef, d2wall, iblank, porI, porJ, porK, fw, dw) &
     !!$OMP& private(ux, uy, uz, vx, vy, vz, wx, wy, wz, qx, qy, qz) &
