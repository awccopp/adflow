module blockette_state_b

    use constants
    ! This temporary module contains all cache-blocked code. It also
    ! contains the statically allocated variables on which the blocked
    ! code operates.

    ! Dummy Block dimensions
    integer(kind=intType), parameter :: BS=8
    integer(kind=intType), parameter :: bbil=BS+1, bbjl=BS+1, bbkl=BS+1
    integer(kind=intType), parameter :: bbie=BS+2, bbje=BS+2, bbke=BS+2
    integer(kind=intType), parameter :: bbib=BS+3, bbjb=BS+3, bbkb=BS+3

    ! Actual dimensions to execute
    integer(kind=intType) :: nx, ny, nz, il, jl, kl, ie, je, ke, ib, jb, kb

    ! Variables to track transferring variables between blockettes
    integer(kind=intType) :: singleHaloStart, doubleHaloStart, nodeStart

    ! Current indices into the original block
    integer(kind=intType) :: ii, jj, kk

    ! Double halos
    real(kind=realType), dimension(0:bbib, 0:bbjb, 0:bbkb, 1:6) :: w
    real(kind=realType), dimension(0:bbib, 0:bbjb, 0:bbkb, 1:6) :: wd
    real(kind=realType), dimension(0:bbib, 0:bbjb, 0:bbkb) :: P, gamma
    ! real(kind=realType), dimension(0:bbib, 0:bbjb, 0:bbkb) :: ss ! Entropy

    ! Double halo derivative seeds
    real(kind=realType), dimension(0:bbib, 0:bbjb, 0:bbkb) :: pd

    ! Single halos
    real(kind=realType), dimension(0:bbie, 0:bbje, 0:bbke, 3) :: x
    real(kind=realType), dimension(1:bbie, 1:bbje, 1:bbke):: rlv, rev, vol, aa
    real(kind=realType), dimension(1:bbie, 1:bbje, 1:bbke) :: radI, radJ, radK, dtl
    ! real(kind=realType),dimension(1:bbie, 1:bbje, 1:bbke, 3) :: dss ! Shock sensor

    ! Single halo derivative seeds
    real(kind=realType), dimension(1:bbie, 1:bbje, 1:bbke):: rlvd, revd
    real(kind=realType), dimension(1:bbie, 1:bbje, 1:bbke) :: radId, radJd, radKd

    ! Single halo access-only variables ? Check if this is access only
    real(kind=realType), dimension(1:bbie, 1:bbje, 1:bbke):: aad

    ! No halos
    real(kind=realType), dimension(2:bbil, 2:bbjl, 2:bbkl) :: volRef, d2wall
    integer(kind=intType), dimension(2:bbil, 2:bbjl, 2:bbkl) :: iblank

    ! Single halo read derivative seed
    real(kind=realType), dimension(1:bbie, 1:bbje, 1:bbke, 1:6) :: dwd

    ! No halo access-only variables
    real(kind=realType), dimension(2:bbil, 2:bbjl, 2:bbkl, 1:1) :: scratchd

    ! Face Porosities
    integer(kind=porType), dimension(1:bbil, 2:bbjl, 2:bbkl) :: porI
    integer(kind=porType), dimension(2:bbil, 1:bbjl, 2:bbkl) :: porJ
    integer(kind=porType), dimension(2:bbil, 2:bbjl, 1:bbkl) :: porK

    ! Single halos (only owned cells significant)
    real(kind=realType), dimension(1:bbie, 1:bbje, 1:bbke, 1:5) :: fw
    real(kind=realType), dimension(1:bbie, 1:bbje, 1:bbke, 1:6) :: dw

    ! Single halos (only owned cells significant), access-only variables
    real(kind=realType), dimension(1:bbie, 1:bbje, 1:bbke, 1:5) :: fwd

    ! Face projected areas
    real(kind=realType), dimension(0:bbie, 1:bbje, 1:bbke, 3) :: sI
    real(kind=realType), dimension(1:bbie, 0:bbje, 1:bbke, 3) :: sJ
    real(kind=realType), dimension(1:bbie, 1:bbje, 0:bbke, 3) :: sK

    ! Face velocities
    real(kind=realType), dimension(0:bbie, 1:bbje, 1:bbke) :: sFaceI
    real(kind=realType), dimension(1:bbie, 0:bbje, 1:bbke) :: sFaceJ
    real(kind=realType), dimension(1:bbie, 1:bbje, 0:bbke) :: sFaceK

    ! Nodal gradients
    real(kind=realType), dimension(1:bbil, 1:bbjl, 1:bbkl) :: ux, uy, uz
    real(kind=realType), dimension(1:bbil, 1:bbjl, 1:bbkl) :: vx, vy, vz
    real(kind=realType), dimension(1:bbil, 1:bbjl, 1:bbkl) :: wx, wy, wz
    real(kind=realType), dimension(1:bbil, 1:bbjl, 1:bbkl) :: qx, qy, qz

    ! And their derivative seeds
    real(kind=realType), dimension(1:bbil, 1:bbjl, 1:bbkl) :: uxd, uyd, uzd
    real(kind=realType), dimension(1:bbil, 1:bbjl, 1:bbkl) :: vxd, vyd, vzd
    real(kind=realType), dimension(1:bbil, 1:bbjl, 1:bbkl) :: wxd, wyd, wzd
    real(kind=realType), dimension(1:bbil, 1:bbjl, 1:bbkl) :: qxd, qyd, qzd

    ! Make *all* of these variables tread-private
    !$OMP THREADPRIVATE(nx, ny, nz, il, jl, kl, ie, je, ke, ib, jb, kb)
    !$OMP THREADPRIVATE(w, p, gamma, ss, x, rlv, rev, vol, aa, radI, radJ, radK)
    !$OMP THREADPRIVATE(dss, volRef, d2wall, iblank, porI, porJ, porK, fw, dw)
    !$OMP THREADPRIVATE(sI, sJ, sK, ux, uy, uz, vx, vy, vz, wx, wy, wz, qx, qy, qz)
contains

! not sure if we need this, check later.
#ifndef USE_COMPLEX
    subroutine blockette_fast_b

        ! This is specialized form of master that *ONLY* computes drdw
        ! products. It uses a few specialzed routines that are
        ! differentiated without including spatial dependencies. This
        ! results in slightly faster code. This specialization is
        ! justififed since this routine is needed for the transpose
        ! matrix-vector products in solving the adjoint system and thus
        ! this routine is called several orders of magniutde more than
        ! master_b. This routine has to be fast!

        use constants
        use iteration, only : currentLevel
        use flowVarRefState, only : nw, viscous
        ! use blockPointers, only : nDom, il, jl, kl, wd, dwd, iblank
        use blockPointers, only : &
             bnx=>nx, bny=>ny, bnz=>nz, &
             bil=>il, bjl=>jl, bkl=>kl, &
             bie=>ie, bje=>je, bke=>ke, &
             bib=>ib, bjb=>jb, bkb=>kb, &
             bw=>w, bp=>p, bgamma=>gamma, baa=>aa, baad=>aad, &
             bx=>x, brlv=>rlv, brev=>rev, bvol=>vol, bVolRef=>volRef, bd2wall=>d2wall, &
             biblank=>iblank, bPorI=>porI, bPorJ=>porJ, bPorK=>porK, bdw=>dw, bfw=>fw, &
             bShockSensor=>shockSensor, &
             bsi=>si, bsj=>sj, bsk=>sk, &
             bsFaceI=>sFaceI, bsFaceJ=>sFaceJ, bsFaceK=>sFaceK , &
             bdtl=>dtl,  &
             addGridVelocities, &
             bwd=>wd, brlvd=>rlvd, brevd=>revd, &
             bradi=>radi, bradj=>radj, bradk=>radk, &
             bradid=>radid, bradjd=>radjd, bradkd=>radkd, &
             bpd=>pd, bdwd=>dwd, &
             bux=>ux, buy=>uy, buz=>uz, &
             bvx=>vx, bvy=>vy, bvz=>vz, &
             bwx=>wx, bwy=>wy, bwz=>wz, &
             bqx=>qx, bqy=>qy, bqz=>qz
        use inputPhysics, only : equationMode, turbModel, equations
        use inputDiscretization, only : lowSpeedPreconditioner, spaceDiscr
        use inputTimeSpectral, only : nTimeIntervalsSpectral
        use sa_fast_b, only : qq

        ! NEED ALL OF THESE HERE
        ! use adjointextra_b, only : resscale_B, sumdwandfw_b
        ! use flowutils_b, only : computeSpeedOfSoundSquared_b
        ! use sa_fast_b, only : saresscale_fast_b, saviscous_fast_b, sasource_fast_b
        ! use turbutils_fast_b, only : turbAdvection_fast_b
        ! use fluxes_fast_b, only :inviscidUpwindFlux_fast_b, inviscidDissFluxScalar_fast_b, &
        ! inviscidDissFluxMatrix_fast_b, viscousFlux_fast_b, inviscidCentralFlux_fast_b
        ! use solverutils_fast_b, only : timeStep_block_fast_b
        ! use flowutils_fast_b, only : allnodalgradients_fast_b

        implicit none

        ! Working Variables
        integer(kind=intType) :: i, j, k, l

        ! here is an exhaustive list of the variables that we are interested in
        ! use blockPointers, only : nDom, il, jl, kl, wd, dwd, iblank
        ! use blockpointers, only : il, jl, kl, nx, ny, nz, volref, dw, dwd
        ! use blockpointers, only : il, jl, kl, dw, dwd, fw, fwd, iblank
        ! use blockpointers IN VISCOUSFLUX_FAST_B, ALLNODALGRADIENTS_FAST_B
        ! use blockpointers, only : ie, je, ke, w, wd, p, pd, aa, aad, gamma
        ! use blockpointers, only : nx, ny, nz, il, jl, kl, ie, je, ke, ib, &
        ! &   jb, kb, w, wd, p, pd, pori, porj, pork, fw, fwd, radi, radid, radj, &
        ! &   radjd, radk, radkd, gamma
        ! use blockpointers, only : nx, ny, nz, il, jl, kl, ie, je, ke, ib, &
        ! &   jb, kb, w, wd, p, pd, pori, porj, pork, fw, fwd, gamma, si, sj, sk, &
        ! &   indfamilyi, indfamilyj, indfamilyk, spectralsol, addgridvelocities, &
        ! &   sfacei, sfacej, sfacek, factfamilyi, factfamilyj, factfamilyk
        ! use blockpointers, only : il, jl, kl, ie, je, ke, ib, jb, kb, w, &
        ! &   wd, p, pd, pori, porj, pork, fw, fwd, gamma, si, sj, sk, indfamilyi,&
        ! &   indfamilyj, indfamilyk, spectralsol, addgridvelocities, sfacei, &
        ! &   sfacej, sfacek, rotmatrixi, rotmatrixj, rotmatrixk, factfamilyi, &
        ! &   factfamilyj, factfamilyk
        ! use blockpointers, only : nx, il, ie, ny, jl, je, nz, kl, ke, &
        ! &   spectralsol, w, wd, si, sj, sk, dw, dwd, pori, porj, pork, &
        ! &   indfamilyi, indfamilyj, indfamilyk, p, pd, sfacei, sfacej, sfacek, &
        ! &   nbkglobal, addgridvelocities, blockismoving, vol, factfamilyi, &
        ! &   factfamilyj, factfamilyk
        ! use blockpointers, IN SARESSCALE_FAST_B, SAVISCOUS_FAST,
        !     use blockpointers, only : nx, ny, nz, il, jl, kl, vol, sfacei, &
        ! &   sfacej, sfacek, w, wd, si, sj, sk, addgridvelocities, bmti1, bmti2, &
        ! &   bmtj1, bmtj2, bmtk1, bmtk2, scratch, scratchd
        ! use blockpointers, IN SASOURCE_FAST_B
        !     use blockpointers, only : ie, je, ke, il, jl, kl, w, wd, p, pd, &
        ! &   rlv, rlvd, rev, revd, radi, radid, radj, radjd, radk, radkd, si, sj,&
        ! &   sk, sfacei, sfacej, sfacek, dtl, gamma, vol, addgridvelocities, &
        ! &   sectionid



        ! We first need to copy the relevant variables from the block to the blockettes
        ! Block loop over the owned cells
        !$OMP parallel do private(i,j,k,l) collapse(2)
        do kk=2, bkl, BS
            do jj=2, bjl, BS
                do ii=2, bil, BS

                    ! Determine the actual size this block will be and set
                    ! the sizes in the blockette module for each of the
                    ! subroutines.

                    nx = min(ii+BS-1, bil) - ii + 1
                    ny = min(jj+BS-1, bjl) - jj + 1
                    nz = min(kk+BS-1, bkl) - kk + 1

                    il = nx + 1; jl = ny + 1; kl = nz + 1
                    ie = nx + 2; je = ny + 2; ke = nz + 2
                    ib = nx + 3; jb = ny + 3; kb = nz + 3

                    ! firstBlockette: if (ii==2) then
                    !
                    !     ! First loop. Need to compute the extra stuff. Set
                    !     ! the generic starts and copy the extra
                    !     ! variables in to the starting slots
                    !     singleHaloStart = 1
                    !     doubleHaloStart = 0
                    !     nodeStart = 1
                    !
                    !     ! Double halos
                    !     do k=0, kb
                    !         do j=0, jb
                    !             do i=0, 3
                    !                 w(i,j,k,1:nw) = bw(i+ii-2, j+jj-2, k+kk-2, 1:nw)
                    !                 p(i,j,k) = bP(i+ii-2, j+jj-2, k+kk-2)
                    !                 gamma(i,j,k) = bgamma(i+ii-2, j+jj-2, k+kk-2)
                    !                 ! if (currentLevel == 1) then
                    !                 !     ss(i,j,k) = bShockSensor(i+ii-2, j+jj-2,k+kk-2)
                    !                 ! end if
                    !             end do
                    !         end do
                    !     end do
                    !
                    !     ! Single halos
                    !     do k=1, ke
                    !         do j=1, je
                    !             do i=1, 2
                    !                 rlv(i,j,k) = brlv(i+ii-2, j+jj-2, k+kk-2)
                    !                 rev(i,j,k) = brev(i+ii-2, j+jj-2, k+kk-2)
                    !                 vol(i,j,k) = bvol(i+ii-2, j+jj-2, k+kk-2)
                    !             end do
                    !         end do
                    !     end do
                    !
                    !     ! X
                    !     do k=0, ke
                    !         do j=0, je
                    !             do i=0, 1
                    !                 x(i,j,k,:) = bx(i+ii-2, j+jj-2, k+kk-2, :)
                    !             end do
                    !         end do
                    !     end do
                    ! else
                    !
                    !     ! Subsequent loop. We can save a bunch of work by
                    !     ! copying some of the pre-computed values from the
                    !     ! previous blockette to this blockette. Basically the
                    !     ! values that are at the "I end" get shuffled back to
                    !     ! the I-start. We *also* do this for some of the
                    !     ! intermediate variables that are costly to compute
                    !     ! like the nodal gradients, and spectral radius which
                    !     ! helps cut back on the amount of data duplication.
                    !
                    !     ! Important Note: This cell is not the first cell. If
                    !     ! this code is being executed, the previous blockette
                    !     ! was copied fully in the i direction.
                    !     ! Therefore, we can just copy the values from
                    !     ! the end of the blockette as it is allocated.
                    !     ! To do this, we ignore the dimensions of the "current"
                    !     ! blockette, and just take the baseline BS dimensions
                    !     ! as the current blockette might be partially filled
                    !     ! in the i direction.
                    !
                    !     singleHaloStart = 3
                    !     doubleHaloStart = 4
                    !     nodeStart = 2
                    !
                    !     ! Double halos
                    !     do k=0, kb
                    !         do j=0, jb
                    !             do i=0, 3
                    !                 w(i,j,k,1:nw) = w(BS+i, j, k, 1:nw)
                    !                 p(i,j,k) = p(BS+i, j, k)
                    !                 gamma(i,j,k) = gamma(BS+i, j, k)
                    !                 ! ss(i,j,k) = ss(BS+i, j, k)
                    !             end do
                    !         end do
                    !     end do
                    !
                    !     ! Single halos
                    !     do k=1, ke
                    !         do j=1, je
                    !             do i=1, 2
                    !                 rlv(i,j,k) = rlv(BS+i, j, k)
                    !                 rev(i,j,k) = rev(BS+i, j, k)
                    !                 vol(i,j,k) = vol(BS+i, j, k)
                    !
                    !                 ! Computed variables
                    !
                    !                 ! DONT Copy the spectral-radii. The loop that calculates
                    !                 ! spectral radii also calculates portion of the time step,
                    !                 ! so we don't want to mess with its boundaries to keep
                    !                 ! it simple.
                    !                 aa(i,j,k) = aa(BS+i, j, k)
                    !                 ! dss(i,j,k,:) = dss(BS+i, j, k, :)
                    !             end do
                    !         end do
                    !     end do
                    !
                    !     ! X
                    !     do k=0, ke
                    !         do j=0, je
                    !             do i=0, 1
                    !                 x(i,j,k,:) = x(BS+i, j, k, :)
                    !             end do
                    !         end do
                    !     end do
                    !
                    !     ! ! Nodal gradients
                    !     ! do k=1, kl
                    !     !     do j=1, jl
                    !     !         ux(1, j, k) = ux(BS+1, j, k)
                    !     !         uy(1, j, k) = uy(BS+1, j, k)
                    !     !         uz(1, j, k) = uz(BS+1, j, k)
                    !     !
                    !     !         vx(1, j, k) = vx(BS+1, j, k)
                    !     !         vy(1, j, k) = vy(BS+1, j, k)
                    !     !         vz(1, j, k) = vz(BS+1, j, k)
                    !     !
                    !     !         wx(1, j, k) = wx(BS+1, j, k)
                    !     !         wy(1, j, k) = wy(BS+1, j, k)
                    !     !         wz(1, j, k) = wz(BS+1, j, k)
                    !     !
                    !     !         qx(1, j, k) = qx(BS+1, j, k)
                    !     !         qy(1, j, k) = qy(BS+1, j, k)
                    !     !         qz(1, j, k) = qz(BS+1, j, k)
                    !     !     end do
                    !     ! end do
                    ! end if firstBlockette

                    ! -------------------------------------
                    !      Fill in the remaining values
                    ! -------------------------------------
                    ! to test, fill in all values

                    ! Double halos
                    do k=0, kb
                        do j=0, jb
                            do i=0, ib
                                w(i,j,k,1:nw) = bw(i+ii-2, j+jj-2, k+kk-2, 1:nw)
                                p(i,j,k) = bP(i+ii-2, j+jj-2, k+kk-2)
                                gamma(i,j,k) = bgamma(i+ii-2, j+jj-2, k+kk-2)
                                ! if (currentLevel == 1) then
                                !     ss(i,j,k) = bShockSensor(i+ii-2, j+jj-2,k+kk-2)
                                ! end if
                            end do
                        end do
                    end do

                    ! Single halos
                    do k=1, ke
                        do j=1, je
                            do i=1, ie
                                rlv(i,j,k) = brlv(i+ii-2, j+jj-2, k+kk-2)
                                rev(i,j,k) = brev(i+ii-2, j+jj-2, k+kk-2)
                                vol(i,j,k) = bvol(i+ii-2, j+jj-2, k+kk-2)
                                aa(i,j,k) = baa(i+ii-2, j+jj-2, k+kk-2)
                                radi(i,j,k) = bradi(i+ii-2, j+jj-2, k+kk-2)
                                radj(i,j,k) = bradj(i+ii-2, j+jj-2, k+kk-2)
                                radk(i,j,k) = bradk(i+ii-2, j+jj-2, k+kk-2)
                            end do
                        end do
                    end do

                    ! X
                    do k=0, ke
                        do j=0, je
                            do i=0, ie
                                x(i,j,k,:) = bx(i+ii-2, j+jj-2, k+kk-2, :)
                            end do
                        end do
                    end do

                    ! No Halos (no change)
                    do k=2, kl
                        do j=2, jl
                            do i=2, il
                                iblank(i,j,k) = biblank(i+ii-2,j+jj-2,k+kk-2)
                                if (equations .eq. ransequations) &
                                d2wall(i,j,k) = bd2wall(i+ii-2,j+jj-2,k+kk-2)
                                volRef(i,j,k) = bvolRef(i+ii-2,j+jj-2,k+kk-2)
                            end do
                        end do
                    end do

                    ! AD seeds for the residual
                    dwd = zero
                    do l=1,6
                      do k=2, kl
                        do j=2, jl
                          do i=2, il
                            dwd(i,j,k,l) = bdwd(i+ii-2,j+jj-2,k+kk-2,l)
                          end do
                        end do
                      end do
                    end do

                    ! Porosities (no change)
                    do k=2, kl
                        do j=2, jl
                            do i=1, il
                                porI(i,j,k) = bporI(i+ii-2,j+jj-2,k+kk-2)
                            end do
                        end do
                    end do

                    do k=2, kl
                        do j=1, jl
                            do i=2, il
                                PorJ(i,j,k) = bporJ(i+ii-2,j+jj-2,k+kk-2)
                            end do
                        end do
                    end do

                    do k=1, kl
                        do j=2, jl
                            do i=2, il
                                PorK(i,j,k) = bporK(i+ii-2,j+jj-2,k+kk-2)
                            end do
                        end do
                    end do

                    ! Face velocities if necessary
                    ! we had a bug here
                    if (addGridVelocities) then
                        do k=1, ke
                            do j=1, je
                                do i=0, ie
                                    sFaceI(i, j, k) = bsFaceI(i+ii-2, j+jj-2, k+kk-2)
                                end do
                            end do
                        end do

                        do k=1, ke
                            do j=0, je
                                do i=1, ie
                                    sFaceJ(i, j, k) = bsFaceJ(i+ii-2, j+jj-2, k+kk-2)
                                end do
                            end do
                        end do

                        do k=0, ke
                            do j=1, je
                                do i=1, ie
                                    sFaceK(i, j, k) = bsFaceK(i+ii-2, j+jj-2, k+kk-2)
                                end do
                            end do
                        end do
                    else
                        sFaceI = zero
                        sFaceJ = zero
                        sFaceK = zero
                    end if

                    ! Nodal gradients
                    do k=1, kl
                        do j=1, jl
                            do i=1, il
                                ! ux(i, j, k) = bsFaceK(ii+ii-2, j+jj-2, k+kk-2)

                                ux(i, j, k) = bux(i+ii-2, j+jj-2, k+kk-2)
                                uy(i, j, k) = buy(i+ii-2, j+jj-2, k+kk-2)
                                uz(i, j, k) = buz(i+ii-2, j+jj-2, k+kk-2)

                                vx(i, j, k) = bvx(i+ii-2, j+jj-2, k+kk-2)
                                vy(i, j, k) = bvy(i+ii-2, j+jj-2, k+kk-2)
                                vz(i, j, k) = bvz(i+ii-2, j+jj-2, k+kk-2)

                                wx(i, j, k) = bwx(i+ii-2, j+jj-2, k+kk-2)
                                wy(i, j, k) = bwy(i+ii-2, j+jj-2, k+kk-2)
                                wz(i, j, k) = bwz(i+ii-2, j+jj-2, k+kk-2)

                                qx(i, j, k) = bqx(i+ii-2, j+jj-2, k+kk-2)
                                qy(i, j, k) = bqy(i+ii-2, j+jj-2, k+kk-2)
                                qz(i, j, k) = bqz(i+ii-2, j+jj-2, k+kk-2)
                            end do
                        end do
                    end do

                    ! Clear the viscous flux before we start.
                    !fw = zero
                    ! transfer the visc. fluxes
                    do l=1, 5
                       do k=1, ke
                          do j=1, je
                             do i=1, ie
                                fw(i, j, k, l) = bfw(i+ii-2,j+jj-2,k+kk-2,l)
                                ! print*,fw(i, j, k, l)
                             end do
                          end do
                       end do
                   end do

                    ! Also clear out the access-only variables as these are all derivative seeds
                    scratchd = zero
                    wd       = zero
                    aad      = zero
                    fwd      = zero
                    pd       = zero
                    rlvd = zero
                    revd = zero
                    radid = zero
                    radjd = zero
                    radkd = zero

                    uxd = zero
                    uyd = zero
                    uzd = zero
                    vxd = zero
                    vyd = zero
                    vzd = zero
                    wxd = zero
                    wyd = zero
                    wzd = zero
                    qxd = zero
                    qyd = zero
                    qzd = zero

                call metrics

                ! Now we start running back through the main residual code:
                call resScale_b
                call sumDwAndFw_b

                ! if (lowSpeedPreconditioner) then
                !    call applyLowSpeedPreconditioner_b
                ! end if

                ! Note that master_b does not include the approximation codes
                ! as those are never needed in reverse.
                if (viscous) then
                    call viscousFlux_fast_b
                    call allNodalGradients_fast_b
                    call computeSpeedOfSoundSquared_b
                end if

                select case (spaceDiscr)
                case (dissScalar)
                    call inviscidDissFluxScalar_fast_b
                case (dissMatrix)
                    call inviscidDissFluxMatrix_fast_b
                case (upwind)
                    call inviscidUpwindFlux_fast_b(.True.)
                end select

                call inviscidCentralFlux_fast_b

                ! Compute turbulence residual for RANS equations
                if( equations == RANSEquations) then
                    select case (turbModel)
                    case (spalartAllmaras)
                        call saResScale_fast_b
                        call saViscous_fast_b
                        !call unsteadyTurbTerm_b(1_intType, 1_intType, itu1-1, qq)
                        call turbAdvection_fast_b(1_intType, 1_intType, itu1-1, qq)
                        call saSource_fast_b
                    end select

                    !call unsteadyTurbSpectral_block_b(itu1, itu1, nn, sps)
                end if

                call timeStep_block_fast_b(.false.)

                ! we now need to copy the blockette variables to the blocks and we are done
                ! Note: we need to write out all the AD seeds that we have modified, turns out this list is pretty extensive, so be careful here. Furthermore, we need to accumulate AD seeds for overlapping variables (i.e. single and double halos)
                ! Now we can just set the part of dw we computed
                ! (owned cells only) and we're done!

                ! Write out the modified AD seeds for the state
                do l=1, 6
                   ! do k=2, kl
                   !    do j=2, jl
                   !       do i=2, il
                    do k=0, kb
                       do j=0, jb
                          do i=0, ib
                            bwd(i+ii-2,j+jj-2,k+kk-2,l) = bwd(i+ii-2,j+jj-2,k+kk-2,l) + wd(i,j,k,l)
                            ! if (wd(i,j,k,l).ne.zero) print *,"state seed not zero!"
                         end do
                      end do
                   end do
               end do

               ! Single halos
               do k=1, ke
                   do j=1, je
                       do i=1, ie
                          brlvd(i+ii-2, j+jj-2, k+kk-2) = brlvd(i+ii-2, j+jj-2, k+kk-2) + rlvd(i,j,k)
                          brevd(i+ii-2, j+jj-2, k+kk-2) = brevd(i+ii-2, j+jj-2, k+kk-2) + revd(i,j,k)
                          baad(i+ii-2, j+jj-2, k+kk-2) = baad(i+ii-2, j+jj-2, k+kk-2) + aad(i,j,k)
                          bradid(i+ii-2, j+jj-2, k+kk-2) = bradid(i+ii-2, j+jj-2, k+kk-2) + radid(i,j,k)
                          bradjd(i+ii-2, j+jj-2, k+kk-2) = bradjd(i+ii-2, j+jj-2, k+kk-2) + radjd(i,j,k)
                          bradkd(i+ii-2, j+jj-2, k+kk-2) = bradkd(i+ii-2, j+jj-2, k+kk-2) + radkd(i,j,k)
                       end do
                   end do
               end do

               ! Double halos
               do k=0, kb
                   do j=0, jb
                       do i=0, ib
                           bPd(i+ii-2, j+jj-2, k+kk-2) = bPd(i+ii-2, j+jj-2, k+kk-2) + pd(i,j,k)
                       end do
                   end do
               end do



               ! ! Also copy out the dtl if we were asked for it
               ! if (updateDt) then
               !     do k=2, kl
               !         do j=2, jl
               !             do i=2, il
               !                 bdtl(i+ii-2, j+jj-2, k+kk-2) = dtl(i, j, k)
               !             end do
               !         end do
               !     end do
               ! end if

           end do
       end do
   end do
   !$OMP END PARALLEL DO


    end subroutine blockette_fast_b
#endif

  subroutine metrics
    ! ---------------------------------------------
    !              Metric computation
    ! ---------------------------------------------

    use constants
    use blockPointers, only : rightHanded
    implicit none

    integer(kind=intType) :: i, j, k, l, m, n
    real(kind=realType), dimension(3) :: v1, v2
    real(kind=realType) :: fact

    ! Projected areas of cell faces in the i direction.
    if (rightHanded) then
       fact = half
    else
       fact = -half
    end if
    do k=1,ke
       n = k -1
       do j=1,je
          m = j -1
          do i=0,ie

             ! Determine the two diagonal vectors of the face.

             v1(1) = x(i,j,n,1) - x(i,m,k,1)
             v1(2) = x(i,j,n,2) - x(i,m,k,2)
             v1(3) = x(i,j,n,3) - x(i,m,k,3)

             v2(1) = x(i,j,k,1) - x(i,m,n,1)
             v2(2) = x(i,j,k,2) - x(i,m,n,2)
             v2(3) = x(i,j,k,3) - x(i,m,n,3)

             ! The face normal, which is the cross product of the two
             ! diagonal vectors times fact; remember that fact is
             ! either -0.5 or 0.5.

             si(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
             si(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
             si(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

          enddo
       enddo
    enddo

    ! Projected areas of cell faces in the j direction.

    do k=1,ke
       n = k -1
       do j=0,je
          do i=1,ie
             l = i -1

             ! Determine the two diagonal vectors of the face.

             v1(1) = x(i,j,n,1) - x(l,j,k,1)
             v1(2) = x(i,j,n,2) - x(l,j,k,2)
             v1(3) = x(i,j,n,3) - x(l,j,k,3)

             v2(1) = x(l,j,n,1) - x(i,j,k,1)
             v2(2) = x(l,j,n,2) - x(i,j,k,2)
             v2(3) = x(l,j,n,3) - x(i,j,k,3)

             ! The face normal, which is the cross product of the two
             ! diagonal vectors times fact; remember that fact is
             ! either -0.5 or 0.5.

             sj(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
             sj(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
             sj(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

          enddo
       enddo
    enddo

    ! Projected areas of cell faces in the k direction.

    do k=0,ke
       do j=1,je
          m = j -1
          do i=1,ie
             l = i -1

             ! Determine the two diagonal vectors of the face.

             v1(1) = x(i,j,k,1) - x(l,m,k,1)
             v1(2) = x(i,j,k,2) - x(l,m,k,2)
             v1(3) = x(i,j,k,3) - x(l,m,k,3)

             v2(1) = x(l,j,k,1) - x(i,m,k,1)
             v2(2) = x(l,j,k,2) - x(i,m,k,2)
             v2(3) = x(l,j,k,3) - x(i,m,k,3)

             ! The face normal, which is the cross product of the two
             ! diagonal vectors times fact; remember that fact is
             ! either -0.5 or 0.5.

             sk(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
             sk(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
             sk(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

          enddo
       enddo
    enddo
  end subroutine metrics

!  differentiation of resscale in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: *dw
!   with respect to varying inputs: *dw
!   rw status of diff variables: *dw:in-out
!   plus diff mem management of: dw:in
subroutine resscale_b()
    use constants
    ! use blockpointers, only : il, jl, kl, nx, ny, nz, volref, dw, dwd
    use flowvarrefstate, only : nwf, nt1, nt2
    use inputiteration, only : turbresscale
    implicit none
    ! local variables
    integer(kind=inttype) :: i, j, k, ii, nturb
    real(kind=realtype) :: ovol
    intrinsic mod
    ! divide through by the reference volume
    nturb = nt2 - nt1 + 1
    do ii=0,nx*ny*nz-1
        i = mod(ii, nx) + 2
        j = mod(ii/nx, ny) + 2
        k = ii/(nx*ny) + 2
        ovol = one/volref(i, j, k)
        dwd(i, j, k, nt1:nt2) = ovol*turbresscale(1:nturb)*dwd(i, j, k, &
        &       nt1:nt2)
        dwd(i, j, k, 1:nwf) = ovol*dwd(i, j, k, 1:nwf)
    end do
end subroutine resscale_b

!  differentiation of sumdwandfw in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: *dw
!   with respect to varying inputs: *dw *fw
!   rw status of diff variables: *dw:in-out *fw:out
!   plus diff mem management of: dw:in fw:in
subroutine sumdwandfw_b()
    use constants
    ! use blockpointers, only : il, jl, kl, dw, dwd, fw, fwd, iblank
    use flowvarrefstate, only : nwf
    implicit none
    ! local variables
    integer(kind=inttype) :: i, j, k, l
    intrinsic real
    intrinsic max
    integer :: branch
    real(kind=realtype) :: x1
    real(kind=realtype) :: max1
    do l=1,nwf
        do k=2,kl
            do j=2,jl
                do i=2,il
                    x1 = real(iblank(i, j, k), realtype)
                    if (x1 .lt. zero) then
                        call pushreal8(max1)
                        max1 = zero
                        call pushcontrol1b(0)
                    else
                        call pushreal8(max1)
                        max1 = x1
                        call pushcontrol1b(1)
                    end if
                end do
            end do
        end do
    end do
    fwd = 0.0_8
    do l=nwf,1,-1
        do k=kl,2,-1
            do j=jl,2,-1
                do i=il,2,-1
                    fwd(i, j, k, l) = fwd(i, j, k, l) + max1*dwd(i, j, k, l)
                    dwd(i, j, k, l) = max1*dwd(i, j, k, l)
                    call popcontrol1b(branch)
                    if (branch .eq. 0) then
                        call popreal8(max1)
                    else
                        call popreal8(max1)
                    end if
                end do
            end do
        end do
    end do
end subroutine sumdwandfw_b

!  differentiation of viscousflux in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: *aa *w *fw
!   with respect to varying inputs: *rev *aa *wx *wy *wz *w *rlv
!                *qx *qy *qz *ux *uy *uz *vx *vy *vz *fw
!   rw status of diff variables: *rev:out *aa:incr *wx:out *wy:out
!                *wz:out *w:incr *rlv:out *qx:out *qy:out *qz:out
!                *ux:out *uy:out *uz:out *vx:out *vy:out *vz:out
!                *fw:in-out
!   plus diff mem management of: rev:in aa:in wx:in wy:in wz:in
!                w:in rlv:in qx:in qy:in qz:in ux:in uy:in uz:in
!                vx:in vy:in vz:in fw:in
subroutine viscousflux_fast_b()
    !
    !       viscousflux computes the viscous fluxes using a central
    !       difference scheme for a block.
    !       it is assumed that the pointers in block pointer already point
    !       to the correct block.
    !
    use constants
    ! use blockpointers
    use flowvarrefstate
    use inputphysics
    use iteration
    implicit none
    ! possibly correct the wall shear stress.
    ! wall function is not aded
    !
    !      local parameter.
    !
    real(kind=realtype), parameter :: twothird=two*third
    real(kind=realtype), parameter :: xminn=1.e-14_realtype
    !
    !      local variables.
    !
    integer(kind=inttype) :: i, j, k, ii
    real(kind=realtype) :: rfilv, por, mul, mue, mut, heatcoef
    real(kind=realtype) :: muld, mued, mutd, heatcoefd
    real(kind=realtype) :: gm1, factlamheat, factturbheat
    real(kind=realtype) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z
    real(kind=realtype) :: u_xd, u_yd, u_zd, v_xd, v_yd, v_zd, w_xd, &
    &   w_yd, w_zd
    real(kind=realtype) :: q_x, q_y, q_z, ubar, vbar, wbar
    real(kind=realtype) :: q_xd, q_yd, q_zd, ubard, vbard, wbard
    real(kind=realtype) :: corr, ssx, ssy, ssz, ss, fracdiv
    real(kind=realtype) :: corrd, fracdivd
    real(kind=realtype) :: tauxx, tauyy, tauzz
    real(kind=realtype) :: tauxxd, tauyyd, tauzzd
    real(kind=realtype) :: tauxy, tauxz, tauyz
    real(kind=realtype) :: tauxyd, tauxzd, tauyzd
    real(kind=realtype) :: tauxxs, tauyys, tauzzs
    real(kind=realtype) :: tauxxsd, tauyysd, tauzzsd
    real(kind=realtype) :: tauxys, tauxzs, tauyzs
    real(kind=realtype) :: tauxysd, tauxzsd, tauyzsd
    real(kind=realtype) :: exx, eyy, ezz
    real(kind=realtype) :: exxd, eyyd, ezzd
    real(kind=realtype) :: exy, exz, eyz
    real(kind=realtype) :: exyd, exzd, eyzd
    real(kind=realtype) :: wxy, wxz, wyz, wyx, wzx, wzy
    real(kind=realtype) :: wxyd, wxzd, wyzd, wyxd, wzxd, wzyd
    real(kind=realtype) :: den, ccr1, fact
    real(kind=realtype) :: dend, factd
    real(kind=realtype) :: fmx, fmy, fmz, frhoe
    real(kind=realtype) :: fmxd, fmyd, fmzd, frhoed
    logical :: correctfork, storewalltensor
    intrinsic abs
    intrinsic mod
    intrinsic sqrt
    intrinsic max
    integer :: branch
    real(kind=realtype) :: tempd14
    real(kind=realtype) :: tempd13
    real(kind=realtype) :: tempd12
    real(kind=realtype) :: tempd49
    real(kind=realtype) :: tempd11
    real(kind=realtype) :: tempd48
    real(kind=realtype) :: tempd10
    real(kind=realtype) :: tempd47
    real(kind=realtype) :: tempd46
    real(kind=realtype) :: tempd45
    real(kind=realtype) :: tempd44
    real(kind=realtype) :: tempd43
    real(kind=realtype) :: tempd42
    real(kind=realtype) :: tempd41
    real(kind=realtype) :: tempd40
    real(kind=realtype) :: tempd70
    real(kind=realtype) :: tempd39
    real(kind=realtype) :: tempd38
    real(kind=realtype) :: tempd37
    real(kind=realtype) :: tempd36
    real(kind=realtype) :: tempd35
    real(kind=realtype) :: tempd34
    real(kind=realtype) :: tempd33
    real(kind=realtype) :: tempd32
    real(kind=realtype) :: tempd69
    real(kind=realtype) :: tempd31
    real(kind=realtype) :: tempd68
    real(kind=realtype) :: tempd30
    real(kind=realtype) :: tempd67
    real(kind=realtype) :: tempd66
    real(kind=realtype) :: tempd65
    real(kind=realtype) :: tempd64
    real(kind=realtype) :: tempd63
    real(kind=realtype) :: tempd62
    real(kind=realtype) :: tempd61
    real(kind=realtype) :: tempd60
    real(kind=realtype) :: tempd9
    real(kind=realtype) :: tempd
    real(kind=realtype) :: tempd8
    real(kind=realtype) :: tempd7
    real(kind=realtype) :: tempd6
    real(kind=realtype) :: tempd5
    real(kind=realtype) :: tempd4
    real(kind=realtype) :: tempd3
    real(kind=realtype) :: tempd2
    real(kind=realtype) :: tempd1
    real(kind=realtype) :: tempd0
    real(kind=realtype) :: tempd29
    real(kind=realtype) :: tempd28
    real(kind=realtype) :: tempd27
    real(kind=realtype) :: tempd26
    real(kind=realtype) :: tempd25
    real(kind=realtype) :: tempd24
    real(kind=realtype) :: tempd23
    real(kind=realtype) :: tempd22
    real(kind=realtype) :: tempd59
    real(kind=realtype) :: tempd21
    real(kind=realtype) :: tempd58
    real(kind=realtype) :: tempd20
    real(kind=realtype) :: tempd57
    real(kind=realtype) :: tempd56
    real(kind=realtype) :: tempd55
    real(kind=realtype) :: tempd54
    real(kind=realtype) :: tempd53
    real(kind=realtype) :: tempd52
    real(kind=realtype) :: tempd51
    real(kind=realtype) :: tempd50
    real(kind=realtype) :: abs0
    real(kind=realtype) :: tempd19
    real(kind=realtype) :: tempd18
    real(kind=realtype) :: tempd17
    real(kind=realtype) :: tempd16
    real(kind=realtype) :: tempd15
    ! set qcr parameters
    ccr1 = 0.3_realtype
    ! set rfilv to rfil to indicate that this is the viscous part.
    ! if rfilv == 0 the viscous residuals need not to be computed
    ! and a return can be made.
    rfilv = rfil
    if (rfilv .ge. 0.) then
        abs0 = rfilv
    else
        abs0 = -rfilv
    end if
    if (abs0 .lt. thresholdreal) then
        revd = 0.0_8
        wxd = 0.0_8
        wyd = 0.0_8
        wzd = 0.0_8
        rlvd = 0.0_8
        qxd = 0.0_8
        qyd = 0.0_8
        qzd = 0.0_8
        uxd = 0.0_8
        uyd = 0.0_8
        uzd = 0.0_8
        vxd = 0.0_8
        vyd = 0.0_8
        vzd = 0.0_8
    else
        revd = 0.0_8
        wxd = 0.0_8
        wyd = 0.0_8
        wzd = 0.0_8
        rlvd = 0.0_8
        qxd = 0.0_8
        qyd = 0.0_8
        qzd = 0.0_8
        uxd = 0.0_8
        uyd = 0.0_8
        uzd = 0.0_8
        vxd = 0.0_8
        vyd = 0.0_8
        vzd = 0.0_8
        mued = 0.0_8
        mue = zero
        revd = 0.0_8
        wxd = 0.0_8
        wyd = 0.0_8
        wzd = 0.0_8
        rlvd = 0.0_8
        qxd = 0.0_8
        qyd = 0.0_8
        qzd = 0.0_8
        uxd = 0.0_8
        uyd = 0.0_8
        uzd = 0.0_8
        vxd = 0.0_8
        vyd = 0.0_8
        vzd = 0.0_8
        mued = 0.0_8
        ! do k=2, kl
        !    do j=2, jl
        !       do i=1, il
        do ii=0,il*ny*nz-1
            i = mod(ii, il) + 1
            j = mod(ii/il, ny) + 2
            k = ii/(il*ny) + 2
            ! set the value of the porosity. if not zero, it is set
            ! to average the eddy-viscosity and to take the factor
            ! rfilv into account.
            por = half*rfilv
            if (pori(i, j, k) .eq. noflux) por = zero
            ! compute the laminar and (if present) the eddy viscosities
            ! multiplied the porosity. compute the factor in front of
            ! the gradients of the speed of sound squared for the heat
            ! flux.
            mul = por*(rlv(i, j, k)+rlv(i+1, j, k))
            if (eddymodel) then
                mue = por*(rev(i, j, k)+rev(i+1, j, k))
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            mut = mul + mue
            gm1 = half*(gamma(i, j, k)+gamma(i+1, j, k)) - one
            factlamheat = one/(prandtl*gm1)
            factturbheat = one/(prandtlturb*gm1)
            heatcoef = mul*factlamheat + mue*factturbheat
            ! compute the gradients at the face by averaging the four
            ! nodal values.
            u_x = fourth*(ux(i, j-1, k-1)+ux(i, j, k-1)+ux(i, j-1, k)+ux(i, &
            &         j, k))
            u_y = fourth*(uy(i, j-1, k-1)+uy(i, j, k-1)+uy(i, j-1, k)+uy(i, &
            &         j, k))
            u_z = fourth*(uz(i, j-1, k-1)+uz(i, j, k-1)+uz(i, j-1, k)+uz(i, &
            &         j, k))
            v_x = fourth*(vx(i, j-1, k-1)+vx(i, j, k-1)+vx(i, j-1, k)+vx(i, &
            &         j, k))
            v_y = fourth*(vy(i, j-1, k-1)+vy(i, j, k-1)+vy(i, j-1, k)+vy(i, &
            &         j, k))
            v_z = fourth*(vz(i, j-1, k-1)+vz(i, j, k-1)+vz(i, j-1, k)+vz(i, &
            &         j, k))
            w_x = fourth*(wx(i, j-1, k-1)+wx(i, j, k-1)+wx(i, j-1, k)+wx(i, &
            &         j, k))
            w_y = fourth*(wy(i, j-1, k-1)+wy(i, j, k-1)+wy(i, j-1, k)+wy(i, &
            &         j, k))
            w_z = fourth*(wz(i, j-1, k-1)+wz(i, j, k-1)+wz(i, j-1, k)+wz(i, &
            &         j, k))
            q_x = fourth*(qx(i, j-1, k-1)+qx(i, j, k-1)+qx(i, j-1, k)+qx(i, &
            &         j, k))
            q_y = fourth*(qy(i, j-1, k-1)+qy(i, j, k-1)+qy(i, j-1, k)+qy(i, &
            &         j, k))
            q_z = fourth*(qz(i, j-1, k-1)+qz(i, j, k-1)+qz(i, j-1, k)+qz(i, &
            &         j, k))
            ! the gradients in the normal direction are corrected, such
            ! that no averaging takes places here.
            ! first determine the vector in the direction from the
            ! cell center i to cell center i+1.
            ssx = eighth*(x(i+1, j-1, k-1, 1)-x(i-1, j-1, k-1, 1)+x(i+1, j-1&
            &         , k, 1)-x(i-1, j-1, k, 1)+x(i+1, j, k-1, 1)-x(i-1, j, k-1, 1)+&
            &         x(i+1, j, k, 1)-x(i-1, j, k, 1))
            ssy = eighth*(x(i+1, j-1, k-1, 2)-x(i-1, j-1, k-1, 2)+x(i+1, j-1&
            &         , k, 2)-x(i-1, j-1, k, 2)+x(i+1, j, k-1, 2)-x(i-1, j, k-1, 2)+&
            &         x(i+1, j, k, 2)-x(i-1, j, k, 2))
            ssz = eighth*(x(i+1, j-1, k-1, 3)-x(i-1, j-1, k-1, 3)+x(i+1, j-1&
            &         , k, 3)-x(i-1, j-1, k, 3)+x(i+1, j, k-1, 3)-x(i-1, j, k-1, 3)+&
            &         x(i+1, j, k, 3)-x(i-1, j, k, 3))
            ! determine the length of this vector and create the
            ! unit normal.
            ss = one/sqrt(ssx*ssx+ssy*ssy+ssz*ssz)
            ssx = ss*ssx
            ssy = ss*ssy
            ssz = ss*ssz
            ! correct the gradients.
            corr = u_x*ssx + u_y*ssy + u_z*ssz - (w(i+1, j, k, ivx)-w(i, j, &
            &         k, ivx))*ss
            u_x = u_x - corr*ssx
            u_y = u_y - corr*ssy
            u_z = u_z - corr*ssz
            corr = v_x*ssx + v_y*ssy + v_z*ssz - (w(i+1, j, k, ivy)-w(i, j, &
            &         k, ivy))*ss
            v_x = v_x - corr*ssx
            v_y = v_y - corr*ssy
            v_z = v_z - corr*ssz
            corr = w_x*ssx + w_y*ssy + w_z*ssz - (w(i+1, j, k, ivz)-w(i, j, &
            &         k, ivz))*ss
            w_x = w_x - corr*ssx
            w_y = w_y - corr*ssy
            w_z = w_z - corr*ssz
            corr = q_x*ssx + q_y*ssy + q_z*ssz + (aa(i+1, j, k)-aa(i, j, k))&
            &         *ss
            q_x = q_x - corr*ssx
            q_y = q_y - corr*ssy
            q_z = q_z - corr*ssz
            ! compute the stress tensor and the heat flux vector.
            ! we remove the viscosity from the stress tensor (tau)
            ! to define taus since we still need to separate between
            ! laminar and turbulent stress for qcr.
            ! therefore, laminar tau = mue*taus, turbulent
            ! tau = mue*taus, and total tau = mut*taus.
            fracdiv = twothird*(u_x+v_y+w_z)
            tauxxs = two*u_x - fracdiv
            tauyys = two*v_y - fracdiv
            tauzzs = two*w_z - fracdiv
            tauxys = u_y + v_x
            tauxzs = u_z + w_x
            tauyzs = v_z + w_y
            ! add qcr corrections if necessary
            if (useqcr) then
                ! in the qcr formulation, we add an extra term to the turbulent stress tensor:
                !
                ! tau_ij,qcr = tau_ij - e_ij
                !
                ! where, according to tmr website (http://turbmodels.larc.nasa.gov/spalart.html):
                !
                ! e_ij = ccr1*(o_ik*tau_jk + o_jk*tau_ik)
                !
                ! we are computing o_ik as follows:
                !
                ! o_ik = 2*w_ik/den
                !
                ! remember that the tau_ij in e_ij should use only the eddy viscosity!
                ! compute denominator
                den = sqrt(u_x*u_x + u_y*u_y + u_z*u_z + v_x*v_x + v_y*v_y + &
                &           v_z*v_z + w_x*w_x + w_y*w_y + w_z*w_z)
                if (den .lt. xminn) then
                    den = xminn
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                    den = den
                end if
                ! compute factor that will multiply all tensor components.
                ! here we add the eddy viscosity that should multiply the stress tensor (tau)
                ! components as well.
                fact = mue*ccr1/den
                ! compute off-diagonal terms of vorticity tensor (we will ommit the 1/2)
                ! the diagonals of the vorticity tensor components are always zero
                wxy = u_y - v_x
                wxz = u_z - w_x
                wyz = v_z - w_y
                wyx = -wxy
                wzx = -wxz
                wzy = -wyz
                ! compute the extra terms of the boussinesq relation
                exx = fact*(wxy*tauxys+wxz*tauxzs)*two
                eyy = fact*(wyx*tauxys+wyz*tauyzs)*two
                ezz = fact*(wzx*tauxzs+wzy*tauyzs)*two
                exy = fact*(wxy*tauyys+wxz*tauyzs+wyx*tauxxs+wyz*tauxzs)
                exz = fact*(wxy*tauyzs+wxz*tauzzs+wzx*tauxxs+wzy*tauxys)
                ! apply the total viscosity to the stress tensor and add extra terms
                eyz = fact*(wyx*tauxzs+wyz*tauzzs+wzx*tauxys+wzy*tauyys)
                tauxx = mut*tauxxs - exx
                tauyy = mut*tauyys - eyy
                tauzz = mut*tauzzs - ezz
                tauxy = mut*tauxys - exy
                tauxz = mut*tauxzs - exz
                tauyz = mut*tauyzs - eyz
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                ! just apply the total viscosity to the stress tensor
                tauxx = mut*tauxxs
                tauyy = mut*tauyys
                tauzz = mut*tauzzs
                tauxy = mut*tauxys
                tauxz = mut*tauxzs
                tauyz = mut*tauyzs
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            ! compute the average velocities for the face. remember that
            ! the velocities are stored and not the momentum.
            ubar = half*(w(i, j, k, ivx)+w(i+1, j, k, ivx))
            vbar = half*(w(i, j, k, ivy)+w(i+1, j, k, ivy))
            wbar = half*(w(i, j, k, ivz)+w(i+1, j, k, ivz))
            ! compute the viscous fluxes for this i-face.
            ! update the residuals of cell i and i+1.
            ! store the stress tensor and the heat flux vector if this
            ! face is part of a viscous subface. both the cases i == 1
            ! and i == il must be tested.
            frhoed = fwd(i+1, j, k, irhoe) - fwd(i, j, k, irhoe)
            fmzd = fwd(i+1, j, k, imz) - fwd(i, j, k, imz)
            fmyd = fwd(i+1, j, k, imy) - fwd(i, j, k, imy)
            fmxd = fwd(i+1, j, k, imx) - fwd(i, j, k, imx)
            tempd68 = si(i, j, k, 1)*frhoed
            tempd69 = si(i, j, k, 2)*frhoed
            tempd70 = si(i, j, k, 3)*frhoed
            ubard = tauxz*tempd70 + tauxy*tempd69 + tauxx*tempd68
            tauxxd = si(i, j, k, 1)*fmxd + ubar*tempd68
            vbard = tauyz*tempd70 + tauyy*tempd69 + tauxy*tempd68
            tauxyd = si(i, j, k, 1)*fmyd + si(i, j, k, 2)*fmxd + ubar*&
            &         tempd69 + vbar*tempd68
            wbard = tauzz*tempd70 + tauyz*tempd69 + tauxz*tempd68
            tauxzd = si(i, j, k, 1)*fmzd + si(i, j, k, 3)*fmxd + ubar*&
            &         tempd70 + wbar*tempd68
            tauyyd = si(i, j, k, 2)*fmyd + vbar*tempd69
            tauyzd = si(i, j, k, 2)*fmzd + si(i, j, k, 3)*fmyd + vbar*&
            &         tempd70 + wbar*tempd69
            tauzzd = si(i, j, k, 3)*fmzd + wbar*tempd70
            q_xd = -(si(i, j, k, 1)*frhoed)
            q_yd = -(si(i, j, k, 2)*frhoed)
            q_zd = -(si(i, j, k, 3)*frhoed)
            wd(i, j, k, ivz) = wd(i, j, k, ivz) + half*wbard
            wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) + half*wbard
            wd(i, j, k, ivy) = wd(i, j, k, ivy) + half*vbard
            wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) + half*vbard
            wd(i, j, k, ivx) = wd(i, j, k, ivx) + half*ubard
            wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) + half*ubard
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                exzd = -tauxzd
                exyd = -tauxyd
                ezzd = -tauzzd
                eyyd = -tauyyd
                tempd61 = fact*exzd
                tempd64 = fact*exyd
                tempd62 = two*fact*ezzd
                tempd63 = two*fact*eyyd
                mutd = tauxzs*tauxzd + tauzzs*tauzzd + tauxxs*tauxxd + tauyys*&
                &           tauyyd + tauxys*tauxyd + tauyzs*tauyzd
                tauyzsd = wxy*tempd61 + wzy*tempd62 + wyz*tempd63 + wxz*&
                &           tempd64 + mut*tauyzd
                eyzd = -tauyzd
                tauxxsd = wzx*tempd61 + wyx*tempd64 + mut*tauxxd
                exxd = -tauxxd
                tempd65 = fact*eyzd
                tauzzsd = wyz*tempd65 + wxz*tempd61 + mut*tauzzd
                tauyysd = wzy*tempd65 + wxy*tempd64 + mut*tauyyd
                factd = (wxy*tauyzs+wxz*tauzzs+wzx*tauxxs+wzy*tauxys)*exzd + &
                &           two*(wzx*tauxzs+wzy*tauyzs)*ezzd + two*(wxy*tauxys+wxz*&
                &           tauxzs)*exxd + two*(wyx*tauxys+wyz*tauyzs)*eyyd + (wxy*&
                &           tauyys+wxz*tauyzs+wyx*tauxxs+wyz*tauxzs)*exyd + (wyx*tauxzs+&
                &           wyz*tauzzs+wzx*tauxys+wzy*tauyys)*eyzd
                wyxd = tauxxs*tempd64 + tauxys*tempd63 + tauxzs*tempd65
                wzxd = tauxxs*tempd61 + tauxzs*tempd62 + tauxys*tempd65
                wzyd = tauxys*tempd61 + tauyzs*tempd62 + tauyys*tempd65
                wyzd = tauxzs*tempd64 - wzyd + tauyzs*tempd63 + tauzzs*tempd65
                tempd66 = two*fact*exxd
                tauxzsd = wyx*tempd65 + wzx*tempd62 + wxz*tempd66 + wyz*&
                &           tempd64 + mut*tauxzd
                tauxysd = wzx*tempd65 + wyx*tempd63 + wxy*tempd66 + wzy*&
                &           tempd61 + mut*tauxyd
                wxyd = tauyys*tempd64 - wyxd + tauxys*tempd66 + tauyzs*tempd61
                wxzd = tauyzs*tempd64 - wzxd + tauxzs*tempd66 + tauzzs*tempd61
                v_zd = wyzd
                w_yd = -wyzd
                u_zd = wxzd
                w_xd = -wxzd
                u_yd = wxyd
                v_xd = -wxyd
                tempd67 = ccr1*factd/den
                mued = mued + tempd67
                dend = -(mue*tempd67/den)
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) dend = 0.0_8
                if (u_x**2 + u_y**2 + u_z**2 + v_x**2 + v_y**2 + v_z**2 + w_x&
                &             **2 + w_y**2 + w_z**2 .eq. 0.0_8) then
                tempd60 = 0.0
            else
                tempd60 = dend/(2.0*sqrt(u_x**2+u_y**2+u_z**2+v_x**2+v_y**2+&
                &             v_z**2+w_x**2+w_y**2+w_z**2))
            end if
            u_xd = 2*u_x*tempd60
            u_yd = u_yd + 2*u_y*tempd60
            u_zd = u_zd + 2*u_z*tempd60
            v_xd = v_xd + 2*v_x*tempd60
            v_yd = 2*v_y*tempd60
            v_zd = v_zd + 2*v_z*tempd60
            w_xd = w_xd + 2*w_x*tempd60
            w_yd = w_yd + 2*w_y*tempd60
            w_zd = 2*w_z*tempd60
        else
            mutd = tauxzs*tauxzd + tauzzs*tauzzd + tauxxs*tauxxd + tauyys*&
            &           tauyyd + tauxys*tauxyd + tauyzs*tauyzd
            tauyzsd = mut*tauyzd
            tauxzsd = mut*tauxzd
            tauxysd = mut*tauxyd
            tauzzsd = mut*tauzzd
            tauyysd = mut*tauyyd
            tauxxsd = mut*tauxxd
            u_xd = 0.0_8
            u_yd = 0.0_8
            u_zd = 0.0_8
            w_xd = 0.0_8
            w_yd = 0.0_8
            w_zd = 0.0_8
            v_xd = 0.0_8
            v_yd = 0.0_8
            v_zd = 0.0_8
        end if
        fracdivd = -tauyysd - tauxxsd - tauzzsd
        tempd47 = twothird*fracdivd
        heatcoefd = q_y*q_yd + q_x*q_xd + q_z*q_zd
        q_zd = heatcoef*q_zd
        q_yd = heatcoef*q_yd
        q_xd = heatcoef*q_xd
        v_zd = v_zd + tauyzsd
        w_yd = w_yd + tauyzsd
        u_zd = u_zd + tauxzsd
        w_xd = w_xd + tauxzsd
        u_yd = u_yd + tauxysd
        v_xd = v_xd + tauxysd
        w_zd = w_zd + tempd47 + two*tauzzsd
        v_yd = v_yd + tempd47 + two*tauyysd
        u_xd = u_xd + tempd47 + two*tauxxsd
        corrd = -(ssy*q_yd) - ssx*q_xd - ssz*q_zd
        q_xd = q_xd + ssx*corrd
        q_yd = q_yd + ssy*corrd
        q_zd = q_zd + ssz*corrd
        aad(i+1, j, k) = aad(i+1, j, k) + ss*corrd
        aad(i, j, k) = aad(i, j, k) - ss*corrd
        corrd = -(ssy*w_yd) - ssx*w_xd - ssz*w_zd
        w_xd = w_xd + ssx*corrd
        w_yd = w_yd + ssy*corrd
        w_zd = w_zd + ssz*corrd
        wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) - ss*corrd
        wd(i, j, k, ivz) = wd(i, j, k, ivz) + ss*corrd
        corrd = -(ssy*v_yd) - ssx*v_xd - ssz*v_zd
        v_xd = v_xd + ssx*corrd
        v_yd = v_yd + ssy*corrd
        v_zd = v_zd + ssz*corrd
        wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) - ss*corrd
        wd(i, j, k, ivy) = wd(i, j, k, ivy) + ss*corrd
        corrd = -(ssy*u_yd) - ssx*u_xd - ssz*u_zd
        u_xd = u_xd + ssx*corrd
        u_yd = u_yd + ssy*corrd
        u_zd = u_zd + ssz*corrd
        wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) - ss*corrd
        wd(i, j, k, ivx) = wd(i, j, k, ivx) + ss*corrd
        tempd48 = fourth*q_zd
        qzd(i, j-1, k-1) = qzd(i, j-1, k-1) + tempd48
        qzd(i, j, k-1) = qzd(i, j, k-1) + tempd48
        qzd(i, j-1, k) = qzd(i, j-1, k) + tempd48
        qzd(i, j, k) = qzd(i, j, k) + tempd48
        tempd49 = fourth*q_yd
        qyd(i, j-1, k-1) = qyd(i, j-1, k-1) + tempd49
        qyd(i, j, k-1) = qyd(i, j, k-1) + tempd49
        qyd(i, j-1, k) = qyd(i, j-1, k) + tempd49
        qyd(i, j, k) = qyd(i, j, k) + tempd49
        tempd50 = fourth*q_xd
        qxd(i, j-1, k-1) = qxd(i, j-1, k-1) + tempd50
        qxd(i, j, k-1) = qxd(i, j, k-1) + tempd50
        qxd(i, j-1, k) = qxd(i, j-1, k) + tempd50
        qxd(i, j, k) = qxd(i, j, k) + tempd50
        tempd51 = fourth*w_zd
        wzd(i, j-1, k-1) = wzd(i, j-1, k-1) + tempd51
        wzd(i, j, k-1) = wzd(i, j, k-1) + tempd51
        wzd(i, j-1, k) = wzd(i, j-1, k) + tempd51
        wzd(i, j, k) = wzd(i, j, k) + tempd51
        tempd52 = fourth*w_yd
        wyd(i, j-1, k-1) = wyd(i, j-1, k-1) + tempd52
        wyd(i, j, k-1) = wyd(i, j, k-1) + tempd52
        wyd(i, j-1, k) = wyd(i, j-1, k) + tempd52
        wyd(i, j, k) = wyd(i, j, k) + tempd52
        tempd53 = fourth*w_xd
        wxd(i, j-1, k-1) = wxd(i, j-1, k-1) + tempd53
        wxd(i, j, k-1) = wxd(i, j, k-1) + tempd53
        wxd(i, j-1, k) = wxd(i, j-1, k) + tempd53
        wxd(i, j, k) = wxd(i, j, k) + tempd53
        tempd54 = fourth*v_zd
        vzd(i, j-1, k-1) = vzd(i, j-1, k-1) + tempd54
        vzd(i, j, k-1) = vzd(i, j, k-1) + tempd54
        vzd(i, j-1, k) = vzd(i, j-1, k) + tempd54
        vzd(i, j, k) = vzd(i, j, k) + tempd54
        tempd55 = fourth*v_yd
        vyd(i, j-1, k-1) = vyd(i, j-1, k-1) + tempd55
        vyd(i, j, k-1) = vyd(i, j, k-1) + tempd55
        vyd(i, j-1, k) = vyd(i, j-1, k) + tempd55
        vyd(i, j, k) = vyd(i, j, k) + tempd55
        tempd56 = fourth*v_xd
        vxd(i, j-1, k-1) = vxd(i, j-1, k-1) + tempd56
        vxd(i, j, k-1) = vxd(i, j, k-1) + tempd56
        vxd(i, j-1, k) = vxd(i, j-1, k) + tempd56
        vxd(i, j, k) = vxd(i, j, k) + tempd56
        tempd57 = fourth*u_zd
        uzd(i, j-1, k-1) = uzd(i, j-1, k-1) + tempd57
        uzd(i, j, k-1) = uzd(i, j, k-1) + tempd57
        uzd(i, j-1, k) = uzd(i, j-1, k) + tempd57
        uzd(i, j, k) = uzd(i, j, k) + tempd57
        tempd58 = fourth*u_yd
        uyd(i, j-1, k-1) = uyd(i, j-1, k-1) + tempd58
        uyd(i, j, k-1) = uyd(i, j, k-1) + tempd58
        uyd(i, j-1, k) = uyd(i, j-1, k) + tempd58
        uyd(i, j, k) = uyd(i, j, k) + tempd58
        tempd59 = fourth*u_xd
        uxd(i, j-1, k-1) = uxd(i, j-1, k-1) + tempd59
        uxd(i, j, k-1) = uxd(i, j, k-1) + tempd59
        uxd(i, j-1, k) = uxd(i, j-1, k) + tempd59
        uxd(i, j, k) = uxd(i, j, k) + tempd59
        muld = mutd + factlamheat*heatcoefd
        mued = mued + mutd + factturbheat*heatcoefd
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) then
            revd(i, j, k) = revd(i, j, k) + por*mued
            revd(i+1, j, k) = revd(i+1, j, k) + por*mued
            mued = 0.0_8
        end if
        rlvd(i, j, k) = rlvd(i, j, k) + por*muld
        rlvd(i+1, j, k) = rlvd(i+1, j, k) + por*muld
        !     end do
        ! end do
    end do
    mued = 0.0_8
    mue = zero
    mued = 0.0_8
    ! do k=2,kl
    !    do j=1,jl
    !       do i=2,il
    do ii=0,nx*jl*nz-1
        i = mod(ii, nx) + 2
        j = mod(ii/nx, jl) + 1
        k = ii/(nx*jl) + 2
        ! set the value of the porosity. if not zero, it is set
        ! to average the eddy-viscosity and to take the factor
        ! rfilv into account.
        por = half*rfilv
        if (porj(i, j, k) .eq. noflux) por = zero
        ! compute the laminar and (if present) the eddy viscosities
        ! multiplied by the porosity. compute the factor in front of
        ! the gradients of the speed of sound squared for the heat
        ! flux.
        mul = por*(rlv(i, j, k)+rlv(i, j+1, k))
        if (eddymodel) then
            mue = por*(rev(i, j, k)+rev(i, j+1, k))
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 0
        else
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 1
        end if
        mut = mul + mue
        gm1 = half*(gamma(i, j, k)+gamma(i, j+1, k)) - one
        factlamheat = one/(prandtl*gm1)
        factturbheat = one/(prandtlturb*gm1)
        heatcoef = mul*factlamheat + mue*factturbheat
        ! compute the gradients at the face by averaging the four
        ! nodal values.
        u_x = fourth*(ux(i-1, j, k-1)+ux(i, j, k-1)+ux(i-1, j, k)+ux(i, &
        &         j, k))
        u_y = fourth*(uy(i-1, j, k-1)+uy(i, j, k-1)+uy(i-1, j, k)+uy(i, &
        &         j, k))
        u_z = fourth*(uz(i-1, j, k-1)+uz(i, j, k-1)+uz(i-1, j, k)+uz(i, &
        &         j, k))
        v_x = fourth*(vx(i-1, j, k-1)+vx(i, j, k-1)+vx(i-1, j, k)+vx(i, &
        &         j, k))
        v_y = fourth*(vy(i-1, j, k-1)+vy(i, j, k-1)+vy(i-1, j, k)+vy(i, &
        &         j, k))
        v_z = fourth*(vz(i-1, j, k-1)+vz(i, j, k-1)+vz(i-1, j, k)+vz(i, &
        &         j, k))
        w_x = fourth*(wx(i-1, j, k-1)+wx(i, j, k-1)+wx(i-1, j, k)+wx(i, &
        &         j, k))
        w_y = fourth*(wy(i-1, j, k-1)+wy(i, j, k-1)+wy(i-1, j, k)+wy(i, &
        &         j, k))
        w_z = fourth*(wz(i-1, j, k-1)+wz(i, j, k-1)+wz(i-1, j, k)+wz(i, &
        &         j, k))
        q_x = fourth*(qx(i-1, j, k-1)+qx(i, j, k-1)+qx(i-1, j, k)+qx(i, &
        &         j, k))
        q_y = fourth*(qy(i-1, j, k-1)+qy(i, j, k-1)+qy(i-1, j, k)+qy(i, &
        &         j, k))
        q_z = fourth*(qz(i-1, j, k-1)+qz(i, j, k-1)+qz(i-1, j, k)+qz(i, &
        &         j, k))
        ! the gradients in the normal direction are corrected, such
        ! that no averaging takes places here.
        ! first determine the vector in the direction from the
        ! cell center j to cell center j+1.
        ssx = eighth*(x(i-1, j+1, k-1, 1)-x(i-1, j-1, k-1, 1)+x(i-1, j+1&
        &         , k, 1)-x(i-1, j-1, k, 1)+x(i, j+1, k-1, 1)-x(i, j-1, k-1, 1)+&
        &         x(i, j+1, k, 1)-x(i, j-1, k, 1))
        ssy = eighth*(x(i-1, j+1, k-1, 2)-x(i-1, j-1, k-1, 2)+x(i-1, j+1&
        &         , k, 2)-x(i-1, j-1, k, 2)+x(i, j+1, k-1, 2)-x(i, j-1, k-1, 2)+&
        &         x(i, j+1, k, 2)-x(i, j-1, k, 2))
        ssz = eighth*(x(i-1, j+1, k-1, 3)-x(i-1, j-1, k-1, 3)+x(i-1, j+1&
        &         , k, 3)-x(i-1, j-1, k, 3)+x(i, j+1, k-1, 3)-x(i, j-1, k-1, 3)+&
        &         x(i, j+1, k, 3)-x(i, j-1, k, 3))
        ! determine the length of this vector and create the
        ! unit normal.
        ss = one/sqrt(ssx*ssx+ssy*ssy+ssz*ssz)
        ssx = ss*ssx
        ssy = ss*ssy
        ssz = ss*ssz
        ! correct the gradients.
        corr = u_x*ssx + u_y*ssy + u_z*ssz - (w(i, j+1, k, ivx)-w(i, j, &
        &         k, ivx))*ss
        u_x = u_x - corr*ssx
        u_y = u_y - corr*ssy
        u_z = u_z - corr*ssz
        corr = v_x*ssx + v_y*ssy + v_z*ssz - (w(i, j+1, k, ivy)-w(i, j, &
        &         k, ivy))*ss
        v_x = v_x - corr*ssx
        v_y = v_y - corr*ssy
        v_z = v_z - corr*ssz
        corr = w_x*ssx + w_y*ssy + w_z*ssz - (w(i, j+1, k, ivz)-w(i, j, &
        &         k, ivz))*ss
        w_x = w_x - corr*ssx
        w_y = w_y - corr*ssy
        w_z = w_z - corr*ssz
        corr = q_x*ssx + q_y*ssy + q_z*ssz + (aa(i, j+1, k)-aa(i, j, k))&
        &         *ss
        q_x = q_x - corr*ssx
        q_y = q_y - corr*ssy
        q_z = q_z - corr*ssz
        ! compute the stress tensor and the heat flux vector.
        ! we remove the viscosity from the stress tensor (tau)
        ! to define taus since we still need to separate between
        ! laminar and turbulent stress for qcr.
        ! therefore, laminar tau = mue*taus, turbulent
        ! tau = mue*taus, and total tau = mut*taus.
        fracdiv = twothird*(u_x+v_y+w_z)
        tauxxs = two*u_x - fracdiv
        tauyys = two*v_y - fracdiv
        tauzzs = two*w_z - fracdiv
        tauxys = u_y + v_x
        tauxzs = u_z + w_x
        tauyzs = v_z + w_y
        ! add qcr corrections if necessary
        if (useqcr) then
            ! in the qcr formulation, we add an extra term to the turbulent stress tensor:
            !
            ! tau_ij,qcr = tau_ij - e_ij
            !
            ! where, according to tmr website (http://turbmodels.larc.nasa.gov/spalart.html):
            !
            ! e_ij = ccr1*(o_ik*tau_jk + o_jk*tau_ik)
            !
            ! we are computing o_ik as follows:
            !
            ! o_ik = 2*w_ik/den
            !
            ! remember that the tau_ij in e_ij should use only the eddy viscosity!
            ! compute denominator
            den = sqrt(u_x*u_x + u_y*u_y + u_z*u_z + v_x*v_x + v_y*v_y + &
            &           v_z*v_z + w_x*w_x + w_y*w_y + w_z*w_z)
            if (den .lt. xminn) then
                den = xminn
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
                den = den
            end if
            ! compute factor that will multiply all tensor components.
            ! here we add the eddy viscosity that should multiply the stress tensor (tau)
            ! components as well.
            fact = mue*ccr1/den
            ! compute off-diagonal terms of vorticity tensor (we will ommit the 1/2)
            ! the diagonals of the vorticity tensor components are always zero
            wxy = u_y - v_x
            wxz = u_z - w_x
            wyz = v_z - w_y
            wyx = -wxy
            wzx = -wxz
            wzy = -wyz
            ! compute the extra terms of the boussinesq relation
            exx = fact*(wxy*tauxys+wxz*tauxzs)*two
            eyy = fact*(wyx*tauxys+wyz*tauyzs)*two
            ezz = fact*(wzx*tauxzs+wzy*tauyzs)*two
            exy = fact*(wxy*tauyys+wxz*tauyzs+wyx*tauxxs+wyz*tauxzs)
            exz = fact*(wxy*tauyzs+wxz*tauzzs+wzx*tauxxs+wzy*tauxys)
            eyz = fact*(wyx*tauxzs+wyz*tauzzs+wzx*tauxys+wzy*tauyys)
            ! apply the total viscosity to the stress tensor and add extra terms
            tauxx = mut*tauxxs - exx
            tauyy = mut*tauyys - eyy
            tauzz = mut*tauzzs - ezz
            tauxy = mut*tauxys - exy
            tauxz = mut*tauxzs - exz
            tauyz = mut*tauyzs - eyz
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 0
        else
            ! just apply the total viscosity to the stress tensor
            tauxx = mut*tauxxs
            tauyy = mut*tauyys
            tauzz = mut*tauzzs
            tauxy = mut*tauxys
            tauxz = mut*tauxzs
            tauyz = mut*tauyzs
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 1
        end if
        ! compute the average velocities for the face. remember that
        ! the velocities are stored and not the momentum.
        ubar = half*(w(i, j, k, ivx)+w(i, j+1, k, ivx))
        vbar = half*(w(i, j, k, ivy)+w(i, j+1, k, ivy))
        wbar = half*(w(i, j, k, ivz)+w(i, j+1, k, ivz))
        ! compute the viscous fluxes for this j-face.
        ! update the residuals of cell j and j+1.
        ! store the stress tensor and the heat flux vector if this
        ! face is part of a viscous subface. both the cases j == 1
        ! and j == jl must be tested.
        frhoed = fwd(i, j+1, k, irhoe) - fwd(i, j, k, irhoe)
        fmzd = fwd(i, j+1, k, imz) - fwd(i, j, k, imz)
        fmyd = fwd(i, j+1, k, imy) - fwd(i, j, k, imy)
        fmxd = fwd(i, j+1, k, imx) - fwd(i, j, k, imx)
        tempd44 = sj(i, j, k, 1)*frhoed
        tempd45 = sj(i, j, k, 2)*frhoed
        tempd46 = sj(i, j, k, 3)*frhoed
        ubard = tauxz*tempd46 + tauxy*tempd45 + tauxx*tempd44
        tauxxd = sj(i, j, k, 1)*fmxd + ubar*tempd44
        vbard = tauyz*tempd46 + tauyy*tempd45 + tauxy*tempd44
        tauxyd = sj(i, j, k, 1)*fmyd + sj(i, j, k, 2)*fmxd + ubar*&
        &         tempd45 + vbar*tempd44
        wbard = tauzz*tempd46 + tauyz*tempd45 + tauxz*tempd44
        tauxzd = sj(i, j, k, 1)*fmzd + sj(i, j, k, 3)*fmxd + ubar*&
        &         tempd46 + wbar*tempd44
        tauyyd = sj(i, j, k, 2)*fmyd + vbar*tempd45
        tauyzd = sj(i, j, k, 2)*fmzd + sj(i, j, k, 3)*fmyd + vbar*&
        &         tempd46 + wbar*tempd45
        tauzzd = sj(i, j, k, 3)*fmzd + wbar*tempd46
        q_xd = -(sj(i, j, k, 1)*frhoed)
        q_yd = -(sj(i, j, k, 2)*frhoed)
        q_zd = -(sj(i, j, k, 3)*frhoed)
        wd(i, j, k, ivz) = wd(i, j, k, ivz) + half*wbard
        wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) + half*wbard
        wd(i, j, k, ivy) = wd(i, j, k, ivy) + half*vbard
        wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) + half*vbard
        wd(i, j, k, ivx) = wd(i, j, k, ivx) + half*ubard
        wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) + half*ubard
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) then
            exzd = -tauxzd
            exyd = -tauxyd
            ezzd = -tauzzd
            eyyd = -tauyyd
            tempd37 = fact*exzd
            tempd40 = fact*exyd
            tempd38 = two*fact*ezzd
            tempd39 = two*fact*eyyd
            mutd = tauxzs*tauxzd + tauzzs*tauzzd + tauxxs*tauxxd + tauyys*&
            &           tauyyd + tauxys*tauxyd + tauyzs*tauyzd
            tauyzsd = wxy*tempd37 + wzy*tempd38 + wyz*tempd39 + wxz*&
            &           tempd40 + mut*tauyzd
            eyzd = -tauyzd
            tauxxsd = wzx*tempd37 + wyx*tempd40 + mut*tauxxd
            exxd = -tauxxd
            tempd41 = fact*eyzd
            tauzzsd = wyz*tempd41 + wxz*tempd37 + mut*tauzzd
            tauyysd = wzy*tempd41 + wxy*tempd40 + mut*tauyyd
            factd = (wxy*tauyzs+wxz*tauzzs+wzx*tauxxs+wzy*tauxys)*exzd + &
            &           two*(wzx*tauxzs+wzy*tauyzs)*ezzd + two*(wxy*tauxys+wxz*&
            &           tauxzs)*exxd + two*(wyx*tauxys+wyz*tauyzs)*eyyd + (wxy*&
            &           tauyys+wxz*tauyzs+wyx*tauxxs+wyz*tauxzs)*exyd + (wyx*tauxzs+&
            &           wyz*tauzzs+wzx*tauxys+wzy*tauyys)*eyzd
            wyxd = tauxxs*tempd40 + tauxys*tempd39 + tauxzs*tempd41
            wzxd = tauxxs*tempd37 + tauxzs*tempd38 + tauxys*tempd41
            wzyd = tauxys*tempd37 + tauyzs*tempd38 + tauyys*tempd41
            wyzd = tauxzs*tempd40 - wzyd + tauyzs*tempd39 + tauzzs*tempd41
            tempd42 = two*fact*exxd
            tauxzsd = wyx*tempd41 + wzx*tempd38 + wxz*tempd42 + wyz*&
            &           tempd40 + mut*tauxzd
            tauxysd = wzx*tempd41 + wyx*tempd39 + wxy*tempd42 + wzy*&
            &           tempd37 + mut*tauxyd
            wxyd = tauyys*tempd40 - wyxd + tauxys*tempd42 + tauyzs*tempd37
            wxzd = tauyzs*tempd40 - wzxd + tauxzs*tempd42 + tauzzs*tempd37
            v_zd = wyzd
            w_yd = -wyzd
            u_zd = wxzd
            w_xd = -wxzd
            u_yd = wxyd
            v_xd = -wxyd
            tempd43 = ccr1*factd/den
            mued = mued + tempd43
            dend = -(mue*tempd43/den)
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) dend = 0.0_8
            if (u_x**2 + u_y**2 + u_z**2 + v_x**2 + v_y**2 + v_z**2 + w_x&
            &             **2 + w_y**2 + w_z**2 .eq. 0.0_8) then
            tempd36 = 0.0
        else
            tempd36 = dend/(2.0*sqrt(u_x**2+u_y**2+u_z**2+v_x**2+v_y**2+&
            &             v_z**2+w_x**2+w_y**2+w_z**2))
        end if
        u_xd = 2*u_x*tempd36
        u_yd = u_yd + 2*u_y*tempd36
        u_zd = u_zd + 2*u_z*tempd36
        v_xd = v_xd + 2*v_x*tempd36
        v_yd = 2*v_y*tempd36
        v_zd = v_zd + 2*v_z*tempd36
        w_xd = w_xd + 2*w_x*tempd36
        w_yd = w_yd + 2*w_y*tempd36
        w_zd = 2*w_z*tempd36
    else
        mutd = tauxzs*tauxzd + tauzzs*tauzzd + tauxxs*tauxxd + tauyys*&
        &           tauyyd + tauxys*tauxyd + tauyzs*tauyzd
        tauyzsd = mut*tauyzd
        tauxzsd = mut*tauxzd
        tauxysd = mut*tauxyd
        tauzzsd = mut*tauzzd
        tauyysd = mut*tauyyd
        tauxxsd = mut*tauxxd
        u_xd = 0.0_8
        u_yd = 0.0_8
        u_zd = 0.0_8
        w_xd = 0.0_8
        w_yd = 0.0_8
        w_zd = 0.0_8
        v_xd = 0.0_8
        v_yd = 0.0_8
        v_zd = 0.0_8
    end if
    fracdivd = -tauyysd - tauxxsd - tauzzsd
    tempd23 = twothird*fracdivd
    heatcoefd = q_y*q_yd + q_x*q_xd + q_z*q_zd
    q_zd = heatcoef*q_zd
    q_yd = heatcoef*q_yd
    q_xd = heatcoef*q_xd
    v_zd = v_zd + tauyzsd
    w_yd = w_yd + tauyzsd
    u_zd = u_zd + tauxzsd
    w_xd = w_xd + tauxzsd
    u_yd = u_yd + tauxysd
    v_xd = v_xd + tauxysd
    w_zd = w_zd + tempd23 + two*tauzzsd
    v_yd = v_yd + tempd23 + two*tauyysd
    u_xd = u_xd + tempd23 + two*tauxxsd
    corrd = -(ssy*q_yd) - ssx*q_xd - ssz*q_zd
    q_xd = q_xd + ssx*corrd
    q_yd = q_yd + ssy*corrd
    q_zd = q_zd + ssz*corrd
    aad(i, j+1, k) = aad(i, j+1, k) + ss*corrd
    aad(i, j, k) = aad(i, j, k) - ss*corrd
    corrd = -(ssy*w_yd) - ssx*w_xd - ssz*w_zd
    w_xd = w_xd + ssx*corrd
    w_yd = w_yd + ssy*corrd
    w_zd = w_zd + ssz*corrd
    wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) - ss*corrd
    wd(i, j, k, ivz) = wd(i, j, k, ivz) + ss*corrd
    corrd = -(ssy*v_yd) - ssx*v_xd - ssz*v_zd
    v_xd = v_xd + ssx*corrd
    v_yd = v_yd + ssy*corrd
    v_zd = v_zd + ssz*corrd
    wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) - ss*corrd
    wd(i, j, k, ivy) = wd(i, j, k, ivy) + ss*corrd
    corrd = -(ssy*u_yd) - ssx*u_xd - ssz*u_zd
    u_xd = u_xd + ssx*corrd
    u_yd = u_yd + ssy*corrd
    u_zd = u_zd + ssz*corrd
    wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) - ss*corrd
    wd(i, j, k, ivx) = wd(i, j, k, ivx) + ss*corrd
    tempd24 = fourth*q_zd
    qzd(i-1, j, k-1) = qzd(i-1, j, k-1) + tempd24
    qzd(i, j, k-1) = qzd(i, j, k-1) + tempd24
    qzd(i-1, j, k) = qzd(i-1, j, k) + tempd24
    qzd(i, j, k) = qzd(i, j, k) + tempd24
    tempd25 = fourth*q_yd
    qyd(i-1, j, k-1) = qyd(i-1, j, k-1) + tempd25
    qyd(i, j, k-1) = qyd(i, j, k-1) + tempd25
    qyd(i-1, j, k) = qyd(i-1, j, k) + tempd25
    qyd(i, j, k) = qyd(i, j, k) + tempd25
    tempd26 = fourth*q_xd
    qxd(i-1, j, k-1) = qxd(i-1, j, k-1) + tempd26
    qxd(i, j, k-1) = qxd(i, j, k-1) + tempd26
    qxd(i-1, j, k) = qxd(i-1, j, k) + tempd26
    qxd(i, j, k) = qxd(i, j, k) + tempd26
    tempd27 = fourth*w_zd
    wzd(i-1, j, k-1) = wzd(i-1, j, k-1) + tempd27
    wzd(i, j, k-1) = wzd(i, j, k-1) + tempd27
    wzd(i-1, j, k) = wzd(i-1, j, k) + tempd27
    wzd(i, j, k) = wzd(i, j, k) + tempd27
    tempd28 = fourth*w_yd
    wyd(i-1, j, k-1) = wyd(i-1, j, k-1) + tempd28
    wyd(i, j, k-1) = wyd(i, j, k-1) + tempd28
    wyd(i-1, j, k) = wyd(i-1, j, k) + tempd28
    wyd(i, j, k) = wyd(i, j, k) + tempd28
    tempd29 = fourth*w_xd
    wxd(i-1, j, k-1) = wxd(i-1, j, k-1) + tempd29
    wxd(i, j, k-1) = wxd(i, j, k-1) + tempd29
    wxd(i-1, j, k) = wxd(i-1, j, k) + tempd29
    wxd(i, j, k) = wxd(i, j, k) + tempd29
    tempd30 = fourth*v_zd
    vzd(i-1, j, k-1) = vzd(i-1, j, k-1) + tempd30
    vzd(i, j, k-1) = vzd(i, j, k-1) + tempd30
    vzd(i-1, j, k) = vzd(i-1, j, k) + tempd30
    vzd(i, j, k) = vzd(i, j, k) + tempd30
    tempd31 = fourth*v_yd
    vyd(i-1, j, k-1) = vyd(i-1, j, k-1) + tempd31
    vyd(i, j, k-1) = vyd(i, j, k-1) + tempd31
    vyd(i-1, j, k) = vyd(i-1, j, k) + tempd31
    vyd(i, j, k) = vyd(i, j, k) + tempd31
    tempd32 = fourth*v_xd
    vxd(i-1, j, k-1) = vxd(i-1, j, k-1) + tempd32
    vxd(i, j, k-1) = vxd(i, j, k-1) + tempd32
    vxd(i-1, j, k) = vxd(i-1, j, k) + tempd32
    vxd(i, j, k) = vxd(i, j, k) + tempd32
    tempd33 = fourth*u_zd
    uzd(i-1, j, k-1) = uzd(i-1, j, k-1) + tempd33
    uzd(i, j, k-1) = uzd(i, j, k-1) + tempd33
    uzd(i-1, j, k) = uzd(i-1, j, k) + tempd33
    uzd(i, j, k) = uzd(i, j, k) + tempd33
    tempd34 = fourth*u_yd
    uyd(i-1, j, k-1) = uyd(i-1, j, k-1) + tempd34
    uyd(i, j, k-1) = uyd(i, j, k-1) + tempd34
    uyd(i-1, j, k) = uyd(i-1, j, k) + tempd34
    uyd(i, j, k) = uyd(i, j, k) + tempd34
    tempd35 = fourth*u_xd
    uxd(i-1, j, k-1) = uxd(i-1, j, k-1) + tempd35
    uxd(i, j, k-1) = uxd(i, j, k-1) + tempd35
    uxd(i-1, j, k) = uxd(i-1, j, k) + tempd35
    uxd(i, j, k) = uxd(i, j, k) + tempd35
    muld = mutd + factlamheat*heatcoefd
    mued = mued + mutd + factturbheat*heatcoefd
    branch = myIntStack(myIntPtr)
    myIntPtr = myIntPtr - 1
    if (branch .eq. 0) then
        revd(i, j, k) = revd(i, j, k) + por*mued
        revd(i, j+1, k) = revd(i, j+1, k) + por*mued
        mued = 0.0_8
    end if
    rlvd(i, j, k) = rlvd(i, j, k) + por*muld
    rlvd(i, j+1, k) = rlvd(i, j+1, k) + por*muld
end do
! end do
! end do
! end do
mued = 0.0_8
!
!         viscous fluxes in the k-direction.
!
mue = zero
mued = 0.0_8
! do k=1,kl
!    do j=2,jl
!       do i=2,il
do ii=0,nx*ny*kl-1
    i = mod(ii, nx) + 2
    j = mod(ii/nx, ny) + 2
    k = ii/(nx*ny) + 1
    ! set the value of the porosity. if not zero, it is set
    ! to average the eddy-viscosity and to take the factor
    ! rfilv into account.
    por = half*rfilv
    if (pork(i, j, k) .eq. noflux) por = zero
    ! compute the laminar and (if present) the eddy viscosities
    ! multiplied by the porosity. compute the factor in front of
    ! the gradients of the speed of sound squared for the heat
    ! flux.
    mul = por*(rlv(i, j, k)+rlv(i, j, k+1))
    if (eddymodel) then
        mue = por*(rev(i, j, k)+rev(i, j, k+1))
        myIntPtr = myIntPtr + 1
        myIntStack(myIntPtr) = 0
    else
        myIntPtr = myIntPtr + 1
        myIntStack(myIntPtr) = 1
    end if
    mut = mul + mue
    gm1 = half*(gamma(i, j, k)+gamma(i, j, k+1)) - one
    factlamheat = one/(prandtl*gm1)
    factturbheat = one/(prandtlturb*gm1)
    heatcoef = mul*factlamheat + mue*factturbheat
    ! compute the gradients at the face by averaging the four
    ! nodal values.
    u_x = fourth*(ux(i-1, j-1, k)+ux(i, j-1, k)+ux(i-1, j, k)+ux(i, &
    &         j, k))
    u_y = fourth*(uy(i-1, j-1, k)+uy(i, j-1, k)+uy(i-1, j, k)+uy(i, &
    &         j, k))
    u_z = fourth*(uz(i-1, j-1, k)+uz(i, j-1, k)+uz(i-1, j, k)+uz(i, &
    &         j, k))
    v_x = fourth*(vx(i-1, j-1, k)+vx(i, j-1, k)+vx(i-1, j, k)+vx(i, &
    &         j, k))
    v_y = fourth*(vy(i-1, j-1, k)+vy(i, j-1, k)+vy(i-1, j, k)+vy(i, &
    &         j, k))
    v_z = fourth*(vz(i-1, j-1, k)+vz(i, j-1, k)+vz(i-1, j, k)+vz(i, &
    &         j, k))
    w_x = fourth*(wx(i-1, j-1, k)+wx(i, j-1, k)+wx(i-1, j, k)+wx(i, &
    &         j, k))
    w_y = fourth*(wy(i-1, j-1, k)+wy(i, j-1, k)+wy(i-1, j, k)+wy(i, &
    &         j, k))
    w_z = fourth*(wz(i-1, j-1, k)+wz(i, j-1, k)+wz(i-1, j, k)+wz(i, &
    &         j, k))
    q_x = fourth*(qx(i-1, j-1, k)+qx(i, j-1, k)+qx(i-1, j, k)+qx(i, &
    &         j, k))
    q_y = fourth*(qy(i-1, j-1, k)+qy(i, j-1, k)+qy(i-1, j, k)+qy(i, &
    &         j, k))
    q_z = fourth*(qz(i-1, j-1, k)+qz(i, j-1, k)+qz(i-1, j, k)+qz(i, &
    &         j, k))
    ! the gradients in the normal direction are corrected, such
    ! that no averaging takes places here.
    ! first determine the vector in the direction from the
    ! cell center k to cell center k+1.
    ssx = eighth*(x(i-1, j-1, k+1, 1)-x(i-1, j-1, k-1, 1)+x(i-1, j, &
    &         k+1, 1)-x(i-1, j, k-1, 1)+x(i, j-1, k+1, 1)-x(i, j-1, k-1, 1)+&
    &         x(i, j, k+1, 1)-x(i, j, k-1, 1))
    ssy = eighth*(x(i-1, j-1, k+1, 2)-x(i-1, j-1, k-1, 2)+x(i-1, j, &
    &         k+1, 2)-x(i-1, j, k-1, 2)+x(i, j-1, k+1, 2)-x(i, j-1, k-1, 2)+&
    &         x(i, j, k+1, 2)-x(i, j, k-1, 2))
    ssz = eighth*(x(i-1, j-1, k+1, 3)-x(i-1, j-1, k-1, 3)+x(i-1, j, &
    &         k+1, 3)-x(i-1, j, k-1, 3)+x(i, j-1, k+1, 3)-x(i, j-1, k-1, 3)+&
    &         x(i, j, k+1, 3)-x(i, j, k-1, 3))
    ! determine the length of this vector and create the
    ! unit normal.
    ss = one/sqrt(ssx*ssx+ssy*ssy+ssz*ssz)
    ssx = ss*ssx
    ssy = ss*ssy
    ssz = ss*ssz
    ! correct the gradients.
    corr = u_x*ssx + u_y*ssy + u_z*ssz - (w(i, j, k+1, ivx)-w(i, j, &
    &         k, ivx))*ss
    u_x = u_x - corr*ssx
    u_y = u_y - corr*ssy
    u_z = u_z - corr*ssz
    corr = v_x*ssx + v_y*ssy + v_z*ssz - (w(i, j, k+1, ivy)-w(i, j, &
    &         k, ivy))*ss
    v_x = v_x - corr*ssx
    v_y = v_y - corr*ssy
    v_z = v_z - corr*ssz
    corr = w_x*ssx + w_y*ssy + w_z*ssz - (w(i, j, k+1, ivz)-w(i, j, &
    &         k, ivz))*ss
    w_x = w_x - corr*ssx
    w_y = w_y - corr*ssy
    w_z = w_z - corr*ssz
    corr = q_x*ssx + q_y*ssy + q_z*ssz + (aa(i, j, k+1)-aa(i, j, k))&
    &         *ss
    q_x = q_x - corr*ssx
    q_y = q_y - corr*ssy
    q_z = q_z - corr*ssz
    ! compute the stress tensor and the heat flux vector.
    ! we remove the viscosity from the stress tensor (tau)
    ! to define taus since we still need to separate between
    ! laminar and turbulent stress for qcr.
    ! therefore, laminar tau = mue*taus, turbulent
    ! tau = mue*taus, and total tau = mut*taus.
    fracdiv = twothird*(u_x+v_y+w_z)
    tauxxs = two*u_x - fracdiv
    tauyys = two*v_y - fracdiv
    tauzzs = two*w_z - fracdiv
    tauxys = u_y + v_x
    tauxzs = u_z + w_x
    tauyzs = v_z + w_y
    ! add qcr corrections if necessary
    if (useqcr) then
        ! in the qcr formulation, we add an extra term to the turbulent stress tensor:
        !
        ! tau_ij,qcr = tau_ij - e_ij
        !
        ! where, according to tmr website (http://turbmodels.larc.nasa.gov/spalart.html):
        !
        ! e_ij = ccr1*(o_ik*tau_jk + o_jk*tau_ik)
        !
        ! we are computing o_ik as follows:
        !
        ! o_ik = 2*w_ik/den
        !
        ! remember that the tau_ij in e_ij should use only the eddy viscosity!
        ! compute denominator
        den = sqrt(u_x*u_x + u_y*u_y + u_z*u_z + v_x*v_x + v_y*v_y + &
        &           v_z*v_z + w_x*w_x + w_y*w_y + w_z*w_z)
        if (den .lt. xminn) then
            den = xminn
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 0
        else
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 1
            den = den
        end if
        ! compute factor that will multiply all tensor components.
        ! here we add the eddy viscosity that should multiply the stress tensor (tau)
        ! components as well.
        fact = mue*ccr1/den
        ! compute off-diagonal terms of vorticity tensor (we will ommit the 1/2)
        ! the diagonals of the vorticity tensor components are always zero
        wxy = u_y - v_x
        wxz = u_z - w_x
        wyz = v_z - w_y
        wyx = -wxy
        wzx = -wxz
        wzy = -wyz
        ! compute the extra terms of the boussinesq relation
        exx = fact*(wxy*tauxys+wxz*tauxzs)*two
        eyy = fact*(wyx*tauxys+wyz*tauyzs)*two
        ezz = fact*(wzx*tauxzs+wzy*tauyzs)*two
        exy = fact*(wxy*tauyys+wxz*tauyzs+wyx*tauxxs+wyz*tauxzs)
        exz = fact*(wxy*tauyzs+wxz*tauzzs+wzx*tauxxs+wzy*tauxys)
        eyz = fact*(wyx*tauxzs+wyz*tauzzs+wzx*tauxys+wzy*tauyys)
        ! apply the total viscosity to the stress tensor and add extra terms
        tauxx = mut*tauxxs - exx
        tauyy = mut*tauyys - eyy
        tauzz = mut*tauzzs - ezz
        tauxy = mut*tauxys - exy
        tauxz = mut*tauxzs - exz
        tauyz = mut*tauyzs - eyz
        myIntPtr = myIntPtr + 1
        myIntStack(myIntPtr) = 0
    else
        ! just apply the total viscosity to the stress tensor
        tauxx = mut*tauxxs
        tauyy = mut*tauyys
        tauzz = mut*tauzzs
        tauxy = mut*tauxys
        tauxz = mut*tauxzs
        tauyz = mut*tauyzs
        myIntPtr = myIntPtr + 1
        myIntStack(myIntPtr) = 1
    end if
    ! compute the average velocities for the face. remember that
    ! the velocities are stored and not the momentum.
    ubar = half*(w(i, j, k, ivx)+w(i, j, k+1, ivx))
    vbar = half*(w(i, j, k, ivy)+w(i, j, k+1, ivy))
    wbar = half*(w(i, j, k, ivz)+w(i, j, k+1, ivz))
    ! compute the viscous fluxes for this k-face.
    ! update the residuals of cell k and k+1.
    ! store the stress tensor and the heat flux vector if this
    ! face is part of a viscous subface. both the cases k == 1
    ! and k == kl must be tested.
    frhoed = fwd(i, j, k+1, irhoe) - fwd(i, j, k, irhoe)
    fmzd = fwd(i, j, k+1, imz) - fwd(i, j, k, imz)
    fmyd = fwd(i, j, k+1, imy) - fwd(i, j, k, imy)
    fmxd = fwd(i, j, k+1, imx) - fwd(i, j, k, imx)
    q_xd = -(sk(i, j, k, 1)*frhoed)
    q_yd = -(sk(i, j, k, 2)*frhoed)
    q_zd = -(sk(i, j, k, 3)*frhoed)
    tempd20 = sk(i, j, k, 3)*frhoed
    tauzzd = sk(i, j, k, 3)*fmzd + wbar*tempd20
    tempd21 = sk(i, j, k, 2)*frhoed
    tauyzd = wbar*tempd21 + sk(i, j, k, 3)*fmyd + sk(i, j, k, 2)*&
    &         fmzd + vbar*tempd20
    tauyyd = sk(i, j, k, 2)*fmyd + vbar*tempd21
    tempd22 = sk(i, j, k, 1)*frhoed
    ubard = tauxy*tempd21 + tauxx*tempd22 + tauxz*tempd20
    tauxzd = wbar*tempd22 + sk(i, j, k, 3)*fmxd + sk(i, j, k, 1)*&
    &         fmzd + ubar*tempd20
    vbard = tauyy*tempd21 + tauxy*tempd22 + tauyz*tempd20
    wbard = tauyz*tempd21 + tauxz*tempd22 + tauzz*tempd20
    tauxyd = vbar*tempd22 + sk(i, j, k, 2)*fmxd + sk(i, j, k, 1)*&
    &         fmyd + ubar*tempd21
    tauxxd = sk(i, j, k, 1)*fmxd + ubar*tempd22
    wd(i, j, k, ivz) = wd(i, j, k, ivz) + half*wbard
    wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) + half*wbard
    wd(i, j, k, ivy) = wd(i, j, k, ivy) + half*vbard
    wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) + half*vbard
    wd(i, j, k, ivx) = wd(i, j, k, ivx) + half*ubard
    wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) + half*ubard
    branch = myIntStack(myIntPtr)
    myIntPtr = myIntPtr - 1
    if (branch .eq. 0) then
        exzd = -tauxzd
        exyd = -tauxyd
        ezzd = -tauzzd
        eyyd = -tauyyd
        tempd13 = fact*exzd
        tempd16 = fact*exyd
        tempd14 = two*fact*ezzd
        tempd15 = two*fact*eyyd
        mutd = tauxzs*tauxzd + tauzzs*tauzzd + tauxxs*tauxxd + tauyys*&
        &           tauyyd + tauxys*tauxyd + tauyzs*tauyzd
        tauyzsd = wxy*tempd13 + wzy*tempd14 + wyz*tempd15 + wxz*&
        &           tempd16 + mut*tauyzd
        eyzd = -tauyzd
        tauxxsd = wzx*tempd13 + wyx*tempd16 + mut*tauxxd
        exxd = -tauxxd
        tempd17 = fact*eyzd
        tauzzsd = wyz*tempd17 + wxz*tempd13 + mut*tauzzd
        tauyysd = wzy*tempd17 + wxy*tempd16 + mut*tauyyd
        factd = (wxy*tauyzs+wxz*tauzzs+wzx*tauxxs+wzy*tauxys)*exzd + &
        &           two*(wzx*tauxzs+wzy*tauyzs)*ezzd + two*(wxy*tauxys+wxz*&
        &           tauxzs)*exxd + two*(wyx*tauxys+wyz*tauyzs)*eyyd + (wxy*&
        &           tauyys+wxz*tauyzs+wyx*tauxxs+wyz*tauxzs)*exyd + (wyx*tauxzs+&
        &           wyz*tauzzs+wzx*tauxys+wzy*tauyys)*eyzd
        wyxd = tauxxs*tempd16 + tauxys*tempd15 + tauxzs*tempd17
        wzxd = tauxxs*tempd13 + tauxzs*tempd14 + tauxys*tempd17
        wzyd = tauxys*tempd13 + tauyzs*tempd14 + tauyys*tempd17
        wyzd = tauxzs*tempd16 - wzyd + tauyzs*tempd15 + tauzzs*tempd17
        tempd18 = two*fact*exxd
        tauxzsd = wyx*tempd17 + wzx*tempd14 + wxz*tempd18 + wyz*&
        &           tempd16 + mut*tauxzd
        tauxysd = wzx*tempd17 + wyx*tempd15 + wxy*tempd18 + wzy*&
        &           tempd13 + mut*tauxyd
        wxyd = tauyys*tempd16 - wyxd + tauxys*tempd18 + tauyzs*tempd13
        wxzd = tauyzs*tempd16 - wzxd + tauxzs*tempd18 + tauzzs*tempd13
        v_zd = wyzd
        w_yd = -wyzd
        u_zd = wxzd
        w_xd = -wxzd
        u_yd = wxyd
        v_xd = -wxyd
        tempd19 = ccr1*factd/den
        mued = mued + tempd19
        dend = -(mue*tempd19/den)
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) dend = 0.0_8
        if (u_x**2 + u_y**2 + u_z**2 + v_x**2 + v_y**2 + v_z**2 + w_x&
        &             **2 + w_y**2 + w_z**2 .eq. 0.0_8) then
        tempd12 = 0.0
    else
        tempd12 = dend/(2.0*sqrt(u_x**2+u_y**2+u_z**2+v_x**2+v_y**2+&
        &             v_z**2+w_x**2+w_y**2+w_z**2))
    end if
    u_xd = 2*u_x*tempd12
    u_yd = u_yd + 2*u_y*tempd12
    u_zd = u_zd + 2*u_z*tempd12
    v_xd = v_xd + 2*v_x*tempd12
    v_yd = 2*v_y*tempd12
    v_zd = v_zd + 2*v_z*tempd12
    w_xd = w_xd + 2*w_x*tempd12
    w_yd = w_yd + 2*w_y*tempd12
    w_zd = 2*w_z*tempd12
else
    mutd = tauxzs*tauxzd + tauzzs*tauzzd + tauxxs*tauxxd + tauyys*&
    &           tauyyd + tauxys*tauxyd + tauyzs*tauyzd
    tauyzsd = mut*tauyzd
    tauxzsd = mut*tauxzd
    tauxysd = mut*tauxyd
    tauzzsd = mut*tauzzd
    tauyysd = mut*tauyyd
    tauxxsd = mut*tauxxd
    u_xd = 0.0_8
    u_yd = 0.0_8
    u_zd = 0.0_8
    w_xd = 0.0_8
    w_yd = 0.0_8
    w_zd = 0.0_8
    v_xd = 0.0_8
    v_yd = 0.0_8
    v_zd = 0.0_8
end if
fracdivd = -tauyysd - tauxxsd - tauzzsd
tempd = twothird*fracdivd
heatcoefd = q_y*q_yd + q_x*q_xd + q_z*q_zd
q_zd = heatcoef*q_zd
q_yd = heatcoef*q_yd
q_xd = heatcoef*q_xd
v_zd = v_zd + tauyzsd
w_yd = w_yd + tauyzsd
u_zd = u_zd + tauxzsd
w_xd = w_xd + tauxzsd
u_yd = u_yd + tauxysd
v_xd = v_xd + tauxysd
w_zd = w_zd + tempd + two*tauzzsd
v_yd = v_yd + tempd + two*tauyysd
u_xd = u_xd + tempd + two*tauxxsd
corrd = -(ssy*q_yd) - ssx*q_xd - ssz*q_zd
q_xd = q_xd + ssx*corrd
q_yd = q_yd + ssy*corrd
q_zd = q_zd + ssz*corrd
aad(i, j, k+1) = aad(i, j, k+1) + ss*corrd
aad(i, j, k) = aad(i, j, k) - ss*corrd
corrd = -(ssy*w_yd) - ssx*w_xd - ssz*w_zd
w_xd = w_xd + ssx*corrd
w_yd = w_yd + ssy*corrd
w_zd = w_zd + ssz*corrd
wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) - ss*corrd
wd(i, j, k, ivz) = wd(i, j, k, ivz) + ss*corrd
corrd = -(ssy*v_yd) - ssx*v_xd - ssz*v_zd
v_xd = v_xd + ssx*corrd
v_yd = v_yd + ssy*corrd
v_zd = v_zd + ssz*corrd
wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) - ss*corrd
wd(i, j, k, ivy) = wd(i, j, k, ivy) + ss*corrd
corrd = -(ssy*u_yd) - ssx*u_xd - ssz*u_zd
u_xd = u_xd + ssx*corrd
u_yd = u_yd + ssy*corrd
u_zd = u_zd + ssz*corrd
wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) - ss*corrd
wd(i, j, k, ivx) = wd(i, j, k, ivx) + ss*corrd
tempd0 = fourth*q_zd
qzd(i-1, j-1, k) = qzd(i-1, j-1, k) + tempd0
qzd(i, j-1, k) = qzd(i, j-1, k) + tempd0
qzd(i-1, j, k) = qzd(i-1, j, k) + tempd0
qzd(i, j, k) = qzd(i, j, k) + tempd0
tempd1 = fourth*q_yd
qyd(i-1, j-1, k) = qyd(i-1, j-1, k) + tempd1
qyd(i, j-1, k) = qyd(i, j-1, k) + tempd1
qyd(i-1, j, k) = qyd(i-1, j, k) + tempd1
qyd(i, j, k) = qyd(i, j, k) + tempd1
tempd2 = fourth*q_xd
qxd(i-1, j-1, k) = qxd(i-1, j-1, k) + tempd2
qxd(i, j-1, k) = qxd(i, j-1, k) + tempd2
qxd(i-1, j, k) = qxd(i-1, j, k) + tempd2
qxd(i, j, k) = qxd(i, j, k) + tempd2
tempd3 = fourth*w_zd
wzd(i-1, j-1, k) = wzd(i-1, j-1, k) + tempd3
wzd(i, j-1, k) = wzd(i, j-1, k) + tempd3
wzd(i-1, j, k) = wzd(i-1, j, k) + tempd3
wzd(i, j, k) = wzd(i, j, k) + tempd3
tempd4 = fourth*w_yd
wyd(i-1, j-1, k) = wyd(i-1, j-1, k) + tempd4
wyd(i, j-1, k) = wyd(i, j-1, k) + tempd4
wyd(i-1, j, k) = wyd(i-1, j, k) + tempd4
wyd(i, j, k) = wyd(i, j, k) + tempd4
tempd5 = fourth*w_xd
wxd(i-1, j-1, k) = wxd(i-1, j-1, k) + tempd5
wxd(i, j-1, k) = wxd(i, j-1, k) + tempd5
wxd(i-1, j, k) = wxd(i-1, j, k) + tempd5
wxd(i, j, k) = wxd(i, j, k) + tempd5
tempd6 = fourth*v_zd
vzd(i-1, j-1, k) = vzd(i-1, j-1, k) + tempd6
vzd(i, j-1, k) = vzd(i, j-1, k) + tempd6
vzd(i-1, j, k) = vzd(i-1, j, k) + tempd6
vzd(i, j, k) = vzd(i, j, k) + tempd6
tempd7 = fourth*v_yd
vyd(i-1, j-1, k) = vyd(i-1, j-1, k) + tempd7
vyd(i, j-1, k) = vyd(i, j-1, k) + tempd7
vyd(i-1, j, k) = vyd(i-1, j, k) + tempd7
vyd(i, j, k) = vyd(i, j, k) + tempd7
tempd8 = fourth*v_xd
vxd(i-1, j-1, k) = vxd(i-1, j-1, k) + tempd8
vxd(i, j-1, k) = vxd(i, j-1, k) + tempd8
vxd(i-1, j, k) = vxd(i-1, j, k) + tempd8
vxd(i, j, k) = vxd(i, j, k) + tempd8
tempd9 = fourth*u_zd
uzd(i-1, j-1, k) = uzd(i-1, j-1, k) + tempd9
uzd(i, j-1, k) = uzd(i, j-1, k) + tempd9
uzd(i-1, j, k) = uzd(i-1, j, k) + tempd9
uzd(i, j, k) = uzd(i, j, k) + tempd9
tempd10 = fourth*u_yd
uyd(i-1, j-1, k) = uyd(i-1, j-1, k) + tempd10
uyd(i, j-1, k) = uyd(i, j-1, k) + tempd10
uyd(i-1, j, k) = uyd(i-1, j, k) + tempd10
uyd(i, j, k) = uyd(i, j, k) + tempd10
tempd11 = fourth*u_xd
uxd(i-1, j-1, k) = uxd(i-1, j-1, k) + tempd11
uxd(i, j-1, k) = uxd(i, j-1, k) + tempd11
uxd(i-1, j, k) = uxd(i-1, j, k) + tempd11
uxd(i, j, k) = uxd(i, j, k) + tempd11
muld = mutd + factlamheat*heatcoefd
mued = mued + mutd + factturbheat*heatcoefd
branch = myIntStack(myIntPtr)
myIntPtr = myIntPtr - 1
if (branch .eq. 0) then
    revd(i, j, k) = revd(i, j, k) + por*mued
    revd(i, j, k+1) = revd(i, j, k+1) + por*mued
    mued = 0.0_8
end if
rlvd(i, j, k) = rlvd(i, j, k) + por*muld
rlvd(i, j, k+1) = rlvd(i, j, k+1) + por*muld
end do
! end do
! end do
! end do
end if
end subroutine viscousflux_fast_b

!  differentiation of allnodalgradients in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: *aa *wx *wy *wz *w *qx *qy
!                *qz *ux *uy *uz *vx *vy *vz
!   with respect to varying inputs: *aa *wx *wy *wz *w *qx *qy
!                *qz *ux *uy *uz *vx *vy *vz
!   rw status of diff variables: *aa:incr *wx:in-zero *wy:in-zero
!                *wz:in-zero *w:incr *qx:in-zero *qy:in-zero *qz:in-zero
!                *ux:in-zero *uy:in-zero *uz:in-zero *vx:in-zero
!                *vy:in-zero *vz:in-zero
!   plus diff mem management of: aa:in wx:in wy:in wz:in w:in qx:in
!                qy:in qz:in ux:in uy:in uz:in vx:in vy:in vz:in
subroutine allnodalgradients_fast_b()
    !
    !         nodalgradients computes the nodal velocity gradients and
    !         minus the gradient of the speed of sound squared. the minus
    !         sign is present, because this is the definition of the heat
    !         flux. these gradients are computed for all nodes.
    !
    use constants
    ! use blockpointers
    implicit none
    !        local variables.
    integer(kind=inttype) :: i, j, k
    integer(kind=inttype) :: k1, kk
    integer(kind=inttype) :: istart, iend, isize, ii
    integer(kind=inttype) :: jstart, jend, jsize
    integer(kind=inttype) :: kstart, kend, ksize
    real(kind=realtype) :: oneoverv, ubar, vbar, wbar, a2
    real(kind=realtype) :: ubard, vbard, wbard, a2d
    real(kind=realtype) :: sx, sx1, sy, sy1, sz, sz1
    intrinsic mod
    ! zero all nodeal gradients:
    integer :: branch
    real(kind=realtype) :: tempd10
    real(kind=realtype) :: tempd9
    real(kind=realtype) :: tempd
    real(kind=realtype) :: tempd8
    real(kind=realtype) :: tempd7
    real(kind=realtype) :: tempd6
    real(kind=realtype) :: tempd5
    real(kind=realtype) :: tempd4
    real(kind=realtype) :: tempd3
    real(kind=realtype) :: tempd2
    real(kind=realtype) :: tempd1
    real(kind=realtype) :: tempd0
    do ii=0,il*jl*kl-1
        i = mod(ii, il) + 1
        j = mod(ii/il, jl) + 1
        k = ii/(il*jl) + 1
        ! compute the inverse of 8 times the volume for this node.
        oneoverv = one/(vol(i, j, k)+vol(i, j, k+1)+vol(i+1, j, k)+vol(i+1&
        &       , j, k+1)+vol(i, j+1, k)+vol(i, j+1, k+1)+vol(i+1, j+1, k)+vol(i&
        &       +1, j+1, k+1))
        ! compute the correct velocity gradients and "unit" heat
        ! fluxes. the velocity gradients are stored in ux, etc.
        qzd(i, j, k) = oneoverv*qzd(i, j, k)
        qyd(i, j, k) = oneoverv*qyd(i, j, k)
        qxd(i, j, k) = oneoverv*qxd(i, j, k)
        wzd(i, j, k) = oneoverv*wzd(i, j, k)
        wyd(i, j, k) = oneoverv*wyd(i, j, k)
        wxd(i, j, k) = oneoverv*wxd(i, j, k)
        vzd(i, j, k) = oneoverv*vzd(i, j, k)
        vyd(i, j, k) = oneoverv*vyd(i, j, k)
        vxd(i, j, k) = oneoverv*vxd(i, j, k)
        uzd(i, j, k) = oneoverv*uzd(i, j, k)
        uyd(i, j, k) = oneoverv*uyd(i, j, k)
        uxd(i, j, k) = oneoverv*uxd(i, j, k)
    end do
    do ii=0,ie*jl*kl-1
        i = mod(ii, ie) + 1
        j = mod(ii/ie, jl) + 1
        k = ii/(ie*jl) + 1
        ! compute 8 times the average normal for this part of
        ! the control volume. the factor 8 is taken care of later
        ! on when the division by the volume takes place.
        sx = si(i-1, j, k, 1) + si(i-1, j+1, k, 1) + si(i-1, j, k+1, 1) + &
        &       si(i-1, j+1, k+1, 1) + si(i, j, k, 1) + si(i, j+1, k, 1) + si(i&
        &       , j, k+1, 1) + si(i, j+1, k+1, 1)
        sy = si(i-1, j, k, 2) + si(i-1, j+1, k, 2) + si(i-1, j, k+1, 2) + &
        &       si(i-1, j+1, k+1, 2) + si(i, j, k, 2) + si(i, j+1, k, 2) + si(i&
        &       , j, k+1, 2) + si(i, j+1, k+1, 2)
        sz = si(i-1, j, k, 3) + si(i-1, j+1, k, 3) + si(i-1, j, k+1, 3) + &
        &       si(i-1, j+1, k+1, 3) + si(i, j, k, 3) + si(i, j+1, k, 3) + si(i&
        &       , j, k+1, 3) + si(i, j+1, k+1, 3)
        ! compute the average velocities and speed of sound squared
        ! for this integration point. node that these variables are
        ! stored in w(ivx), w(ivy), w(ivz) and p.
        ! add the contributions to the surface integral to the node
        ! j-1 and substract it from the node j. for the heat flux it
        ! is reversed, because the negative of the gradient of the
        ! speed of sound must be computed.
        if (i .gt. 1) then
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 0
        else
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 1
        end if
        if (i .lt. ie) then
            a2d = sy*qyd(i, j, k) + sx*qxd(i, j, k) + sz*qzd(i, j, k)
            wbard = -(sy*wyd(i, j, k)) - sx*wxd(i, j, k) - sz*wzd(i, j, k)
            vbard = -(sy*vyd(i, j, k)) - sx*vxd(i, j, k) - sz*vzd(i, j, k)
            ubard = -(sy*uyd(i, j, k)) - sx*uxd(i, j, k) - sz*uzd(i, j, k)
        else
            wbard = 0.0_8
            vbard = 0.0_8
            ubard = 0.0_8
            a2d = 0.0_8
        end if
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) then
            a2d = a2d - sy*qyd(i-1, j, k) - sx*qxd(i-1, j, k) - sz*qzd(i-1, &
            &         j, k)
            wbard = wbard + sy*wyd(i-1, j, k) + sx*wxd(i-1, j, k) + sz*wzd(i&
            &         -1, j, k)
            vbard = vbard + sy*vyd(i-1, j, k) + sx*vxd(i-1, j, k) + sz*vzd(i&
            &         -1, j, k)
            ubard = ubard + sy*uyd(i-1, j, k) + sx*uxd(i-1, j, k) + sz*uzd(i&
            &         -1, j, k)
        end if
        tempd7 = fourth*a2d
        aad(i, j, k) = aad(i, j, k) + tempd7
        aad(i, j+1, k) = aad(i, j+1, k) + tempd7
        aad(i, j, k+1) = aad(i, j, k+1) + tempd7
        aad(i, j+1, k+1) = aad(i, j+1, k+1) + tempd7
        tempd8 = fourth*wbard
        wd(i, j, k, ivz) = wd(i, j, k, ivz) + tempd8
        wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) + tempd8
        wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) + tempd8
        wd(i, j+1, k+1, ivz) = wd(i, j+1, k+1, ivz) + tempd8
        tempd9 = fourth*vbard
        wd(i, j, k, ivy) = wd(i, j, k, ivy) + tempd9
        wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) + tempd9
        wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) + tempd9
        wd(i, j+1, k+1, ivy) = wd(i, j+1, k+1, ivy) + tempd9
        tempd10 = fourth*ubard
        wd(i, j, k, ivx) = wd(i, j, k, ivx) + tempd10
        wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) + tempd10
        wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) + tempd10
        wd(i, j+1, k+1, ivx) = wd(i, j+1, k+1, ivx) + tempd10
    end do
    do ii=0,il*je*kl-1
        i = mod(ii, il) + 1
        j = mod(ii/il, je) + 1
        k = ii/(il*je) + 1
        ! compute 8 times the average normal for this part of
        ! the control volume. the factor 8 is taken care of later
        ! on when the division by the volume takes place.
        sx = sj(i, j-1, k, 1) + sj(i+1, j-1, k, 1) + sj(i, j-1, k+1, 1) + &
        &       sj(i+1, j-1, k+1, 1) + sj(i, j, k, 1) + sj(i+1, j, k, 1) + sj(i&
        &       , j, k+1, 1) + sj(i+1, j, k+1, 1)
        sy = sj(i, j-1, k, 2) + sj(i+1, j-1, k, 2) + sj(i, j-1, k+1, 2) + &
        &       sj(i+1, j-1, k+1, 2) + sj(i, j, k, 2) + sj(i+1, j, k, 2) + sj(i&
        &       , j, k+1, 2) + sj(i+1, j, k+1, 2)
        sz = sj(i, j-1, k, 3) + sj(i+1, j-1, k, 3) + sj(i, j-1, k+1, 3) + &
        &       sj(i+1, j-1, k+1, 3) + sj(i, j, k, 3) + sj(i+1, j, k, 3) + sj(i&
        &       , j, k+1, 3) + sj(i+1, j, k+1, 3)
        ! compute the average velocities and speed of sound squared
        ! for this integration point. node that these variables are
        ! stored in w(ivx), w(ivy), w(ivz) and p.
        ! add the contributions to the surface integral to the node
        ! j-1 and substract it from the node j. for the heat flux it
        ! is reversed, because the negative of the gradient of the
        ! speed of sound must be computed.
        if (j .gt. 1) then
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 0
        else
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 1
        end if
        if (j .lt. je) then
            a2d = sy*qyd(i, j, k) + sx*qxd(i, j, k) + sz*qzd(i, j, k)
            wbard = -(sy*wyd(i, j, k)) - sx*wxd(i, j, k) - sz*wzd(i, j, k)
            vbard = -(sy*vyd(i, j, k)) - sx*vxd(i, j, k) - sz*vzd(i, j, k)
            ubard = -(sy*uyd(i, j, k)) - sx*uxd(i, j, k) - sz*uzd(i, j, k)
        else
            wbard = 0.0_8
            vbard = 0.0_8
            ubard = 0.0_8
            a2d = 0.0_8
        end if
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) then
            a2d = a2d - sy*qyd(i, j-1, k) - sx*qxd(i, j-1, k) - sz*qzd(i, j-&
            &         1, k)
            wbard = wbard + sy*wyd(i, j-1, k) + sx*wxd(i, j-1, k) + sz*wzd(i&
            &         , j-1, k)
            vbard = vbard + sy*vyd(i, j-1, k) + sx*vxd(i, j-1, k) + sz*vzd(i&
            &         , j-1, k)
            ubard = ubard + sy*uyd(i, j-1, k) + sx*uxd(i, j-1, k) + sz*uzd(i&
            &         , j-1, k)
        end if
        tempd3 = fourth*a2d
        aad(i, j, k) = aad(i, j, k) + tempd3
        aad(i+1, j, k) = aad(i+1, j, k) + tempd3
        aad(i, j, k+1) = aad(i, j, k+1) + tempd3
        aad(i+1, j, k+1) = aad(i+1, j, k+1) + tempd3
        tempd4 = fourth*wbard
        wd(i, j, k, ivz) = wd(i, j, k, ivz) + tempd4
        wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) + tempd4
        wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) + tempd4
        wd(i+1, j, k+1, ivz) = wd(i+1, j, k+1, ivz) + tempd4
        tempd5 = fourth*vbard
        wd(i, j, k, ivy) = wd(i, j, k, ivy) + tempd5
        wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) + tempd5
        wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) + tempd5
        wd(i+1, j, k+1, ivy) = wd(i+1, j, k+1, ivy) + tempd5
        tempd6 = fourth*ubard
        wd(i, j, k, ivx) = wd(i, j, k, ivx) + tempd6
        wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) + tempd6
        wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) + tempd6
        wd(i+1, j, k+1, ivx) = wd(i+1, j, k+1, ivx) + tempd6
    end do
    do ii=0,il*jl*ke-1
        i = mod(ii, il) + 1
        j = mod(ii/il, jl) + 1
        k = ii/(il*jl) + 1
        ! compute 8 times the average normal for this part of
        ! the control volume. the factor 8 is taken care of later
        ! on when the division by the volume takes place.
        sx = sk(i, j, k-1, 1) + sk(i+1, j, k-1, 1) + sk(i, j+1, k-1, 1) + &
        &       sk(i+1, j+1, k-1, 1) + sk(i, j, k, 1) + sk(i+1, j, k, 1) + sk(i&
        &       , j+1, k, 1) + sk(i+1, j+1, k, 1)
        sy = sk(i, j, k-1, 2) + sk(i+1, j, k-1, 2) + sk(i, j+1, k-1, 2) + &
        &       sk(i+1, j+1, k-1, 2) + sk(i, j, k, 2) + sk(i+1, j, k, 2) + sk(i&
        &       , j+1, k, 2) + sk(i+1, j+1, k, 2)
        sz = sk(i, j, k-1, 3) + sk(i+1, j, k-1, 3) + sk(i, j+1, k-1, 3) + &
        &       sk(i+1, j+1, k-1, 3) + sk(i, j, k, 3) + sk(i+1, j, k, 3) + sk(i&
        &       , j+1, k, 3) + sk(i+1, j+1, k, 3)
        ! compute the average velocities and speed of sound squared
        ! for this integration point. node that these variables are
        ! stored in w(ivx), w(ivy), w(ivz) and p.
        ! add the contributions to the surface integral to the node
        ! j-1 and substract it from the node j. for the heat flux it
        ! is reversed, because the negative of the gradient of the
        ! speed of sound must be computed.
        if (k .gt. 1) then
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 0
        else
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 1
        end if
        if (k .lt. ke) then
            a2d = sy*qyd(i, j, k) + sx*qxd(i, j, k) + sz*qzd(i, j, k)
            wbard = -(sy*wyd(i, j, k)) - sx*wxd(i, j, k) - sz*wzd(i, j, k)
            vbard = -(sy*vyd(i, j, k)) - sx*vxd(i, j, k) - sz*vzd(i, j, k)
            ubard = -(sy*uyd(i, j, k)) - sx*uxd(i, j, k) - sz*uzd(i, j, k)
        else
            wbard = 0.0_8
            vbard = 0.0_8
            ubard = 0.0_8
            a2d = 0.0_8
        end if
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) then
            a2d = a2d - sy*qyd(i, j, k-1) - sx*qxd(i, j, k-1) - sz*qzd(i, j&
            &         , k-1)
            wbard = wbard + sy*wyd(i, j, k-1) + sx*wxd(i, j, k-1) + sz*wzd(i&
            &         , j, k-1)
            vbard = vbard + sy*vyd(i, j, k-1) + sx*vxd(i, j, k-1) + sz*vzd(i&
            &         , j, k-1)
            ubard = ubard + sy*uyd(i, j, k-1) + sx*uxd(i, j, k-1) + sz*uzd(i&
            &         , j, k-1)
        end if
        tempd = fourth*a2d
        aad(i, j, k) = aad(i, j, k) + tempd
        aad(i+1, j, k) = aad(i+1, j, k) + tempd
        aad(i, j+1, k) = aad(i, j+1, k) + tempd
        aad(i+1, j+1, k) = aad(i+1, j+1, k) + tempd
        tempd0 = fourth*wbard
        wd(i, j, k, ivz) = wd(i, j, k, ivz) + tempd0
        wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) + tempd0
        wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) + tempd0
        wd(i+1, j+1, k, ivz) = wd(i+1, j+1, k, ivz) + tempd0
        tempd1 = fourth*vbard
        wd(i, j, k, ivy) = wd(i, j, k, ivy) + tempd1
        wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) + tempd1
        wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) + tempd1
        wd(i+1, j+1, k, ivy) = wd(i+1, j+1, k, ivy) + tempd1
        tempd2 = fourth*ubard
        wd(i, j, k, ivx) = wd(i, j, k, ivx) + tempd2
        wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) + tempd2
        wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) + tempd2
        wd(i+1, j+1, k, ivx) = wd(i+1, j+1, k, ivx) + tempd2
    end do
    wxd = 0.0_8
    wyd = 0.0_8
    wzd = 0.0_8
    qxd = 0.0_8
    qyd = 0.0_8
    qzd = 0.0_8
    uxd = 0.0_8
    uyd = 0.0_8
    uzd = 0.0_8
    vxd = 0.0_8
    vyd = 0.0_8
    vzd = 0.0_8
end subroutine allnodalgradients_fast_b

!  differentiation of computespeedofsoundsquared in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: *aa *p *w
!   with respect to varying inputs: *aa *p *w
!   rw status of diff variables: *aa:in-out *p:incr *w:incr
!   plus diff mem management of: aa:in p:in w:in
subroutine computespeedofsoundsquared_b()
    !
    !       computespeedofsoundsquared does what it says.
    !
    use constants
    ! use blockpointers, only : ie, je, ke, w, wd, p, pd, aa, aad, gamma
    use utils_b, only : getcorrectfork
    implicit none
    !
    !      local variables.
    !
    real(kind=realtype), parameter :: twothird=two*third
    integer(kind=inttype) :: i, j, k, ii
    real(kind=realtype) :: pp
    real(kind=realtype) :: ppd
    logical :: correctfork
    intrinsic mod
    real(kind=realtype) :: temp0
    real(kind=realtype) :: tempd
    real(kind=realtype) :: tempd0
    real(kind=realtype) :: temp
    ! determine if we need to correct for k
    correctfork = getcorrectfork()
    if (correctfork) then
        do ii=0,ie*je*ke-1
            i = mod(ii, ie) + 1
            j = mod(ii/ie, je) + 1
            k = ii/(ie*je) + 1
            pp = p(i, j, k) - twothird*w(i, j, k, irho)*w(i, j, k, itu1)
            temp = w(i, j, k, irho)
            tempd = gamma(i, j, k)*aad(i, j, k)/temp
            ppd = tempd
            wd(i, j, k, irho) = wd(i, j, k, irho) - pp*tempd/temp
            aad(i, j, k) = 0.0_8
            pd(i, j, k) = pd(i, j, k) + ppd
            wd(i, j, k, irho) = wd(i, j, k, irho) - twothird*w(i, j, k, itu1&
            &         )*ppd
            wd(i, j, k, itu1) = wd(i, j, k, itu1) - twothird*w(i, j, k, irho&
            &         )*ppd
        end do
    else
        do ii=0,ie*je*ke-1
            i = mod(ii, ie) + 1
            j = mod(ii/ie, je) + 1
            k = ii/(ie*je) + 1
            temp0 = w(i, j, k, irho)
            tempd0 = gamma(i, j, k)*aad(i, j, k)/temp0
            pd(i, j, k) = pd(i, j, k) + tempd0
            wd(i, j, k, irho) = wd(i, j, k, irho) - p(i, j, k)*tempd0/temp0
            aad(i, j, k) = 0.0_8
        end do
    end if
end subroutine computespeedofsoundsquared_b

!  differentiation of invisciddissfluxscalar in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: *p *w *fw
!   with respect to varying inputs: *p *w *fw *radi *radj *radk
!   rw status of diff variables: *p:incr *w:incr *fw:in-out *radi:out
!                *radj:out *radk:out
!   plus diff mem management of: p:in w:in fw:in radi:in radj:in
!                radk:in
subroutine invisciddissfluxscalar_fast_b()
    !
    !       invisciddissfluxscalar computes the scalar artificial
    !       dissipation, see aiaa paper 81-1259, for a given block.
    !       therefore it is assumed that the pointers in  blockpointers
    !       already point to the correct block.
    !
    use constants
    ! use blockpointers, only : nx, ny, nz, il, jl, kl, ie, je, ke, ib, &
    ! &   jb, kb, w, wd, p, pd, pori, porj, pork, fw, fwd, radi, radid, radj, &
    ! &   radjd, radk, radkd, gamma
    use flowvarrefstate, only : gammainf, pinfcorr, rhoinf
    use inputdiscretization, only : vis2, vis4
    use inputphysics, only : equations
    use iteration, only : rfil
    use utils_fast_b, only : mydim, mydim_fast_b
    implicit none
    !
    !      local parameter.
    !
    real(kind=realtype), parameter :: dssmax=0.25_realtype
    !
    !      local variables.
    !
    integer(kind=inttype) :: i, j, k, ind, ii
    real(kind=realtype) :: sslim, rhoi
    real(kind=realtype) :: sfil, fis2, fis4
    real(kind=realtype) :: ppor, rrad, dis2, dis4
    real(kind=realtype) :: rradd, dis2d, dis4d
    real(kind=realtype) :: ddw1, ddw2, ddw3, ddw4, ddw5, fs
    real(kind=realtype) :: ddw1d, ddw2d, ddw3d, ddw4d, ddw5d, fsd
    real(kind=realtype), dimension(ie, je, ke, 3) :: dss
    real(kind=realtype), dimension(ie, je, ke, 3) :: dssd
    real(kind=realtype), dimension(0:ib, 0:jb, 0:kb) :: ss
    real(kind=realtype), dimension(0:ib, 0:jb, 0:kb) :: ssd
    intrinsic abs
    intrinsic mod
    intrinsic max
    intrinsic min
    real(kind=realtype) :: arg1
    real(kind=realtype) :: arg1d
    integer :: branch
    real(kind=realtype) :: temp3
    real(kind=realtype) :: tempd14
    real(kind=realtype) :: temp29
    real(kind=realtype) :: temp2
    real(kind=realtype) :: temp28
    real(kind=realtype) :: tempd13
    real(kind=realtype) :: temp1
    real(kind=realtype) :: temp27
    real(kind=realtype) :: tempd12
    real(kind=realtype) :: temp0
    real(kind=realtype) :: temp26
    real(kind=realtype) :: tempd11
    real(kind=realtype) :: temp25
    real(kind=realtype) :: tempd10
    real(kind=realtype) :: temp24
    real(kind=realtype) :: temp23
    real(kind=realtype) :: temp22
    real(kind=realtype) :: temp21
    real(kind=realtype) :: temp20
    real(kind=realtype) :: min3
    real(kind=realtype) :: min2
    real(kind=realtype) :: min1
    real(kind=realtype) :: min1d
    real(kind=realtype) :: x3
    real(kind=realtype) :: x2
    real(kind=realtype) :: x2d
    real(kind=realtype) :: x1
    real(kind=realtype) :: temp19
    real(kind=realtype) :: temp18
    real(kind=realtype) :: temp17
    real(kind=realtype) :: temp16
    real(kind=realtype) :: temp15
    real(kind=realtype) :: temp14
    real(kind=realtype) :: temp13
    real(kind=realtype) :: y3d
    real(kind=realtype) :: temp12
    real(kind=realtype) :: temp11
    real(kind=realtype) :: temp10
    real(kind=realtype) :: temp40
    real(kind=realtype) :: tempd9
    real(kind=realtype) :: tempd
    real(kind=realtype) :: tempd8
    real(kind=realtype) :: tempd7
    real(kind=realtype) :: tempd6
    real(kind=realtype) :: tempd5
    real(kind=realtype) :: tempd4
    real(kind=realtype) :: tempd3
    real(kind=realtype) :: tempd2
    real(kind=realtype) :: tempd1
    real(kind=realtype) :: tempd0
    real(kind=realtype) :: x1d
    real(kind=realtype) :: min3d
    real(kind=realtype) :: y2d
    real(kind=realtype) :: temp39
    real(kind=realtype) :: temp38
    real(kind=realtype) :: temp37
    real(kind=realtype) :: temp36
    real(kind=realtype) :: temp35
    real(kind=realtype) :: temp34
    real(kind=realtype) :: temp33
    real(kind=realtype) :: temp32
    real(kind=realtype) :: temp31
    real(kind=realtype) :: temp30
    real(kind=realtype) :: abs0
    real(kind=realtype) :: temp
    real(kind=realtype) :: temp9
    real(kind=realtype) :: temp8
    real(kind=realtype) :: min2d
    real(kind=realtype) :: tempd19
    real(kind=realtype) :: temp7
    real(kind=realtype) :: tempd18
    real(kind=realtype) :: y3
    real(kind=realtype) :: temp6
    real(kind=realtype) :: tempd17
    real(kind=realtype) :: y2
    real(kind=realtype) :: x3d
    real(kind=realtype) :: temp5
    real(kind=realtype) :: tempd16
    real(kind=realtype) :: y1
    real(kind=realtype) :: temp4
    real(kind=realtype) :: y1d
    real(kind=realtype) :: tempd15
    if (rfil .ge. 0.) then
        abs0 = rfil
    else
        abs0 = -rfil
    end if
    ! check if rfil == 0. if so, the dissipative flux needs not to
    ! be computed.
    if (abs0 .lt. thresholdreal) then
        radid = 0.0_8
        radjd = 0.0_8
        radkd = 0.0_8
    else
        ! determine the variables used to compute the switch.
        ! for the inviscid case this is the pressure; for the viscous
        ! case it is the entropy.
        select case  (equations)
        case (eulerequations)
            ! inviscid case. pressure switch is based on the pressure.
            ! also set the value of sslim. to be fully consistent this
            ! must have the dimension of pressure and it is therefore
            ! set to a fraction of the free stream value.
            sslim = 0.001_realtype*pinfcorr
            ! copy the pressure in ss. only need the entries used in the
            ! discretization, i.e. not including the corner halo's, but we'll
            ! just copy all anyway.
            ss = p
            call pushcontrol2b(1)
        case (nsequations, ransequations)
            !===============================================================
            ! viscous case. pressure switch is based on the entropy.
            ! also set the value of sslim. to be fully consistent this
            ! must have the dimension of entropy and it is therefore
            ! set to a fraction of the free stream value.
            sslim = 0.001_realtype*pinfcorr/rhoinf**gammainf
            ! store the entropy in ss. see above.
            do ii=0,(ib+1)*(jb+1)*(kb+1)-1
                i = mod(ii, ib + 1)
                j = mod(ii/(ib+1), jb + 1)
                k = ii/((ib+1)*(jb+1))
                ss(i, j, k) = p(i, j, k)/w(i, j, k, irho)**gamma(i, j, k)
            end do
            call pushcontrol2b(0)
        case default
            call pushcontrol2b(2)
        end select
        ! compute the pressure sensor for each cell, in each direction:
        do ii=0,ie*je*ke-1
            i = mod(ii, ie) + 1
            j = mod(ii/ie, je) + 1
            k = ii/(ie*je) + 1
            x1 = (ss(i+1, j, k)-two*ss(i, j, k)+ss(i-1, j, k))/(ss(i+1, j, k&
            &         )+two*ss(i, j, k)+ss(i-1, j, k)+sslim)
            if (x1 .ge. 0.) then
                dss(i, j, k, 1) = x1
            else
                dss(i, j, k, 1) = -x1
            end if
            x2 = (ss(i, j+1, k)-two*ss(i, j, k)+ss(i, j-1, k))/(ss(i, j+1, k&
            &         )+two*ss(i, j, k)+ss(i, j-1, k)+sslim)
            if (x2 .ge. 0.) then
                dss(i, j, k, 2) = x2
            else
                dss(i, j, k, 2) = -x2
            end if
            x3 = (ss(i, j, k+1)-two*ss(i, j, k)+ss(i, j, k-1))/(ss(i, j, k+1&
            &         )+two*ss(i, j, k)+ss(i, j, k-1)+sslim)
            if (x3 .ge. 0.) then
                dss(i, j, k, 3) = x3
            else
                dss(i, j, k, 3) = -x3
            end if
        end do
        ! set a couple of constants for the scheme.
        fis2 = rfil*vis2
        fis4 = rfil*vis4
        sfil = one - rfil
        ! initialize the dissipative residual to a certain times,
        ! possibly zero, the previously stored value. owned cells
        ! only, because the halo values do not matter.
        radkd = 0.0_8
        dssd = 0.0_8
        do ii=0,nx*ny*kl-1
            i = mod(ii, nx) + 2
            j = mod(ii/nx, ny) + 2
            k = ii/(nx*ny) + 1
            ! compute the dissipation coefficients for this face.
            ppor = zero
            if (pork(i, j, k) .eq. normalflux) ppor = half
            rrad = ppor*(radk(i, j, k)+radk(i, j, k+1))
            if (dss(i, j, k, 3) .lt. dss(i, j, k+1, 3)) then
                y3 = dss(i, j, k+1, 3)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                y3 = dss(i, j, k, 3)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            if (dssmax .gt. y3) then
                min3 = y3
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                min3 = dssmax
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            dis2 = fis2*rrad*min3
            arg1 = fis4*rrad
            dis4 = mydim(arg1, dis2)
            ! compute and scatter the dissipative flux.
            ! density. store it in the mass flow of the
            ! appropriate sliding mesh interface.
            ddw1 = w(i, j, k+1, irho) - w(i, j, k, irho)
            ! x-momentum.
            ddw2 = w(i, j, k+1, ivx)*w(i, j, k+1, irho) - w(i, j, k, ivx)*w(&
            &         i, j, k, irho)
            ! y-momentum.
            ddw3 = w(i, j, k+1, ivy)*w(i, j, k+1, irho) - w(i, j, k, ivy)*w(&
            &         i, j, k, irho)
            ! z-momentum.
            ddw4 = w(i, j, k+1, ivz)*w(i, j, k+1, irho) - w(i, j, k, ivz)*w(&
            &         i, j, k, irho)
            ! energy.
            ddw5 = w(i, j, k+1, irhoe) + p(i, j, k+1) - (w(i, j, k, irhoe)+p&
            &         (i, j, k))
            fsd = fwd(i, j, k+1, irhoe) - fwd(i, j, k, irhoe)
            tempd15 = -(dis4*fsd)
            dis2d = ddw5*fsd
            ddw5d = dis2*fsd - three*tempd15
            dis4d = -((w(i, j, k+2, irhoe)+p(i, j, k+2)-w(i, j, k-1, irhoe)-&
            &         p(i, j, k-1)-three*ddw5)*fsd)
            wd(i, j, k+2, irhoe) = wd(i, j, k+2, irhoe) + tempd15
            pd(i, j, k+2) = pd(i, j, k+2) + tempd15
            wd(i, j, k-1, irhoe) = wd(i, j, k-1, irhoe) - tempd15
            pd(i, j, k-1) = pd(i, j, k-1) - tempd15
            wd(i, j, k+1, irhoe) = wd(i, j, k+1, irhoe) + ddw5d
            pd(i, j, k+1) = pd(i, j, k+1) + ddw5d
            wd(i, j, k, irhoe) = wd(i, j, k, irhoe) - ddw5d
            pd(i, j, k) = pd(i, j, k) - ddw5d
            fsd = fwd(i, j, k+1, imz) - fwd(i, j, k, imz)
            temp40 = w(i, j, k-1, irho)
            temp39 = w(i, j, k-1, ivz)
            temp38 = w(i, j, k+2, irho)
            temp37 = w(i, j, k+2, ivz)
            tempd16 = -(dis4*fsd)
            dis2d = dis2d + ddw4*fsd
            ddw4d = dis2*fsd - three*tempd16
            dis4d = dis4d - (temp37*temp38-temp39*temp40-three*ddw4)*fsd
            wd(i, j, k+2, ivz) = wd(i, j, k+2, ivz) + temp38*tempd16
            wd(i, j, k+2, irho) = wd(i, j, k+2, irho) + temp37*tempd16
            wd(i, j, k-1, ivz) = wd(i, j, k-1, ivz) - temp40*tempd16
            wd(i, j, k-1, irho) = wd(i, j, k-1, irho) - temp39*tempd16
            wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) + w(i, j, k+1, irho)*&
            &         ddw4d
            wd(i, j, k+1, irho) = wd(i, j, k+1, irho) + w(i, j, k+1, ivz)*&
            &         ddw4d
            wd(i, j, k, ivz) = wd(i, j, k, ivz) - w(i, j, k, irho)*ddw4d
            wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, ivz)*ddw4d
            fsd = fwd(i, j, k+1, imy) - fwd(i, j, k, imy)
            temp36 = w(i, j, k-1, irho)
            temp35 = w(i, j, k-1, ivy)
            temp34 = w(i, j, k+2, irho)
            temp33 = w(i, j, k+2, ivy)
            tempd17 = -(dis4*fsd)
            dis2d = dis2d + ddw3*fsd
            ddw3d = dis2*fsd - three*tempd17
            dis4d = dis4d - (temp33*temp34-temp35*temp36-three*ddw3)*fsd
            wd(i, j, k+2, ivy) = wd(i, j, k+2, ivy) + temp34*tempd17
            wd(i, j, k+2, irho) = wd(i, j, k+2, irho) + temp33*tempd17
            wd(i, j, k-1, ivy) = wd(i, j, k-1, ivy) - temp36*tempd17
            wd(i, j, k-1, irho) = wd(i, j, k-1, irho) - temp35*tempd17
            wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) + w(i, j, k+1, irho)*&
            &         ddw3d
            wd(i, j, k+1, irho) = wd(i, j, k+1, irho) + w(i, j, k+1, ivy)*&
            &         ddw3d
            wd(i, j, k, ivy) = wd(i, j, k, ivy) - w(i, j, k, irho)*ddw3d
            wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, ivy)*ddw3d
            fsd = fwd(i, j, k+1, imx) - fwd(i, j, k, imx)
            temp32 = w(i, j, k-1, irho)
            temp31 = w(i, j, k-1, ivx)
            temp30 = w(i, j, k+2, irho)
            temp29 = w(i, j, k+2, ivx)
            tempd18 = -(dis4*fsd)
            dis2d = dis2d + ddw2*fsd
            ddw2d = dis2*fsd - three*tempd18
            dis4d = dis4d - (temp29*temp30-temp31*temp32-three*ddw2)*fsd
            wd(i, j, k+2, ivx) = wd(i, j, k+2, ivx) + temp30*tempd18
            wd(i, j, k+2, irho) = wd(i, j, k+2, irho) + temp29*tempd18
            wd(i, j, k-1, ivx) = wd(i, j, k-1, ivx) - temp32*tempd18
            wd(i, j, k-1, irho) = wd(i, j, k-1, irho) - temp31*tempd18
            wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) + w(i, j, k+1, irho)*&
            &         ddw2d
            wd(i, j, k+1, irho) = wd(i, j, k+1, irho) + w(i, j, k+1, ivx)*&
            &         ddw2d
            wd(i, j, k, ivx) = wd(i, j, k, ivx) - w(i, j, k, irho)*ddw2d
            wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, ivx)*ddw2d
            fsd = fwd(i, j, k+1, irho) - fwd(i, j, k, irho)
            tempd19 = -(dis4*fsd)
            dis2d = dis2d + ddw1*fsd
            ddw1d = dis2*fsd - three*tempd19
            dis4d = dis4d - (w(i, j, k+2, irho)-w(i, j, k-1, irho)-three*&
            &         ddw1)*fsd
            wd(i, j, k+2, irho) = wd(i, j, k+2, irho) + tempd19
            wd(i, j, k-1, irho) = wd(i, j, k-1, irho) - tempd19
            wd(i, j, k+1, irho) = wd(i, j, k+1, irho) + ddw1d
            wd(i, j, k, irho) = wd(i, j, k, irho) - ddw1d
            call mydim_fast_b(arg1, arg1d, dis2, dis2d, dis4d)
            rradd = fis2*min3*dis2d + fis4*arg1d
            min3d = fis2*rrad*dis2d
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                y3d = min3d
            else
                y3d = 0.0_8
            end if
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                dssd(i, j, k+1, 3) = dssd(i, j, k+1, 3) + y3d
            else
                dssd(i, j, k, 3) = dssd(i, j, k, 3) + y3d
            end if
            radkd(i, j, k) = radkd(i, j, k) + ppor*rradd
            radkd(i, j, k+1) = radkd(i, j, k+1) + ppor*rradd
        end do
        radjd = 0.0_8
        do ii=0,nx*jl*nz-1
            i = mod(ii, nx) + 2
            j = mod(ii/nx, jl) + 1
            k = ii/(nx*jl) + 2
            ! compute the dissipation coefficients for this face.
            ppor = zero
            if (porj(i, j, k) .eq. normalflux) ppor = half
            rrad = ppor*(radj(i, j, k)+radj(i, j+1, k))
            if (dss(i, j, k, 2) .lt. dss(i, j+1, k, 2)) then
                y2 = dss(i, j+1, k, 2)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                y2 = dss(i, j, k, 2)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            if (dssmax .gt. y2) then
                min2 = y2
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                min2 = dssmax
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            dis2 = fis2*rrad*min2
            arg1 = fis4*rrad
            dis4 = mydim(arg1, dis2)
            ! compute and scatter the dissipative flux.
            ! density. store it in the mass flow of the
            ! appropriate sliding mesh interface.
            ddw1 = w(i, j+1, k, irho) - w(i, j, k, irho)
            ! x-momentum.
            ddw2 = w(i, j+1, k, ivx)*w(i, j+1, k, irho) - w(i, j, k, ivx)*w(&
            &         i, j, k, irho)
            ! y-momentum.
            ddw3 = w(i, j+1, k, ivy)*w(i, j+1, k, irho) - w(i, j, k, ivy)*w(&
            &         i, j, k, irho)
            ! z-momentum.
            ddw4 = w(i, j+1, k, ivz)*w(i, j+1, k, irho) - w(i, j, k, ivz)*w(&
            &         i, j, k, irho)
            ! energy.
            ddw5 = w(i, j+1, k, irhoe) + p(i, j+1, k) - (w(i, j, k, irhoe)+p&
            &         (i, j, k))
            fsd = fwd(i, j+1, k, irhoe) - fwd(i, j, k, irhoe)
            tempd10 = -(dis4*fsd)
            dis2d = ddw5*fsd
            ddw5d = dis2*fsd - three*tempd10
            dis4d = -((w(i, j+2, k, irhoe)+p(i, j+2, k)-w(i, j-1, k, irhoe)-&
            &         p(i, j-1, k)-three*ddw5)*fsd)
            wd(i, j+2, k, irhoe) = wd(i, j+2, k, irhoe) + tempd10
            pd(i, j+2, k) = pd(i, j+2, k) + tempd10
            wd(i, j-1, k, irhoe) = wd(i, j-1, k, irhoe) - tempd10
            pd(i, j-1, k) = pd(i, j-1, k) - tempd10
            wd(i, j+1, k, irhoe) = wd(i, j+1, k, irhoe) + ddw5d
            pd(i, j+1, k) = pd(i, j+1, k) + ddw5d
            wd(i, j, k, irhoe) = wd(i, j, k, irhoe) - ddw5d
            pd(i, j, k) = pd(i, j, k) - ddw5d
            fsd = fwd(i, j+1, k, imz) - fwd(i, j, k, imz)
            temp28 = w(i, j-1, k, irho)
            temp27 = w(i, j-1, k, ivz)
            temp26 = w(i, j+2, k, irho)
            temp25 = w(i, j+2, k, ivz)
            tempd11 = -(dis4*fsd)
            dis2d = dis2d + ddw4*fsd
            ddw4d = dis2*fsd - three*tempd11
            dis4d = dis4d - (temp25*temp26-temp27*temp28-three*ddw4)*fsd
            wd(i, j+2, k, ivz) = wd(i, j+2, k, ivz) + temp26*tempd11
            wd(i, j+2, k, irho) = wd(i, j+2, k, irho) + temp25*tempd11
            wd(i, j-1, k, ivz) = wd(i, j-1, k, ivz) - temp28*tempd11
            wd(i, j-1, k, irho) = wd(i, j-1, k, irho) - temp27*tempd11
            wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) + w(i, j+1, k, irho)*&
            &         ddw4d
            wd(i, j+1, k, irho) = wd(i, j+1, k, irho) + w(i, j+1, k, ivz)*&
            &         ddw4d
            wd(i, j, k, ivz) = wd(i, j, k, ivz) - w(i, j, k, irho)*ddw4d
            wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, ivz)*ddw4d
            fsd = fwd(i, j+1, k, imy) - fwd(i, j, k, imy)
            temp24 = w(i, j-1, k, irho)
            temp23 = w(i, j-1, k, ivy)
            temp22 = w(i, j+2, k, irho)
            temp21 = w(i, j+2, k, ivy)
            tempd12 = -(dis4*fsd)
            dis2d = dis2d + ddw3*fsd
            ddw3d = dis2*fsd - three*tempd12
            dis4d = dis4d - (temp21*temp22-temp23*temp24-three*ddw3)*fsd
            wd(i, j+2, k, ivy) = wd(i, j+2, k, ivy) + temp22*tempd12
            wd(i, j+2, k, irho) = wd(i, j+2, k, irho) + temp21*tempd12
            wd(i, j-1, k, ivy) = wd(i, j-1, k, ivy) - temp24*tempd12
            wd(i, j-1, k, irho) = wd(i, j-1, k, irho) - temp23*tempd12
            wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) + w(i, j+1, k, irho)*&
            &         ddw3d
            wd(i, j+1, k, irho) = wd(i, j+1, k, irho) + w(i, j+1, k, ivy)*&
            &         ddw3d
            wd(i, j, k, ivy) = wd(i, j, k, ivy) - w(i, j, k, irho)*ddw3d
            wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, ivy)*ddw3d
            fsd = fwd(i, j+1, k, imx) - fwd(i, j, k, imx)
            temp20 = w(i, j-1, k, irho)
            temp19 = w(i, j-1, k, ivx)
            temp18 = w(i, j+2, k, irho)
            temp17 = w(i, j+2, k, ivx)
            tempd13 = -(dis4*fsd)
            dis2d = dis2d + ddw2*fsd
            ddw2d = dis2*fsd - three*tempd13
            dis4d = dis4d - (temp17*temp18-temp19*temp20-three*ddw2)*fsd
            wd(i, j+2, k, ivx) = wd(i, j+2, k, ivx) + temp18*tempd13
            wd(i, j+2, k, irho) = wd(i, j+2, k, irho) + temp17*tempd13
            wd(i, j-1, k, ivx) = wd(i, j-1, k, ivx) - temp20*tempd13
            wd(i, j-1, k, irho) = wd(i, j-1, k, irho) - temp19*tempd13
            wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) + w(i, j+1, k, irho)*&
            &         ddw2d
            wd(i, j+1, k, irho) = wd(i, j+1, k, irho) + w(i, j+1, k, ivx)*&
            &         ddw2d
            wd(i, j, k, ivx) = wd(i, j, k, ivx) - w(i, j, k, irho)*ddw2d
            wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, ivx)*ddw2d
            fsd = fwd(i, j+1, k, irho) - fwd(i, j, k, irho)
            tempd14 = -(dis4*fsd)
            dis2d = dis2d + ddw1*fsd
            ddw1d = dis2*fsd - three*tempd14
            dis4d = dis4d - (w(i, j+2, k, irho)-w(i, j-1, k, irho)-three*&
            &         ddw1)*fsd
            wd(i, j+2, k, irho) = wd(i, j+2, k, irho) + tempd14
            wd(i, j-1, k, irho) = wd(i, j-1, k, irho) - tempd14
            wd(i, j+1, k, irho) = wd(i, j+1, k, irho) + ddw1d
            wd(i, j, k, irho) = wd(i, j, k, irho) - ddw1d
            call mydim_fast_b(arg1, arg1d, dis2, dis2d, dis4d)
            rradd = fis2*min2*dis2d + fis4*arg1d
            min2d = fis2*rrad*dis2d
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                y2d = min2d
            else
                y2d = 0.0_8
            end if
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                dssd(i, j+1, k, 2) = dssd(i, j+1, k, 2) + y2d
            else
                dssd(i, j, k, 2) = dssd(i, j, k, 2) + y2d
            end if
            radjd(i, j, k) = radjd(i, j, k) + ppor*rradd
            radjd(i, j+1, k) = radjd(i, j+1, k) + ppor*rradd
        end do
        radid = 0.0_8
        do ii=0,il*ny*nz-1
            i = mod(ii, il) + 1
            j = mod(ii/il, ny) + 2
            k = ii/(il*ny) + 2
            ! compute the dissipation coefficients for this face.
            ppor = zero
            if (pori(i, j, k) .eq. normalflux) ppor = half
            rrad = ppor*(radi(i, j, k)+radi(i+1, j, k))
            if (dss(i, j, k, 1) .lt. dss(i+1, j, k, 1)) then
                y1 = dss(i+1, j, k, 1)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                y1 = dss(i, j, k, 1)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            if (dssmax .gt. y1) then
                min1 = y1
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                min1 = dssmax
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            dis2 = fis2*rrad*min1
            arg1 = fis4*rrad
            dis4 = mydim(arg1, dis2)
            ! compute and scatter the dissipative flux.
            ! density. store it in the mass flow of the
            ! appropriate sliding mesh interface.
            ddw1 = w(i+1, j, k, irho) - w(i, j, k, irho)
            ! x-momentum.
            ddw2 = w(i+1, j, k, ivx)*w(i+1, j, k, irho) - w(i, j, k, ivx)*w(&
            &         i, j, k, irho)
            ! y-momentum.
            ddw3 = w(i+1, j, k, ivy)*w(i+1, j, k, irho) - w(i, j, k, ivy)*w(&
            &         i, j, k, irho)
            ! z-momentum.
            ddw4 = w(i+1, j, k, ivz)*w(i+1, j, k, irho) - w(i, j, k, ivz)*w(&
            &         i, j, k, irho)
            ! energy.
            ddw5 = w(i+1, j, k, irhoe) + p(i+1, j, k) - (w(i, j, k, irhoe)+p&
            &         (i, j, k))
            fsd = fwd(i+1, j, k, irhoe) - fwd(i, j, k, irhoe)
            tempd5 = -(dis4*fsd)
            dis2d = ddw5*fsd
            ddw5d = dis2*fsd - three*tempd5
            dis4d = -((w(i+2, j, k, irhoe)+p(i+2, j, k)-w(i-1, j, k, irhoe)-&
            &         p(i-1, j, k)-three*ddw5)*fsd)
            wd(i+2, j, k, irhoe) = wd(i+2, j, k, irhoe) + tempd5
            pd(i+2, j, k) = pd(i+2, j, k) + tempd5
            wd(i-1, j, k, irhoe) = wd(i-1, j, k, irhoe) - tempd5
            pd(i-1, j, k) = pd(i-1, j, k) - tempd5
            wd(i+1, j, k, irhoe) = wd(i+1, j, k, irhoe) + ddw5d
            pd(i+1, j, k) = pd(i+1, j, k) + ddw5d
            wd(i, j, k, irhoe) = wd(i, j, k, irhoe) - ddw5d
            pd(i, j, k) = pd(i, j, k) - ddw5d
            fsd = fwd(i+1, j, k, imz) - fwd(i, j, k, imz)
            temp16 = w(i-1, j, k, irho)
            temp15 = w(i-1, j, k, ivz)
            temp14 = w(i+2, j, k, irho)
            temp13 = w(i+2, j, k, ivz)
            tempd6 = -(dis4*fsd)
            dis2d = dis2d + ddw4*fsd
            ddw4d = dis2*fsd - three*tempd6
            dis4d = dis4d - (temp13*temp14-temp15*temp16-three*ddw4)*fsd
            wd(i+2, j, k, ivz) = wd(i+2, j, k, ivz) + temp14*tempd6
            wd(i+2, j, k, irho) = wd(i+2, j, k, irho) + temp13*tempd6
            wd(i-1, j, k, ivz) = wd(i-1, j, k, ivz) - temp16*tempd6
            wd(i-1, j, k, irho) = wd(i-1, j, k, irho) - temp15*tempd6
            wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) + w(i+1, j, k, irho)*&
            &         ddw4d
            wd(i+1, j, k, irho) = wd(i+1, j, k, irho) + w(i+1, j, k, ivz)*&
            &         ddw4d
            wd(i, j, k, ivz) = wd(i, j, k, ivz) - w(i, j, k, irho)*ddw4d
            wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, ivz)*ddw4d
            fsd = fwd(i+1, j, k, imy) - fwd(i, j, k, imy)
            temp12 = w(i-1, j, k, irho)
            temp11 = w(i-1, j, k, ivy)
            temp10 = w(i+2, j, k, irho)
            temp9 = w(i+2, j, k, ivy)
            tempd7 = -(dis4*fsd)
            dis2d = dis2d + ddw3*fsd
            ddw3d = dis2*fsd - three*tempd7
            dis4d = dis4d - (temp9*temp10-temp11*temp12-three*ddw3)*fsd
            wd(i+2, j, k, ivy) = wd(i+2, j, k, ivy) + temp10*tempd7
            wd(i+2, j, k, irho) = wd(i+2, j, k, irho) + temp9*tempd7
            wd(i-1, j, k, ivy) = wd(i-1, j, k, ivy) - temp12*tempd7
            wd(i-1, j, k, irho) = wd(i-1, j, k, irho) - temp11*tempd7
            wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) + w(i+1, j, k, irho)*&
            &         ddw3d
            wd(i+1, j, k, irho) = wd(i+1, j, k, irho) + w(i+1, j, k, ivy)*&
            &         ddw3d
            wd(i, j, k, ivy) = wd(i, j, k, ivy) - w(i, j, k, irho)*ddw3d
            wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, ivy)*ddw3d
            fsd = fwd(i+1, j, k, imx) - fwd(i, j, k, imx)
            temp8 = w(i-1, j, k, irho)
            temp7 = w(i-1, j, k, ivx)
            temp6 = w(i+2, j, k, irho)
            temp5 = w(i+2, j, k, ivx)
            tempd8 = -(dis4*fsd)
            dis2d = dis2d + ddw2*fsd
            ddw2d = dis2*fsd - three*tempd8
            dis4d = dis4d - (temp5*temp6-temp7*temp8-three*ddw2)*fsd
            wd(i+2, j, k, ivx) = wd(i+2, j, k, ivx) + temp6*tempd8
            wd(i+2, j, k, irho) = wd(i+2, j, k, irho) + temp5*tempd8
            wd(i-1, j, k, ivx) = wd(i-1, j, k, ivx) - temp8*tempd8
            wd(i-1, j, k, irho) = wd(i-1, j, k, irho) - temp7*tempd8
            wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) + w(i+1, j, k, irho)*&
            &         ddw2d
            wd(i+1, j, k, irho) = wd(i+1, j, k, irho) + w(i+1, j, k, ivx)*&
            &         ddw2d
            wd(i, j, k, ivx) = wd(i, j, k, ivx) - w(i, j, k, irho)*ddw2d
            wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, ivx)*ddw2d
            fsd = fwd(i+1, j, k, irho) - fwd(i, j, k, irho)
            tempd9 = -(dis4*fsd)
            dis2d = dis2d + ddw1*fsd
            ddw1d = dis2*fsd - three*tempd9
            dis4d = dis4d - (w(i+2, j, k, irho)-w(i-1, j, k, irho)-three*&
            &         ddw1)*fsd
            wd(i+2, j, k, irho) = wd(i+2, j, k, irho) + tempd9
            wd(i-1, j, k, irho) = wd(i-1, j, k, irho) - tempd9
            wd(i+1, j, k, irho) = wd(i+1, j, k, irho) + ddw1d
            wd(i, j, k, irho) = wd(i, j, k, irho) - ddw1d
            call mydim_fast_b(arg1, arg1d, dis2, dis2d, dis4d)
            rradd = fis2*min1*dis2d + fis4*arg1d
            min1d = fis2*rrad*dis2d
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                y1d = min1d
            else
                y1d = 0.0_8
            end if
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                dssd(i+1, j, k, 1) = dssd(i+1, j, k, 1) + y1d
            else
                dssd(i, j, k, 1) = dssd(i, j, k, 1) + y1d
            end if
            radid(i, j, k) = radid(i, j, k) + ppor*rradd
            radid(i+1, j, k) = radid(i+1, j, k) + ppor*rradd
        end do
        fwd = sfil*fwd
        ssd = 0.0_8
        do ii=0,ie*je*ke-1
            i = mod(ii, ie) + 1
            j = mod(ii/ie, je) + 1
            k = ii/(ie*je) + 1
            x1 = (ss(i+1, j, k)-two*ss(i, j, k)+ss(i-1, j, k))/(ss(i+1, j, k&
            &         )+two*ss(i, j, k)+ss(i-1, j, k)+sslim)
            if (x1 .ge. 0.) then
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            x2 = (ss(i, j+1, k)-two*ss(i, j, k)+ss(i, j-1, k))/(ss(i, j+1, k&
            &         )+two*ss(i, j, k)+ss(i, j-1, k)+sslim)
            if (x2 .ge. 0.) then
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            x3 = (ss(i, j, k+1)-two*ss(i, j, k)+ss(i, j, k-1))/(ss(i, j, k+1&
            &         )+two*ss(i, j, k)+ss(i, j, k-1)+sslim)
            if (x3 .ge. 0.) then
                x3d = dssd(i, j, k, 3)
                dssd(i, j, k, 3) = 0.0_8
            else
                x3d = -dssd(i, j, k, 3)
                dssd(i, j, k, 3) = 0.0_8
            end if
            temp4 = sslim + ss(i, j, k+1) + two*ss(i, j, k) + ss(i, j, k-1)
            tempd3 = x3d/temp4
            tempd4 = -((ss(i, j, k+1)-two*ss(i, j, k)+ss(i, j, k-1))*tempd3/&
            &         temp4)
            ssd(i, j, k+1) = ssd(i, j, k+1) + tempd4 + tempd3
            ssd(i, j, k) = ssd(i, j, k) + two*tempd4 - two*tempd3
            ssd(i, j, k-1) = ssd(i, j, k-1) + tempd4 + tempd3
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                x2d = dssd(i, j, k, 2)
                dssd(i, j, k, 2) = 0.0_8
            else
                x2d = -dssd(i, j, k, 2)
                dssd(i, j, k, 2) = 0.0_8
            end if
            temp3 = sslim + ss(i, j+1, k) + two*ss(i, j, k) + ss(i, j-1, k)
            tempd1 = x2d/temp3
            tempd2 = -((ss(i, j+1, k)-two*ss(i, j, k)+ss(i, j-1, k))*tempd1/&
            &         temp3)
            ssd(i, j+1, k) = ssd(i, j+1, k) + tempd2 + tempd1
            ssd(i, j, k) = ssd(i, j, k) + two*tempd2 - two*tempd1
            ssd(i, j-1, k) = ssd(i, j-1, k) + tempd2 + tempd1
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                x1d = dssd(i, j, k, 1)
                dssd(i, j, k, 1) = 0.0_8
            else
                x1d = -dssd(i, j, k, 1)
                dssd(i, j, k, 1) = 0.0_8
            end if
            temp2 = sslim + ss(i+1, j, k) + two*ss(i, j, k) + ss(i-1, j, k)
            tempd = x1d/temp2
            tempd0 = -((ss(i+1, j, k)-two*ss(i, j, k)+ss(i-1, j, k))*tempd/&
            &         temp2)
            ssd(i+1, j, k) = ssd(i+1, j, k) + tempd0 + tempd
            ssd(i, j, k) = ssd(i, j, k) + two*tempd0 - two*tempd
            ssd(i-1, j, k) = ssd(i-1, j, k) + tempd0 + tempd
        end do
        call popcontrol2b(branch)
        if (branch .eq. 0) then
            do ii=0,(ib+1)*(jb+1)*(kb+1)-1
                i = mod(ii, ib + 1)
                j = mod(ii/(ib+1), jb + 1)
                k = ii/((ib+1)*(jb+1))
                temp1 = gamma(i, j, k)
                temp0 = w(i, j, k, irho)
                temp = temp0**temp1
                pd(i, j, k) = pd(i, j, k) + ssd(i, j, k)/temp
                if (.not.(temp0 .le. 0.0_8 .and. (temp1 .eq. 0.0_8 .or. temp1 &
                &             .ne. int(temp1)))) wd(i, j, k, irho) = wd(i, j, k, irho) -&
                &             p(i, j, k)*temp1*temp0**(temp1-1)*ssd(i, j, k)/temp**2
                ssd(i, j, k) = 0.0_8
            end do
        else if (branch .eq. 1) then
            pd = pd + ssd
        end if
    end if
end subroutine invisciddissfluxscalar_fast_b

!  differentiation of invisciddissfluxmatrix in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: *p *w *fw
!   with respect to varying inputs: *p *w *fw
!   rw status of diff variables: *p:incr *w:incr *fw:in-out
!   plus diff mem management of: p:in w:in fw:in
subroutine invisciddissfluxmatrix_fast_b()
    !
    !       invisciddissfluxmatrix computes the matrix artificial
    !       dissipation term. instead of the spectral radius, as used in
    !       the scalar dissipation scheme, the absolute value of the flux
    !       jacobian is used. this leads to a less diffusive and
    !       consequently more accurate scheme. it is assumed that the
    !       pointers in blockpointers already point to the correct block.
    !
    use constants
    ! use blockpointers, only : nx, ny, nz, il, jl, kl, ie, je, ke, ib, &
    ! &   jb, kb, w, wd, p, pd, pori, porj, pork, fw, fwd, gamma, si, sj, sk, &
    ! &   indfamilyi, indfamilyj, indfamilyk, spectralsol, addgridvelocities, &
    ! &   sfacei, sfacej, sfacek, factfamilyi, factfamilyj, factfamilyk
    use blockpointers, only : addgridvelocities
    use flowvarrefstate, only : pinfcorr
    use inputdiscretization, only : vis2, vis4
    use inputphysics, only : equations
    use iteration, only : rfil
    use cgnsgrid, only : massflowfamilydiss
    use utils_fast_b, only : getcorrectfork, mydim, mydim_fast_b
    implicit none
    !
    !      local parameters.
    !
    real(kind=realtype), parameter :: dpmax=0.25_realtype
    real(kind=realtype), parameter :: epsacoustic=0.25_realtype
    real(kind=realtype), parameter :: epsshear=0.025_realtype
    real(kind=realtype), parameter :: omega=0.5_realtype
    real(kind=realtype), parameter :: oneminomega=one-omega
    !
    !      local variables.
    !
    integer(kind=inttype) :: i, j, k, ind, ii
    real(kind=realtype) :: plim, sface
    real(kind=realtype) :: sfil, fis2, fis4
    real(kind=realtype) :: gammaavg, gm1, ovgm1, gm53
    real(kind=realtype) :: ppor, rrad, dis2, dis4
    real(kind=realtype) :: rradd, dis2d, dis4d
    real(kind=realtype) :: dp1, dp2, tmp, fs
    real(kind=realtype) :: fsd
    real(kind=realtype) :: ddw1, ddw2, ddw3, ddw4, ddw5, ddw6
    real(kind=realtype) :: ddw1d, ddw2d, ddw3d, ddw4d, ddw5d, ddw6d
    real(kind=realtype) :: dr, dru, drv, drw, dre, drk, sx, sy, sz
    real(kind=realtype) :: drd, drud, drvd, drwd, dred, drkd
    real(kind=realtype) :: uavg, vavg, wavg, a2avg, aavg, havg
    real(kind=realtype) :: uavgd, vavgd, wavgd, a2avgd, aavgd, havgd
    real(kind=realtype) :: alphaavg, unavg, ovaavg, ova2avg
    real(kind=realtype) :: alphaavgd, unavgd, ovaavgd, ova2avgd
    real(kind=realtype) :: kavg, lam1, lam2, lam3, area
    real(kind=realtype) :: kavgd, lam1d, lam2d, lam3d
    real(kind=realtype) :: abv1, abv2, abv3, abv4, abv5, abv6, abv7
    real(kind=realtype) :: abv1d, abv2d, abv3d, abv4d, abv5d, abv6d, &
    &   abv7d
    real(kind=realtype), dimension(ie, je, ke, 3) :: dss
    real(kind=realtype), dimension(ie, je, ke, 3) :: dssd
    logical :: correctfork
    intrinsic abs
    intrinsic mod
    intrinsic max
    intrinsic min
    intrinsic sqrt
    real(kind=realtype) :: arg1
    real(kind=realtype) :: arg1d
    integer :: branch
    real(kind=realtype) :: temp3
    real(kind=realtype) :: tempd14
    real(kind=realtype) :: temp29
    real(kind=realtype) :: temp2
    real(kind=realtype) :: tempd13
    real(kind=realtype) :: temp28
    real(kind=realtype) :: temp1
    real(kind=realtype) :: tempd12
    real(kind=realtype) :: temp27
    real(kind=realtype) :: max10d
    real(kind=realtype) :: temp0
    real(kind=realtype) :: tempd11
    real(kind=realtype) :: temp26
    real(kind=realtype) :: tempd10
    real(kind=realtype) :: temp25
    real(kind=realtype) :: temp24
    real(kind=realtype) :: temp23
    real(kind=realtype) :: temp22
    real(kind=realtype) :: temp21
    real(kind=realtype) :: temp20
    real(kind=realtype) :: abs1d
    real(kind=realtype) :: temp55
    real(kind=realtype) :: temp54
    real(kind=realtype) :: max2d
    real(kind=realtype) :: temp53
    real(kind=realtype) :: min3
    real(kind=realtype) :: temp52
    real(kind=realtype) :: min2
    real(kind=realtype) :: temp51
    real(kind=realtype) :: min1
    real(kind=realtype) :: temp50
    real(kind=realtype) :: abs4d
    real(kind=realtype) :: min1d
    real(kind=realtype) :: max8d
    real(kind=realtype) :: x3
    real(kind=realtype) :: x2
    real(kind=realtype) :: x2d
    real(kind=realtype) :: x1
    real(kind=realtype) :: temp19
    real(kind=realtype) :: temp18
    real(kind=realtype) :: temp17
    real(kind=realtype) :: temp16
    real(kind=realtype) :: temp15
    real(kind=realtype) :: tempd37
    real(kind=realtype) :: temp14
    real(kind=realtype) :: tempd36
    real(kind=realtype) :: temp13
    real(kind=realtype) :: y3d
    real(kind=realtype) :: tempd35
    real(kind=realtype) :: temp12
    real(kind=realtype) :: temp49
    real(kind=realtype) :: tempd34
    real(kind=realtype) :: temp11
    real(kind=realtype) :: temp48
    real(kind=realtype) :: tempd33
    real(kind=realtype) :: temp10
    real(kind=realtype) :: temp47
    real(kind=realtype) :: tempd32
    real(kind=realtype) :: max12d
    real(kind=realtype) :: temp46
    real(kind=realtype) :: tempd31
    real(kind=realtype) :: temp45
    real(kind=realtype) :: tempd30
    real(kind=realtype) :: temp44
    real(kind=realtype) :: temp43
    real(kind=realtype) :: temp42
    real(kind=realtype) :: temp41
    real(kind=realtype) :: temp40
    real(kind=realtype) :: abs3d
    real(kind=realtype) :: tempd9
    real(kind=realtype) :: tempd
    real(kind=realtype) :: tempd8
    real(kind=realtype) :: max4d
    real(kind=realtype) :: tempd7
    real(kind=realtype) :: tempd6
    real(kind=realtype) :: tempd5
    real(kind=realtype) :: tempd4
    real(kind=realtype) :: abs6d
    real(kind=realtype) :: tempd3
    real(kind=realtype) :: tempd2
    real(kind=realtype) :: tempd1
    real(kind=realtype) :: max7d
    real(kind=realtype) :: tempd0
    real(kind=realtype) :: x1d
    real(kind=realtype) :: min3d
    real(kind=realtype) :: tempd29
    real(kind=realtype) :: tempd28
    real(kind=realtype) :: tempd27
    real(kind=realtype) :: tempd26
    real(kind=realtype) :: y2d
    real(kind=realtype) :: tempd25
    real(kind=realtype) :: tempd24
    real(kind=realtype) :: temp39
    real(kind=realtype) :: tempd23
    real(kind=realtype) :: temp38
    real(kind=realtype) :: tempd22
    real(kind=realtype) :: temp37
    real(kind=realtype) :: max11d
    real(kind=realtype) :: tempd21
    real(kind=realtype) :: temp36
    real(kind=realtype) :: tempd20
    real(kind=realtype) :: temp35
    real(kind=realtype) :: temp34
    real(kind=realtype) :: abs6
    real(kind=realtype) :: temp33
    real(kind=realtype) :: abs5
    real(kind=realtype) :: temp32
    real(kind=realtype) :: abs4
    real(kind=realtype) :: temp31
    real(kind=realtype) :: abs3
    real(kind=realtype) :: temp30
    real(kind=realtype) :: abs2
    real(kind=realtype) :: abs2d
    real(kind=realtype) :: abs1
    real(kind=realtype) :: abs0
    real(kind=realtype) :: max3d
    real(kind=realtype) :: max9
    real(kind=realtype) :: abs5d
    real(kind=realtype) :: max8
    real(kind=realtype) :: max7
    real(kind=realtype) :: max6
    real(kind=realtype) :: max6d
    real(kind=realtype) :: max5
    real(kind=realtype) :: max4
    real(kind=realtype) :: temp
    real(kind=realtype) :: max3
    real(kind=realtype) :: max2
    real(kind=realtype) :: max1
    real(kind=realtype) :: max12
    real(kind=realtype) :: temp9
    real(kind=realtype) :: max11
    real(kind=realtype) :: temp8
    real(kind=realtype) :: min2d
    real(kind=realtype) :: tempd19
    real(kind=realtype) :: max10
    real(kind=realtype) :: temp7
    real(kind=realtype) :: tempd18
    real(kind=realtype) :: y3
    real(kind=realtype) :: temp6
    real(kind=realtype) :: tempd17
    real(kind=realtype) :: y2
    real(kind=realtype) :: x3d
    real(kind=realtype) :: temp5
    real(kind=realtype) :: tempd16
    real(kind=realtype) :: y1
    real(kind=realtype) :: y1d
    real(kind=realtype) :: temp4
    real(kind=realtype) :: tempd15
    if (rfil .ge. 0.) then
        abs0 = rfil
    else
        abs0 = -rfil
    end if
    ! check if rfil == 0. if so, the dissipative flux needs not to
    ! be computed.
    if (abs0 .ge. thresholdreal) then
        ! set the value of plim. to be fully consistent this must have
        ! the dimension of a pressure. therefore a fraction of pinfcorr
        ! is used.
        plim = 0.001_realtype*pinfcorr
        ! determine whether or not the total energy must be corrected
        ! for the presence of the turbulent kinetic energy.
        correctfork = getcorrectfork()
        ! initialize sface to zero. this value will be used if the
        ! block is not moving.
        sface = zero
        ! set a couple of constants for the scheme.
        fis2 = rfil*vis2
        fis4 = rfil*vis4
        sfil = one - rfil
        ! initialize the dissipative residual to a certain times,
        ! possibly zero, the previously stored value.
        ! compute the pressure sensor for each cell, in each direction:
        do ii=0,ie*je*ke-1
            i = mod(ii, ie) + 1
            j = mod(ii/ie, je) + 1
            k = ii/(ie*je) + 1
            if (p(i+1, j, k) - p(i, j, k) .ge. 0.) then
                abs1 = p(i+1, j, k) - p(i, j, k)
            else
                abs1 = -(p(i+1, j, k)-p(i, j, k))
            end if
            if (p(i, j, k) - p(i-1, j, k) .ge. 0.) then
                abs4 = p(i, j, k) - p(i-1, j, k)
            else
                abs4 = -(p(i, j, k)-p(i-1, j, k))
            end if
            x1 = (p(i+1, j, k)-two*p(i, j, k)+p(i-1, j, k))/(omega*(p(i+1, j&
            &         , k)+two*p(i, j, k)+p(i-1, j, k))+oneminomega*(abs1+abs4)+plim&
            &         )
            if (x1 .ge. 0.) then
                dss(i, j, k, 1) = x1
            else
                dss(i, j, k, 1) = -x1
            end if
            if (p(i, j+1, k) - p(i, j, k) .ge. 0.) then
                abs2 = p(i, j+1, k) - p(i, j, k)
            else
                abs2 = -(p(i, j+1, k)-p(i, j, k))
            end if
            if (p(i, j, k) - p(i, j-1, k) .ge. 0.) then
                abs5 = p(i, j, k) - p(i, j-1, k)
            else
                abs5 = -(p(i, j, k)-p(i, j-1, k))
            end if
            x2 = (p(i, j+1, k)-two*p(i, j, k)+p(i, j-1, k))/(omega*(p(i, j+1&
            &         , k)+two*p(i, j, k)+p(i, j-1, k))+oneminomega*(abs2+abs5)+plim&
            &         )
            if (x2 .ge. 0.) then
                dss(i, j, k, 2) = x2
            else
                dss(i, j, k, 2) = -x2
            end if
            if (p(i, j, k+1) - p(i, j, k) .ge. 0.) then
                abs3 = p(i, j, k+1) - p(i, j, k)
            else
                abs3 = -(p(i, j, k+1)-p(i, j, k))
            end if
            if (p(i, j, k) - p(i, j, k-1) .ge. 0.) then
                abs6 = p(i, j, k) - p(i, j, k-1)
            else
                abs6 = -(p(i, j, k)-p(i, j, k-1))
            end if
            x3 = (p(i, j, k+1)-two*p(i, j, k)+p(i, j, k-1))/(omega*(p(i, j, &
            &         k+1)+two*p(i, j, k)+p(i, j, k-1))+oneminomega*(abs3+abs6)+plim&
            &         )
            if (x3 .ge. 0.) then
                dss(i, j, k, 3) = x3
            else
                dss(i, j, k, 3) = -x3
            end if
        end do
        !
        !       dissipative fluxes in the i-direction.
        !
        do ii=0,il*ny*nz-1
            i = mod(ii, il) + 1
            j = mod(ii/il, ny) + 2
            k = ii/(il*ny) + 2
            ! compute the dissipation coefficients for this face.
            ppor = zero
            if (pori(i, j, k) .eq. normalflux) ppor = one
            if (dss(i, j, k, 1) .lt. dss(i+1, j, k, 1)) then
                y1 = dss(i+1, j, k, 1)
            else
                y1 = dss(i, j, k, 1)
            end if
            if (dpmax .gt. y1) then
                min1 = y1
            else
                min1 = dpmax
            end if
            dis2 = ppor*fis2*min1
            arg1 = ppor*fis4
            dis4 = mydim(arg1, dis2)
            ! construct the vector of the first and third differences
            ! multiplied by the appropriate constants.
            ddw1 = w(i+1, j, k, irho) - w(i, j, k, irho)
            dr = dis2*ddw1 - dis4*(w(i+2, j, k, irho)-w(i-1, j, k, irho)-&
            &         three*ddw1)
            ddw2 = w(i+1, j, k, irho)*w(i+1, j, k, ivx) - w(i, j, k, irho)*w&
            &         (i, j, k, ivx)
            dru = dis2*ddw2 - dis4*(w(i+2, j, k, irho)*w(i+2, j, k, ivx)-w(i&
            &         -1, j, k, irho)*w(i-1, j, k, ivx)-three*ddw2)
            ddw3 = w(i+1, j, k, irho)*w(i+1, j, k, ivy) - w(i, j, k, irho)*w&
            &         (i, j, k, ivy)
            drv = dis2*ddw3 - dis4*(w(i+2, j, k, irho)*w(i+2, j, k, ivy)-w(i&
            &         -1, j, k, irho)*w(i-1, j, k, ivy)-three*ddw3)
            ddw4 = w(i+1, j, k, irho)*w(i+1, j, k, ivz) - w(i, j, k, irho)*w&
            &         (i, j, k, ivz)
            drw = dis2*ddw4 - dis4*(w(i+2, j, k, irho)*w(i+2, j, k, ivz)-w(i&
            &         -1, j, k, irho)*w(i-1, j, k, ivz)-three*ddw4)
            ddw5 = w(i+1, j, k, irhoe) - w(i, j, k, irhoe)
            dre = dis2*ddw5 - dis4*(w(i+2, j, k, irhoe)-w(i-1, j, k, irhoe)-&
            &         three*ddw5)
            ! in case a k-equation is present, compute the difference
            ! of rhok and store the average value of k. if not present,
            ! set both these values to zero, such that later on no
            ! decision needs to be made anymore.
            if (correctfork) then
                ddw6 = w(i+1, j, k, irho)*w(i+1, j, k, itu1) - w(i, j, k, irho&
                &           )*w(i, j, k, itu1)
                drk = dis2*ddw6 - dis4*(w(i+2, j, k, irho)*w(i+2, j, k, itu1)-&
                &           w(i-1, j, k, irho)*w(i-1, j, k, itu1)-three*ddw6)
                kavg = half*(w(i, j, k, itu1)+w(i+1, j, k, itu1))
            else
                drk = zero
                kavg = zero
            end if
            ! compute the average value of gamma and compute some
            ! expressions in which it occurs.
            gammaavg = half*(gamma(i+1, j, k)+gamma(i, j, k))
            gm1 = gammaavg - one
            ovgm1 = one/gm1
            gm53 = gammaavg - five*third
            ! compute the average state at the interface.
            uavg = half*(w(i+1, j, k, ivx)+w(i, j, k, ivx))
            vavg = half*(w(i+1, j, k, ivy)+w(i, j, k, ivy))
            wavg = half*(w(i+1, j, k, ivz)+w(i, j, k, ivz))
            a2avg = half*(gamma(i+1, j, k)*p(i+1, j, k)/w(i+1, j, k, irho)+&
            &         gamma(i, j, k)*p(i, j, k)/w(i, j, k, irho))
            area = sqrt(si(i, j, k, 1)**2 + si(i, j, k, 2)**2 + si(i, j, k, &
            &         3)**2)
            if (1.e-25_realtype .lt. area) then
                max1 = area
            else
                max1 = 1.e-25_realtype
            end if
            tmp = one/max1
            sx = si(i, j, k, 1)*tmp
            sy = si(i, j, k, 2)*tmp
            sz = si(i, j, k, 3)*tmp
            alphaavg = half*(uavg**2+vavg**2+wavg**2)
            havg = alphaavg + ovgm1*(a2avg-gm53*kavg)
            aavg = sqrt(a2avg)
            unavg = uavg*sx + vavg*sy + wavg*sz
            ovaavg = one/aavg
            ova2avg = one/a2avg
            ! the mesh velocity if the face is moving. it must be
            ! divided by the area to obtain a true velocity.
            if (addgridvelocities) sface = sfacei(i, j, k)*tmp
            if (unavg - sface + aavg .ge. 0.) then
                lam1 = unavg - sface + aavg
            else
                lam1 = -(unavg-sface+aavg)
            end if
            if (unavg - sface - aavg .ge. 0.) then
                lam2 = unavg - sface - aavg
            else
                lam2 = -(unavg-sface-aavg)
            end if
            if (unavg - sface .ge. 0.) then
                lam3 = unavg - sface
            else
                lam3 = -(unavg-sface)
            end if
            rrad = lam3 + aavg
            if (lam1 .lt. epsacoustic*rrad) then
                max2 = epsacoustic*rrad
            else
                max2 = lam1
            end if
            ! multiply the eigenvalues by the area to obtain
            ! the correct values for the dissipation term.
            lam1 = max2*area
            if (lam2 .lt. epsacoustic*rrad) then
                max3 = epsacoustic*rrad
            else
                max3 = lam2
            end if
            lam2 = max3*area
            if (lam3 .lt. epsshear*rrad) then
                max4 = epsshear*rrad
            else
                max4 = lam3
            end if
            lam3 = max4*area
            ! some abbreviations, which occur quite often in the
            ! dissipation terms.
            abv1 = half*(lam1+lam2)
            abv2 = half*(lam1-lam2)
            abv3 = abv1 - lam3
            abv4 = gm1*(alphaavg*dr-uavg*dru-vavg*drv-wavg*drw+dre) - gm53*&
            &         drk
            abv5 = sx*dru + sy*drv + sz*drw - unavg*dr
            abv6 = abv3*abv4*ova2avg + abv2*abv5*ovaavg
            abv7 = abv2*abv4*ovaavg + abv3*abv5
            ! compute and scatter the dissipative flux.
            ! density.
            fs = lam3*dr + abv6
            fw(i+1, j, k, irho) = fw(i+1, j, k, irho) + fs
            fw(i, j, k, irho) = fw(i, j, k, irho) - fs
            ! x-momentum.
            fs = lam3*dru + uavg*abv6 + sx*abv7
            fw(i+1, j, k, imx) = fw(i+1, j, k, imx) + fs
            fw(i, j, k, imx) = fw(i, j, k, imx) - fs
            ! y-momentum.
            fs = lam3*drv + vavg*abv6 + sy*abv7
            fw(i+1, j, k, imy) = fw(i+1, j, k, imy) + fs
            fw(i, j, k, imy) = fw(i, j, k, imy) - fs
            ! z-momentum.
            fs = lam3*drw + wavg*abv6 + sz*abv7
            fw(i+1, j, k, imz) = fw(i+1, j, k, imz) + fs
            fw(i, j, k, imz) = fw(i, j, k, imz) - fs
            ! energy.
            fs = lam3*dre + havg*abv6 + unavg*abv7
            fw(i+1, j, k, irhoe) = fw(i+1, j, k, irhoe) + fs
            fw(i, j, k, irhoe) = fw(i, j, k, irhoe) - fs
        end do
        !
        !       dissipative fluxes in the j-direction.
        !
        do ii=0,nx*jl*nz-1
            i = mod(ii, nx) + 2
            j = mod(ii/nx, jl) + 1
            k = ii/(nx*jl) + 2
            ! compute the dissipation coefficients for this face.
            ppor = zero
            if (porj(i, j, k) .eq. normalflux) ppor = one
            if (dss(i, j, k, 2) .lt. dss(i, j+1, k, 2)) then
                y2 = dss(i, j+1, k, 2)
            else
                y2 = dss(i, j, k, 2)
            end if
            if (dpmax .gt. y2) then
                min2 = y2
            else
                min2 = dpmax
            end if
            dis2 = ppor*fis2*min2
            arg1 = ppor*fis4
            dis4 = mydim(arg1, dis2)
            ! construct the vector of the first and third differences
            ! multiplied by the appropriate constants.
            ddw1 = w(i, j+1, k, irho) - w(i, j, k, irho)
            dr = dis2*ddw1 - dis4*(w(i, j+2, k, irho)-w(i, j-1, k, irho)-&
            &         three*ddw1)
            ddw2 = w(i, j+1, k, irho)*w(i, j+1, k, ivx) - w(i, j, k, irho)*w&
            &         (i, j, k, ivx)
            dru = dis2*ddw2 - dis4*(w(i, j+2, k, irho)*w(i, j+2, k, ivx)-w(i&
            &         , j-1, k, irho)*w(i, j-1, k, ivx)-three*ddw2)
            ddw3 = w(i, j+1, k, irho)*w(i, j+1, k, ivy) - w(i, j, k, irho)*w&
            &         (i, j, k, ivy)
            drv = dis2*ddw3 - dis4*(w(i, j+2, k, irho)*w(i, j+2, k, ivy)-w(i&
            &         , j-1, k, irho)*w(i, j-1, k, ivy)-three*ddw3)
            ddw4 = w(i, j+1, k, irho)*w(i, j+1, k, ivz) - w(i, j, k, irho)*w&
            &         (i, j, k, ivz)
            drw = dis2*ddw4 - dis4*(w(i, j+2, k, irho)*w(i, j+2, k, ivz)-w(i&
            &         , j-1, k, irho)*w(i, j-1, k, ivz)-three*ddw4)
            ddw5 = w(i, j+1, k, irhoe) - w(i, j, k, irhoe)
            dre = dis2*ddw5 - dis4*(w(i, j+2, k, irhoe)-w(i, j-1, k, irhoe)-&
            &         three*ddw5)
            ! in case a k-equation is present, compute the difference
            ! of rhok and store the average value of k. if not present,
            ! set both these values to zero, such that later on no
            ! decision needs to be made anymore.
            if (correctfork) then
                ddw6 = w(i, j+1, k, irho)*w(i, j+1, k, itu1) - w(i, j, k, irho&
                &           )*w(i, j, k, itu1)
                drk = dis2*ddw6 - dis4*(w(i, j+2, k, irho)*w(i, j+2, k, itu1)-&
                &           w(i, j-1, k, irho)*w(i, j-1, k, itu1)-three*ddw6)
                kavg = half*(w(i, j, k, itu1)+w(i, j+1, k, itu1))
            else
                drk = zero
                kavg = zero
            end if
            ! compute the average value of gamma and compute some
            ! expressions in which it occurs.
            gammaavg = half*(gamma(i, j+1, k)+gamma(i, j, k))
            gm1 = gammaavg - one
            ovgm1 = one/gm1
            gm53 = gammaavg - five*third
            ! compute the average state at the interface.
            uavg = half*(w(i, j+1, k, ivx)+w(i, j, k, ivx))
            vavg = half*(w(i, j+1, k, ivy)+w(i, j, k, ivy))
            wavg = half*(w(i, j+1, k, ivz)+w(i, j, k, ivz))
            a2avg = half*(gamma(i, j+1, k)*p(i, j+1, k)/w(i, j+1, k, irho)+&
            &         gamma(i, j, k)*p(i, j, k)/w(i, j, k, irho))
            area = sqrt(sj(i, j, k, 1)**2 + sj(i, j, k, 2)**2 + sj(i, j, k, &
            &         3)**2)
            if (1.e-25_realtype .lt. area) then
                max5 = area
            else
                max5 = 1.e-25_realtype
            end if
            tmp = one/max5
            sx = sj(i, j, k, 1)*tmp
            sy = sj(i, j, k, 2)*tmp
            sz = sj(i, j, k, 3)*tmp
            alphaavg = half*(uavg**2+vavg**2+wavg**2)
            havg = alphaavg + ovgm1*(a2avg-gm53*kavg)
            aavg = sqrt(a2avg)
            unavg = uavg*sx + vavg*sy + wavg*sz
            ovaavg = one/aavg
            ova2avg = one/a2avg
            ! the mesh velocity if the face is moving. it must be
            ! divided by the area to obtain a true velocity.
            if (addgridvelocities) sface = sfacej(i, j, k)*tmp
            if (unavg - sface + aavg .ge. 0.) then
                lam1 = unavg - sface + aavg
            else
                lam1 = -(unavg-sface+aavg)
            end if
            if (unavg - sface - aavg .ge. 0.) then
                lam2 = unavg - sface - aavg
            else
                lam2 = -(unavg-sface-aavg)
            end if
            if (unavg - sface .ge. 0.) then
                lam3 = unavg - sface
            else
                lam3 = -(unavg-sface)
            end if
            rrad = lam3 + aavg
            if (lam1 .lt. epsacoustic*rrad) then
                max6 = epsacoustic*rrad
            else
                max6 = lam1
            end if
            ! multiply the eigenvalues by the area to obtain
            ! the correct values for the dissipation term.
            lam1 = max6*area
            if (lam2 .lt. epsacoustic*rrad) then
                max7 = epsacoustic*rrad
            else
                max7 = lam2
            end if
            lam2 = max7*area
            if (lam3 .lt. epsshear*rrad) then
                max8 = epsshear*rrad
            else
                max8 = lam3
            end if
            lam3 = max8*area
            ! some abbreviations, which occur quite often in the
            ! dissipation terms.
            abv1 = half*(lam1+lam2)
            abv2 = half*(lam1-lam2)
            abv3 = abv1 - lam3
            abv4 = gm1*(alphaavg*dr-uavg*dru-vavg*drv-wavg*drw+dre) - gm53*&
            &         drk
            abv5 = sx*dru + sy*drv + sz*drw - unavg*dr
            abv6 = abv3*abv4*ova2avg + abv2*abv5*ovaavg
            abv7 = abv2*abv4*ovaavg + abv3*abv5
            ! compute and scatter the dissipative flux.
            ! density.
            fs = lam3*dr + abv6
            fw(i, j+1, k, irho) = fw(i, j+1, k, irho) + fs
            fw(i, j, k, irho) = fw(i, j, k, irho) - fs
            ! x-momentum.
            fs = lam3*dru + uavg*abv6 + sx*abv7
            fw(i, j+1, k, imx) = fw(i, j+1, k, imx) + fs
            fw(i, j, k, imx) = fw(i, j, k, imx) - fs
            ! y-momentum.
            fs = lam3*drv + vavg*abv6 + sy*abv7
            fw(i, j+1, k, imy) = fw(i, j+1, k, imy) + fs
            fw(i, j, k, imy) = fw(i, j, k, imy) - fs
            ! z-momentum.
            fs = lam3*drw + wavg*abv6 + sz*abv7
            fw(i, j+1, k, imz) = fw(i, j+1, k, imz) + fs
            fw(i, j, k, imz) = fw(i, j, k, imz) - fs
            ! energy.
            fs = lam3*dre + havg*abv6 + unavg*abv7
            fw(i, j+1, k, irhoe) = fw(i, j+1, k, irhoe) + fs
            fw(i, j, k, irhoe) = fw(i, j, k, irhoe) - fs
        end do
        dssd = 0.0_8
        do ii=0,nx*ny*kl-1
            i = mod(ii, nx) + 2
            j = mod(ii/nx, ny) + 2
            k = ii/(nx*ny) + 1
            ! compute the dissipation coefficients for this face.
            ppor = zero
            if (pork(i, j, k) .eq. normalflux) ppor = one
            if (dss(i, j, k, 3) .lt. dss(i, j, k+1, 3)) then
                y3 = dss(i, j, k+1, 3)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                y3 = dss(i, j, k, 3)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            if (dpmax .gt. y3) then
                min3 = y3
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                min3 = dpmax
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            dis2 = ppor*fis2*min3
            arg1 = ppor*fis4
            dis4 = mydim(arg1, dis2)
            ! construct the vector of the first and third differences
            ! multiplied by the appropriate constants.
            ddw1 = w(i, j, k+1, irho) - w(i, j, k, irho)
            dr = dis2*ddw1 - dis4*(w(i, j, k+2, irho)-w(i, j, k-1, irho)-&
            &         three*ddw1)
            ddw2 = w(i, j, k+1, irho)*w(i, j, k+1, ivx) - w(i, j, k, irho)*w&
            &         (i, j, k, ivx)
            dru = dis2*ddw2 - dis4*(w(i, j, k+2, irho)*w(i, j, k+2, ivx)-w(i&
            &         , j, k-1, irho)*w(i, j, k-1, ivx)-three*ddw2)
            ddw3 = w(i, j, k+1, irho)*w(i, j, k+1, ivy) - w(i, j, k, irho)*w&
            &         (i, j, k, ivy)
            drv = dis2*ddw3 - dis4*(w(i, j, k+2, irho)*w(i, j, k+2, ivy)-w(i&
            &         , j, k-1, irho)*w(i, j, k-1, ivy)-three*ddw3)
            ddw4 = w(i, j, k+1, irho)*w(i, j, k+1, ivz) - w(i, j, k, irho)*w&
            &         (i, j, k, ivz)
            drw = dis2*ddw4 - dis4*(w(i, j, k+2, irho)*w(i, j, k+2, ivz)-w(i&
            &         , j, k-1, irho)*w(i, j, k-1, ivz)-three*ddw4)
            ddw5 = w(i, j, k+1, irhoe) - w(i, j, k, irhoe)
            dre = dis2*ddw5 - dis4*(w(i, j, k+2, irhoe)-w(i, j, k-1, irhoe)-&
            &         three*ddw5)
            ! in case a k-equation is present, compute the difference
            ! of rhok and store the average value of k. if not present,
            ! set both these values to zero, such that later on no
            ! decision needs to be made anymore.
            if (correctfork) then
                ddw6 = w(i, j, k+1, irho)*w(i, j, k+1, itu1) - w(i, j, k, irho&
                &           )*w(i, j, k, itu1)
                drk = dis2*ddw6 - dis4*(w(i, j, k+2, irho)*w(i, j, k+2, itu1)-&
                &           w(i, j, k-1, irho)*w(i, j, k-1, itu1)-three*ddw6)
                kavg = half*(w(i, j, k+1, itu1)+w(i, j, k, itu1))
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            else
                drk = zero
                kavg = zero
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            end if
            ! compute the average value of gamma and compute some
            ! expressions in which it occurs.
            gammaavg = half*(gamma(i, j, k+1)+gamma(i, j, k))
            gm1 = gammaavg - one
            ovgm1 = one/gm1
            gm53 = gammaavg - five*third
            ! compute the average state at the interface.
            uavg = half*(w(i, j, k+1, ivx)+w(i, j, k, ivx))
            vavg = half*(w(i, j, k+1, ivy)+w(i, j, k, ivy))
            wavg = half*(w(i, j, k+1, ivz)+w(i, j, k, ivz))
            a2avg = half*(gamma(i, j, k+1)*p(i, j, k+1)/w(i, j, k+1, irho)+&
            &         gamma(i, j, k)*p(i, j, k)/w(i, j, k, irho))
            area = sqrt(sk(i, j, k, 1)**2 + sk(i, j, k, 2)**2 + sk(i, j, k, &
            &         3)**2)
            if (1.e-25_realtype .lt. area) then
                max9 = area
            else
                max9 = 1.e-25_realtype
            end if
            tmp = one/max9
            sx = sk(i, j, k, 1)*tmp
            sy = sk(i, j, k, 2)*tmp
            sz = sk(i, j, k, 3)*tmp
            alphaavg = half*(uavg**2+vavg**2+wavg**2)
            havg = alphaavg + ovgm1*(a2avg-gm53*kavg)
            aavg = sqrt(a2avg)
            unavg = uavg*sx + vavg*sy + wavg*sz
            ovaavg = one/aavg
            ova2avg = one/a2avg
            ! the mesh velocity if the face is moving. it must be
            ! divided by the area to obtain a true velocity.
            if (addgridvelocities) sface = sfacek(i, j, k)*tmp
            if (unavg - sface + aavg .ge. 0.) then
                lam1 = unavg - sface + aavg
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                lam1 = -(unavg-sface+aavg)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            if (unavg - sface - aavg .ge. 0.) then
                lam2 = unavg - sface - aavg
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                lam2 = -(unavg-sface-aavg)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            if (unavg - sface .ge. 0.) then
                lam3 = unavg - sface
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                lam3 = -(unavg-sface)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            rrad = lam3 + aavg
            if (lam1 .lt. epsacoustic*rrad) then
                max10 = epsacoustic*rrad
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                max10 = lam1
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            ! multiply the eigenvalues by the area to obtain
            ! the correct values for the dissipation term.
            lam1 = max10*area
            if (lam2 .lt. epsacoustic*rrad) then
                max11 = epsacoustic*rrad
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                max11 = lam2
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            lam2 = max11*area
            if (lam3 .lt. epsshear*rrad) then
                max12 = epsshear*rrad
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                max12 = lam3
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            lam3 = max12*area
            ! some abbreviations, which occur quite often in the
            ! dissipation terms.
            abv1 = half*(lam1+lam2)
            abv2 = half*(lam1-lam2)
            abv3 = abv1 - lam3
            abv4 = gm1*(alphaavg*dr-uavg*dru-vavg*drv-wavg*drw+dre) - gm53*&
            &         drk
            abv5 = sx*dru + sy*drv + sz*drw - unavg*dr
            abv6 = abv3*abv4*ova2avg + abv2*abv5*ovaavg
            abv7 = abv2*abv4*ovaavg + abv3*abv5
            ! compute and scatter the dissipative flux.
            ! density.
            ! x-momentum.
            ! y-momentum.
            ! z-momentum.
            ! energy.
            fsd = fwd(i, j, k+1, irhoe) - fwd(i, j, k, irhoe)
            lam3d = dre*fsd
            dred = lam3*fsd
            havgd = abv6*fsd
            abv6d = havg*fsd
            unavgd = abv7*fsd
            abv7d = unavg*fsd
            fsd = fwd(i, j, k+1, imz) - fwd(i, j, k, imz)
            lam3d = lam3d + drw*fsd
            drwd = lam3*fsd
            wavgd = abv6*fsd
            abv6d = abv6d + wavg*fsd
            abv7d = abv7d + sz*fsd
            fsd = fwd(i, j, k+1, imy) - fwd(i, j, k, imy)
            lam3d = lam3d + drv*fsd
            drvd = lam3*fsd
            vavgd = abv6*fsd
            abv6d = abv6d + vavg*fsd
            abv7d = abv7d + sy*fsd
            fsd = fwd(i, j, k+1, imx) - fwd(i, j, k, imx)
            lam3d = lam3d + dru*fsd
            drud = lam3*fsd
            uavgd = abv6*fsd
            abv6d = abv6d + uavg*fsd
            abv7d = abv7d + sx*fsd
            fsd = fwd(i, j, k+1, irho) - fwd(i, j, k, irho)
            abv6d = abv6d + fsd
            abv2d = ovaavg*abv5*abv6d + ovaavg*abv4*abv7d
            abv4d = ova2avg*abv3*abv6d + ovaavg*abv2*abv7d
            ovaavgd = abv2*abv5*abv6d + abv2*abv4*abv7d
            abv3d = ova2avg*abv4*abv6d + abv5*abv7d
            lam3d = lam3d + dr*fsd - abv3d
            abv5d = ovaavg*abv2*abv6d + abv3*abv7d
            ova2avgd = abv3*abv4*abv6d
            unavgd = unavgd - dr*abv5d
            tempd37 = gm1*abv4d
            drd = alphaavg*tempd37 - unavg*abv5d + lam3*fsd
            drud = drud + sx*abv5d - uavg*tempd37
            drvd = drvd + sy*abv5d - vavg*tempd37
            drwd = drwd + sz*abv5d - wavg*tempd37
            alphaavgd = dr*tempd37
            uavgd = uavgd - dru*tempd37
            vavgd = vavgd - drv*tempd37
            dred = dred + tempd37
            wavgd = wavgd - drw*tempd37
            drkd = -(gm53*abv4d)
            abv1d = abv3d
            lam1d = half*abv1d + half*abv2d
            lam2d = half*abv1d - half*abv2d
            max12d = area*lam3d
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                rradd = epsshear*max12d
                lam3d = 0.0_8
            else
                lam3d = max12d
                rradd = 0.0_8
            end if
            max11d = area*lam2d
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                rradd = rradd + epsacoustic*max11d
                lam2d = 0.0_8
            else
                lam2d = max11d
            end if
            max10d = area*lam1d
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                rradd = rradd + epsacoustic*max10d
                lam1d = 0.0_8
            else
                lam1d = max10d
            end if
            lam3d = lam3d + rradd
            aavgd = rradd
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                unavgd = unavgd + lam3d
            else
                unavgd = unavgd - lam3d
            end if
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                unavgd = unavgd + lam2d
                aavgd = aavgd - lam2d
            else
                aavgd = aavgd + lam2d
                unavgd = unavgd - lam2d
            end if
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                unavgd = unavgd + lam1d
                aavgd = aavgd + lam1d
            else
                unavgd = unavgd - lam1d
                aavgd = aavgd - lam1d
            end if
            alphaavgd = alphaavgd + havgd
            tempd36 = half*alphaavgd
            aavgd = aavgd - one*ovaavgd/aavg**2
            if (a2avg .eq. 0.0_8) then
                a2avgd = ovgm1*havgd - one*ova2avgd/a2avg**2
            else
                a2avgd = aavgd/(2.0*sqrt(a2avg)) + ovgm1*havgd - one*ova2avgd/&
                &           a2avg**2
            end if
            uavgd = uavgd + 2*uavg*tempd36 + sx*unavgd
            vavgd = vavgd + 2*vavg*tempd36 + sy*unavgd
            wavgd = wavgd + 2*wavg*tempd36 + sz*unavgd
            kavgd = -(ovgm1*gm53*havgd)
            temp55 = w(i, j, k, irho)
            temp54 = w(i, j, k+1, irho)
            tempd34 = gamma(i, j, k+1)*half*a2avgd/temp54
            tempd35 = gamma(i, j, k)*half*a2avgd/temp55
            pd(i, j, k+1) = pd(i, j, k+1) + tempd34
            wd(i, j, k+1, irho) = wd(i, j, k+1, irho) - p(i, j, k+1)*tempd34&
            &         /temp54
            pd(i, j, k) = pd(i, j, k) + tempd35
            wd(i, j, k, irho) = wd(i, j, k, irho) - p(i, j, k)*tempd35/&
            &         temp55
            wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) + half*wavgd
            wd(i, j, k, ivz) = wd(i, j, k, ivz) + half*wavgd
            wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) + half*vavgd
            wd(i, j, k, ivy) = wd(i, j, k, ivy) + half*vavgd
            wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) + half*uavgd
            wd(i, j, k, ivx) = wd(i, j, k, ivx) + half*uavgd
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                dis2d = 0.0_8
                dis4d = 0.0_8
            else
                wd(i, j, k+1, itu1) = wd(i, j, k+1, itu1) + half*kavgd
                wd(i, j, k, itu1) = wd(i, j, k, itu1) + half*kavgd
                temp53 = w(i, j, k-1, itu1)
                temp52 = w(i, j, k-1, irho)
                temp51 = w(i, j, k+2, itu1)
                temp50 = w(i, j, k+2, irho)
                tempd33 = -(dis4*drkd)
                dis2d = ddw6*drkd
                ddw6d = dis2*drkd - three*tempd33
                dis4d = -((temp50*temp51-temp52*temp53-three*ddw6)*drkd)
                wd(i, j, k+2, irho) = wd(i, j, k+2, irho) + temp51*tempd33
                wd(i, j, k+2, itu1) = wd(i, j, k+2, itu1) + temp50*tempd33
                wd(i, j, k-1, irho) = wd(i, j, k-1, irho) - temp53*tempd33
                wd(i, j, k-1, itu1) = wd(i, j, k-1, itu1) - temp52*tempd33
                wd(i, j, k+1, irho) = wd(i, j, k+1, irho) + w(i, j, k+1, itu1)&
                &           *ddw6d
                wd(i, j, k+1, itu1) = wd(i, j, k+1, itu1) + w(i, j, k+1, irho)&
                &           *ddw6d
                wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, itu1)*ddw6d
                wd(i, j, k, itu1) = wd(i, j, k, itu1) - w(i, j, k, irho)*ddw6d
            end if
            temp38 = w(i, j, k+2, irho)
            temp39 = w(i, j, k+2, ivx)
            temp40 = w(i, j, k-1, irho)
            temp41 = w(i, j, k-1, ivx)
            temp42 = w(i, j, k+2, irho)
            temp43 = w(i, j, k+2, ivy)
            temp44 = w(i, j, k-1, irho)
            temp45 = w(i, j, k-1, ivy)
            temp46 = w(i, j, k+2, irho)
            temp47 = w(i, j, k+2, ivz)
            temp48 = w(i, j, k-1, irho)
            temp49 = w(i, j, k-1, ivz)
            tempd28 = -(dis4*dred)
            dis2d = dis2d + ddw4*drwd + ddw2*drud + ddw1*drd + ddw3*drvd + &
            &         ddw5*dred
            ddw5d = dis2*dred - three*tempd28
            dis4d = dis4d - (temp46*temp47-temp48*temp49-three*ddw4)*drwd - &
            &         (temp38*temp39-temp40*temp41-three*ddw2)*drud - (w(i, j, k+2, &
            &         irho)-w(i, j, k-1, irho)-three*ddw1)*drd - (temp42*temp43-&
            &         temp44*temp45-three*ddw3)*drvd - (w(i, j, k+2, irhoe)-w(i, j, &
            &         k-1, irhoe)-three*ddw5)*dred
            wd(i, j, k+2, irhoe) = wd(i, j, k+2, irhoe) + tempd28
            wd(i, j, k-1, irhoe) = wd(i, j, k-1, irhoe) - tempd28
            wd(i, j, k+1, irhoe) = wd(i, j, k+1, irhoe) + ddw5d
            wd(i, j, k, irhoe) = wd(i, j, k, irhoe) - ddw5d
            tempd29 = -(dis4*drwd)
            ddw4d = dis2*drwd - three*tempd29
            wd(i, j, k+2, irho) = wd(i, j, k+2, irho) + temp47*tempd29
            wd(i, j, k+2, ivz) = wd(i, j, k+2, ivz) + temp46*tempd29
            wd(i, j, k-1, irho) = wd(i, j, k-1, irho) - temp49*tempd29
            wd(i, j, k-1, ivz) = wd(i, j, k-1, ivz) - temp48*tempd29
            wd(i, j, k+1, irho) = wd(i, j, k+1, irho) + w(i, j, k+1, ivz)*&
            &         ddw4d
            wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) + w(i, j, k+1, irho)*&
            &         ddw4d
            wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, ivz)*ddw4d
            wd(i, j, k, ivz) = wd(i, j, k, ivz) - w(i, j, k, irho)*ddw4d
            tempd30 = -(dis4*drvd)
            ddw3d = dis2*drvd - three*tempd30
            wd(i, j, k+2, irho) = wd(i, j, k+2, irho) + temp43*tempd30
            wd(i, j, k+2, ivy) = wd(i, j, k+2, ivy) + temp42*tempd30
            wd(i, j, k-1, irho) = wd(i, j, k-1, irho) - temp45*tempd30
            wd(i, j, k-1, ivy) = wd(i, j, k-1, ivy) - temp44*tempd30
            wd(i, j, k+1, irho) = wd(i, j, k+1, irho) + w(i, j, k+1, ivy)*&
            &         ddw3d
            wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) + w(i, j, k+1, irho)*&
            &         ddw3d
            wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, ivy)*ddw3d
            wd(i, j, k, ivy) = wd(i, j, k, ivy) - w(i, j, k, irho)*ddw3d
            tempd31 = -(dis4*drud)
            ddw2d = dis2*drud - three*tempd31
            wd(i, j, k+2, irho) = wd(i, j, k+2, irho) + temp39*tempd31
            wd(i, j, k+2, ivx) = wd(i, j, k+2, ivx) + temp38*tempd31
            wd(i, j, k-1, irho) = wd(i, j, k-1, irho) - temp41*tempd31
            wd(i, j, k-1, ivx) = wd(i, j, k-1, ivx) - temp40*tempd31
            wd(i, j, k+1, irho) = wd(i, j, k+1, irho) + w(i, j, k+1, ivx)*&
            &         ddw2d
            wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) + w(i, j, k+1, irho)*&
            &         ddw2d
            wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, ivx)*ddw2d
            wd(i, j, k, ivx) = wd(i, j, k, ivx) - w(i, j, k, irho)*ddw2d
            tempd32 = -(dis4*drd)
            ddw1d = dis2*drd - three*tempd32
            wd(i, j, k+2, irho) = wd(i, j, k+2, irho) + tempd32
            wd(i, j, k-1, irho) = wd(i, j, k-1, irho) - tempd32
            wd(i, j, k+1, irho) = wd(i, j, k+1, irho) + ddw1d
            wd(i, j, k, irho) = wd(i, j, k, irho) - ddw1d
            call mydim_fast_b(arg1, arg1d, dis2, dis2d, dis4d)
            min3d = ppor*fis2*dis2d
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                y3d = min3d
            else
                y3d = 0.0_8
            end if
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                dssd(i, j, k+1, 3) = dssd(i, j, k+1, 3) + y3d
            else
                dssd(i, j, k, 3) = dssd(i, j, k, 3) + y3d
            end if
        end do
        do ii=0,nx*jl*nz-1
            i = mod(ii, nx) + 2
            j = mod(ii/nx, jl) + 1
            k = ii/(nx*jl) + 2
            ! compute the dissipation coefficients for this face.
            ppor = zero
            if (porj(i, j, k) .eq. normalflux) ppor = one
            if (dss(i, j, k, 2) .lt. dss(i, j+1, k, 2)) then
                y2 = dss(i, j+1, k, 2)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                y2 = dss(i, j, k, 2)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            if (dpmax .gt. y2) then
                min2 = y2
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                min2 = dpmax
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            dis2 = ppor*fis2*min2
            arg1 = ppor*fis4
            dis4 = mydim(arg1, dis2)
            ! construct the vector of the first and third differences
            ! multiplied by the appropriate constants.
            ddw1 = w(i, j+1, k, irho) - w(i, j, k, irho)
            dr = dis2*ddw1 - dis4*(w(i, j+2, k, irho)-w(i, j-1, k, irho)-&
            &         three*ddw1)
            ddw2 = w(i, j+1, k, irho)*w(i, j+1, k, ivx) - w(i, j, k, irho)*w&
            &         (i, j, k, ivx)
            dru = dis2*ddw2 - dis4*(w(i, j+2, k, irho)*w(i, j+2, k, ivx)-w(i&
            &         , j-1, k, irho)*w(i, j-1, k, ivx)-three*ddw2)
            ddw3 = w(i, j+1, k, irho)*w(i, j+1, k, ivy) - w(i, j, k, irho)*w&
            &         (i, j, k, ivy)
            drv = dis2*ddw3 - dis4*(w(i, j+2, k, irho)*w(i, j+2, k, ivy)-w(i&
            &         , j-1, k, irho)*w(i, j-1, k, ivy)-three*ddw3)
            ddw4 = w(i, j+1, k, irho)*w(i, j+1, k, ivz) - w(i, j, k, irho)*w&
            &         (i, j, k, ivz)
            drw = dis2*ddw4 - dis4*(w(i, j+2, k, irho)*w(i, j+2, k, ivz)-w(i&
            &         , j-1, k, irho)*w(i, j-1, k, ivz)-three*ddw4)
            ddw5 = w(i, j+1, k, irhoe) - w(i, j, k, irhoe)
            dre = dis2*ddw5 - dis4*(w(i, j+2, k, irhoe)-w(i, j-1, k, irhoe)-&
            &         three*ddw5)
            ! in case a k-equation is present, compute the difference
            ! of rhok and store the average value of k. if not present,
            ! set both these values to zero, such that later on no
            ! decision needs to be made anymore.
            if (correctfork) then
                ddw6 = w(i, j+1, k, irho)*w(i, j+1, k, itu1) - w(i, j, k, irho&
                &           )*w(i, j, k, itu1)
                drk = dis2*ddw6 - dis4*(w(i, j+2, k, irho)*w(i, j+2, k, itu1)-&
                &           w(i, j-1, k, irho)*w(i, j-1, k, itu1)-three*ddw6)
                kavg = half*(w(i, j, k, itu1)+w(i, j+1, k, itu1))
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            else
                drk = zero
                kavg = zero
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            end if
            ! compute the average value of gamma and compute some
            ! expressions in which it occurs.
            gammaavg = half*(gamma(i, j+1, k)+gamma(i, j, k))
            gm1 = gammaavg - one
            ovgm1 = one/gm1
            gm53 = gammaavg - five*third
            ! compute the average state at the interface.
            uavg = half*(w(i, j+1, k, ivx)+w(i, j, k, ivx))
            vavg = half*(w(i, j+1, k, ivy)+w(i, j, k, ivy))
            wavg = half*(w(i, j+1, k, ivz)+w(i, j, k, ivz))
            a2avg = half*(gamma(i, j+1, k)*p(i, j+1, k)/w(i, j+1, k, irho)+&
            &         gamma(i, j, k)*p(i, j, k)/w(i, j, k, irho))
            area = sqrt(sj(i, j, k, 1)**2 + sj(i, j, k, 2)**2 + sj(i, j, k, &
            &         3)**2)
            if (1.e-25_realtype .lt. area) then
                max5 = area
            else
                max5 = 1.e-25_realtype
            end if
            tmp = one/max5
            sx = sj(i, j, k, 1)*tmp
            sy = sj(i, j, k, 2)*tmp
            sz = sj(i, j, k, 3)*tmp
            alphaavg = half*(uavg**2+vavg**2+wavg**2)
            havg = alphaavg + ovgm1*(a2avg-gm53*kavg)
            aavg = sqrt(a2avg)
            unavg = uavg*sx + vavg*sy + wavg*sz
            ovaavg = one/aavg
            ova2avg = one/a2avg
            ! the mesh velocity if the face is moving. it must be
            ! divided by the area to obtain a true velocity.
            if (addgridvelocities) sface = sfacej(i, j, k)*tmp
            if (unavg - sface + aavg .ge. 0.) then
                lam1 = unavg - sface + aavg
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                lam1 = -(unavg-sface+aavg)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            if (unavg - sface - aavg .ge. 0.) then
                lam2 = unavg - sface - aavg
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                lam2 = -(unavg-sface-aavg)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            if (unavg - sface .ge. 0.) then
                lam3 = unavg - sface
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                lam3 = -(unavg-sface)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            rrad = lam3 + aavg
            if (lam1 .lt. epsacoustic*rrad) then
                max6 = epsacoustic*rrad
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                max6 = lam1
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            ! multiply the eigenvalues by the area to obtain
            ! the correct values for the dissipation term.
            lam1 = max6*area
            if (lam2 .lt. epsacoustic*rrad) then
                max7 = epsacoustic*rrad
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                max7 = lam2
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            lam2 = max7*area
            if (lam3 .lt. epsshear*rrad) then
                max8 = epsshear*rrad
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                max8 = lam3
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            lam3 = max8*area
            ! some abbreviations, which occur quite often in the
            ! dissipation terms.
            abv1 = half*(lam1+lam2)
            abv2 = half*(lam1-lam2)
            abv3 = abv1 - lam3
            abv4 = gm1*(alphaavg*dr-uavg*dru-vavg*drv-wavg*drw+dre) - gm53*&
            &         drk
            abv5 = sx*dru + sy*drv + sz*drw - unavg*dr
            abv6 = abv3*abv4*ova2avg + abv2*abv5*ovaavg
            abv7 = abv2*abv4*ovaavg + abv3*abv5
            ! compute and scatter the dissipative flux.
            ! density.
            ! x-momentum.
            ! y-momentum.
            ! z-momentum.
            ! energy.
            fsd = fwd(i, j+1, k, irhoe) - fwd(i, j, k, irhoe)
            lam3d = dre*fsd
            dred = lam3*fsd
            havgd = abv6*fsd
            abv6d = havg*fsd
            unavgd = abv7*fsd
            abv7d = unavg*fsd
            fsd = fwd(i, j+1, k, imz) - fwd(i, j, k, imz)
            lam3d = lam3d + drw*fsd
            drwd = lam3*fsd
            wavgd = abv6*fsd
            abv6d = abv6d + wavg*fsd
            abv7d = abv7d + sz*fsd
            fsd = fwd(i, j+1, k, imy) - fwd(i, j, k, imy)
            lam3d = lam3d + drv*fsd
            drvd = lam3*fsd
            vavgd = abv6*fsd
            abv6d = abv6d + vavg*fsd
            abv7d = abv7d + sy*fsd
            fsd = fwd(i, j+1, k, imx) - fwd(i, j, k, imx)
            lam3d = lam3d + dru*fsd
            drud = lam3*fsd
            uavgd = abv6*fsd
            abv6d = abv6d + uavg*fsd
            abv7d = abv7d + sx*fsd
            fsd = fwd(i, j+1, k, irho) - fwd(i, j, k, irho)
            abv6d = abv6d + fsd
            abv2d = ovaavg*abv5*abv6d + ovaavg*abv4*abv7d
            abv4d = ova2avg*abv3*abv6d + ovaavg*abv2*abv7d
            ovaavgd = abv2*abv5*abv6d + abv2*abv4*abv7d
            abv3d = ova2avg*abv4*abv6d + abv5*abv7d
            lam3d = lam3d + dr*fsd - abv3d
            abv5d = ovaavg*abv2*abv6d + abv3*abv7d
            ova2avgd = abv3*abv4*abv6d
            unavgd = unavgd - dr*abv5d
            tempd27 = gm1*abv4d
            drd = alphaavg*tempd27 - unavg*abv5d + lam3*fsd
            drud = drud + sx*abv5d - uavg*tempd27
            drvd = drvd + sy*abv5d - vavg*tempd27
            drwd = drwd + sz*abv5d - wavg*tempd27
            alphaavgd = dr*tempd27
            uavgd = uavgd - dru*tempd27
            vavgd = vavgd - drv*tempd27
            dred = dred + tempd27
            wavgd = wavgd - drw*tempd27
            drkd = -(gm53*abv4d)
            abv1d = abv3d
            lam1d = half*abv1d + half*abv2d
            lam2d = half*abv1d - half*abv2d
            max8d = area*lam3d
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                rradd = epsshear*max8d
                lam3d = 0.0_8
            else
                lam3d = max8d
                rradd = 0.0_8
            end if
            max7d = area*lam2d
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                rradd = rradd + epsacoustic*max7d
                lam2d = 0.0_8
            else
                lam2d = max7d
            end if
            max6d = area*lam1d
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                rradd = rradd + epsacoustic*max6d
                lam1d = 0.0_8
            else
                lam1d = max6d
            end if
            lam3d = lam3d + rradd
            aavgd = rradd
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                unavgd = unavgd + lam3d
            else
                unavgd = unavgd - lam3d
            end if
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                unavgd = unavgd + lam2d
                aavgd = aavgd - lam2d
            else
                aavgd = aavgd + lam2d
                unavgd = unavgd - lam2d
            end if
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                unavgd = unavgd + lam1d
                aavgd = aavgd + lam1d
            else
                unavgd = unavgd - lam1d
                aavgd = aavgd - lam1d
            end if
            alphaavgd = alphaavgd + havgd
            tempd26 = half*alphaavgd
            aavgd = aavgd - one*ovaavgd/aavg**2
            if (a2avg .eq. 0.0_8) then
                a2avgd = ovgm1*havgd - one*ova2avgd/a2avg**2
            else
                a2avgd = aavgd/(2.0*sqrt(a2avg)) + ovgm1*havgd - one*ova2avgd/&
                &           a2avg**2
            end if
            uavgd = uavgd + 2*uavg*tempd26 + sx*unavgd
            vavgd = vavgd + 2*vavg*tempd26 + sy*unavgd
            wavgd = wavgd + 2*wavg*tempd26 + sz*unavgd
            kavgd = -(ovgm1*gm53*havgd)
            temp37 = w(i, j, k, irho)
            temp36 = w(i, j+1, k, irho)
            tempd24 = gamma(i, j+1, k)*half*a2avgd/temp36
            tempd25 = gamma(i, j, k)*half*a2avgd/temp37
            pd(i, j+1, k) = pd(i, j+1, k) + tempd24
            wd(i, j+1, k, irho) = wd(i, j+1, k, irho) - p(i, j+1, k)*tempd24&
            &         /temp36
            pd(i, j, k) = pd(i, j, k) + tempd25
            wd(i, j, k, irho) = wd(i, j, k, irho) - p(i, j, k)*tempd25/&
            &         temp37
            wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) + half*wavgd
            wd(i, j, k, ivz) = wd(i, j, k, ivz) + half*wavgd
            wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) + half*vavgd
            wd(i, j, k, ivy) = wd(i, j, k, ivy) + half*vavgd
            wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) + half*uavgd
            wd(i, j, k, ivx) = wd(i, j, k, ivx) + half*uavgd
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                dis2d = 0.0_8
                dis4d = 0.0_8
            else
                wd(i, j, k, itu1) = wd(i, j, k, itu1) + half*kavgd
                wd(i, j+1, k, itu1) = wd(i, j+1, k, itu1) + half*kavgd
                temp35 = w(i, j-1, k, itu1)
                temp34 = w(i, j-1, k, irho)
                temp33 = w(i, j+2, k, itu1)
                temp32 = w(i, j+2, k, irho)
                tempd23 = -(dis4*drkd)
                dis2d = ddw6*drkd
                ddw6d = dis2*drkd - three*tempd23
                dis4d = -((temp32*temp33-temp34*temp35-three*ddw6)*drkd)
                wd(i, j+2, k, irho) = wd(i, j+2, k, irho) + temp33*tempd23
                wd(i, j+2, k, itu1) = wd(i, j+2, k, itu1) + temp32*tempd23
                wd(i, j-1, k, irho) = wd(i, j-1, k, irho) - temp35*tempd23
                wd(i, j-1, k, itu1) = wd(i, j-1, k, itu1) - temp34*tempd23
                wd(i, j+1, k, irho) = wd(i, j+1, k, irho) + w(i, j+1, k, itu1)&
                &           *ddw6d
                wd(i, j+1, k, itu1) = wd(i, j+1, k, itu1) + w(i, j+1, k, irho)&
                &           *ddw6d
                wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, itu1)*ddw6d
                wd(i, j, k, itu1) = wd(i, j, k, itu1) - w(i, j, k, irho)*ddw6d
            end if
            temp20 = w(i, j+2, k, irho)
            temp21 = w(i, j+2, k, ivx)
            temp22 = w(i, j-1, k, irho)
            temp23 = w(i, j-1, k, ivx)
            temp24 = w(i, j+2, k, irho)
            temp25 = w(i, j+2, k, ivy)
            temp26 = w(i, j-1, k, irho)
            temp27 = w(i, j-1, k, ivy)
            temp28 = w(i, j+2, k, irho)
            temp29 = w(i, j+2, k, ivz)
            temp30 = w(i, j-1, k, irho)
            temp31 = w(i, j-1, k, ivz)
            tempd18 = -(dis4*dred)
            dis2d = dis2d + ddw4*drwd + ddw2*drud + ddw1*drd + ddw3*drvd + &
            &         ddw5*dred
            ddw5d = dis2*dred - three*tempd18
            dis4d = dis4d - (temp28*temp29-temp30*temp31-three*ddw4)*drwd - &
            &         (temp20*temp21-temp22*temp23-three*ddw2)*drud - (w(i, j+2, k, &
            &         irho)-w(i, j-1, k, irho)-three*ddw1)*drd - (temp24*temp25-&
            &         temp26*temp27-three*ddw3)*drvd - (w(i, j+2, k, irhoe)-w(i, j-1&
            &         , k, irhoe)-three*ddw5)*dred
            wd(i, j+2, k, irhoe) = wd(i, j+2, k, irhoe) + tempd18
            wd(i, j-1, k, irhoe) = wd(i, j-1, k, irhoe) - tempd18
            wd(i, j+1, k, irhoe) = wd(i, j+1, k, irhoe) + ddw5d
            wd(i, j, k, irhoe) = wd(i, j, k, irhoe) - ddw5d
            tempd19 = -(dis4*drwd)
            ddw4d = dis2*drwd - three*tempd19
            wd(i, j+2, k, irho) = wd(i, j+2, k, irho) + temp29*tempd19
            wd(i, j+2, k, ivz) = wd(i, j+2, k, ivz) + temp28*tempd19
            wd(i, j-1, k, irho) = wd(i, j-1, k, irho) - temp31*tempd19
            wd(i, j-1, k, ivz) = wd(i, j-1, k, ivz) - temp30*tempd19
            wd(i, j+1, k, irho) = wd(i, j+1, k, irho) + w(i, j+1, k, ivz)*&
            &         ddw4d
            wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) + w(i, j+1, k, irho)*&
            &         ddw4d
            wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, ivz)*ddw4d
            wd(i, j, k, ivz) = wd(i, j, k, ivz) - w(i, j, k, irho)*ddw4d
            tempd20 = -(dis4*drvd)
            ddw3d = dis2*drvd - three*tempd20
            wd(i, j+2, k, irho) = wd(i, j+2, k, irho) + temp25*tempd20
            wd(i, j+2, k, ivy) = wd(i, j+2, k, ivy) + temp24*tempd20
            wd(i, j-1, k, irho) = wd(i, j-1, k, irho) - temp27*tempd20
            wd(i, j-1, k, ivy) = wd(i, j-1, k, ivy) - temp26*tempd20
            wd(i, j+1, k, irho) = wd(i, j+1, k, irho) + w(i, j+1, k, ivy)*&
            &         ddw3d
            wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) + w(i, j+1, k, irho)*&
            &         ddw3d
            wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, ivy)*ddw3d
            wd(i, j, k, ivy) = wd(i, j, k, ivy) - w(i, j, k, irho)*ddw3d
            tempd21 = -(dis4*drud)
            ddw2d = dis2*drud - three*tempd21
            wd(i, j+2, k, irho) = wd(i, j+2, k, irho) + temp21*tempd21
            wd(i, j+2, k, ivx) = wd(i, j+2, k, ivx) + temp20*tempd21
            wd(i, j-1, k, irho) = wd(i, j-1, k, irho) - temp23*tempd21
            wd(i, j-1, k, ivx) = wd(i, j-1, k, ivx) - temp22*tempd21
            wd(i, j+1, k, irho) = wd(i, j+1, k, irho) + w(i, j+1, k, ivx)*&
            &         ddw2d
            wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) + w(i, j+1, k, irho)*&
            &         ddw2d
            wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, ivx)*ddw2d
            wd(i, j, k, ivx) = wd(i, j, k, ivx) - w(i, j, k, irho)*ddw2d
            tempd22 = -(dis4*drd)
            ddw1d = dis2*drd - three*tempd22
            wd(i, j+2, k, irho) = wd(i, j+2, k, irho) + tempd22
            wd(i, j-1, k, irho) = wd(i, j-1, k, irho) - tempd22
            wd(i, j+1, k, irho) = wd(i, j+1, k, irho) + ddw1d
            wd(i, j, k, irho) = wd(i, j, k, irho) - ddw1d
            call mydim_fast_b(arg1, arg1d, dis2, dis2d, dis4d)
            min2d = ppor*fis2*dis2d
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                y2d = min2d
            else
                y2d = 0.0_8
            end if
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                dssd(i, j+1, k, 2) = dssd(i, j+1, k, 2) + y2d
            else
                dssd(i, j, k, 2) = dssd(i, j, k, 2) + y2d
            end if
        end do
        do ii=0,il*ny*nz-1
            i = mod(ii, il) + 1
            j = mod(ii/il, ny) + 2
            k = ii/(il*ny) + 2
            ! compute the dissipation coefficients for this face.
            ppor = zero
            if (pori(i, j, k) .eq. normalflux) ppor = one
            if (dss(i, j, k, 1) .lt. dss(i+1, j, k, 1)) then
                y1 = dss(i+1, j, k, 1)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                y1 = dss(i, j, k, 1)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            if (dpmax .gt. y1) then
                min1 = y1
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                min1 = dpmax
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            dis2 = ppor*fis2*min1
            arg1 = ppor*fis4
            dis4 = mydim(arg1, dis2)
            ! construct the vector of the first and third differences
            ! multiplied by the appropriate constants.
            ddw1 = w(i+1, j, k, irho) - w(i, j, k, irho)
            dr = dis2*ddw1 - dis4*(w(i+2, j, k, irho)-w(i-1, j, k, irho)-&
            &         three*ddw1)
            ddw2 = w(i+1, j, k, irho)*w(i+1, j, k, ivx) - w(i, j, k, irho)*w&
            &         (i, j, k, ivx)
            dru = dis2*ddw2 - dis4*(w(i+2, j, k, irho)*w(i+2, j, k, ivx)-w(i&
            &         -1, j, k, irho)*w(i-1, j, k, ivx)-three*ddw2)
            ddw3 = w(i+1, j, k, irho)*w(i+1, j, k, ivy) - w(i, j, k, irho)*w&
            &         (i, j, k, ivy)
            drv = dis2*ddw3 - dis4*(w(i+2, j, k, irho)*w(i+2, j, k, ivy)-w(i&
            &         -1, j, k, irho)*w(i-1, j, k, ivy)-three*ddw3)
            ddw4 = w(i+1, j, k, irho)*w(i+1, j, k, ivz) - w(i, j, k, irho)*w&
            &         (i, j, k, ivz)
            drw = dis2*ddw4 - dis4*(w(i+2, j, k, irho)*w(i+2, j, k, ivz)-w(i&
            &         -1, j, k, irho)*w(i-1, j, k, ivz)-three*ddw4)
            ddw5 = w(i+1, j, k, irhoe) - w(i, j, k, irhoe)
            dre = dis2*ddw5 - dis4*(w(i+2, j, k, irhoe)-w(i-1, j, k, irhoe)-&
            &         three*ddw5)
            ! in case a k-equation is present, compute the difference
            ! of rhok and store the average value of k. if not present,
            ! set both these values to zero, such that later on no
            ! decision needs to be made anymore.
            if (correctfork) then
                ddw6 = w(i+1, j, k, irho)*w(i+1, j, k, itu1) - w(i, j, k, irho&
                &           )*w(i, j, k, itu1)
                drk = dis2*ddw6 - dis4*(w(i+2, j, k, irho)*w(i+2, j, k, itu1)-&
                &           w(i-1, j, k, irho)*w(i-1, j, k, itu1)-three*ddw6)
                kavg = half*(w(i, j, k, itu1)+w(i+1, j, k, itu1))
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            else
                drk = zero
                kavg = zero
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            end if
            ! compute the average value of gamma and compute some
            ! expressions in which it occurs.
            gammaavg = half*(gamma(i+1, j, k)+gamma(i, j, k))
            gm1 = gammaavg - one
            ovgm1 = one/gm1
            gm53 = gammaavg - five*third
            ! compute the average state at the interface.
            uavg = half*(w(i+1, j, k, ivx)+w(i, j, k, ivx))
            vavg = half*(w(i+1, j, k, ivy)+w(i, j, k, ivy))
            wavg = half*(w(i+1, j, k, ivz)+w(i, j, k, ivz))
            a2avg = half*(gamma(i+1, j, k)*p(i+1, j, k)/w(i+1, j, k, irho)+&
            &         gamma(i, j, k)*p(i, j, k)/w(i, j, k, irho))
            area = sqrt(si(i, j, k, 1)**2 + si(i, j, k, 2)**2 + si(i, j, k, &
            &         3)**2)
            if (1.e-25_realtype .lt. area) then
                max1 = area
            else
                max1 = 1.e-25_realtype
            end if
            tmp = one/max1
            sx = si(i, j, k, 1)*tmp
            sy = si(i, j, k, 2)*tmp
            sz = si(i, j, k, 3)*tmp
            alphaavg = half*(uavg**2+vavg**2+wavg**2)
            havg = alphaavg + ovgm1*(a2avg-gm53*kavg)
            aavg = sqrt(a2avg)
            unavg = uavg*sx + vavg*sy + wavg*sz
            ovaavg = one/aavg
            ova2avg = one/a2avg
            ! the mesh velocity if the face is moving. it must be
            ! divided by the area to obtain a true velocity.
            if (addgridvelocities) sface = sfacei(i, j, k)*tmp
            if (unavg - sface + aavg .ge. 0.) then
                lam1 = unavg - sface + aavg
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                lam1 = -(unavg-sface+aavg)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            if (unavg - sface - aavg .ge. 0.) then
                lam2 = unavg - sface - aavg
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                lam2 = -(unavg-sface-aavg)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            if (unavg - sface .ge. 0.) then
                lam3 = unavg - sface
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                lam3 = -(unavg-sface)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            rrad = lam3 + aavg
            if (lam1 .lt. epsacoustic*rrad) then
                max2 = epsacoustic*rrad
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                max2 = lam1
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            ! multiply the eigenvalues by the area to obtain
            ! the correct values for the dissipation term.
            lam1 = max2*area
            if (lam2 .lt. epsacoustic*rrad) then
                max3 = epsacoustic*rrad
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                max3 = lam2
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            lam2 = max3*area
            if (lam3 .lt. epsshear*rrad) then
                max4 = epsshear*rrad
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                max4 = lam3
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            lam3 = max4*area
            ! some abbreviations, which occur quite often in the
            ! dissipation terms.
            abv1 = half*(lam1+lam2)
            abv2 = half*(lam1-lam2)
            abv3 = abv1 - lam3
            abv4 = gm1*(alphaavg*dr-uavg*dru-vavg*drv-wavg*drw+dre) - gm53*&
            &         drk
            abv5 = sx*dru + sy*drv + sz*drw - unavg*dr
            abv6 = abv3*abv4*ova2avg + abv2*abv5*ovaavg
            abv7 = abv2*abv4*ovaavg + abv3*abv5
            ! compute and scatter the dissipative flux.
            ! density.
            ! x-momentum.
            ! y-momentum.
            ! z-momentum.
            ! energy.
            fsd = fwd(i+1, j, k, irhoe) - fwd(i, j, k, irhoe)
            lam3d = dre*fsd
            dred = lam3*fsd
            havgd = abv6*fsd
            abv6d = havg*fsd
            unavgd = abv7*fsd
            abv7d = unavg*fsd
            fsd = fwd(i+1, j, k, imz) - fwd(i, j, k, imz)
            lam3d = lam3d + drw*fsd
            drwd = lam3*fsd
            wavgd = abv6*fsd
            abv6d = abv6d + wavg*fsd
            abv7d = abv7d + sz*fsd
            fsd = fwd(i+1, j, k, imy) - fwd(i, j, k, imy)
            lam3d = lam3d + drv*fsd
            drvd = lam3*fsd
            vavgd = abv6*fsd
            abv6d = abv6d + vavg*fsd
            abv7d = abv7d + sy*fsd
            fsd = fwd(i+1, j, k, imx) - fwd(i, j, k, imx)
            lam3d = lam3d + dru*fsd
            drud = lam3*fsd
            uavgd = abv6*fsd
            abv6d = abv6d + uavg*fsd
            abv7d = abv7d + sx*fsd
            fsd = fwd(i+1, j, k, irho) - fwd(i, j, k, irho)
            abv6d = abv6d + fsd
            abv2d = ovaavg*abv5*abv6d + ovaavg*abv4*abv7d
            abv4d = ova2avg*abv3*abv6d + ovaavg*abv2*abv7d
            ovaavgd = abv2*abv5*abv6d + abv2*abv4*abv7d
            abv3d = ova2avg*abv4*abv6d + abv5*abv7d
            lam3d = lam3d + dr*fsd - abv3d
            abv5d = ovaavg*abv2*abv6d + abv3*abv7d
            ova2avgd = abv3*abv4*abv6d
            unavgd = unavgd - dr*abv5d
            tempd17 = gm1*abv4d
            drd = alphaavg*tempd17 - unavg*abv5d + lam3*fsd
            drud = drud + sx*abv5d - uavg*tempd17
            drvd = drvd + sy*abv5d - vavg*tempd17
            drwd = drwd + sz*abv5d - wavg*tempd17
            alphaavgd = dr*tempd17
            uavgd = uavgd - dru*tempd17
            vavgd = vavgd - drv*tempd17
            dred = dred + tempd17
            wavgd = wavgd - drw*tempd17
            drkd = -(gm53*abv4d)
            abv1d = abv3d
            lam1d = half*abv1d + half*abv2d
            lam2d = half*abv1d - half*abv2d
            max4d = area*lam3d
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                rradd = epsshear*max4d
                lam3d = 0.0_8
            else
                lam3d = max4d
                rradd = 0.0_8
            end if
            max3d = area*lam2d
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                rradd = rradd + epsacoustic*max3d
                lam2d = 0.0_8
            else
                lam2d = max3d
            end if
            max2d = area*lam1d
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                rradd = rradd + epsacoustic*max2d
                lam1d = 0.0_8
            else
                lam1d = max2d
            end if
            lam3d = lam3d + rradd
            aavgd = rradd
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                unavgd = unavgd + lam3d
            else
                unavgd = unavgd - lam3d
            end if
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                unavgd = unavgd + lam2d
                aavgd = aavgd - lam2d
            else
                aavgd = aavgd + lam2d
                unavgd = unavgd - lam2d
            end if
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                unavgd = unavgd + lam1d
                aavgd = aavgd + lam1d
            else
                unavgd = unavgd - lam1d
                aavgd = aavgd - lam1d
            end if
            alphaavgd = alphaavgd + havgd
            tempd16 = half*alphaavgd
            aavgd = aavgd - one*ovaavgd/aavg**2
            if (a2avg .eq. 0.0_8) then
                a2avgd = ovgm1*havgd - one*ova2avgd/a2avg**2
            else
                a2avgd = aavgd/(2.0*sqrt(a2avg)) + ovgm1*havgd - one*ova2avgd/&
                &           a2avg**2
            end if
            uavgd = uavgd + 2*uavg*tempd16 + sx*unavgd
            vavgd = vavgd + 2*vavg*tempd16 + sy*unavgd
            wavgd = wavgd + 2*wavg*tempd16 + sz*unavgd
            kavgd = -(ovgm1*gm53*havgd)
            temp19 = w(i, j, k, irho)
            temp18 = w(i+1, j, k, irho)
            tempd14 = gamma(i+1, j, k)*half*a2avgd/temp18
            tempd15 = gamma(i, j, k)*half*a2avgd/temp19
            pd(i+1, j, k) = pd(i+1, j, k) + tempd14
            wd(i+1, j, k, irho) = wd(i+1, j, k, irho) - p(i+1, j, k)*tempd14&
            &         /temp18
            pd(i, j, k) = pd(i, j, k) + tempd15
            wd(i, j, k, irho) = wd(i, j, k, irho) - p(i, j, k)*tempd15/&
            &         temp19
            wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) + half*wavgd
            wd(i, j, k, ivz) = wd(i, j, k, ivz) + half*wavgd
            wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) + half*vavgd
            wd(i, j, k, ivy) = wd(i, j, k, ivy) + half*vavgd
            wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) + half*uavgd
            wd(i, j, k, ivx) = wd(i, j, k, ivx) + half*uavgd
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                dis2d = 0.0_8
                dis4d = 0.0_8
            else
                wd(i, j, k, itu1) = wd(i, j, k, itu1) + half*kavgd
                wd(i+1, j, k, itu1) = wd(i+1, j, k, itu1) + half*kavgd
                temp17 = w(i-1, j, k, itu1)
                temp16 = w(i-1, j, k, irho)
                temp15 = w(i+2, j, k, itu1)
                temp14 = w(i+2, j, k, irho)
                tempd13 = -(dis4*drkd)
                dis2d = ddw6*drkd
                ddw6d = dis2*drkd - three*tempd13
                dis4d = -((temp14*temp15-temp16*temp17-three*ddw6)*drkd)
                wd(i+2, j, k, irho) = wd(i+2, j, k, irho) + temp15*tempd13
                wd(i+2, j, k, itu1) = wd(i+2, j, k, itu1) + temp14*tempd13
                wd(i-1, j, k, irho) = wd(i-1, j, k, irho) - temp17*tempd13
                wd(i-1, j, k, itu1) = wd(i-1, j, k, itu1) - temp16*tempd13
                wd(i+1, j, k, irho) = wd(i+1, j, k, irho) + w(i+1, j, k, itu1)&
                &           *ddw6d
                wd(i+1, j, k, itu1) = wd(i+1, j, k, itu1) + w(i+1, j, k, irho)&
                &           *ddw6d
                wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, itu1)*ddw6d
                wd(i, j, k, itu1) = wd(i, j, k, itu1) - w(i, j, k, irho)*ddw6d
            end if
            temp2 = w(i+2, j, k, irho)
            temp3 = w(i+2, j, k, ivx)
            temp4 = w(i-1, j, k, irho)
            temp5 = w(i-1, j, k, ivx)
            temp6 = w(i+2, j, k, irho)
            temp7 = w(i+2, j, k, ivy)
            temp8 = w(i-1, j, k, irho)
            temp9 = w(i-1, j, k, ivy)
            temp10 = w(i+2, j, k, irho)
            temp11 = w(i+2, j, k, ivz)
            temp12 = w(i-1, j, k, irho)
            temp13 = w(i-1, j, k, ivz)
            tempd8 = -(dis4*dred)
            dis2d = dis2d + ddw4*drwd + ddw2*drud + ddw1*drd + ddw3*drvd + &
            &         ddw5*dred
            ddw5d = dis2*dred - three*tempd8
            dis4d = dis4d - (temp10*temp11-temp12*temp13-three*ddw4)*drwd - &
            &         (temp2*temp3-temp4*temp5-three*ddw2)*drud - (w(i+2, j, k, irho&
            &         )-w(i-1, j, k, irho)-three*ddw1)*drd - (temp6*temp7-temp8*&
            &         temp9-three*ddw3)*drvd - (w(i+2, j, k, irhoe)-w(i-1, j, k, &
            &         irhoe)-three*ddw5)*dred
            wd(i+2, j, k, irhoe) = wd(i+2, j, k, irhoe) + tempd8
            wd(i-1, j, k, irhoe) = wd(i-1, j, k, irhoe) - tempd8
            wd(i+1, j, k, irhoe) = wd(i+1, j, k, irhoe) + ddw5d
            wd(i, j, k, irhoe) = wd(i, j, k, irhoe) - ddw5d
            tempd9 = -(dis4*drwd)
            ddw4d = dis2*drwd - three*tempd9
            wd(i+2, j, k, irho) = wd(i+2, j, k, irho) + temp11*tempd9
            wd(i+2, j, k, ivz) = wd(i+2, j, k, ivz) + temp10*tempd9
            wd(i-1, j, k, irho) = wd(i-1, j, k, irho) - temp13*tempd9
            wd(i-1, j, k, ivz) = wd(i-1, j, k, ivz) - temp12*tempd9
            wd(i+1, j, k, irho) = wd(i+1, j, k, irho) + w(i+1, j, k, ivz)*&
            &         ddw4d
            wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) + w(i+1, j, k, irho)*&
            &         ddw4d
            wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, ivz)*ddw4d
            wd(i, j, k, ivz) = wd(i, j, k, ivz) - w(i, j, k, irho)*ddw4d
            tempd10 = -(dis4*drvd)
            ddw3d = dis2*drvd - three*tempd10
            wd(i+2, j, k, irho) = wd(i+2, j, k, irho) + temp7*tempd10
            wd(i+2, j, k, ivy) = wd(i+2, j, k, ivy) + temp6*tempd10
            wd(i-1, j, k, irho) = wd(i-1, j, k, irho) - temp9*tempd10
            wd(i-1, j, k, ivy) = wd(i-1, j, k, ivy) - temp8*tempd10
            wd(i+1, j, k, irho) = wd(i+1, j, k, irho) + w(i+1, j, k, ivy)*&
            &         ddw3d
            wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) + w(i+1, j, k, irho)*&
            &         ddw3d
            wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, ivy)*ddw3d
            wd(i, j, k, ivy) = wd(i, j, k, ivy) - w(i, j, k, irho)*ddw3d
            tempd11 = -(dis4*drud)
            ddw2d = dis2*drud - three*tempd11
            wd(i+2, j, k, irho) = wd(i+2, j, k, irho) + temp3*tempd11
            wd(i+2, j, k, ivx) = wd(i+2, j, k, ivx) + temp2*tempd11
            wd(i-1, j, k, irho) = wd(i-1, j, k, irho) - temp5*tempd11
            wd(i-1, j, k, ivx) = wd(i-1, j, k, ivx) - temp4*tempd11
            wd(i+1, j, k, irho) = wd(i+1, j, k, irho) + w(i+1, j, k, ivx)*&
            &         ddw2d
            wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) + w(i+1, j, k, irho)*&
            &         ddw2d
            wd(i, j, k, irho) = wd(i, j, k, irho) - w(i, j, k, ivx)*ddw2d
            wd(i, j, k, ivx) = wd(i, j, k, ivx) - w(i, j, k, irho)*ddw2d
            tempd12 = -(dis4*drd)
            ddw1d = dis2*drd - three*tempd12
            wd(i+2, j, k, irho) = wd(i+2, j, k, irho) + tempd12
            wd(i-1, j, k, irho) = wd(i-1, j, k, irho) - tempd12
            wd(i+1, j, k, irho) = wd(i+1, j, k, irho) + ddw1d
            wd(i, j, k, irho) = wd(i, j, k, irho) - ddw1d
            call mydim_fast_b(arg1, arg1d, dis2, dis2d, dis4d)
            min1d = ppor*fis2*dis2d
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                y1d = min1d
            else
                y1d = 0.0_8
            end if
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                dssd(i+1, j, k, 1) = dssd(i+1, j, k, 1) + y1d
            else
                dssd(i, j, k, 1) = dssd(i, j, k, 1) + y1d
            end if
        end do
        do ii=0,ie*je*ke-1
            i = mod(ii, ie) + 1
            j = mod(ii/ie, je) + 1
            k = ii/(ie*je) + 1
            if (p(i+1, j, k) - p(i, j, k) .ge. 0.) then
                abs1 = p(i+1, j, k) - p(i, j, k)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            else
                abs1 = -(p(i+1, j, k)-p(i, j, k))
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            end if
            if (p(i, j, k) - p(i-1, j, k) .ge. 0.) then
                abs4 = p(i, j, k) - p(i-1, j, k)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                abs4 = -(p(i, j, k)-p(i-1, j, k))
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            x1 = (p(i+1, j, k)-two*p(i, j, k)+p(i-1, j, k))/(omega*(p(i+1, j&
            &         , k)+two*p(i, j, k)+p(i-1, j, k))+oneminomega*(abs1+abs4)+plim&
            &         )
            if (x1 .ge. 0.) then
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            if (p(i, j+1, k) - p(i, j, k) .ge. 0.) then
                abs2 = p(i, j+1, k) - p(i, j, k)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            else
                abs2 = -(p(i, j+1, k)-p(i, j, k))
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            end if
            if (p(i, j, k) - p(i, j-1, k) .ge. 0.) then
                abs5 = p(i, j, k) - p(i, j-1, k)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                abs5 = -(p(i, j, k)-p(i, j-1, k))
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            x2 = (p(i, j+1, k)-two*p(i, j, k)+p(i, j-1, k))/(omega*(p(i, j+1&
            &         , k)+two*p(i, j, k)+p(i, j-1, k))+oneminomega*(abs2+abs5)+plim&
            &         )
            if (x2 .ge. 0.) then
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            if (p(i, j, k+1) - p(i, j, k) .ge. 0.) then
                abs3 = p(i, j, k+1) - p(i, j, k)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            else
                abs3 = -(p(i, j, k+1)-p(i, j, k))
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            end if
            if (p(i, j, k) - p(i, j, k-1) .ge. 0.) then
                abs6 = p(i, j, k) - p(i, j, k-1)
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                abs6 = -(p(i, j, k)-p(i, j, k-1))
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            x3 = (p(i, j, k+1)-two*p(i, j, k)+p(i, j, k-1))/(omega*(p(i, j, &
            &         k+1)+two*p(i, j, k)+p(i, j, k-1))+oneminomega*(abs3+abs6)+plim&
            &         )
            if (x3 .ge. 0.) then
                x3d = dssd(i, j, k, 3)
                dssd(i, j, k, 3) = 0.0_8
            else
                x3d = -dssd(i, j, k, 3)
                dssd(i, j, k, 3) = 0.0_8
            end if
            temp1 = plim + omega*(p(i, j, k+1)+two*p(i, j, k)+p(i, j, k-1)) &
            &         + oneminomega*(abs3+abs6)
            tempd5 = x3d/temp1
            tempd6 = -((p(i, j, k+1)-two*p(i, j, k)+p(i, j, k-1))*tempd5/&
            &         temp1)
            tempd7 = omega*tempd6
            pd(i, j, k+1) = pd(i, j, k+1) + tempd7 + tempd5
            pd(i, j, k) = pd(i, j, k) + two*tempd7 - two*tempd5
            pd(i, j, k-1) = pd(i, j, k-1) + tempd7 + tempd5
            abs3d = oneminomega*tempd6
            abs6d = oneminomega*tempd6
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                pd(i, j, k) = pd(i, j, k) + abs6d
                pd(i, j, k-1) = pd(i, j, k-1) - abs6d
            else
                pd(i, j, k-1) = pd(i, j, k-1) + abs6d
                pd(i, j, k) = pd(i, j, k) - abs6d
            end if
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                pd(i, j, k) = pd(i, j, k) + abs3d
                pd(i, j, k+1) = pd(i, j, k+1) - abs3d
            else
                pd(i, j, k+1) = pd(i, j, k+1) + abs3d
                pd(i, j, k) = pd(i, j, k) - abs3d
            end if
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                x2d = dssd(i, j, k, 2)
                dssd(i, j, k, 2) = 0.0_8
            else
                x2d = -dssd(i, j, k, 2)
                dssd(i, j, k, 2) = 0.0_8
            end if
            temp0 = plim + omega*(p(i, j+1, k)+two*p(i, j, k)+p(i, j-1, k)) &
            &         + oneminomega*(abs2+abs5)
            tempd2 = x2d/temp0
            tempd3 = -((p(i, j+1, k)-two*p(i, j, k)+p(i, j-1, k))*tempd2/&
            &         temp0)
            tempd4 = omega*tempd3
            pd(i, j+1, k) = pd(i, j+1, k) + tempd4 + tempd2
            pd(i, j, k) = pd(i, j, k) + two*tempd4 - two*tempd2
            pd(i, j-1, k) = pd(i, j-1, k) + tempd4 + tempd2
            abs2d = oneminomega*tempd3
            abs5d = oneminomega*tempd3
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                pd(i, j, k) = pd(i, j, k) + abs5d
                pd(i, j-1, k) = pd(i, j-1, k) - abs5d
            else
                pd(i, j-1, k) = pd(i, j-1, k) + abs5d
                pd(i, j, k) = pd(i, j, k) - abs5d
            end if
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                pd(i, j, k) = pd(i, j, k) + abs2d
                pd(i, j+1, k) = pd(i, j+1, k) - abs2d
            else
                pd(i, j+1, k) = pd(i, j+1, k) + abs2d
                pd(i, j, k) = pd(i, j, k) - abs2d
            end if
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                x1d = dssd(i, j, k, 1)
                dssd(i, j, k, 1) = 0.0_8
            else
                x1d = -dssd(i, j, k, 1)
                dssd(i, j, k, 1) = 0.0_8
            end if
            temp = plim + omega*(p(i+1, j, k)+two*p(i, j, k)+p(i-1, j, k)) +&
            &         oneminomega*(abs1+abs4)
            tempd = x1d/temp
            tempd0 = -((p(i+1, j, k)-two*p(i, j, k)+p(i-1, j, k))*tempd/temp&
            &         )
            tempd1 = omega*tempd0
            pd(i+1, j, k) = pd(i+1, j, k) + tempd1 + tempd
            pd(i, j, k) = pd(i, j, k) + two*tempd1 - two*tempd
            pd(i-1, j, k) = pd(i-1, j, k) + tempd1 + tempd
            abs1d = oneminomega*tempd0
            abs4d = oneminomega*tempd0
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                pd(i, j, k) = pd(i, j, k) + abs4d
                pd(i-1, j, k) = pd(i-1, j, k) - abs4d
            else
                pd(i-1, j, k) = pd(i-1, j, k) + abs4d
                pd(i, j, k) = pd(i, j, k) - abs4d
            end if
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                pd(i, j, k) = pd(i, j, k) + abs1d
                pd(i+1, j, k) = pd(i+1, j, k) - abs1d
            else
                pd(i+1, j, k) = pd(i+1, j, k) + abs1d
                pd(i, j, k) = pd(i, j, k) - abs1d
            end if
        end do
        fwd = sfil*fwd
    end if
end subroutine invisciddissfluxmatrix_fast_b

!  differentiation of inviscidupwindflux in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: *p *w *fw
!   with respect to varying inputs: *p *w *fw
!   rw status of diff variables: *p:incr *w:incr *fw:in-out
!   plus diff mem management of: p:in w:in fw:in
subroutine inviscidupwindflux_fast_b(finegrid)
    !
    !       inviscidupwindflux computes the artificial dissipation part of
    !       the euler fluxes by means of an approximate solution of the 1d
    !       riemann problem on the face. for first order schemes,
    !       finegrid == .false., the states in the cells are assumed to
    !       be constant; for the second order schemes on the fine grid a
    !       nonlinear reconstruction of the left and right state is done
    !       for which several options exist.
    !       it is assumed that the pointers in blockpointers already
    !       point to the correct block.
    !
    use constants
    ! use blockpointers, only : il, jl, kl, ie, je, ke, ib, jb, kb, w, &
    ! &   wd, p, pd, pori, porj, pork, fw, fwd, gamma, si, sj, sk, indfamilyi,&
    ! &   indfamilyj, indfamilyk, spectralsol, addgridvelocities, sfacei, &
    ! &   sfacej, sfacek, rotmatrixi, rotmatrixj, rotmatrixk, factfamilyi, &
    ! &   factfamilyj, factfamilyk
    use blockpointers, only : addgridvelocities, rotmatrixi, rotmatrixj, &
    & rotmatrixk, factfamilyi, factfamilyj, factfamilyk
    use flowvarrefstate, only : kpresent, nw, nwf, rgas, tref
    use inputdiscretization, only : limiter, lumpeddiss, precond, &
    &   riemann, riemanncoarse, orderturb, kappacoef
    use inputphysics, only : equations
    use iteration, only : rfil, currentlevel, groundlevel
    use cgnsgrid, only : massflowfamilydiss
    use utils_fast_b, only : getcorrectfork, terminate
    use flowutils_fast_b, only : etot, etot_fast_b
    implicit none
    !
    !      subroutine arguments.
    !
    logical, intent(in) :: finegrid
    !
    !      local variables.
    !
    integer(kind=portype) :: por
    integer(kind=inttype) :: nwint
    integer(kind=inttype) :: i, j, k, ind
    integer(kind=inttype) :: limused, riemannused
    real(kind=realtype) :: sx, sy, sz, omk, opk, sfil, gammaface
    real(kind=realtype) :: factminmod, sface
    real(kind=realtype), dimension(nw) :: left, right
    real(kind=realtype), dimension(nw) :: leftd, rightd
    real(kind=realtype), dimension(nw) :: du1, du2, du3
    real(kind=realtype), dimension(nw) :: du1d, du2d, du3d
    real(kind=realtype), dimension(nwf) :: flux
    real(kind=realtype), dimension(nwf) :: fluxd
    logical :: firstorderk, correctfork, rotationalperiodic
    intrinsic abs
    intrinsic associated
    intrinsic max
    integer :: branch
    real(kind=realtype) :: abs0
    real(kind=realtype) :: max1
    if (rfil .ge. 0.) then
        abs0 = rfil
    else
        abs0 = -rfil
    end if
    !
    ! check if rfil == 0. if so, the dissipative flux needs not to
    ! be computed.
    if (abs0 .ge. thresholdreal) then
        ! check if the formulation for rotational periodic problems
        ! must be used.
        if (associated(rotmatrixi)) then
            rotationalperiodic = .true.
        else
            rotationalperiodic = .false.
        end if
        ! initialize the dissipative residual to a certain times,
        ! possibly zero, the previously stored value. owned cells
        ! only, because the halo values do not matter.
        sfil = one - rfil
        ! determine whether or not the total energy must be corrected
        ! for the presence of the turbulent kinetic energy.
        correctfork = getcorrectfork()
        if (1.e-10_realtype .lt. one - kappacoef) then
            max1 = one - kappacoef
        else
            max1 = 1.e-10_realtype
        end if
        ! compute the factor used in the minmod limiter.
        factminmod = (three-kappacoef)/max1
        ! determine the limiter scheme to be used. on the fine grid the
        ! user specified scheme is used; on the coarse grid a first order
        ! scheme is computed.
        limused = firstorder
        if (finegrid) limused = limiter
        ! lumped diss is true for doing approx pc
        if (lumpeddiss) limused = firstorder
        ! determine the riemann solver which must be used.
        riemannused = riemanncoarse
        if (finegrid) riemannused = riemann
        ! store 1-kappa and 1+kappa a bit easier and multiply it by 0.25.
        omk = fourth*(one-kappacoef)
        opk = fourth*(one+kappacoef)
        ! initialize sface to zero. this value will be used if the
        ! block is not moving.
        sface = zero
        ! set the number of variables to be interpolated depending
        ! whether or not a k-equation is present. if a k-equation is
        ! present also set the logical firstorderk. this indicates
        ! whether or not only a first order approximation is to be used
        ! for the turbulent kinetic energy.
        if (correctfork) then
            if (orderturb .eq. firstorder) then
                nwint = nwf
                firstorderk = .true.
            else
                nwint = itu1
                firstorderk = .false.
            end if
        else
            nwint = nwf
            firstorderk = .false.
        end if
        !
        !       flux computation. a distinction is made between first and
        !       second order schemes to avoid the overhead for the first order
        !       scheme.
        !
        if (limused .eq. firstorder) then
            !
            !         first order reconstruction. the states in the cells are
            !         constant. the left and right states are constructed easily.
            !
            ! fluxes in the i-direction.
            do k=2,kl
                do j=2,jl
                    do i=1,il
                        ! store the normal vector, the porosity and the
                        ! mesh velocity if present.
                        sx = si(i, j, k, 1)
                        sy = si(i, j, k, 2)
                        sz = si(i, j, k, 3)
                        if (addgridvelocities) then
                            sface = sfacei(i, j, k)
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 0
                        else
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 1
                        end if
                        ! determine the left and right state.
                        left(irho) = w(i, j, k, irho)
                        left(ivx) = w(i, j, k, ivx)
                        left(ivy) = w(i, j, k, ivy)
                        left(ivz) = w(i, j, k, ivz)
                        left(irhoe) = p(i, j, k)
                        if (correctfork) then
                            left(itu1) = w(i, j, k, itu1)
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 0
                        else
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 1
                        end if
                        right(irho) = w(i+1, j, k, irho)
                        right(ivx) = w(i+1, j, k, ivx)
                        right(ivy) = w(i+1, j, k, ivy)
                        right(ivz) = w(i+1, j, k, ivz)
                        right(irhoe) = p(i+1, j, k)
                        if (correctfork) then
                            right(itu1) = w(i+1, j, k, itu1)
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 0
                        else
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 1
                        end if
                    end do
                end do
            end do
            ! store the density flux in the mass flow of the
            ! appropriate sliding mesh interface.
            ! fluxes in j-direction.
            do k=2,kl
                do j=1,jl
                    do i=2,il
                        ! store the normal vector, the porosity and the
                        ! mesh velocity if present.
                        sx = sj(i, j, k, 1)
                        sy = sj(i, j, k, 2)
                        sz = sj(i, j, k, 3)
                        if (addgridvelocities) then
                            sface = sfacej(i, j, k)
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 0
                        else
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 1
                        end if
                        ! determine the left and right state.
                        left(irho) = w(i, j, k, irho)
                        left(ivx) = w(i, j, k, ivx)
                        left(ivy) = w(i, j, k, ivy)
                        left(ivz) = w(i, j, k, ivz)
                        left(irhoe) = p(i, j, k)
                        if (correctfork) then
                            left(itu1) = w(i, j, k, itu1)
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 0
                        else
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 1
                        end if
                        right(irho) = w(i, j+1, k, irho)
                        right(ivx) = w(i, j+1, k, ivx)
                        right(ivy) = w(i, j+1, k, ivy)
                        right(ivz) = w(i, j+1, k, ivz)
                        right(irhoe) = p(i, j+1, k)
                        if (correctfork) then
                            right(itu1) = w(i, j+1, k, itu1)
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 0
                        else
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 1
                        end if
                    end do
                end do
            end do
            ! store the density flux in the mass flow of the
            ! appropriate sliding mesh interface.
            ! fluxes in k-direction.
            do k=1,kl
                do j=2,jl
                    do i=2,il
                        ! store the normal vector, the porosity and the
                        ! mesh velocity if present.
                        sx = sk(i, j, k, 1)
                        sy = sk(i, j, k, 2)
                        sz = sk(i, j, k, 3)
                        if (addgridvelocities) then
                            sface = sfacek(i, j, k)
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 0
                        else
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 1
                        end if
                        ! determine the left and right state.
                        left(irho) = w(i, j, k, irho)
                        left(ivx) = w(i, j, k, ivx)
                        left(ivy) = w(i, j, k, ivy)
                        left(ivz) = w(i, j, k, ivz)
                        left(irhoe) = p(i, j, k)
                        if (correctfork) then
                            left(itu1) = w(i, j, k, itu1)
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 0
                        else
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 1
                        end if
                        right(irho) = w(i, j, k+1, irho)
                        right(ivx) = w(i, j, k+1, ivx)
                        right(ivy) = w(i, j, k+1, ivy)
                        right(ivz) = w(i, j, k+1, ivz)
                        right(irhoe) = p(i, j, k+1)
                        if (correctfork) then
                            right(itu1) = w(i, j, k+1, itu1)
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 0
                        else
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 1
                        end if
                    end do
                end do
            end do
            fluxd = 0.0_8
            leftd = 0.0_8
            rightd = 0.0_8
            do k=kl,1,-1
                do j=jl,2,-1
                    do i=il,2,-1
                        fluxd(irhoe) = fluxd(irhoe) - fwd(i, j, k+1, irhoe)
                        fluxd(imz) = fluxd(imz) - fwd(i, j, k+1, imz)
                        fluxd(imy) = fluxd(imy) - fwd(i, j, k+1, imy)
                        fluxd(imx) = fluxd(imx) - fwd(i, j, k+1, imx)
                        fluxd(irho) = fluxd(irho) - fwd(i, j, k+1, irho)
                        fluxd(irhoe) = fluxd(irhoe) + fwd(i, j, k, irhoe)
                        fluxd(imz) = fluxd(imz) + fwd(i, j, k, imz)
                        fluxd(imy) = fluxd(imy) + fwd(i, j, k, imy)
                        fluxd(imx) = fluxd(imx) + fwd(i, j, k, imx)
                        fluxd(irho) = fluxd(irho) + fwd(i, j, k, irho)
                        gammaface = half*(gamma(i, j, k)+gamma(i, j, k+1))
                        por = pork(i, j, k)
                        call riemannflux_fast_b(left, leftd, right, rightd, flux, &
                        &                               fluxd)
                        branch = myIntStack(myIntPtr)
                        myIntPtr = myIntPtr - 1
                        if (branch .eq. 0) then
                            wd(i, j, k+1, itu1) = wd(i, j, k+1, itu1) + rightd(itu1)
                            rightd(itu1) = 0.0_8
                        end if
                        pd(i, j, k+1) = pd(i, j, k+1) + rightd(irhoe)
                        rightd(irhoe) = 0.0_8
                        wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) + rightd(ivz)
                        rightd(ivz) = 0.0_8
                        wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) + rightd(ivy)
                        rightd(ivy) = 0.0_8
                        wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) + rightd(ivx)
                        rightd(ivx) = 0.0_8
                        wd(i, j, k+1, irho) = wd(i, j, k+1, irho) + rightd(irho)
                        rightd(irho) = 0.0_8
                        branch = myIntStack(myIntPtr)
                        myIntPtr = myIntPtr - 1
                        if (branch .eq. 0) then
                            wd(i, j, k, itu1) = wd(i, j, k, itu1) + leftd(itu1)
                            leftd(itu1) = 0.0_8
                        end if
                        pd(i, j, k) = pd(i, j, k) + leftd(irhoe)
                        leftd(irhoe) = 0.0_8
                        wd(i, j, k, ivz) = wd(i, j, k, ivz) + leftd(ivz)
                        leftd(ivz) = 0.0_8
                        wd(i, j, k, ivy) = wd(i, j, k, ivy) + leftd(ivy)
                        leftd(ivy) = 0.0_8
                        wd(i, j, k, ivx) = wd(i, j, k, ivx) + leftd(ivx)
                        leftd(ivx) = 0.0_8
                        wd(i, j, k, irho) = wd(i, j, k, irho) + leftd(irho)
                        leftd(irho) = 0.0_8
                        branch = myIntStack(myIntPtr)
                        myIntPtr = myIntPtr - 1
                        if (branch .eq. 0) call popreal8(sface)
                    end do
                end do
            end do
            do k=kl,2,-1
                do j=jl,1,-1
                    do i=il,2,-1
                        fluxd(irhoe) = fluxd(irhoe) - fwd(i, j+1, k, irhoe)
                        fluxd(imz) = fluxd(imz) - fwd(i, j+1, k, imz)
                        fluxd(imy) = fluxd(imy) - fwd(i, j+1, k, imy)
                        fluxd(imx) = fluxd(imx) - fwd(i, j+1, k, imx)
                        fluxd(irho) = fluxd(irho) - fwd(i, j+1, k, irho)
                        fluxd(irhoe) = fluxd(irhoe) + fwd(i, j, k, irhoe)
                        fluxd(imz) = fluxd(imz) + fwd(i, j, k, imz)
                        fluxd(imy) = fluxd(imy) + fwd(i, j, k, imy)
                        fluxd(imx) = fluxd(imx) + fwd(i, j, k, imx)
                        fluxd(irho) = fluxd(irho) + fwd(i, j, k, irho)
                        gammaface = half*(gamma(i, j, k)+gamma(i, j+1, k))
                        por = porj(i, j, k)
                        call riemannflux_fast_b(left, leftd, right, rightd, flux, &
                        &                               fluxd)
                        branch = myIntStack(myIntPtr)
                        myIntPtr = myIntPtr - 1
                        if (branch .eq. 0) then
                            wd(i, j+1, k, itu1) = wd(i, j+1, k, itu1) + rightd(itu1)
                            rightd(itu1) = 0.0_8
                        end if
                        pd(i, j+1, k) = pd(i, j+1, k) + rightd(irhoe)
                        rightd(irhoe) = 0.0_8
                        wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) + rightd(ivz)
                        rightd(ivz) = 0.0_8
                        wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) + rightd(ivy)
                        rightd(ivy) = 0.0_8
                        wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) + rightd(ivx)
                        rightd(ivx) = 0.0_8
                        wd(i, j+1, k, irho) = wd(i, j+1, k, irho) + rightd(irho)
                        rightd(irho) = 0.0_8
                        branch = myIntStack(myIntPtr)
                        myIntPtr = myIntPtr - 1
                        if (branch .eq. 0) then
                            wd(i, j, k, itu1) = wd(i, j, k, itu1) + leftd(itu1)
                            leftd(itu1) = 0.0_8
                        end if
                        pd(i, j, k) = pd(i, j, k) + leftd(irhoe)
                        leftd(irhoe) = 0.0_8
                        wd(i, j, k, ivz) = wd(i, j, k, ivz) + leftd(ivz)
                        leftd(ivz) = 0.0_8
                        wd(i, j, k, ivy) = wd(i, j, k, ivy) + leftd(ivy)
                        leftd(ivy) = 0.0_8
                        wd(i, j, k, ivx) = wd(i, j, k, ivx) + leftd(ivx)
                        leftd(ivx) = 0.0_8
                        wd(i, j, k, irho) = wd(i, j, k, irho) + leftd(irho)
                        leftd(irho) = 0.0_8
                        branch = myIntStack(myIntPtr)
                        myIntPtr = myIntPtr - 1
                        if (branch .eq. 0) call popreal8(sface)
                    end do
                end do
            end do
            do k=kl,2,-1
                do j=jl,2,-1
                    do i=il,1,-1
                        fluxd(irhoe) = fluxd(irhoe) - fwd(i+1, j, k, irhoe)
                        fluxd(imz) = fluxd(imz) - fwd(i+1, j, k, imz)
                        fluxd(imy) = fluxd(imy) - fwd(i+1, j, k, imy)
                        fluxd(imx) = fluxd(imx) - fwd(i+1, j, k, imx)
                        fluxd(irho) = fluxd(irho) - fwd(i+1, j, k, irho)
                        fluxd(irhoe) = fluxd(irhoe) + fwd(i, j, k, irhoe)
                        fluxd(imz) = fluxd(imz) + fwd(i, j, k, imz)
                        fluxd(imy) = fluxd(imy) + fwd(i, j, k, imy)
                        fluxd(imx) = fluxd(imx) + fwd(i, j, k, imx)
                        fluxd(irho) = fluxd(irho) + fwd(i, j, k, irho)
                        gammaface = half*(gamma(i, j, k)+gamma(i+1, j, k))
                        por = pori(i, j, k)
                        call riemannflux_fast_b(left, leftd, right, rightd, flux, &
                        &                               fluxd)
                        branch = myIntStack(myIntPtr)
                        myIntPtr = myIntPtr - 1
                        if (branch .eq. 0) then
                            wd(i+1, j, k, itu1) = wd(i+1, j, k, itu1) + rightd(itu1)
                            rightd(itu1) = 0.0_8
                        end if
                        pd(i+1, j, k) = pd(i+1, j, k) + rightd(irhoe)
                        rightd(irhoe) = 0.0_8
                        wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) + rightd(ivz)
                        rightd(ivz) = 0.0_8
                        wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) + rightd(ivy)
                        rightd(ivy) = 0.0_8
                        wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) + rightd(ivx)
                        rightd(ivx) = 0.0_8
                        wd(i+1, j, k, irho) = wd(i+1, j, k, irho) + rightd(irho)
                        rightd(irho) = 0.0_8
                        branch = myIntStack(myIntPtr)
                        myIntPtr = myIntPtr - 1
                        if (branch .eq. 0) then
                            wd(i, j, k, itu1) = wd(i, j, k, itu1) + leftd(itu1)
                            leftd(itu1) = 0.0_8
                        end if
                        pd(i, j, k) = pd(i, j, k) + leftd(irhoe)
                        leftd(irhoe) = 0.0_8
                        wd(i, j, k, ivz) = wd(i, j, k, ivz) + leftd(ivz)
                        leftd(ivz) = 0.0_8
                        wd(i, j, k, ivy) = wd(i, j, k, ivy) + leftd(ivy)
                        leftd(ivy) = 0.0_8
                        wd(i, j, k, ivx) = wd(i, j, k, ivx) + leftd(ivx)
                        leftd(ivx) = 0.0_8
                        wd(i, j, k, irho) = wd(i, j, k, irho) + leftd(irho)
                        leftd(irho) = 0.0_8
                        branch = myIntStack(myIntPtr)
                        myIntPtr = myIntPtr - 1
                        if (branch .eq. 0) call popreal8(sface)
                    end do
                end do
            end do
        else
            ! store the density flux in the mass flow of the
            ! appropriate sliding mesh interface.
            !      ==================================================================
            !      ==================================================================
            !
            !         second order reconstruction of the left and right state.
            !         the three differences used in the, possibly nonlinear,
            !         interpolation are constructed here; the actual left and
            !         right states, or at least the differences from the first
            !         order interpolation, are computed in the subroutine
            !         leftrightstate.
            !
            ! fluxes in the i-direction.
            do k=2,kl
                do j=2,jl
                    do i=1,il
                        ! store the three differences used in the interpolation
                        ! in du1, du2, du3.
                        du1(irho) = w(i, j, k, irho) - w(i-1, j, k, irho)
                        du2(irho) = w(i+1, j, k, irho) - w(i, j, k, irho)
                        du3(irho) = w(i+2, j, k, irho) - w(i+1, j, k, irho)
                        du1(ivx) = w(i, j, k, ivx) - w(i-1, j, k, ivx)
                        du2(ivx) = w(i+1, j, k, ivx) - w(i, j, k, ivx)
                        du3(ivx) = w(i+2, j, k, ivx) - w(i+1, j, k, ivx)
                        du1(ivy) = w(i, j, k, ivy) - w(i-1, j, k, ivy)
                        du2(ivy) = w(i+1, j, k, ivy) - w(i, j, k, ivy)
                        du3(ivy) = w(i+2, j, k, ivy) - w(i+1, j, k, ivy)
                        du1(ivz) = w(i, j, k, ivz) - w(i-1, j, k, ivz)
                        du2(ivz) = w(i+1, j, k, ivz) - w(i, j, k, ivz)
                        du3(ivz) = w(i+2, j, k, ivz) - w(i+1, j, k, ivz)
                        du1(irhoe) = p(i, j, k) - p(i-1, j, k)
                        du2(irhoe) = p(i+1, j, k) - p(i, j, k)
                        du3(irhoe) = p(i+2, j, k) - p(i+1, j, k)
                        if (correctfork) then
                            du1(itu1) = w(i, j, k, itu1) - w(i-1, j, k, itu1)
                            du2(itu1) = w(i+1, j, k, itu1) - w(i, j, k, itu1)
                            du3(itu1) = w(i+2, j, k, itu1) - w(i+1, j, k, itu1)
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 0
                        else
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 1
                        end if
                        ! compute the differences from the first order scheme.
                        call leftrightstate(du1, du2, du3, rotmatrixi, left, right&
                        &                          )
                        ! add the first order part to the currently stored
                        ! differences, such that the correct state vector
                        ! is stored.
                        left(irho) = left(irho) + w(i, j, k, irho)
                        left(ivx) = left(ivx) + w(i, j, k, ivx)
                        left(ivy) = left(ivy) + w(i, j, k, ivy)
                        left(ivz) = left(ivz) + w(i, j, k, ivz)
                        left(irhoe) = left(irhoe) + p(i, j, k)
                        right(irho) = right(irho) + w(i+1, j, k, irho)
                        right(ivx) = right(ivx) + w(i+1, j, k, ivx)
                        right(ivy) = right(ivy) + w(i+1, j, k, ivy)
                        right(ivz) = right(ivz) + w(i+1, j, k, ivz)
                        right(irhoe) = right(irhoe) + p(i+1, j, k)
                        if (correctfork) then
                            left(itu1) = left(itu1) + w(i, j, k, itu1)
                            right(itu1) = right(itu1) + w(i+1, j, k, itu1)
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 0
                        else
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 1
                        end if
                        ! store the normal vector, the porosity and the
                        ! mesh velocity if present.
                        sx = si(i, j, k, 1)
                        sy = si(i, j, k, 2)
                        sz = si(i, j, k, 3)
                        if (addgridvelocities) then
                            sface = sfacei(i, j, k)
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 0
                        else
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 1
                        end if
                    end do
                end do
            end do
            ! store the density flux in the mass flow of the
            ! appropriate sliding mesh interface.
            ! fluxes in the j-direction.
            do k=2,kl
                do j=1,jl
                    do i=2,il
                        ! store the three differences used in the interpolation
                        ! in du1, du2, du3.
                        du1(irho) = w(i, j, k, irho) - w(i, j-1, k, irho)
                        du2(irho) = w(i, j+1, k, irho) - w(i, j, k, irho)
                        du3(irho) = w(i, j+2, k, irho) - w(i, j+1, k, irho)
                        du1(ivx) = w(i, j, k, ivx) - w(i, j-1, k, ivx)
                        du2(ivx) = w(i, j+1, k, ivx) - w(i, j, k, ivx)
                        du3(ivx) = w(i, j+2, k, ivx) - w(i, j+1, k, ivx)
                        du1(ivy) = w(i, j, k, ivy) - w(i, j-1, k, ivy)
                        du2(ivy) = w(i, j+1, k, ivy) - w(i, j, k, ivy)
                        du3(ivy) = w(i, j+2, k, ivy) - w(i, j+1, k, ivy)
                        du1(ivz) = w(i, j, k, ivz) - w(i, j-1, k, ivz)
                        du2(ivz) = w(i, j+1, k, ivz) - w(i, j, k, ivz)
                        du3(ivz) = w(i, j+2, k, ivz) - w(i, j+1, k, ivz)
                        du1(irhoe) = p(i, j, k) - p(i, j-1, k)
                        du2(irhoe) = p(i, j+1, k) - p(i, j, k)
                        du3(irhoe) = p(i, j+2, k) - p(i, j+1, k)
                        if (correctfork) then
                            du1(itu1) = w(i, j, k, itu1) - w(i, j-1, k, itu1)
                            du2(itu1) = w(i, j+1, k, itu1) - w(i, j, k, itu1)
                            du3(itu1) = w(i, j+2, k, itu1) - w(i, j+1, k, itu1)
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 0
                        else
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 1
                        end if
                        ! compute the differences from the first order scheme.
                        call leftrightstate(du1, du2, du3, rotmatrixj, left, right&
                        &                          )
                        ! add the first order part to the currently stored
                        ! differences, such that the correct state vector
                        ! is stored.
                        left(irho) = left(irho) + w(i, j, k, irho)
                        left(ivx) = left(ivx) + w(i, j, k, ivx)
                        left(ivy) = left(ivy) + w(i, j, k, ivy)
                        left(ivz) = left(ivz) + w(i, j, k, ivz)
                        left(irhoe) = left(irhoe) + p(i, j, k)
                        right(irho) = right(irho) + w(i, j+1, k, irho)
                        right(ivx) = right(ivx) + w(i, j+1, k, ivx)
                        right(ivy) = right(ivy) + w(i, j+1, k, ivy)
                        right(ivz) = right(ivz) + w(i, j+1, k, ivz)
                        right(irhoe) = right(irhoe) + p(i, j+1, k)
                        if (correctfork) then
                            left(itu1) = left(itu1) + w(i, j, k, itu1)
                            right(itu1) = right(itu1) + w(i, j+1, k, itu1)
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 0
                        else
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 1
                        end if
                        ! store the normal vector, the porosity and the
                        ! mesh velocity if present.
                        sx = sj(i, j, k, 1)
                        sy = sj(i, j, k, 2)
                        sz = sj(i, j, k, 3)
                        if (addgridvelocities) then
                            sface = sfacej(i, j, k)
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 0
                        else
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 1
                        end if
                    end do
                end do
            end do
            ! store the density flux in the mass flow of the
            ! appropriate sliding mesh interface.
            ! fluxes in the k-direction.
            do k=1,kl
                do j=2,jl
                    do i=2,il
                        ! store the three differences used in the interpolation
                        ! in du1, du2, du3.
                        du1(irho) = w(i, j, k, irho) - w(i, j, k-1, irho)
                        du2(irho) = w(i, j, k+1, irho) - w(i, j, k, irho)
                        du3(irho) = w(i, j, k+2, irho) - w(i, j, k+1, irho)
                        du1(ivx) = w(i, j, k, ivx) - w(i, j, k-1, ivx)
                        du2(ivx) = w(i, j, k+1, ivx) - w(i, j, k, ivx)
                        du3(ivx) = w(i, j, k+2, ivx) - w(i, j, k+1, ivx)
                        du1(ivy) = w(i, j, k, ivy) - w(i, j, k-1, ivy)
                        du2(ivy) = w(i, j, k+1, ivy) - w(i, j, k, ivy)
                        du3(ivy) = w(i, j, k+2, ivy) - w(i, j, k+1, ivy)
                        du1(ivz) = w(i, j, k, ivz) - w(i, j, k-1, ivz)
                        du2(ivz) = w(i, j, k+1, ivz) - w(i, j, k, ivz)
                        du3(ivz) = w(i, j, k+2, ivz) - w(i, j, k+1, ivz)
                        du1(irhoe) = p(i, j, k) - p(i, j, k-1)
                        du2(irhoe) = p(i, j, k+1) - p(i, j, k)
                        du3(irhoe) = p(i, j, k+2) - p(i, j, k+1)
                        if (correctfork) then
                            du1(itu1) = w(i, j, k, itu1) - w(i, j, k-1, itu1)
                            du2(itu1) = w(i, j, k+1, itu1) - w(i, j, k, itu1)
                            du3(itu1) = w(i, j, k+2, itu1) - w(i, j, k+1, itu1)
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 0
                        else
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 1
                        end if
                        ! compute the differences from the first order scheme.
                        call leftrightstate(du1, du2, du3, rotmatrixk, left, right&
                        &                          )
                        ! add the first order part to the currently stored
                        ! differences, such that the correct state vector
                        ! is stored.
                        left(irho) = left(irho) + w(i, j, k, irho)
                        left(ivx) = left(ivx) + w(i, j, k, ivx)
                        left(ivy) = left(ivy) + w(i, j, k, ivy)
                        left(ivz) = left(ivz) + w(i, j, k, ivz)
                        left(irhoe) = left(irhoe) + p(i, j, k)
                        right(irho) = right(irho) + w(i, j, k+1, irho)
                        right(ivx) = right(ivx) + w(i, j, k+1, ivx)
                        right(ivy) = right(ivy) + w(i, j, k+1, ivy)
                        right(ivz) = right(ivz) + w(i, j, k+1, ivz)
                        right(irhoe) = right(irhoe) + p(i, j, k+1)
                        if (correctfork) then
                            left(itu1) = left(itu1) + w(i, j, k, itu1)
                            right(itu1) = right(itu1) + w(i, j, k+1, itu1)
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 0
                        else
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 1
                        end if
                        ! store the normal vector, the porosity and the
                        ! mesh velocity if present.
                        sx = sk(i, j, k, 1)
                        sy = sk(i, j, k, 2)
                        sz = sk(i, j, k, 3)
                        if (addgridvelocities) then
                            sface = sfacek(i, j, k)
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 0
                        else
                            myIntPtr = myIntPtr + 1
                            myIntStack(myIntPtr) = 1
                        end if
                    end do
                end do
            end do
            fluxd = 0.0_8
            leftd = 0.0_8
            rightd = 0.0_8
            du1d = 0.0_8
            du2d = 0.0_8
            du3d = 0.0_8
            do k=kl,1,-1
                do j=jl,2,-1
                    do i=il,2,-1
                        fluxd(irhoe) = fluxd(irhoe) - fwd(i, j, k+1, irhoe)
                        fluxd(imz) = fluxd(imz) - fwd(i, j, k+1, imz)
                        fluxd(imy) = fluxd(imy) - fwd(i, j, k+1, imy)
                        fluxd(imx) = fluxd(imx) - fwd(i, j, k+1, imx)
                        fluxd(irho) = fluxd(irho) - fwd(i, j, k+1, irho)
                        fluxd(irhoe) = fluxd(irhoe) + fwd(i, j, k, irhoe)
                        fluxd(imz) = fluxd(imz) + fwd(i, j, k, imz)
                        fluxd(imy) = fluxd(imy) + fwd(i, j, k, imy)
                        fluxd(imx) = fluxd(imx) + fwd(i, j, k, imx)
                        fluxd(irho) = fluxd(irho) + fwd(i, j, k, irho)
                        gammaface = half*(gamma(i, j, k)+gamma(i, j, k+1))
                        por = pork(i, j, k)
                        call riemannflux_fast_b(left, leftd, right, rightd, flux, &
                        &                               fluxd)
                        branch = myIntStack(myIntPtr)
                        myIntPtr = myIntPtr - 1
                        if (branch .eq. 0) call popreal8(sface)
                        branch = myIntStack(myIntPtr)
                        myIntPtr = myIntPtr - 1
                        if (branch .eq. 0) then
                            wd(i, j, k+1, itu1) = wd(i, j, k+1, itu1) + rightd(itu1)
                            wd(i, j, k, itu1) = wd(i, j, k, itu1) + leftd(itu1)
                        end if
                        pd(i, j, k+1) = pd(i, j, k+1) + rightd(irhoe)
                        wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) + rightd(ivz)
                        wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) + rightd(ivy)
                        wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) + rightd(ivx)
                        wd(i, j, k+1, irho) = wd(i, j, k+1, irho) + rightd(irho)
                        pd(i, j, k) = pd(i, j, k) + leftd(irhoe)
                        wd(i, j, k, ivz) = wd(i, j, k, ivz) + leftd(ivz)
                        wd(i, j, k, ivy) = wd(i, j, k, ivy) + leftd(ivy)
                        wd(i, j, k, ivx) = wd(i, j, k, ivx) + leftd(ivx)
                        wd(i, j, k, irho) = wd(i, j, k, irho) + leftd(irho)
                        call leftrightstate_fast_b(du1, du1d, du2, du2d, du3, du3d&
                        &                                  , rotmatrixk, left, leftd, right, &
                        &                                  rightd)
                        branch = myIntStack(myIntPtr)
                        myIntPtr = myIntPtr - 1
                        if (branch .eq. 0) then
                            wd(i, j, k+2, itu1) = wd(i, j, k+2, itu1) + du3d(itu1)
                            wd(i, j, k+1, itu1) = wd(i, j, k+1, itu1) - du3d(itu1)
                            du3d(itu1) = 0.0_8
                            wd(i, j, k+1, itu1) = wd(i, j, k+1, itu1) + du2d(itu1)
                            wd(i, j, k, itu1) = wd(i, j, k, itu1) - du2d(itu1)
                            du2d(itu1) = 0.0_8
                            wd(i, j, k, itu1) = wd(i, j, k, itu1) + du1d(itu1)
                            wd(i, j, k-1, itu1) = wd(i, j, k-1, itu1) - du1d(itu1)
                            du1d(itu1) = 0.0_8
                        end if
                        pd(i, j, k+2) = pd(i, j, k+2) + du3d(irhoe)
                        pd(i, j, k+1) = pd(i, j, k+1) - du3d(irhoe)
                        du3d(irhoe) = 0.0_8
                        pd(i, j, k+1) = pd(i, j, k+1) + du2d(irhoe)
                        pd(i, j, k) = pd(i, j, k) - du2d(irhoe)
                        du2d(irhoe) = 0.0_8
                        pd(i, j, k) = pd(i, j, k) + du1d(irhoe)
                        pd(i, j, k-1) = pd(i, j, k-1) - du1d(irhoe)
                        du1d(irhoe) = 0.0_8
                        wd(i, j, k+2, ivz) = wd(i, j, k+2, ivz) + du3d(ivz)
                        wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) - du3d(ivz)
                        du3d(ivz) = 0.0_8
                        wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) + du2d(ivz)
                        wd(i, j, k, ivz) = wd(i, j, k, ivz) - du2d(ivz)
                        du2d(ivz) = 0.0_8
                        wd(i, j, k, ivz) = wd(i, j, k, ivz) + du1d(ivz)
                        wd(i, j, k-1, ivz) = wd(i, j, k-1, ivz) - du1d(ivz)
                        du1d(ivz) = 0.0_8
                        wd(i, j, k+2, ivy) = wd(i, j, k+2, ivy) + du3d(ivy)
                        wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) - du3d(ivy)
                        du3d(ivy) = 0.0_8
                        wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) + du2d(ivy)
                        wd(i, j, k, ivy) = wd(i, j, k, ivy) - du2d(ivy)
                        du2d(ivy) = 0.0_8
                        wd(i, j, k, ivy) = wd(i, j, k, ivy) + du1d(ivy)
                        wd(i, j, k-1, ivy) = wd(i, j, k-1, ivy) - du1d(ivy)
                        du1d(ivy) = 0.0_8
                        wd(i, j, k+2, ivx) = wd(i, j, k+2, ivx) + du3d(ivx)
                        wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) - du3d(ivx)
                        du3d(ivx) = 0.0_8
                        wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) + du2d(ivx)
                        wd(i, j, k, ivx) = wd(i, j, k, ivx) - du2d(ivx)
                        du2d(ivx) = 0.0_8
                        wd(i, j, k, ivx) = wd(i, j, k, ivx) + du1d(ivx)
                        wd(i, j, k-1, ivx) = wd(i, j, k-1, ivx) - du1d(ivx)
                        du1d(ivx) = 0.0_8
                        wd(i, j, k+2, irho) = wd(i, j, k+2, irho) + du3d(irho)
                        wd(i, j, k+1, irho) = wd(i, j, k+1, irho) - du3d(irho)
                        du3d(irho) = 0.0_8
                        wd(i, j, k+1, irho) = wd(i, j, k+1, irho) + du2d(irho)
                        wd(i, j, k, irho) = wd(i, j, k, irho) - du2d(irho)
                        du2d(irho) = 0.0_8
                        wd(i, j, k, irho) = wd(i, j, k, irho) + du1d(irho)
                        wd(i, j, k-1, irho) = wd(i, j, k-1, irho) - du1d(irho)
                        du1d(irho) = 0.0_8
                    end do
                end do
            end do
            do k=kl,2,-1
                do j=jl,1,-1
                    do i=il,2,-1
                        fluxd(irhoe) = fluxd(irhoe) - fwd(i, j+1, k, irhoe)
                        fluxd(imz) = fluxd(imz) - fwd(i, j+1, k, imz)
                        fluxd(imy) = fluxd(imy) - fwd(i, j+1, k, imy)
                        fluxd(imx) = fluxd(imx) - fwd(i, j+1, k, imx)
                        fluxd(irho) = fluxd(irho) - fwd(i, j+1, k, irho)
                        fluxd(irhoe) = fluxd(irhoe) + fwd(i, j, k, irhoe)
                        fluxd(imz) = fluxd(imz) + fwd(i, j, k, imz)
                        fluxd(imy) = fluxd(imy) + fwd(i, j, k, imy)
                        fluxd(imx) = fluxd(imx) + fwd(i, j, k, imx)
                        fluxd(irho) = fluxd(irho) + fwd(i, j, k, irho)
                        gammaface = half*(gamma(i, j, k)+gamma(i, j+1, k))
                        por = porj(i, j, k)
                        call riemannflux_fast_b(left, leftd, right, rightd, flux, &
                        &                               fluxd)
                        branch = myIntStack(myIntPtr)
                        myIntPtr = myIntPtr - 1
                        if (branch .eq. 0) call popreal8(sface)
                        branch = myIntStack(myIntPtr)
                        myIntPtr = myIntPtr - 1
                        if (branch .eq. 0) then
                            wd(i, j+1, k, itu1) = wd(i, j+1, k, itu1) + rightd(itu1)
                            wd(i, j, k, itu1) = wd(i, j, k, itu1) + leftd(itu1)
                        end if
                        pd(i, j+1, k) = pd(i, j+1, k) + rightd(irhoe)
                        wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) + rightd(ivz)
                        wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) + rightd(ivy)
                        wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) + rightd(ivx)
                        wd(i, j+1, k, irho) = wd(i, j+1, k, irho) + rightd(irho)
                        pd(i, j, k) = pd(i, j, k) + leftd(irhoe)
                        wd(i, j, k, ivz) = wd(i, j, k, ivz) + leftd(ivz)
                        wd(i, j, k, ivy) = wd(i, j, k, ivy) + leftd(ivy)
                        wd(i, j, k, ivx) = wd(i, j, k, ivx) + leftd(ivx)
                        wd(i, j, k, irho) = wd(i, j, k, irho) + leftd(irho)
                        call leftrightstate_fast_b(du1, du1d, du2, du2d, du3, du3d&
                        &                                  , rotmatrixj, left, leftd, right, &
                        &                                  rightd)
                        branch = myIntStack(myIntPtr)
                        myIntPtr = myIntPtr - 1
                        if (branch .eq. 0) then
                            wd(i, j+2, k, itu1) = wd(i, j+2, k, itu1) + du3d(itu1)
                            wd(i, j+1, k, itu1) = wd(i, j+1, k, itu1) - du3d(itu1)
                            du3d(itu1) = 0.0_8
                            wd(i, j+1, k, itu1) = wd(i, j+1, k, itu1) + du2d(itu1)
                            wd(i, j, k, itu1) = wd(i, j, k, itu1) - du2d(itu1)
                            du2d(itu1) = 0.0_8
                            wd(i, j, k, itu1) = wd(i, j, k, itu1) + du1d(itu1)
                            wd(i, j-1, k, itu1) = wd(i, j-1, k, itu1) - du1d(itu1)
                            du1d(itu1) = 0.0_8
                        end if
                        pd(i, j+2, k) = pd(i, j+2, k) + du3d(irhoe)
                        pd(i, j+1, k) = pd(i, j+1, k) - du3d(irhoe)
                        du3d(irhoe) = 0.0_8
                        pd(i, j+1, k) = pd(i, j+1, k) + du2d(irhoe)
                        pd(i, j, k) = pd(i, j, k) - du2d(irhoe)
                        du2d(irhoe) = 0.0_8
                        pd(i, j, k) = pd(i, j, k) + du1d(irhoe)
                        pd(i, j-1, k) = pd(i, j-1, k) - du1d(irhoe)
                        du1d(irhoe) = 0.0_8
                        wd(i, j+2, k, ivz) = wd(i, j+2, k, ivz) + du3d(ivz)
                        wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) - du3d(ivz)
                        du3d(ivz) = 0.0_8
                        wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) + du2d(ivz)
                        wd(i, j, k, ivz) = wd(i, j, k, ivz) - du2d(ivz)
                        du2d(ivz) = 0.0_8
                        wd(i, j, k, ivz) = wd(i, j, k, ivz) + du1d(ivz)
                        wd(i, j-1, k, ivz) = wd(i, j-1, k, ivz) - du1d(ivz)
                        du1d(ivz) = 0.0_8
                        wd(i, j+2, k, ivy) = wd(i, j+2, k, ivy) + du3d(ivy)
                        wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) - du3d(ivy)
                        du3d(ivy) = 0.0_8
                        wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) + du2d(ivy)
                        wd(i, j, k, ivy) = wd(i, j, k, ivy) - du2d(ivy)
                        du2d(ivy) = 0.0_8
                        wd(i, j, k, ivy) = wd(i, j, k, ivy) + du1d(ivy)
                        wd(i, j-1, k, ivy) = wd(i, j-1, k, ivy) - du1d(ivy)
                        du1d(ivy) = 0.0_8
                        wd(i, j+2, k, ivx) = wd(i, j+2, k, ivx) + du3d(ivx)
                        wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) - du3d(ivx)
                        du3d(ivx) = 0.0_8
                        wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) + du2d(ivx)
                        wd(i, j, k, ivx) = wd(i, j, k, ivx) - du2d(ivx)
                        du2d(ivx) = 0.0_8
                        wd(i, j, k, ivx) = wd(i, j, k, ivx) + du1d(ivx)
                        wd(i, j-1, k, ivx) = wd(i, j-1, k, ivx) - du1d(ivx)
                        du1d(ivx) = 0.0_8
                        wd(i, j+2, k, irho) = wd(i, j+2, k, irho) + du3d(irho)
                        wd(i, j+1, k, irho) = wd(i, j+1, k, irho) - du3d(irho)
                        du3d(irho) = 0.0_8
                        wd(i, j+1, k, irho) = wd(i, j+1, k, irho) + du2d(irho)
                        wd(i, j, k, irho) = wd(i, j, k, irho) - du2d(irho)
                        du2d(irho) = 0.0_8
                        wd(i, j, k, irho) = wd(i, j, k, irho) + du1d(irho)
                        wd(i, j-1, k, irho) = wd(i, j-1, k, irho) - du1d(irho)
                        du1d(irho) = 0.0_8
                    end do
                end do
            end do
            do k=kl,2,-1
                do j=jl,2,-1
                    do i=il,1,-1
                        fluxd(irhoe) = fluxd(irhoe) - fwd(i+1, j, k, irhoe)
                        fluxd(imz) = fluxd(imz) - fwd(i+1, j, k, imz)
                        fluxd(imy) = fluxd(imy) - fwd(i+1, j, k, imy)
                        fluxd(imx) = fluxd(imx) - fwd(i+1, j, k, imx)
                        fluxd(irho) = fluxd(irho) - fwd(i+1, j, k, irho)
                        fluxd(irhoe) = fluxd(irhoe) + fwd(i, j, k, irhoe)
                        fluxd(imz) = fluxd(imz) + fwd(i, j, k, imz)
                        fluxd(imy) = fluxd(imy) + fwd(i, j, k, imy)
                        fluxd(imx) = fluxd(imx) + fwd(i, j, k, imx)
                        fluxd(irho) = fluxd(irho) + fwd(i, j, k, irho)
                        gammaface = half*(gamma(i, j, k)+gamma(i+1, j, k))
                        por = pori(i, j, k)
                        call riemannflux_fast_b(left, leftd, right, rightd, flux, &
                        &                               fluxd)
                        branch = myIntStack(myIntPtr)
                        myIntPtr = myIntPtr - 1
                        if (branch .eq. 0) call popreal8(sface)
                        branch = myIntStack(myIntPtr)
                        myIntPtr = myIntPtr - 1
                        if (branch .eq. 0) then
                            wd(i+1, j, k, itu1) = wd(i+1, j, k, itu1) + rightd(itu1)
                            wd(i, j, k, itu1) = wd(i, j, k, itu1) + leftd(itu1)
                        end if
                        pd(i+1, j, k) = pd(i+1, j, k) + rightd(irhoe)
                        wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) + rightd(ivz)
                        wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) + rightd(ivy)
                        wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) + rightd(ivx)
                        wd(i+1, j, k, irho) = wd(i+1, j, k, irho) + rightd(irho)
                        pd(i, j, k) = pd(i, j, k) + leftd(irhoe)
                        wd(i, j, k, ivz) = wd(i, j, k, ivz) + leftd(ivz)
                        wd(i, j, k, ivy) = wd(i, j, k, ivy) + leftd(ivy)
                        wd(i, j, k, ivx) = wd(i, j, k, ivx) + leftd(ivx)
                        wd(i, j, k, irho) = wd(i, j, k, irho) + leftd(irho)
                        call leftrightstate_fast_b(du1, du1d, du2, du2d, du3, du3d&
                        &                                  , rotmatrixi, left, leftd, right, &
                        &                                  rightd)
                        branch = myIntStack(myIntPtr)
                        myIntPtr = myIntPtr - 1
                        if (branch .eq. 0) then
                            wd(i+2, j, k, itu1) = wd(i+2, j, k, itu1) + du3d(itu1)
                            wd(i+1, j, k, itu1) = wd(i+1, j, k, itu1) - du3d(itu1)
                            du3d(itu1) = 0.0_8
                            wd(i+1, j, k, itu1) = wd(i+1, j, k, itu1) + du2d(itu1)
                            wd(i, j, k, itu1) = wd(i, j, k, itu1) - du2d(itu1)
                            du2d(itu1) = 0.0_8
                            wd(i, j, k, itu1) = wd(i, j, k, itu1) + du1d(itu1)
                            wd(i-1, j, k, itu1) = wd(i-1, j, k, itu1) - du1d(itu1)
                            du1d(itu1) = 0.0_8
                        end if
                        pd(i+2, j, k) = pd(i+2, j, k) + du3d(irhoe)
                        pd(i+1, j, k) = pd(i+1, j, k) - du3d(irhoe)
                        du3d(irhoe) = 0.0_8
                        pd(i+1, j, k) = pd(i+1, j, k) + du2d(irhoe)
                        pd(i, j, k) = pd(i, j, k) - du2d(irhoe)
                        du2d(irhoe) = 0.0_8
                        pd(i, j, k) = pd(i, j, k) + du1d(irhoe)
                        pd(i-1, j, k) = pd(i-1, j, k) - du1d(irhoe)
                        du1d(irhoe) = 0.0_8
                        wd(i+2, j, k, ivz) = wd(i+2, j, k, ivz) + du3d(ivz)
                        wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) - du3d(ivz)
                        du3d(ivz) = 0.0_8
                        wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) + du2d(ivz)
                        wd(i, j, k, ivz) = wd(i, j, k, ivz) - du2d(ivz)
                        du2d(ivz) = 0.0_8
                        wd(i, j, k, ivz) = wd(i, j, k, ivz) + du1d(ivz)
                        wd(i-1, j, k, ivz) = wd(i-1, j, k, ivz) - du1d(ivz)
                        du1d(ivz) = 0.0_8
                        wd(i+2, j, k, ivy) = wd(i+2, j, k, ivy) + du3d(ivy)
                        wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) - du3d(ivy)
                        du3d(ivy) = 0.0_8
                        wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) + du2d(ivy)
                        wd(i, j, k, ivy) = wd(i, j, k, ivy) - du2d(ivy)
                        du2d(ivy) = 0.0_8
                        wd(i, j, k, ivy) = wd(i, j, k, ivy) + du1d(ivy)
                        wd(i-1, j, k, ivy) = wd(i-1, j, k, ivy) - du1d(ivy)
                        du1d(ivy) = 0.0_8
                        wd(i+2, j, k, ivx) = wd(i+2, j, k, ivx) + du3d(ivx)
                        wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) - du3d(ivx)
                        du3d(ivx) = 0.0_8
                        wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) + du2d(ivx)
                        wd(i, j, k, ivx) = wd(i, j, k, ivx) - du2d(ivx)
                        du2d(ivx) = 0.0_8
                        wd(i, j, k, ivx) = wd(i, j, k, ivx) + du1d(ivx)
                        wd(i-1, j, k, ivx) = wd(i-1, j, k, ivx) - du1d(ivx)
                        du1d(ivx) = 0.0_8
                        wd(i+2, j, k, irho) = wd(i+2, j, k, irho) + du3d(irho)
                        wd(i+1, j, k, irho) = wd(i+1, j, k, irho) - du3d(irho)
                        du3d(irho) = 0.0_8
                        wd(i+1, j, k, irho) = wd(i+1, j, k, irho) + du2d(irho)
                        wd(i, j, k, irho) = wd(i, j, k, irho) - du2d(irho)
                        du2d(irho) = 0.0_8
                        wd(i, j, k, irho) = wd(i, j, k, irho) + du1d(irho)
                        wd(i-1, j, k, irho) = wd(i-1, j, k, irho) - du1d(irho)
                        du1d(irho) = 0.0_8
                    end do
                end do
            end do
        end if
        do k=kl,2,-1
            do j=jl,2,-1
                do i=il,2,-1
                    fwd(i, j, k, irhoe) = sfil*fwd(i, j, k, irhoe)
                    fwd(i, j, k, imz) = sfil*fwd(i, j, k, imz)
                    fwd(i, j, k, imy) = sfil*fwd(i, j, k, imy)
                    fwd(i, j, k, imx) = sfil*fwd(i, j, k, imx)
                    fwd(i, j, k, irho) = sfil*fwd(i, j, k, irho)
                end do
            end do
        end do
    end if

contains
    !  differentiation of leftrightstate in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
    !   gradient     of useful results: left right du1 du2 du3
    !   with respect to varying inputs: left right du1 du2 du3
    ! store the density flux in the mass flow of the
    ! appropriate sliding mesh interface.
    !      ==================================================================
    subroutine leftrightstate_fast_b(du1, du1d, du2, du2d, du3, du3d, &
        &     rotmatrix, left, leftd, right, rightd)
        implicit none
        !
        !        local parameter.
        !
        real(kind=realtype), parameter :: epslim=1.e-10_realtype
        !
        !        subroutine arguments.
        !
        real(kind=realtype), dimension(:), intent(inout) :: du1, du2, du3
        real(kind=realtype), dimension(:), intent(inout) :: du1d
        real(kind=realtype), dimension(:) :: left, right
        real(kind=realtype), dimension(:) :: leftd, rightd
        real(kind=realtype), dimension(:, :, :, :, :), pointer :: &
        &     rotmatrix
        !
        !        local variables.
        !
        integer(kind=inttype) :: l
        real(kind=realtype) :: rl1, rl2, rr1, rr2, tmp, dvx, dvy, dvz
        real(kind=realtype) :: rl1d, rl2d, rr1d, rr2d, tmpd, dvxd, dvyd, &
        &     dvzd
        real(kind=realtype), dimension(3, 3) :: rot
        intrinsic abs
        intrinsic max
        intrinsic sign
        intrinsic min
        integer :: branch
        real(kind=realtype), dimension(:), intent(inout) :: du3d
        real(kind=realtype), dimension(:), intent(inout) :: du2d
        real(kind=realtype) :: temp3
        real(kind=realtype) :: temp2
        real(kind=realtype) :: temp1
        real(kind=realtype) :: temp0
        real(kind=realtype) :: x6d
        real(kind=realtype) :: y4d
        real(kind=realtype) :: max2d
        real(kind=realtype) :: max5d
        real(kind=realtype) :: x6
        real(kind=realtype) :: x5
        real(kind=realtype) :: x4
        real(kind=realtype) :: x3
        real(kind=realtype) :: x2
        real(kind=realtype) :: x2d
        real(kind=realtype) :: x1
        real(kind=realtype) :: x5d
        real(kind=realtype) :: y3d
        real(kind=realtype) :: tempd
        real(kind=realtype) :: max4d
        real(kind=realtype) :: tempd8
        real(kind=realtype) :: tempd7
        real(kind=realtype) :: tempd6
        real(kind=realtype) :: tempd5
        real(kind=realtype) :: tempd4
        real(kind=realtype) :: tempd3
        real(kind=realtype) :: tempd2
        real(kind=realtype) :: tempd1
        real(kind=realtype) :: max7d
        real(kind=realtype) :: tempd0
        real(kind=realtype) :: x1d
        real(kind=realtype) :: x4d
        real(kind=realtype) :: y2d
        real(kind=realtype) :: max3d
        real(kind=realtype) :: max7
        real(kind=realtype) :: max6
        real(kind=realtype) :: max6d
        real(kind=realtype) :: max5
        real(kind=realtype) :: max4
        real(kind=realtype) :: temp
        real(kind=realtype) :: max3
        real(kind=realtype) :: max2
        real(kind=realtype) :: y4
        real(kind=realtype) :: y3
        real(kind=realtype) :: y2
        real(kind=realtype) :: x3d
        real(kind=realtype) :: y1
        real(kind=realtype) :: y1d
        real(kind=realtype) :: temp4
        ! check if the velocity components should be transformed to
        ! the cylindrical frame.
        if (rotationalperiodic) then
            ! store the rotation matrix a bit easier. note that the i,j,k
            ! come from the main subroutine.
            rot(1, 1) = rotmatrix(i, j, k, 1, 1)
            rot(1, 2) = rotmatrix(i, j, k, 1, 2)
            rot(1, 3) = rotmatrix(i, j, k, 1, 3)
            rot(2, 1) = rotmatrix(i, j, k, 2, 1)
            rot(2, 2) = rotmatrix(i, j, k, 2, 2)
            rot(2, 3) = rotmatrix(i, j, k, 2, 3)
            rot(3, 1) = rotmatrix(i, j, k, 3, 1)
            rot(3, 2) = rotmatrix(i, j, k, 3, 2)
            rot(3, 3) = rotmatrix(i, j, k, 3, 3)
            ! apply the transformation to the velocity components
            ! of du1, du2 and du3.
            dvx = du1(ivx)
            dvy = du1(ivy)
            dvz = du1(ivz)
            du1(ivx) = rot(1, 1)*dvx + rot(1, 2)*dvy + rot(1, 3)*dvz
            du1(ivy) = rot(2, 1)*dvx + rot(2, 2)*dvy + rot(2, 3)*dvz
            du1(ivz) = rot(3, 1)*dvx + rot(3, 2)*dvy + rot(3, 3)*dvz
            dvx = du2(ivx)
            dvy = du2(ivy)
            dvz = du2(ivz)
            du2(ivx) = rot(1, 1)*dvx + rot(1, 2)*dvy + rot(1, 3)*dvz
            du2(ivy) = rot(2, 1)*dvx + rot(2, 2)*dvy + rot(2, 3)*dvz
            du2(ivz) = rot(3, 1)*dvx + rot(3, 2)*dvy + rot(3, 3)*dvz
            dvx = du3(ivx)
            dvy = du3(ivy)
            dvz = du3(ivz)
            du3(ivx) = rot(1, 1)*dvx + rot(1, 2)*dvy + rot(1, 3)*dvz
            du3(ivy) = rot(2, 1)*dvx + rot(2, 2)*dvy + rot(2, 3)*dvz
            du3(ivz) = rot(3, 1)*dvx + rot(3, 2)*dvy + rot(3, 3)*dvz
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 0
        else
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 1
        end if
        ! determine the limiter used.
        select case  (limused)
        case (nolimiter)
            call pushcontrol2b(1)
        case (vanalbeda)
            !          ==============================================================
            ! nonlinear interpolation using the van albeda limiter.
            ! loop over the number of variables to be interpolated.
            do l=1,nwint
                if (du2(l) .ge. 0.) then
                    x1 = du2(l)
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    x1 = -du2(l)
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (x1 .lt. epslim) then
                    max2 = epslim
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    max2 = x1
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                ! compute the limiter argument rl1, rl2, rr1 and rr2.
                ! note the cut off to 0.0.
                tmp = one/sign(max2, du2(l))
                if (du1(l) .ge. 0.) then
                    x3 = du1(l)
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    x3 = -du1(l)
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (x3 .lt. epslim) then
                    max4 = epslim
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    max4 = x3
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                y1 = du2(l)/sign(max4, du1(l))
                if (zero .lt. y1) then
                    rl1 = y1
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    rl1 = zero
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (zero .lt. du1(l)*tmp) then
                    rl2 = du1(l)*tmp
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    rl2 = zero
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (zero .lt. du3(l)*tmp) then
                    rr1 = du3(l)*tmp
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    rr1 = zero
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (du3(l) .ge. 0.) then
                    x4 = du3(l)
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    x4 = -du3(l)
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (x4 .lt. epslim) then
                    max5 = epslim
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    max5 = x4
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                y2 = du2(l)/sign(max5, du3(l))
                if (zero .lt. y2) then
                    rr2 = y2
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    rr2 = zero
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                ! compute the corresponding limiter values.
                rl1 = rl1*(rl1+one)/(rl1*rl1+one)
                rl2 = rl2*(rl2+one)/(rl2*rl2+one)
                rr1 = rr1*(rr1+one)/(rr1*rr1+one)
                rr2 = rr2*(rr2+one)/(rr2*rr2+one)
                ! compute the nonlinear corrections to the first order
                ! scheme.
            end do
            call pushcontrol2b(2)
        case (minmod)
            !          ==============================================================
            ! nonlinear interpolation using the minmod limiter.
            ! loop over the number of variables to be interpolated.
            do l=1,nwint
                if (du2(l) .ge. 0.) then
                    x2 = du2(l)
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    x2 = -du2(l)
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (x2 .lt. epslim) then
                    max3 = epslim
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    max3 = x2
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                ! compute the limiter argument rl1, rl2, rr1 and rr2.
                ! note the cut off to 0.0.
                tmp = one/sign(max3, du2(l))
                if (du1(l) .ge. 0.) then
                    x5 = du1(l)
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    x5 = -du1(l)
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (x5 .lt. epslim) then
                    max6 = epslim
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    max6 = x5
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                y3 = du2(l)/sign(max6, du1(l))
                if (zero .lt. y3) then
                    rl1 = y3
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    rl1 = zero
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (zero .lt. du1(l)*tmp) then
                    rl2 = du1(l)*tmp
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    rl2 = zero
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (zero .lt. du3(l)*tmp) then
                    rr1 = du3(l)*tmp
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    rr1 = zero
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (du3(l) .ge. 0.) then
                    x6 = du3(l)
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    x6 = -du3(l)
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (x6 .lt. epslim) then
                    max7 = epslim
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    max7 = x6
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                y4 = du2(l)/sign(max7, du3(l))
                if (zero .lt. y4) then
                    rr2 = y4
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    rr2 = zero
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (one .gt. factminmod*rl1) then
                    rl1 = factminmod*rl1
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    rl1 = one
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (one .gt. factminmod*rl2) then
                    rl2 = factminmod*rl2
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    rl2 = one
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (one .gt. factminmod*rr1) then
                    rr1 = factminmod*rr1
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    rr1 = one
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (one .gt. factminmod*rr2) then
                    rr2 = factminmod*rr2
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    rr2 = one
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
            end do
            call pushcontrol2b(3)
        case default
            call pushcontrol2b(0)
        end select
        ! in case only a first order scheme must be used for the
        ! turbulent transport equations, set the correction for the
        ! turbulent kinetic energy to 0.
        if (firstorderk) then
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 0
        else
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 1
        end if
        ! for rotational periodic problems transform the velocity
        ! differences back to cartesian again. note that now the
        ! transpose of the rotation matrix must be used.
        if (rotationalperiodic) then
            dvxd = rot(1, 3)*rightd(ivz)
            dvyd = rot(2, 3)*rightd(ivz)
            dvzd = rot(3, 3)*rightd(ivz)
            rightd(ivz) = 0.0_8
            dvxd = dvxd + rot(1, 2)*rightd(ivy)
            dvyd = dvyd + rot(2, 2)*rightd(ivy)
            dvzd = dvzd + rot(3, 2)*rightd(ivy)
            rightd(ivy) = 0.0_8
            dvxd = dvxd + rot(1, 1)*rightd(ivx)
            dvyd = dvyd + rot(2, 1)*rightd(ivx)
            dvzd = dvzd + rot(3, 1)*rightd(ivx)
            rightd(ivx) = 0.0_8
            rightd(ivz) = rightd(ivz) + dvzd
            rightd(ivy) = rightd(ivy) + dvyd
            rightd(ivx) = rightd(ivx) + dvxd
            dvxd = rot(1, 3)*leftd(ivz)
            dvyd = rot(2, 3)*leftd(ivz)
            dvzd = rot(3, 3)*leftd(ivz)
            leftd(ivz) = 0.0_8
            dvxd = dvxd + rot(1, 2)*leftd(ivy)
            dvyd = dvyd + rot(2, 2)*leftd(ivy)
            dvzd = dvzd + rot(3, 2)*leftd(ivy)
            leftd(ivy) = 0.0_8
            dvxd = dvxd + rot(1, 1)*leftd(ivx)
            dvyd = dvyd + rot(2, 1)*leftd(ivx)
            dvzd = dvzd + rot(3, 1)*leftd(ivx)
            leftd(ivx) = 0.0_8
            leftd(ivz) = leftd(ivz) + dvzd
            leftd(ivy) = leftd(ivy) + dvyd
            leftd(ivx) = leftd(ivx) + dvxd
        end if
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) then
            rightd(itu1) = 0.0_8
            leftd(itu1) = 0.0_8
        end if
        call popcontrol2b(branch)
        if (branch .lt. 2) then
            if (branch .ne. 0) then
                do l=nwint,1,-1
                    du3d(l) = du3d(l) - omk*rightd(l)
                    du2d(l) = du2d(l) + opk*leftd(l) - opk*rightd(l)
                    rightd(l) = 0.0_8
                    du1d(l) = du1d(l) + omk*leftd(l)
                    leftd(l) = 0.0_8
                end do
            end if
        else if (branch .eq. 2) then
            do l=nwint,1,-1
                rr1d = -(opk*du2(l)*rightd(l))
                du2d(l) = du2d(l) + opk*rl2*leftd(l) - opk*rr1*rightd(l)
                rr2d = -(omk*du3(l)*rightd(l))
                du3d(l) = du3d(l) - omk*rr2*rightd(l)
                rightd(l) = 0.0_8
                rl1d = omk*du1(l)*leftd(l)
                du1d(l) = du1d(l) + omk*rl1*leftd(l)
                rl2d = opk*du2(l)*leftd(l)
                leftd(l) = 0.0_8
                tempd2 = rr2d/(one+rr2**2)
                rr2d = (2*rr2-rr2**2*(one+rr2)*2/(one+rr2**2)+one)*tempd2
                tempd3 = rr1d/(one+rr1**2)
                rr1d = (2*rr1-rr1**2*(one+rr1)*2/(one+rr1**2)+one)*tempd3
                tempd4 = rl2d/(one+rl2**2)
                rl2d = (2*rl2-rl2**2*(one+rl2)*2/(one+rl2**2)+one)*tempd4
                tempd5 = rl1d/(one+rl1**2)
                rl1d = (2*rl1-rl1**2*(one+rl1)*2/(one+rl1**2)+one)*tempd5
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    y2d = rr2d
                else
                    y2d = 0.0_8
                end if
                temp1 = sign(max5, du3(l))
                tempd1 = -(du2(l)*y2d/temp1**2)
                du2d(l) = du2d(l) + y2d/temp1
                max5d = sign(1.d0, max5*du3(l))*tempd1
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    x4d = 0.0_8
                else
                    x4d = max5d
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    du3d(l) = du3d(l) + x4d
                else
                    du3d(l) = du3d(l) - x4d
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    du3d(l) = du3d(l) + tmp*rr1d
                    tmpd = du3(l)*rr1d
                else
                    tmpd = 0.0_8
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    du1d(l) = du1d(l) + tmp*rl2d
                    tmpd = tmpd + du1(l)*rl2d
                else
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    y1d = rl1d
                else
                    y1d = 0.0_8
                end if
                temp0 = sign(max4, du1(l))
                tempd0 = -(du2(l)*y1d/temp0**2)
                du2d(l) = du2d(l) + y1d/temp0
                max4d = sign(1.d0, max4*du1(l))*tempd0
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    x3d = 0.0_8
                else
                    x3d = max4d
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    du1d(l) = du1d(l) + x3d
                else
                    du1d(l) = du1d(l) - x3d
                end if
                temp = sign(max2, du2(l))
                tempd = -(one*tmpd/temp**2)
                max2d = sign(1.d0, max2*du2(l))*tempd
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    x1d = 0.0_8
                else
                    x1d = max2d
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    du2d(l) = du2d(l) + x1d
                else
                    du2d(l) = du2d(l) - x1d
                end if
            end do
        else
            do l=nwint,1,-1
                rr1d = -(opk*du2(l)*rightd(l))
                du2d(l) = du2d(l) + opk*rl2*leftd(l) - opk*rr1*rightd(l)
                rr2d = -(omk*du3(l)*rightd(l))
                du3d(l) = du3d(l) - omk*rr2*rightd(l)
                rightd(l) = 0.0_8
                rl1d = omk*du1(l)*leftd(l)
                du1d(l) = du1d(l) + omk*rl1*leftd(l)
                rl2d = opk*du2(l)*leftd(l)
                leftd(l) = 0.0_8
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    rr2d = factminmod*rr2d
                else
                    rr2d = 0.0_8
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    rr1d = factminmod*rr1d
                else
                    rr1d = 0.0_8
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    rl2d = factminmod*rl2d
                else
                    rl2d = 0.0_8
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    rl1d = factminmod*rl1d
                else
                    rl1d = 0.0_8
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    y4d = rr2d
                else
                    y4d = 0.0_8
                end if
                temp4 = sign(max7, du3(l))
                tempd8 = -(du2(l)*y4d/temp4**2)
                du2d(l) = du2d(l) + y4d/temp4
                max7d = sign(1.d0, max7*du3(l))*tempd8
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    x6d = 0.0_8
                else
                    x6d = max7d
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    du3d(l) = du3d(l) + x6d
                else
                    du3d(l) = du3d(l) - x6d
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    du3d(l) = du3d(l) + tmp*rr1d
                    tmpd = du3(l)*rr1d
                else
                    tmpd = 0.0_8
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    du1d(l) = du1d(l) + tmp*rl2d
                    tmpd = tmpd + du1(l)*rl2d
                else
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    y3d = rl1d
                else
                    y3d = 0.0_8
                end if
                temp3 = sign(max6, du1(l))
                tempd7 = -(du2(l)*y3d/temp3**2)
                du2d(l) = du2d(l) + y3d/temp3
                max6d = sign(1.d0, max6*du1(l))*tempd7
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    x5d = 0.0_8
                else
                    x5d = max6d
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    du1d(l) = du1d(l) + x5d
                else
                    du1d(l) = du1d(l) - x5d
                end if
                temp2 = sign(max3, du2(l))
                tempd6 = -(one*tmpd/temp2**2)
                max3d = sign(1.d0, max3*du2(l))*tempd6
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    x2d = 0.0_8
                else
                    x2d = max3d
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    du2d(l) = du2d(l) + x2d
                else
                    du2d(l) = du2d(l) - x2d
                end if
            end do
        end if
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) then
            dvxd = rot(3, 1)*du3d(ivz)
            dvyd = rot(3, 2)*du3d(ivz)
            dvzd = rot(3, 3)*du3d(ivz)
            du3d(ivz) = 0.0_8
            dvxd = dvxd + rot(2, 1)*du3d(ivy)
            dvyd = dvyd + rot(2, 2)*du3d(ivy)
            dvzd = dvzd + rot(2, 3)*du3d(ivy)
            du3d(ivy) = 0.0_8
            dvxd = dvxd + rot(1, 1)*du3d(ivx)
            dvyd = dvyd + rot(1, 2)*du3d(ivx)
            dvzd = dvzd + rot(1, 3)*du3d(ivx)
            du3d(ivx) = 0.0_8
            du3d(ivz) = du3d(ivz) + dvzd
            du3d(ivy) = du3d(ivy) + dvyd
            du3d(ivx) = du3d(ivx) + dvxd
            dvxd = rot(3, 1)*du2d(ivz)
            dvyd = rot(3, 2)*du2d(ivz)
            dvzd = rot(3, 3)*du2d(ivz)
            du2d(ivz) = 0.0_8
            dvxd = dvxd + rot(2, 1)*du2d(ivy)
            dvyd = dvyd + rot(2, 2)*du2d(ivy)
            dvzd = dvzd + rot(2, 3)*du2d(ivy)
            du2d(ivy) = 0.0_8
            dvxd = dvxd + rot(1, 1)*du2d(ivx)
            dvyd = dvyd + rot(1, 2)*du2d(ivx)
            dvzd = dvzd + rot(1, 3)*du2d(ivx)
            du2d(ivx) = 0.0_8
            du2d(ivz) = du2d(ivz) + dvzd
            du2d(ivy) = du2d(ivy) + dvyd
            du2d(ivx) = du2d(ivx) + dvxd
            dvxd = rot(3, 1)*du1d(ivz)
            dvyd = rot(3, 2)*du1d(ivz)
            dvzd = rot(3, 3)*du1d(ivz)
            du1d(ivz) = 0.0_8
            dvxd = dvxd + rot(2, 1)*du1d(ivy)
            dvyd = dvyd + rot(2, 2)*du1d(ivy)
            dvzd = dvzd + rot(2, 3)*du1d(ivy)
            du1d(ivy) = 0.0_8
            dvxd = dvxd + rot(1, 1)*du1d(ivx)
            dvyd = dvyd + rot(1, 2)*du1d(ivx)
            dvzd = dvzd + rot(1, 3)*du1d(ivx)
            du1d(ivx) = 0.0_8
            du1d(ivz) = du1d(ivz) + dvzd
            du1d(ivy) = du1d(ivy) + dvyd
            du1d(ivx) = du1d(ivx) + dvxd
        end if
    end subroutine leftrightstate_fast_b
    ! store the density flux in the mass flow of the
    ! appropriate sliding mesh interface.
    !      ==================================================================
    subroutine leftrightstate(du1, du2, du3, rotmatrix, left, right)
        implicit none
        !
        !        local parameter.
        !
        real(kind=realtype), parameter :: epslim=1.e-10_realtype
        !
        !        subroutine arguments.
        !
        real(kind=realtype), dimension(:), intent(inout) :: du1, du2, du3
        real(kind=realtype), dimension(:), intent(out) :: left, right
        real(kind=realtype), dimension(:, :, :, :, :), pointer :: &
        &     rotmatrix
        !
        !        local variables.
        !
        integer(kind=inttype) :: l
        real(kind=realtype) :: rl1, rl2, rr1, rr2, tmp, dvx, dvy, dvz
        real(kind=realtype), dimension(3, 3) :: rot
        intrinsic abs
        intrinsic max
        intrinsic sign
        intrinsic min
        real(kind=realtype) :: x6
        real(kind=realtype) :: x5
        real(kind=realtype) :: x4
        real(kind=realtype) :: x3
        real(kind=realtype) :: x2
        real(kind=realtype) :: x1
        real(kind=realtype) :: max7
        real(kind=realtype) :: max6
        real(kind=realtype) :: max5
        real(kind=realtype) :: max4
        real(kind=realtype) :: max3
        real(kind=realtype) :: max2
        real(kind=realtype) :: y4
        real(kind=realtype) :: y3
        real(kind=realtype) :: y2
        real(kind=realtype) :: y1
        ! check if the velocity components should be transformed to
        ! the cylindrical frame.
        if (rotationalperiodic) then
            ! store the rotation matrix a bit easier. note that the i,j,k
            ! come from the main subroutine.
            rot(1, 1) = rotmatrix(i, j, k, 1, 1)
            rot(1, 2) = rotmatrix(i, j, k, 1, 2)
            rot(1, 3) = rotmatrix(i, j, k, 1, 3)
            rot(2, 1) = rotmatrix(i, j, k, 2, 1)
            rot(2, 2) = rotmatrix(i, j, k, 2, 2)
            rot(2, 3) = rotmatrix(i, j, k, 2, 3)
            rot(3, 1) = rotmatrix(i, j, k, 3, 1)
            rot(3, 2) = rotmatrix(i, j, k, 3, 2)
            rot(3, 3) = rotmatrix(i, j, k, 3, 3)
            ! apply the transformation to the velocity components
            ! of du1, du2 and du3.
            dvx = du1(ivx)
            dvy = du1(ivy)
            dvz = du1(ivz)
            du1(ivx) = rot(1, 1)*dvx + rot(1, 2)*dvy + rot(1, 3)*dvz
            du1(ivy) = rot(2, 1)*dvx + rot(2, 2)*dvy + rot(2, 3)*dvz
            du1(ivz) = rot(3, 1)*dvx + rot(3, 2)*dvy + rot(3, 3)*dvz
            dvx = du2(ivx)
            dvy = du2(ivy)
            dvz = du2(ivz)
            du2(ivx) = rot(1, 1)*dvx + rot(1, 2)*dvy + rot(1, 3)*dvz
            du2(ivy) = rot(2, 1)*dvx + rot(2, 2)*dvy + rot(2, 3)*dvz
            du2(ivz) = rot(3, 1)*dvx + rot(3, 2)*dvy + rot(3, 3)*dvz
            dvx = du3(ivx)
            dvy = du3(ivy)
            dvz = du3(ivz)
            du3(ivx) = rot(1, 1)*dvx + rot(1, 2)*dvy + rot(1, 3)*dvz
            du3(ivy) = rot(2, 1)*dvx + rot(2, 2)*dvy + rot(2, 3)*dvz
            du3(ivz) = rot(3, 1)*dvx + rot(3, 2)*dvy + rot(3, 3)*dvz
        end if
        ! determine the limiter used.
        select case  (limused)
        case (nolimiter)
            ! linear interpolation; no limiter.
            ! loop over the number of variables to be interpolated.
            do l=1,nwint
                left(l) = omk*du1(l) + opk*du2(l)
                right(l) = -(omk*du3(l)) - opk*du2(l)
            end do
        case (vanalbeda)
            !          ==============================================================
            ! nonlinear interpolation using the van albeda limiter.
            ! loop over the number of variables to be interpolated.
            do l=1,nwint
                if (du2(l) .ge. 0.) then
                    x1 = du2(l)
                else
                    x1 = -du2(l)
                end if
                if (x1 .lt. epslim) then
                    max2 = epslim
                else
                    max2 = x1
                end if
                ! compute the limiter argument rl1, rl2, rr1 and rr2.
                ! note the cut off to 0.0.
                tmp = one/sign(max2, du2(l))
                if (du1(l) .ge. 0.) then
                    x3 = du1(l)
                else
                    x3 = -du1(l)
                end if
                if (x3 .lt. epslim) then
                    max4 = epslim
                else
                    max4 = x3
                end if
                y1 = du2(l)/sign(max4, du1(l))
                if (zero .lt. y1) then
                    rl1 = y1
                else
                    rl1 = zero
                end if
                if (zero .lt. du1(l)*tmp) then
                    rl2 = du1(l)*tmp
                else
                    rl2 = zero
                end if
                if (zero .lt. du3(l)*tmp) then
                    rr1 = du3(l)*tmp
                else
                    rr1 = zero
                end if
                if (du3(l) .ge. 0.) then
                    x4 = du3(l)
                else
                    x4 = -du3(l)
                end if
                if (x4 .lt. epslim) then
                    max5 = epslim
                else
                    max5 = x4
                end if
                y2 = du2(l)/sign(max5, du3(l))
                if (zero .lt. y2) then
                    rr2 = y2
                else
                    rr2 = zero
                end if
                ! compute the corresponding limiter values.
                rl1 = rl1*(rl1+one)/(rl1*rl1+one)
                rl2 = rl2*(rl2+one)/(rl2*rl2+one)
                rr1 = rr1*(rr1+one)/(rr1*rr1+one)
                rr2 = rr2*(rr2+one)/(rr2*rr2+one)
                ! compute the nonlinear corrections to the first order
                ! scheme.
                left(l) = omk*rl1*du1(l) + opk*rl2*du2(l)
                right(l) = -(opk*rr1*du2(l)) - omk*rr2*du3(l)
            end do
        case (minmod)
            !          ==============================================================
            ! nonlinear interpolation using the minmod limiter.
            ! loop over the number of variables to be interpolated.
            do l=1,nwint
                if (du2(l) .ge. 0.) then
                    x2 = du2(l)
                else
                    x2 = -du2(l)
                end if
                if (x2 .lt. epslim) then
                    max3 = epslim
                else
                    max3 = x2
                end if
                ! compute the limiter argument rl1, rl2, rr1 and rr2.
                ! note the cut off to 0.0.
                tmp = one/sign(max3, du2(l))
                if (du1(l) .ge. 0.) then
                    x5 = du1(l)
                else
                    x5 = -du1(l)
                end if
                if (x5 .lt. epslim) then
                    max6 = epslim
                else
                    max6 = x5
                end if
                y3 = du2(l)/sign(max6, du1(l))
                if (zero .lt. y3) then
                    rl1 = y3
                else
                    rl1 = zero
                end if
                if (zero .lt. du1(l)*tmp) then
                    rl2 = du1(l)*tmp
                else
                    rl2 = zero
                end if
                if (zero .lt. du3(l)*tmp) then
                    rr1 = du3(l)*tmp
                else
                    rr1 = zero
                end if
                if (du3(l) .ge. 0.) then
                    x6 = du3(l)
                else
                    x6 = -du3(l)
                end if
                if (x6 .lt. epslim) then
                    max7 = epslim
                else
                    max7 = x6
                end if
                y4 = du2(l)/sign(max7, du3(l))
                if (zero .lt. y4) then
                    rr2 = y4
                else
                    rr2 = zero
                end if
                if (one .gt. factminmod*rl1) then
                    rl1 = factminmod*rl1
                else
                    rl1 = one
                end if
                if (one .gt. factminmod*rl2) then
                    rl2 = factminmod*rl2
                else
                    rl2 = one
                end if
                if (one .gt. factminmod*rr1) then
                    rr1 = factminmod*rr1
                else
                    rr1 = one
                end if
                if (one .gt. factminmod*rr2) then
                    rr2 = factminmod*rr2
                else
                    rr2 = one
                end if
                ! compute the nonlinear corrections to the first order
                ! scheme.
                left(l) = omk*rl1*du1(l) + opk*rl2*du2(l)
                right(l) = -(opk*rr1*du2(l)) - omk*rr2*du3(l)
            end do
        end select
        ! in case only a first order scheme must be used for the
        ! turbulent transport equations, set the correction for the
        ! turbulent kinetic energy to 0.
        if (firstorderk) then
            left(itu1) = zero
            right(itu1) = zero
        end if
        ! for rotational periodic problems transform the velocity
        ! differences back to cartesian again. note that now the
        ! transpose of the rotation matrix must be used.
        if (rotationalperiodic) then
            ! left state.
            dvx = left(ivx)
            dvy = left(ivy)
            dvz = left(ivz)
            left(ivx) = rot(1, 1)*dvx + rot(2, 1)*dvy + rot(3, 1)*dvz
            left(ivy) = rot(1, 2)*dvx + rot(2, 2)*dvy + rot(3, 2)*dvz
            left(ivz) = rot(1, 3)*dvx + rot(2, 3)*dvy + rot(3, 3)*dvz
            ! right state.
            dvx = right(ivx)
            dvy = right(ivy)
            dvz = right(ivz)
            right(ivx) = rot(1, 1)*dvx + rot(2, 1)*dvy + rot(3, 1)*dvz
            right(ivy) = rot(1, 2)*dvx + rot(2, 2)*dvy + rot(3, 2)*dvz
            right(ivz) = rot(1, 3)*dvx + rot(2, 3)*dvy + rot(3, 3)*dvz
        end if
    end subroutine leftrightstate
    !  differentiation of riemannflux in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
    !   gradient     of useful results: flux left right
    !   with respect to varying inputs: flux left right
    !        ================================================================
    subroutine riemannflux_fast_b(left, leftd, right, rightd, flux, &
        &     fluxd)
        implicit none
        !
        !        subroutine arguments.
        !
        real(kind=realtype), dimension(*), intent(in) :: left, right
        real(kind=realtype), dimension(*) :: leftd, rightd
        real(kind=realtype), dimension(*) :: flux
        real(kind=realtype), dimension(*) :: fluxd
        !
        !        local variables.
        !
        real(kind=realtype) :: porflux, rface
        real(kind=realtype) :: etl, etr, z1l, z1r, tmp
        real(kind=realtype) :: etld, etrd, z1ld, z1rd, tmpd
        real(kind=realtype) :: dr, dru, drv, drw, dre, drk
        real(kind=realtype) :: drd, drud, drvd, drwd, dred, drkd
        real(kind=realtype) :: ravg, uavg, vavg, wavg, havg, kavg
        real(kind=realtype) :: uavgd, vavgd, wavgd, havgd, kavgd
        real(kind=realtype) :: alphaavg, a2avg, aavg, unavg
        real(kind=realtype) :: alphaavgd, a2avgd, aavgd, unavgd
        real(kind=realtype) :: ovaavg, ova2avg, area, eta
        real(kind=realtype) :: ovaavgd, ova2avgd, etad
        real(kind=realtype) :: gm1, gm53
        real(kind=realtype) :: lam1, lam2, lam3
        real(kind=realtype) :: lam1d, lam2d, lam3d
        real(kind=realtype) :: abv1, abv2, abv3, abv4, abv5, abv6, abv7
        real(kind=realtype) :: abv1d, abv2d, abv3d, abv4d, abv5d, abv6d, &
        &     abv7d
        real(kind=realtype), dimension(2) :: ktmp
        real(kind=realtype), dimension(2) :: ktmpd
        intrinsic sqrt
        intrinsic max
        intrinsic abs
        integer :: branch
        real(kind=realtype) :: tempd14
        real(kind=realtype) :: temp2
        real(kind=realtype) :: tempd13
        real(kind=realtype) :: temp1
        real(kind=realtype) :: tempd12
        real(kind=realtype) :: temp0
        real(kind=realtype) :: tempd11
        real(kind=realtype) :: tempd10
        real(kind=realtype) :: abs1d
        real(kind=realtype) :: x2
        real(kind=realtype) :: x2d
        real(kind=realtype) :: x1
        real(kind=realtype) :: tempd9
        real(kind=realtype) :: tempd
        real(kind=realtype) :: tempd8
        real(kind=realtype) :: tempd7
        real(kind=realtype) :: tempd6
        real(kind=realtype) :: tempd5
        real(kind=realtype) :: tempd4
        real(kind=realtype) :: tempd3
        real(kind=realtype) :: tempd2
        real(kind=realtype) :: tempd1
        real(kind=realtype) :: tempd0
        real(kind=realtype) :: x1d
        real(kind=realtype) :: abs2
        real(kind=realtype) :: abs2d
        real(kind=realtype) :: abs1
        real(kind=realtype) :: temp
        real(kind=realtype) :: max2
        real(kind=realtype) :: tempd18
        real(kind=realtype) :: tempd17
        real(kind=realtype) :: tempd16
        real(kind=realtype) :: tempd15
        ! set the porosity for the flux. the default value, 0.5*rfil, is
        ! a scaling factor where an rfil != 1 is taken into account.
        porflux = half*rfil
        if (por .eq. noflux .or. por .eq. boundflux) porflux = zero
        ! abbreviate some expressions in which gamma occurs.
        gm1 = gammaface - one
        gm53 = gammaface - five*third
        ! determine which riemann solver must be solved.
        select case  (riemannused)
        case (roe)
            ! determine the preconditioner used.
            select case  (precond)
            case (noprecond)
                ! no preconditioner used. use the roe scheme of the
                ! standard equations.
                ! compute the square root of the left and right densities
                ! and the inverse of the sum.
                z1l = sqrt(left(irho))
                z1r = sqrt(right(irho))
                tmp = one/(z1l+z1r)
                ! compute some variables depending whether or not a
                ! k-equation is present.
                if (correctfork) then
                    ! store the left and right kinetic energy in ktmp,
                    ! which is needed to compute the total energy.
                    ktmp(1) = left(itu1)
                    ktmp(2) = right(itu1)
                    ! store the difference of the turbulent kinetic energy
                    ! per unit volume, i.e. the conserved variable.
                    drk = right(irho)*right(itu1) - left(irho)*left(itu1)
                    ! compute the average turbulent energy per unit mass
                    ! using roe averages.
                    kavg = tmp*(z1l*left(itu1)+z1r*right(itu1))
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                else
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                    ! set the difference of the turbulent kinetic energy
                    ! per unit volume and the averaged kinetic energy per
                    ! unit mass to zero.
                    drk = 0.0
                    kavg = 0.0
                end if
                ! compute the total energy of the left and right state.
                call etot(left(irho), left(ivx), left(ivy), left(ivz), left(&
                &             irhoe), ktmp(1), etl, correctfork)
                call etot(right(irho), right(ivx), right(ivy), right(ivz), &
                &             right(irhoe), ktmp(2), etr, correctfork)
                ! compute the difference of the conservative mean
                ! flow variables.
                dr = right(irho) - left(irho)
                dru = right(irho)*right(ivx) - left(irho)*left(ivx)
                drv = right(irho)*right(ivy) - left(irho)*left(ivy)
                drw = right(irho)*right(ivz) - left(irho)*left(ivz)
                dre = etr - etl
                ! compute the roe average variables, which can be
                ! computed directly from the average roe vector.
                uavg = tmp*(z1l*left(ivx)+z1r*right(ivx))
                vavg = tmp*(z1l*left(ivy)+z1r*right(ivy))
                wavg = tmp*(z1l*left(ivz)+z1r*right(ivz))
                havg = tmp*((etl+left(irhoe))/z1l+(etr+right(irhoe))/z1r)
                ! compute the unit vector and store the area of the
                ! normal. also compute the unit normal velocity of the face.
                area = sqrt(sx**2 + sy**2 + sz**2)
                if (1.e-25_realtype .lt. area) then
                    max2 = area
                else
                    max2 = 1.e-25_realtype
                end if
                tmp = one/max2
                sx = sx*tmp
                sy = sy*tmp
                sz = sz*tmp
                rface = sface*tmp
                ! compute some dependent variables at the roe
                ! average state.
                alphaavg = half*(uavg**2+vavg**2+wavg**2)
                if (gm1*(havg-alphaavg) - gm53*kavg .ge. 0.) then
                    a2avg = gm1*(havg-alphaavg) - gm53*kavg
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    a2avg = -(gm1*(havg-alphaavg)-gm53*kavg)
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                aavg = sqrt(a2avg)
                unavg = uavg*sx + vavg*sy + wavg*sz
                ovaavg = one/aavg
                ova2avg = one/a2avg
                ! set for a boundary the normal velocity to rface, the
                ! normal velocity of the boundary.
                if (por .eq. boundflux) then
                    unavg = rface
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                else
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                end if
                x1 = (left(ivx)-right(ivx))*sx + (left(ivy)-right(ivy))*sy + (&
                &           left(ivz)-right(ivz))*sz
                if (x1 .ge. 0.) then
                    abs1 = x1
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                else
                    abs1 = -x1
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                end if
                x2 = sqrt(gammaface*left(irhoe)/left(irho)) - sqrt(gammaface*&
                &           right(irhoe)/right(irho))
                if (x2 .ge. 0.) then
                    abs2 = x2
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    abs2 = -x2
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                ! compute the coefficient eta for the entropy correction.
                ! at the moment a 1d entropy correction is used, which
                ! removes expansion shocks. although it also reduces the
                ! carbuncle phenomenon, it does not remove it completely.
                ! in other to do that a multi-dimensional entropy fix is
                ! needed, see sanders et. al, jcp, vol. 145, 1998,
                ! pp. 511 - 537. although relatively easy to implement,
                ! an efficient implementation requires the storage of
                ! all the left and right states, which is rather
                ! expensive in terms of memory.
                eta = half*(abs1+abs2)
                if (unavg - rface + aavg .ge. 0.) then
                    lam1 = unavg - rface + aavg
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    lam1 = -(unavg-rface+aavg)
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (unavg - rface - aavg .ge. 0.) then
                    lam2 = unavg - rface - aavg
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    lam2 = -(unavg-rface-aavg)
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (unavg - rface .ge. 0.) then
                    lam3 = unavg - rface
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    lam3 = -(unavg-rface)
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                ! apply the entropy correction to the eigenvalues.
                tmp = two*eta
                if (lam1 .lt. tmp) then
                    lam1 = eta + fourth*lam1*lam1/eta
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (lam2 .lt. tmp) then
                    lam2 = eta + fourth*lam2*lam2/eta
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                if (lam3 .lt. tmp) then
                    lam3 = eta + fourth*lam3*lam3/eta
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                ! multiply the eigenvalues by the area to obtain
                ! the correct values for the dissipation term.
                lam1 = lam1*area
                lam2 = lam2*area
                lam3 = lam3*area
                ! some abbreviations, which occur quite often in the
                ! dissipation terms.
                abv1 = half*(lam1+lam2)
                abv2 = half*(lam1-lam2)
                abv3 = abv1 - lam3
                abv4 = gm1*(alphaavg*dr-uavg*dru-vavg*drv-wavg*drw+dre) - gm53&
                &           *drk
                abv5 = sx*dru + sy*drv + sz*drw - unavg*dr
                abv6 = abv3*abv4*ova2avg + abv2*abv5*ovaavg
                abv7 = abv2*abv4*ovaavg + abv3*abv5
                ! compute the dissipation term, -|a| (wr - wl), which is
                ! multiplied by porflux. note that porflux is either
                ! 0.0 or 0.5*rfil.
                tempd13 = -(porflux*fluxd(irhoe))
                havgd = abv6*tempd13
                fluxd(irhoe) = 0.0_8
                tempd14 = -(porflux*fluxd(imz))
                fluxd(imz) = 0.0_8
                tempd17 = -(porflux*fluxd(imy))
                fluxd(imy) = 0.0_8
                tempd15 = -(porflux*fluxd(imx))
                abv7d = sz*tempd14 + sx*tempd15 + sy*tempd17 + unavg*tempd13
                fluxd(imx) = 0.0_8
                tempd16 = -(porflux*fluxd(irho))
                abv6d = wavg*tempd14 + uavg*tempd15 + tempd16 + vavg*tempd17 +&
                &           havg*tempd13
                fluxd(irho) = 0.0_8
                abv2d = ovaavg*abv5*abv6d + ovaavg*abv4*abv7d
                abv4d = ova2avg*abv3*abv6d + ovaavg*abv2*abv7d
                ovaavgd = abv2*abv5*abv6d + abv2*abv4*abv7d
                abv3d = ova2avg*abv4*abv6d + abv5*abv7d
                lam3d = drw*tempd14 + dru*tempd15 - abv3d + dr*tempd16 + drv*&
                &           tempd17 + dre*tempd13
                abv5d = ovaavg*abv2*abv6d + abv3*abv7d
                unavgd = abv7*tempd13 - dr*abv5d
                ova2avgd = abv3*abv4*abv6d
                tempd18 = gm1*abv4d
                dred = tempd18 + lam3*tempd13
                drwd = sz*abv5d - wavg*tempd18 + lam3*tempd14
                wavgd = abv6*tempd14 - drw*tempd18
                drvd = sy*abv5d - vavg*tempd18 + lam3*tempd17
                vavgd = abv6*tempd17 - drv*tempd18
                drud = sx*abv5d - uavg*tempd18 + lam3*tempd15
                uavgd = abv6*tempd15 - dru*tempd18
                drd = alphaavg*tempd18 - unavg*abv5d + lam3*tempd16
                alphaavgd = dr*tempd18
                drkd = -(gm53*abv4d)
                abv1d = abv3d
                lam1d = half*abv1d + half*abv2d
                lam2d = half*abv1d - half*abv2d
                lam3d = area*lam3d
                lam2d = area*lam2d
                lam1d = area*lam1d
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    tempd12 = fourth*lam3d/eta
                    etad = lam3d - lam3**2*tempd12/eta
                    lam3d = 2*lam3*tempd12
                else
                    etad = 0.0_8
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    tempd11 = fourth*lam2d/eta
                    etad = etad + lam2d - lam2**2*tempd11/eta
                    lam2d = 2*lam2*tempd11
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    tempd10 = fourth*lam1d/eta
                    etad = etad + lam1d - lam1**2*tempd10/eta
                    lam1d = 2*lam1*tempd10
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    unavgd = unavgd + lam3d
                else
                    unavgd = unavgd - lam3d
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    unavgd = unavgd + lam2d
                    aavgd = -lam2d
                else
                    aavgd = lam2d
                    unavgd = unavgd - lam2d
                end if
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    unavgd = unavgd + lam1d
                    aavgd = aavgd + lam1d
                else
                    unavgd = unavgd - lam1d
                    aavgd = aavgd - lam1d
                end if
                abs1d = half*etad
                abs2d = half*etad
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    x2d = abs2d
                else
                    x2d = -abs2d
                end if
                temp1 = left(irhoe)/left(irho)
                if (gammaface*temp1 .eq. 0.0_8) then
                    tempd8 = 0.0
                else
                    tempd8 = gammaface*x2d/(2.0*sqrt(gammaface*temp1)*left(irho)&
                    &             )
                end if
                temp2 = right(irhoe)/right(irho)
                if (gammaface*temp2 .eq. 0.0_8) then
                    tempd9 = 0.0
                else
                    tempd9 = -(gammaface*x2d/(2.0*sqrt(gammaface*temp2)*right(&
                    &             irho)))
                end if
                leftd(irhoe) = leftd(irhoe) + tempd8
                leftd(irho) = leftd(irho) - temp1*tempd8
                rightd(irhoe) = rightd(irhoe) + tempd9
                rightd(irho) = rightd(irho) - temp2*tempd9
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    x1d = -abs1d
                else
                    x1d = abs1d
                end if
                leftd(ivx) = leftd(ivx) + sx*x1d
                rightd(ivx) = rightd(ivx) - sx*x1d
                leftd(ivy) = leftd(ivy) + sy*x1d
                rightd(ivy) = rightd(ivy) - sy*x1d
                leftd(ivz) = leftd(ivz) + sz*x1d
                rightd(ivz) = rightd(ivz) - sz*x1d
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .ne. 0) unavgd = 0.0_8
                aavgd = aavgd - one*ovaavgd/aavg**2
                if (a2avg .eq. 0.0_8) then
                    a2avgd = -(one*ova2avgd/a2avg**2)
                else
                    a2avgd = aavgd/(2.0*sqrt(a2avg)) - one*ova2avgd/a2avg**2
                end if
                uavgd = uavgd + sx*unavgd
                vavgd = vavgd + sy*unavgd
                wavgd = wavgd + sz*unavgd
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .eq. 0) then
                    havgd = havgd + gm1*a2avgd
                    alphaavgd = alphaavgd - gm1*a2avgd
                    kavgd = -(gm53*a2avgd)
                else
                    kavgd = gm53*a2avgd
                    havgd = havgd - gm1*a2avgd
                    alphaavgd = alphaavgd + gm1*a2avgd
                end if
                tempd7 = half*alphaavgd
                uavgd = uavgd + 2*uavg*tempd7
                vavgd = vavgd + 2*vavg*tempd7
                wavgd = wavgd + 2*wavg*tempd7
                tmp = one/(z1l+z1r)
                tempd5 = tmp*uavgd
                tempd6 = tmp*vavgd
                tempd4 = tmp*wavgd
                temp0 = (etr+right(irhoe))/z1r
                temp = (etl+left(irhoe))/z1l
                tempd1 = tmp*havgd
                tempd2 = tempd1/z1l
                tempd3 = tempd1/z1r
                tmpd = (z1l*left(ivz)+z1r*right(ivz))*wavgd + (z1l*left(ivx)+&
                &           z1r*right(ivx))*uavgd + (z1l*left(ivy)+z1r*right(ivy))*vavgd&
                &           + (temp+temp0)*havgd
                etld = tempd2 - dred
                leftd(irhoe) = leftd(irhoe) + tempd2
                z1ld = left(ivz)*tempd4 + left(ivx)*tempd5 + left(ivy)*tempd6 &
                &           - temp*tempd2
                etrd = dred + tempd3
                rightd(irhoe) = rightd(irhoe) + tempd3
                z1rd = right(ivz)*tempd4 + right(ivx)*tempd5 + right(ivy)*&
                &           tempd6 - temp0*tempd3
                leftd(ivz) = leftd(ivz) + z1l*tempd4
                rightd(ivz) = rightd(ivz) + z1r*tempd4
                leftd(ivy) = leftd(ivy) + z1l*tempd6
                rightd(ivy) = rightd(ivy) + z1r*tempd6
                leftd(ivx) = leftd(ivx) + z1l*tempd5
                rightd(ivx) = rightd(ivx) + z1r*tempd5
                rightd(irho) = rightd(irho) + right(ivz)*drwd
                rightd(ivz) = rightd(ivz) + right(irho)*drwd
                leftd(irho) = leftd(irho) - left(ivz)*drwd
                leftd(ivz) = leftd(ivz) - left(irho)*drwd
                rightd(irho) = rightd(irho) + right(ivy)*drvd
                rightd(ivy) = rightd(ivy) + right(irho)*drvd
                leftd(irho) = leftd(irho) - left(ivy)*drvd
                leftd(ivy) = leftd(ivy) - left(irho)*drvd
                rightd(irho) = rightd(irho) + right(ivx)*drud
                rightd(ivx) = rightd(ivx) + right(irho)*drud
                leftd(irho) = leftd(irho) - left(ivx)*drud
                leftd(ivx) = leftd(ivx) - left(irho)*drud
                rightd(irho) = rightd(irho) + drd
                leftd(irho) = leftd(irho) - drd
                ktmpd = 0.0_8
                call etot_fast_b(right(irho), rightd(irho), right(ivx), rightd&
                &                    (ivx), right(ivy), rightd(ivy), right(ivz), rightd(&
                &                    ivz), right(irhoe), rightd(irhoe), ktmp(2), ktmpd(2&
                &                    ), etr, etrd, correctfork)
                call etot_fast_b(left(irho), leftd(irho), left(ivx), leftd(ivx&
                &                    ), left(ivy), leftd(ivy), left(ivz), leftd(ivz), &
                &                    left(irhoe), leftd(irhoe), ktmp(1), ktmpd(1), etl, &
                &                    etld, correctfork)
                branch = myIntStack(myIntPtr)
                myIntPtr = myIntPtr - 1
                if (branch .ne. 0) then
                    tempd0 = tmp*kavgd
                    tmpd = tmpd + (z1l*left(itu1)+z1r*right(itu1))*kavgd
                    z1ld = z1ld + left(itu1)*tempd0
                    leftd(itu1) = leftd(itu1) + z1l*tempd0
                    z1rd = z1rd + right(itu1)*tempd0
                    rightd(itu1) = rightd(itu1) + z1r*tempd0
                    rightd(irho) = rightd(irho) + right(itu1)*drkd
                    rightd(itu1) = rightd(itu1) + ktmpd(2) + right(irho)*drkd
                    leftd(irho) = leftd(irho) - left(itu1)*drkd
                    ktmpd(2) = 0.0_8
                    leftd(itu1) = leftd(itu1) + ktmpd(1) - left(irho)*drkd
                end if
                tempd = -(one*tmpd/(z1l+z1r)**2)
                z1ld = z1ld + tempd
                z1rd = z1rd + tempd
                if (.not.right(irho) .eq. 0.0_8) rightd(irho) = rightd(irho) +&
                &             z1rd/(2.0*sqrt(right(irho)))
                if (.not.left(irho) .eq. 0.0_8) leftd(irho) = leftd(irho) + &
                &             z1ld/(2.0*sqrt(left(irho)))
            end select
        end select
    end subroutine riemannflux_fast_b
    !        ================================================================
    subroutine riemannflux(left, right, flux)
        implicit none
        !
        !        subroutine arguments.
        !
        real(kind=realtype), dimension(*), intent(in) :: left, right
        real(kind=realtype), dimension(*), intent(out) :: flux
        !
        !        local variables.
        !
        real(kind=realtype) :: porflux, rface
        real(kind=realtype) :: etl, etr, z1l, z1r, tmp
        real(kind=realtype) :: dr, dru, drv, drw, dre, drk
        real(kind=realtype) :: ravg, uavg, vavg, wavg, havg, kavg
        real(kind=realtype) :: alphaavg, a2avg, aavg, unavg
        real(kind=realtype) :: ovaavg, ova2avg, area, eta
        real(kind=realtype) :: gm1, gm53
        real(kind=realtype) :: lam1, lam2, lam3
        real(kind=realtype) :: abv1, abv2, abv3, abv4, abv5, abv6, abv7
        real(kind=realtype), dimension(2) :: ktmp
        intrinsic sqrt
        intrinsic max
        intrinsic abs
        real(kind=realtype) :: x2
        real(kind=realtype) :: x1
        real(kind=realtype) :: abs2
        real(kind=realtype) :: abs1
        real(kind=realtype) :: max2
        ! set the porosity for the flux. the default value, 0.5*rfil, is
        ! a scaling factor where an rfil != 1 is taken into account.
        porflux = half*rfil
        if (por .eq. noflux .or. por .eq. boundflux) porflux = zero
        ! abbreviate some expressions in which gamma occurs.
        gm1 = gammaface - one
        gm53 = gammaface - five*third
        ! determine which riemann solver must be solved.
        select case  (riemannused)
        case (roe)
            ! determine the preconditioner used.
            select case  (precond)
            case (noprecond)
                ! no preconditioner used. use the roe scheme of the
                ! standard equations.
                ! compute the square root of the left and right densities
                ! and the inverse of the sum.
                z1l = sqrt(left(irho))
                z1r = sqrt(right(irho))
                tmp = one/(z1l+z1r)
                ! compute some variables depending whether or not a
                ! k-equation is present.
                if (correctfork) then
                    ! store the left and right kinetic energy in ktmp,
                    ! which is needed to compute the total energy.
                    ktmp(1) = left(itu1)
                    ktmp(2) = right(itu1)
                    ! store the difference of the turbulent kinetic energy
                    ! per unit volume, i.e. the conserved variable.
                    drk = right(irho)*right(itu1) - left(irho)*left(itu1)
                    ! compute the average turbulent energy per unit mass
                    ! using roe averages.
                    kavg = tmp*(z1l*left(itu1)+z1r*right(itu1))
                else
                    ! set the difference of the turbulent kinetic energy
                    ! per unit volume and the averaged kinetic energy per
                    ! unit mass to zero.
                    drk = 0.0
                    kavg = 0.0
                end if
                ! compute the total energy of the left and right state.
                call etot(left(irho), left(ivx), left(ivy), left(ivz), left(&
                &             irhoe), ktmp(1), etl, correctfork)
                call etot(right(irho), right(ivx), right(ivy), right(ivz), &
                &             right(irhoe), ktmp(2), etr, correctfork)
                ! compute the difference of the conservative mean
                ! flow variables.
                dr = right(irho) - left(irho)
                dru = right(irho)*right(ivx) - left(irho)*left(ivx)
                drv = right(irho)*right(ivy) - left(irho)*left(ivy)
                drw = right(irho)*right(ivz) - left(irho)*left(ivz)
                dre = etr - etl
                ! compute the roe average variables, which can be
                ! computed directly from the average roe vector.
                ravg = fourth*(z1r+z1l)**2
                uavg = tmp*(z1l*left(ivx)+z1r*right(ivx))
                vavg = tmp*(z1l*left(ivy)+z1r*right(ivy))
                wavg = tmp*(z1l*left(ivz)+z1r*right(ivz))
                havg = tmp*((etl+left(irhoe))/z1l+(etr+right(irhoe))/z1r)
                ! compute the unit vector and store the area of the
                ! normal. also compute the unit normal velocity of the face.
                area = sqrt(sx**2 + sy**2 + sz**2)
                if (1.e-25_realtype .lt. area) then
                    max2 = area
                else
                    max2 = 1.e-25_realtype
                end if
                tmp = one/max2
                sx = sx*tmp
                sy = sy*tmp
                sz = sz*tmp
                rface = sface*tmp
                ! compute some dependent variables at the roe
                ! average state.
                alphaavg = half*(uavg**2+vavg**2+wavg**2)
                if (gm1*(havg-alphaavg) - gm53*kavg .ge. 0.) then
                    a2avg = gm1*(havg-alphaavg) - gm53*kavg
                else
                    a2avg = -(gm1*(havg-alphaavg)-gm53*kavg)
                end if
                aavg = sqrt(a2avg)
                unavg = uavg*sx + vavg*sy + wavg*sz
                ovaavg = one/aavg
                ova2avg = one/a2avg
                ! set for a boundary the normal velocity to rface, the
                ! normal velocity of the boundary.
                if (por .eq. boundflux) unavg = rface
                x1 = (left(ivx)-right(ivx))*sx + (left(ivy)-right(ivy))*sy + (&
                &           left(ivz)-right(ivz))*sz
                if (x1 .ge. 0.) then
                    abs1 = x1
                else
                    abs1 = -x1
                end if
                x2 = sqrt(gammaface*left(irhoe)/left(irho)) - sqrt(gammaface*&
                &           right(irhoe)/right(irho))
                if (x2 .ge. 0.) then
                    abs2 = x2
                else
                    abs2 = -x2
                end if
                ! compute the coefficient eta for the entropy correction.
                ! at the moment a 1d entropy correction is used, which
                ! removes expansion shocks. although it also reduces the
                ! carbuncle phenomenon, it does not remove it completely.
                ! in other to do that a multi-dimensional entropy fix is
                ! needed, see sanders et. al, jcp, vol. 145, 1998,
                ! pp. 511 - 537. although relatively easy to implement,
                ! an efficient implementation requires the storage of
                ! all the left and right states, which is rather
                ! expensive in terms of memory.
                eta = half*(abs1+abs2)
                if (unavg - rface + aavg .ge. 0.) then
                    lam1 = unavg - rface + aavg
                else
                    lam1 = -(unavg-rface+aavg)
                end if
                if (unavg - rface - aavg .ge. 0.) then
                    lam2 = unavg - rface - aavg
                else
                    lam2 = -(unavg-rface-aavg)
                end if
                if (unavg - rface .ge. 0.) then
                    lam3 = unavg - rface
                else
                    lam3 = -(unavg-rface)
                end if
                ! apply the entropy correction to the eigenvalues.
                tmp = two*eta
                if (lam1 .lt. tmp) lam1 = eta + fourth*lam1*lam1/eta
                if (lam2 .lt. tmp) lam2 = eta + fourth*lam2*lam2/eta
                if (lam3 .lt. tmp) lam3 = eta + fourth*lam3*lam3/eta
                ! multiply the eigenvalues by the area to obtain
                ! the correct values for the dissipation term.
                lam1 = lam1*area
                lam2 = lam2*area
                lam3 = lam3*area
                ! some abbreviations, which occur quite often in the
                ! dissipation terms.
                abv1 = half*(lam1+lam2)
                abv2 = half*(lam1-lam2)
                abv3 = abv1 - lam3
                abv4 = gm1*(alphaavg*dr-uavg*dru-vavg*drv-wavg*drw+dre) - gm53&
                &           *drk
                abv5 = sx*dru + sy*drv + sz*drw - unavg*dr
                abv6 = abv3*abv4*ova2avg + abv2*abv5*ovaavg
                abv7 = abv2*abv4*ovaavg + abv3*abv5
                ! compute the dissipation term, -|a| (wr - wl), which is
                ! multiplied by porflux. note that porflux is either
                ! 0.0 or 0.5*rfil.
                flux(irho) = -(porflux*(lam3*dr+abv6))
                flux(imx) = -(porflux*(lam3*dru+uavg*abv6+sx*abv7))
                flux(imy) = -(porflux*(lam3*drv+vavg*abv6+sy*abv7))
                flux(imz) = -(porflux*(lam3*drw+wavg*abv6+sz*abv7))
                flux(irhoe) = -(porflux*(lam3*dre+havg*abv6+unavg*abv7))
            case (turkel)
                !          tmp = max(lam1,lam2,lam3)
                !          flux(irho)  = -porflux*(tmp*dr)
                !          flux(imx)   = -porflux*(tmp*dru)
                !          flux(imy)   = -porflux*(tmp*drv)
                !          flux(imz)   = -porflux*(tmp*drw)
                !          flux(irhoe) = -porflux*(tmp*dre)
                call terminate('riemannflux', &
                &                  'turkel preconditioner not implemented yet')
            case (choimerkle)
                call terminate('riemannflux', &
                &                  'choi merkle preconditioner not implemented yet')
            end select
        case (vanleer)
            call terminate('riemannflux', 'van leer fvs not implemented yet'&
            &               )
        case (ausmdv)
            call terminate('riemannflux', 'ausmdv fvs not implemented yet')
        end select
    end subroutine riemannflux
end subroutine inviscidupwindflux_fast_b

!  differentiation of inviscidcentralflux in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: *p *dw *w
!   with respect to varying inputs: *p *dw *w
!   rw status of diff variables: *p:incr *dw:in-out *w:incr
!   plus diff mem management of: p:in dw:in w:in
subroutine inviscidcentralflux_fast_b()
    !
    !       inviscidcentralflux computes the euler fluxes using a central
    !       discretization for a given block. therefore it is assumed that
    !       the pointers in block pointer already point to the correct
    !       block on the correct multigrid level.
    !
    use constants
!     use blockpointers, only : nx, il, ie, ny, jl, je, nz, kl, ke, &
! &   spectralsol, w, wd, si, sj, sk, dw, dwd, pori, porj, pork, &
! &   indfamilyi, indfamilyj, indfamilyk, p, pd, sfacei, sfacej, sfacek, &
! &   nbkglobal, addgridvelocities, blockismoving, vol, factfamilyi, &
! &   factfamilyj, factfamilyk

    use blockpointers, only : nbkglobal, addgridvelocities, blockismoving

    use cgnsgrid, only : cgnsdoms, massflowfamilyinv
    use flowvarrefstate, only : timeref
    use inputphysics, only : equationmode
    implicit none
    !
    !      local variables.
    !
    integer(kind=inttype) :: i, j, k, ind, ii
    real(kind=realtype) :: qsp, qsm, rqsp, rqsm, porvel, porflux
    real(kind=realtype) :: qspd, qsmd, rqspd, rqsmd
    real(kind=realtype) :: pa, fs, sface, vnp, vnm
    real(kind=realtype) :: pad, fsd, vnpd, vnmd
    real(kind=realtype) :: wwx, wwy, wwz, rvol
    real(kind=realtype) :: rvold
    intrinsic mod
    integer :: branch
    real(kind=realtype) :: tempd
    real(kind=realtype) :: tempd1
    real(kind=realtype) :: tempd0
    if (blockismoving .and. equationmode .eq. steady) then
        ! compute the three nondimensional angular velocities.
        wwx = timeref*cgnsdoms(nbkglobal)%rotrate(1)
        wwy = timeref*cgnsdoms(nbkglobal)%rotrate(2)
        wwz = timeref*cgnsdoms(nbkglobal)%rotrate(3)
        do ii=0,nx*ny*nz-1
            i = mod(ii, nx) + 2
            j = mod(ii/nx, ny) + 2
            k = ii/(nx*ny) + 2
            rvol = w(i, j, k, irho)*vol(i, j, k)
            rvold = (wwx*w(i, j, k, ivy)-wwy*w(i, j, k, ivx))*dwd(i, j, k, &
            &         imz)
            wd(i, j, k, ivy) = wd(i, j, k, ivy) + rvol*wwx*dwd(i, j, k, imz)
            wd(i, j, k, ivx) = wd(i, j, k, ivx) - rvol*wwy*dwd(i, j, k, imz)
            rvold = rvold + (wwz*w(i, j, k, ivx)-wwx*w(i, j, k, ivz))*dwd(i&
            &         , j, k, imy)
            wd(i, j, k, ivx) = wd(i, j, k, ivx) + rvol*wwz*dwd(i, j, k, imy)
            wd(i, j, k, ivz) = wd(i, j, k, ivz) - rvol*wwx*dwd(i, j, k, imy)
            rvold = rvold + (wwy*w(i, j, k, ivz)-wwz*w(i, j, k, ivy))*dwd(i&
            &         , j, k, imx)
            wd(i, j, k, ivz) = wd(i, j, k, ivz) + rvol*wwy*dwd(i, j, k, imx)
            wd(i, j, k, ivy) = wd(i, j, k, ivy) - rvol*wwz*dwd(i, j, k, imx)
            wd(i, j, k, irho) = wd(i, j, k, irho) + vol(i, j, k)*rvold
        end do
    end if
    sface = zero
    do ii=0,nx*ny*kl-1
        i = mod(ii, nx) + 2
        j = mod(ii/nx, ny) + 2
        k = ii/(nx*ny) + 1
        ! set the dot product of the grid velocity and the
        ! normal in k-direction for a moving face.
        if (addgridvelocities) sface = sfacek(i, j, k)
        ! compute the normal velocities of the left and right state.
        vnp = w(i, j, k+1, ivx)*sk(i, j, k, 1) + w(i, j, k+1, ivy)*sk(i, j&
        &       , k, 2) + w(i, j, k+1, ivz)*sk(i, j, k, 3)
        vnm = w(i, j, k, ivx)*sk(i, j, k, 1) + w(i, j, k, ivy)*sk(i, j, k&
        &       , 2) + w(i, j, k, ivz)*sk(i, j, k, 3)
        ! set the values of the porosities for this face.
        ! porvel defines the porosity w.r.t. velocity;
        ! porflux defines the porosity w.r.t. the entire flux.
        ! the latter is only zero for a discontinuous block
        ! block boundary that must be treated conservatively.
        ! the default value of porflux is 0.5, such that the
        ! correct central flux is scattered to both cells.
        ! in case of a boundflux the normal velocity is set
        ! to sface.
        porvel = one
        porflux = half
        if (pork(i, j, k) .eq. noflux) porflux = zero
        if (pork(i, j, k) .eq. boundflux) then
            porvel = zero
            vnp = sface
            vnm = sface
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 0
        else
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 1
        end if
        ! incorporate porflux in porvel.
        porvel = porvel*porflux
        ! compute the normal velocities for the face as well as the
        ! mass fluxes.
        qsp = (vnp-sface)*porvel
        qsm = (vnm-sface)*porvel
        rqsp = qsp*w(i, j, k+1, irho)
        rqsm = qsm*w(i, j, k, irho)
        ! compute the sum of the pressure multiplied by porflux.
        ! for the default value of porflux, 0.5, this leads to
        ! the average pressure.
        ! compute the fluxes and scatter them to the cells
        ! i,j,k and i,j,k+1. store the density flux in the
        ! mass flow of the appropriate sliding mesh interface.
        fsd = dwd(i, j, k, irhoe) - dwd(i, j, k+1, irhoe)
        tempd1 = porflux*fsd
        qspd = w(i, j, k+1, irhoe)*fsd
        wd(i, j, k+1, irhoe) = wd(i, j, k+1, irhoe) + qsp*fsd
        qsmd = w(i, j, k, irhoe)*fsd
        wd(i, j, k, irhoe) = wd(i, j, k, irhoe) + qsm*fsd
        pd(i, j, k+1) = pd(i, j, k+1) + vnp*tempd1
        pd(i, j, k) = pd(i, j, k) + vnm*tempd1
        fsd = dwd(i, j, k, imz) - dwd(i, j, k+1, imz)
        rqspd = w(i, j, k+1, ivz)*fsd
        wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) + rqsp*fsd
        rqsmd = w(i, j, k, ivz)*fsd
        wd(i, j, k, ivz) = wd(i, j, k, ivz) + rqsm*fsd
        pad = sk(i, j, k, 3)*fsd
        fsd = dwd(i, j, k, imy) - dwd(i, j, k+1, imy)
        rqspd = rqspd + w(i, j, k+1, ivy)*fsd
        wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) + rqsp*fsd
        rqsmd = rqsmd + w(i, j, k, ivy)*fsd
        wd(i, j, k, ivy) = wd(i, j, k, ivy) + rqsm*fsd
        pad = pad + sk(i, j, k, 2)*fsd
        fsd = dwd(i, j, k, imx) - dwd(i, j, k+1, imx)
        rqspd = rqspd + w(i, j, k+1, ivx)*fsd
        wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) + rqsp*fsd
        rqsmd = rqsmd + w(i, j, k, ivx)*fsd
        wd(i, j, k, ivx) = wd(i, j, k, ivx) + rqsm*fsd
        pad = pad + sk(i, j, k, 1)*fsd
        fsd = dwd(i, j, k, irho) - dwd(i, j, k+1, irho)
        rqspd = rqspd + fsd
        rqsmd = rqsmd + fsd
        pd(i, j, k+1) = pd(i, j, k+1) + porflux*pad
        pd(i, j, k) = pd(i, j, k) + porflux*pad
        qsmd = qsmd + w(i, j, k, irho)*rqsmd
        vnmd = porvel*qsmd + p(i, j, k)*tempd1
        wd(i, j, k, irho) = wd(i, j, k, irho) + qsm*rqsmd
        qspd = qspd + w(i, j, k+1, irho)*rqspd
        vnpd = porvel*qspd + p(i, j, k+1)*tempd1
        wd(i, j, k+1, irho) = wd(i, j, k+1, irho) + qsp*rqspd
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) then
            vnmd = 0.0_8
            vnpd = 0.0_8
        end if
        wd(i, j, k, ivx) = wd(i, j, k, ivx) + sk(i, j, k, 1)*vnmd
        wd(i, j, k, ivy) = wd(i, j, k, ivy) + sk(i, j, k, 2)*vnmd
        wd(i, j, k, ivz) = wd(i, j, k, ivz) + sk(i, j, k, 3)*vnmd
        wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) + sk(i, j, k, 1)*vnpd
        wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) + sk(i, j, k, 2)*vnpd
        wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) + sk(i, j, k, 3)*vnpd
    end do
    sface = zero
    do ii=0,nx*jl*nz-1
        i = mod(ii, nx) + 2
        j = mod(ii/nx, jl) + 1
        k = ii/(nx*jl) + 2
        ! set the dot product of the grid velocity and the
        ! normal in j-direction for a moving face.
        if (addgridvelocities) sface = sfacej(i, j, k)
        ! compute the normal velocities of the left and right state.
        vnp = w(i, j+1, k, ivx)*sj(i, j, k, 1) + w(i, j+1, k, ivy)*sj(i, j&
        &       , k, 2) + w(i, j+1, k, ivz)*sj(i, j, k, 3)
        vnm = w(i, j, k, ivx)*sj(i, j, k, 1) + w(i, j, k, ivy)*sj(i, j, k&
        &       , 2) + w(i, j, k, ivz)*sj(i, j, k, 3)
        ! set the values of the porosities for this face.
        ! porvel defines the porosity w.r.t. velocity;
        ! porflux defines the porosity w.r.t. the entire flux.
        ! the latter is only zero for a discontinuous block
        ! boundary that must be treated conservatively.
        ! the default value of porflux is 0.5, such that the
        ! correct central flux is scattered to both cells.
        ! in case of a boundflux the normal velocity is set
        ! to sface.
        porvel = one
        porflux = half
        if (porj(i, j, k) .eq. noflux) porflux = zero
        if (porj(i, j, k) .eq. boundflux) then
            porvel = zero
            vnp = sface
            vnm = sface
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 0
        else
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 1
        end if
        ! incorporate porflux in porvel.
        porvel = porvel*porflux
        ! compute the normal velocities for the face as well as the
        ! mass fluxes.
        qsp = (vnp-sface)*porvel
        qsm = (vnm-sface)*porvel
        rqsp = qsp*w(i, j+1, k, irho)
        rqsm = qsm*w(i, j, k, irho)
        ! compute the sum of the pressure multiplied by porflux.
        ! for the default value of porflux, 0.5, this leads to
        ! the average pressure.
        ! compute the fluxes and scatter them to the cells
        ! i,j,k and i,j+1,k. store the density flux in the
        ! mass flow of the appropriate sliding mesh interface.
        fsd = dwd(i, j, k, irhoe) - dwd(i, j+1, k, irhoe)
        tempd0 = porflux*fsd
        qspd = w(i, j+1, k, irhoe)*fsd
        wd(i, j+1, k, irhoe) = wd(i, j+1, k, irhoe) + qsp*fsd
        qsmd = w(i, j, k, irhoe)*fsd
        wd(i, j, k, irhoe) = wd(i, j, k, irhoe) + qsm*fsd
        pd(i, j+1, k) = pd(i, j+1, k) + vnp*tempd0
        pd(i, j, k) = pd(i, j, k) + vnm*tempd0
        fsd = dwd(i, j, k, imz) - dwd(i, j+1, k, imz)
        rqspd = w(i, j+1, k, ivz)*fsd
        wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) + rqsp*fsd
        rqsmd = w(i, j, k, ivz)*fsd
        wd(i, j, k, ivz) = wd(i, j, k, ivz) + rqsm*fsd
        pad = sj(i, j, k, 3)*fsd
        fsd = dwd(i, j, k, imy) - dwd(i, j+1, k, imy)
        rqspd = rqspd + w(i, j+1, k, ivy)*fsd
        wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) + rqsp*fsd
        rqsmd = rqsmd + w(i, j, k, ivy)*fsd
        wd(i, j, k, ivy) = wd(i, j, k, ivy) + rqsm*fsd
        pad = pad + sj(i, j, k, 2)*fsd
        fsd = dwd(i, j, k, imx) - dwd(i, j+1, k, imx)
        rqspd = rqspd + w(i, j+1, k, ivx)*fsd
        wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) + rqsp*fsd
        rqsmd = rqsmd + w(i, j, k, ivx)*fsd
        wd(i, j, k, ivx) = wd(i, j, k, ivx) + rqsm*fsd
        pad = pad + sj(i, j, k, 1)*fsd
        fsd = dwd(i, j, k, irho) - dwd(i, j+1, k, irho)
        rqspd = rqspd + fsd
        rqsmd = rqsmd + fsd
        pd(i, j+1, k) = pd(i, j+1, k) + porflux*pad
        pd(i, j, k) = pd(i, j, k) + porflux*pad
        qsmd = qsmd + w(i, j, k, irho)*rqsmd
        vnmd = porvel*qsmd + p(i, j, k)*tempd0
        wd(i, j, k, irho) = wd(i, j, k, irho) + qsm*rqsmd
        qspd = qspd + w(i, j+1, k, irho)*rqspd
        vnpd = porvel*qspd + p(i, j+1, k)*tempd0
        wd(i, j+1, k, irho) = wd(i, j+1, k, irho) + qsp*rqspd
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) then
            vnmd = 0.0_8
            vnpd = 0.0_8
        end if
        wd(i, j, k, ivx) = wd(i, j, k, ivx) + sj(i, j, k, 1)*vnmd
        wd(i, j, k, ivy) = wd(i, j, k, ivy) + sj(i, j, k, 2)*vnmd
        wd(i, j, k, ivz) = wd(i, j, k, ivz) + sj(i, j, k, 3)*vnmd
        wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) + sj(i, j, k, 1)*vnpd
        wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) + sj(i, j, k, 2)*vnpd
        wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) + sj(i, j, k, 3)*vnpd
    end do
    ! initialize sface to zero. this value will be used if the
    ! block is not moving.
    sface = zero
    do ii=0,il*ny*nz-1
        i = mod(ii, il) + 1
        j = mod(ii/il, ny) + 2
        k = ii/(il*ny) + 2
        ! set the dot product of the grid velocity and the
        ! normal in i-direction for a moving face.
        if (addgridvelocities) sface = sfacei(i, j, k)
        ! compute the normal velocities of the left and right state.
        vnp = w(i+1, j, k, ivx)*si(i, j, k, 1) + w(i+1, j, k, ivy)*si(i, j&
        &       , k, 2) + w(i+1, j, k, ivz)*si(i, j, k, 3)
        vnm = w(i, j, k, ivx)*si(i, j, k, 1) + w(i, j, k, ivy)*si(i, j, k&
        &       , 2) + w(i, j, k, ivz)*si(i, j, k, 3)
        ! set the values of the porosities for this face.
        ! porvel defines the porosity w.r.t. velocity;
        ! porflux defines the porosity w.r.t. the entire flux.
        ! the latter is only zero for a discontinuous block
        ! boundary that must be treated conservatively.
        ! the default value of porflux is 0.5, such that the
        ! correct central flux is scattered to both cells.
        ! in case of a boundflux the normal velocity is set
        ! to sface.
        porvel = one
        porflux = half
        if (pori(i, j, k) .eq. noflux) porflux = zero
        if (pori(i, j, k) .eq. boundflux) then
            porvel = zero
            vnp = sface
            vnm = sface
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 0
        else
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 1
        end if
        ! incorporate porflux in porvel.
        porvel = porvel*porflux
        ! compute the normal velocities relative to the grid for
        ! the face as well as the mass fluxes.
        qsp = (vnp-sface)*porvel
        qsm = (vnm-sface)*porvel
        rqsp = qsp*w(i+1, j, k, irho)
        rqsm = qsm*w(i, j, k, irho)
        ! compute the sum of the pressure multiplied by porflux.
        ! for the default value of porflux, 0.5, this leads to
        ! the average pressure.
        ! compute the fluxes and scatter them to the cells
        ! i,j,k and i+1,j,k. store the density flux in the
        ! mass flow of the appropriate sliding mesh interface.
        fsd = dwd(i, j, k, irhoe) - dwd(i+1, j, k, irhoe)
        tempd = porflux*fsd
        qspd = w(i+1, j, k, irhoe)*fsd
        wd(i+1, j, k, irhoe) = wd(i+1, j, k, irhoe) + qsp*fsd
        qsmd = w(i, j, k, irhoe)*fsd
        wd(i, j, k, irhoe) = wd(i, j, k, irhoe) + qsm*fsd
        pd(i+1, j, k) = pd(i+1, j, k) + vnp*tempd
        pd(i, j, k) = pd(i, j, k) + vnm*tempd
        fsd = dwd(i, j, k, imz) - dwd(i+1, j, k, imz)
        rqspd = w(i+1, j, k, ivz)*fsd
        wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) + rqsp*fsd
        rqsmd = w(i, j, k, ivz)*fsd
        wd(i, j, k, ivz) = wd(i, j, k, ivz) + rqsm*fsd
        pad = si(i, j, k, 3)*fsd
        fsd = dwd(i, j, k, imy) - dwd(i+1, j, k, imy)
        rqspd = rqspd + w(i+1, j, k, ivy)*fsd
        wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) + rqsp*fsd
        rqsmd = rqsmd + w(i, j, k, ivy)*fsd
        wd(i, j, k, ivy) = wd(i, j, k, ivy) + rqsm*fsd
        pad = pad + si(i, j, k, 2)*fsd
        fsd = dwd(i, j, k, imx) - dwd(i+1, j, k, imx)
        rqspd = rqspd + w(i+1, j, k, ivx)*fsd
        wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) + rqsp*fsd
        rqsmd = rqsmd + w(i, j, k, ivx)*fsd
        wd(i, j, k, ivx) = wd(i, j, k, ivx) + rqsm*fsd
        pad = pad + si(i, j, k, 1)*fsd
        fsd = dwd(i, j, k, irho) - dwd(i+1, j, k, irho)
        rqspd = rqspd + fsd
        rqsmd = rqsmd + fsd
        pd(i+1, j, k) = pd(i+1, j, k) + porflux*pad
        pd(i, j, k) = pd(i, j, k) + porflux*pad
        qsmd = qsmd + w(i, j, k, irho)*rqsmd
        vnmd = porvel*qsmd + p(i, j, k)*tempd
        wd(i, j, k, irho) = wd(i, j, k, irho) + qsm*rqsmd
        qspd = qspd + w(i+1, j, k, irho)*rqspd
        vnpd = porvel*qspd + p(i+1, j, k)*tempd
        wd(i+1, j, k, irho) = wd(i+1, j, k, irho) + qsp*rqspd
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) then
            vnmd = 0.0_8
            vnpd = 0.0_8
        end if
        wd(i, j, k, ivx) = wd(i, j, k, ivx) + si(i, j, k, 1)*vnmd
        wd(i, j, k, ivy) = wd(i, j, k, ivy) + si(i, j, k, 2)*vnmd
        wd(i, j, k, ivz) = wd(i, j, k, ivz) + si(i, j, k, 3)*vnmd
        wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) + si(i, j, k, 1)*vnpd
        wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) + si(i, j, k, 2)*vnpd
        wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) + si(i, j, k, 3)*vnpd
    end do
end subroutine inviscidcentralflux_fast_b

!  differentiation of saresscale in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: *dw
!   with respect to varying inputs: *dw *scratch
!   rw status of diff variables: *dw:in-out *scratch:out
!   plus diff mem management of: dw:in scratch:in
subroutine saresscale_fast_b()
    !
    !  multiply the residual by the volume and store this in dw; this
    ! * is done for monitoring reasons only. the multiplication with the
    ! * volume is present to be consistent with the flow residuals; also
    !  the negative value is taken, again to be consistent with the
    ! * flow equations. also multiply by iblank so that no updates occur
    !  in holes or the overset boundary.
    ! use blockpointers
    implicit none
    ! local variables
    integer(kind=inttype) :: i, j, k, ii
    real(kind=realtype) :: rblank
    intrinsic mod
    intrinsic real
    intrinsic max
    real(kind=realtype) :: x1
    scratchd = 0.0_8
    do ii=0,nx*ny*nz-1
        i = mod(ii, nx) + 2
        j = mod(ii/nx, ny) + 2
        k = ii/(nx*ny) + 2
        x1 = real(iblank(i, j, k), realtype)
        if (x1 .lt. zero) then
            rblank = zero
        else
            rblank = x1
        end if
        scratchd(i, j, k, idvt) = scratchd(i, j, k, idvt) - volref(i, j, k&
        &       )*rblank*dwd(i, j, k, itu1)
        dwd(i, j, k, itu1) = 0.0_8
    end do
end subroutine saresscale_fast_b

!  differentiation of saviscous in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: *w *rlv *scratch
!   with respect to varying inputs: *w *rlv *scratch
!   rw status of diff variables: *w:incr *rlv:incr *scratch:in-out
!   plus diff mem management of: w:in rlv:in scratch:in
subroutine saviscous_fast_b()
    !
    !  viscous term.
    !  determine the viscous contribution to the residual
    !  for all internal cells of the block.
    ! use blockpointers
    use paramturb
    use sa_fast_b, only : cv13, kar2Inv, cw36, cb3Inv
    implicit none
    ! local variables.
    integer(kind=inttype) :: i, j, k, nn, ii
    real(kind=realtype) :: nu
    real(kind=realtype) :: nud
    real(kind=realtype) :: fv1, fv2, ft2
    real(kind=realtype) :: voli, volmi, volpi, xm, ym, zm, xp, yp, zp
    real(kind=realtype) :: xa, ya, za, ttm, ttp, cnud, cam, cap
    real(kind=realtype) :: cnudd, camd, capd
    real(kind=realtype) :: nutm, nutp, num, nup, cdm, cdp
    real(kind=realtype) :: nutmd, nutpd, numd, nupd, cdmd, cdpd
    real(kind=realtype) :: c1m, c1p, c10, b1, c1, d1, qs
    real(kind=realtype) :: c1md, c1pd, c10d
    intrinsic mod
    intrinsic max
    integer :: branch
    real(kind=realtype) :: temp3
    real(kind=realtype) :: temp2
    real(kind=realtype) :: temp1
    real(kind=realtype) :: temp0
    real(kind=realtype) :: tempd10
    real(kind=realtype) :: tempd9
    real(kind=realtype) :: tempd
    real(kind=realtype) :: tempd8
    real(kind=realtype) :: tempd7
    real(kind=realtype) :: tempd6
    real(kind=realtype) :: tempd5
    real(kind=realtype) :: tempd4
    real(kind=realtype) :: tempd3
    real(kind=realtype) :: tempd2
    real(kind=realtype) :: tempd1
    real(kind=realtype) :: tempd0
    real(kind=realtype) :: temp
    real(kind=realtype) :: temp7
    real(kind=realtype) :: temp6
    real(kind=realtype) :: temp5
    real(kind=realtype) :: temp4
    ! set model constants
    cb3inv = one/rsacb3
    do ii=0,nx*ny*nz-1
        i = mod(ii, nx) + 2
        j = mod(ii/nx, ny) + 2
        k = ii/(nx*ny) + 2
        ! compute the metrics in xi-direction, i.e. along the
        ! line i = constant.
        voli = one/vol(i, j, k)
        volmi = two/(vol(i, j, k)+vol(i-1, j, k))
        volpi = two/(vol(i, j, k)+vol(i+1, j, k))
        xm = si(i-1, j, k, 1)*volmi
        ym = si(i-1, j, k, 2)*volmi
        zm = si(i-1, j, k, 3)*volmi
        xp = si(i, j, k, 1)*volpi
        yp = si(i, j, k, 2)*volpi
        zp = si(i, j, k, 3)*volpi
        xa = half*(si(i, j, k, 1)+si(i-1, j, k, 1))*voli
        ya = half*(si(i, j, k, 2)+si(i-1, j, k, 2))*voli
        za = half*(si(i, j, k, 3)+si(i-1, j, k, 3))*voli
        ttm = xm*xa + ym*ya + zm*za
        ttp = xp*xa + yp*ya + zp*za
        ! computation of the viscous terms in xi-direction; note
        ! that cross-derivatives are neglected, i.e. the mesh is
        ! assumed to be orthogonal.
        ! furthermore, the grad(nu)**2 has been rewritten as
        ! div(nu grad(nu)) - nu div(grad nu) to enhance stability.
        ! the second derivative in xi-direction is constructed as
        ! the central difference of the first order derivatives, i.e.
        ! d^2/dxi^2 = d/dxi (d/dxi i+1/2 - d/dxi i-1/2).
        ! in this way the metric can be taken into account.
        ! compute the diffusion coefficients multiplying the nodes
        ! i+1, i and i-1 in the second derivative. make sure that
        ! these coefficients are nonnegative.
        cnud = -(rsacb2*w(i, j, k, itu1)*cb3inv)
        cam = ttm*cnud
        cap = ttp*cnud
        nutm = half*(w(i-1, j, k, itu1)+w(i, j, k, itu1))
        nutp = half*(w(i+1, j, k, itu1)+w(i, j, k, itu1))
        nu = rlv(i, j, k)/w(i, j, k, irho)
        num = half*(rlv(i-1, j, k)/w(i-1, j, k, irho)+nu)
        nup = half*(rlv(i+1, j, k)/w(i+1, j, k, irho)+nu)
        cdm = (num+(one+rsacb2)*nutm)*ttm*cb3inv
        cdp = (nup+(one+rsacb2)*nutp)*ttp*cb3inv
        if (cdm + cam .lt. zero) then
            c1m = zero
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 0
        else
            c1m = cdm + cam
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 1
        end if
        if (cdp + cap .lt. zero) then
            c1p = zero
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 0
        else
            c1p = cdp + cap
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 1
        end if
        c10 = c1m + c1p
        ! update the residual for this cell and store the possible
        ! coefficients for the matrix in b1, c1 and d1.
        c1md = w(i-1, j, k, itu1)*scratchd(i, j, k, idvt)
        wd(i-1, j, k, itu1) = wd(i-1, j, k, itu1) + c1m*scratchd(i, j, k, &
        &       idvt)
        c1pd = w(i+1, j, k, itu1)*scratchd(i, j, k, idvt)
        wd(i+1, j, k, itu1) = wd(i+1, j, k, itu1) + c1p*scratchd(i, j, k, &
        &       idvt)
        c10d = -(w(i, j, k, itu1)*scratchd(i, j, k, idvt))
        wd(i, j, k, itu1) = wd(i, j, k, itu1) - c10*scratchd(i, j, k, idvt&
        &       )
        c1md = c1md + c10d
        c1pd = c1pd + c10d
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) then
            cdpd = 0.0_8
            capd = 0.0_8
        else
            cdpd = c1pd
            capd = c1pd
        end if
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) then
            cdmd = 0.0_8
            camd = 0.0_8
        else
            cdmd = c1md
            camd = c1md
        end if
        tempd7 = ttp*cb3inv*cdpd
        nupd = tempd7
        nutpd = (one+rsacb2)*tempd7
        tempd8 = ttm*cb3inv*cdmd
        numd = tempd8
        nutmd = (one+rsacb2)*tempd8
        temp7 = w(i+1, j, k, irho)
        tempd9 = half*nupd/temp7
        rlvd(i+1, j, k) = rlvd(i+1, j, k) + tempd9
        wd(i+1, j, k, irho) = wd(i+1, j, k, irho) - rlv(i+1, j, k)*tempd9/&
        &       temp7
        nud = half*numd + half*nupd
        temp6 = w(i-1, j, k, irho)
        tempd10 = half*numd/temp6
        rlvd(i-1, j, k) = rlvd(i-1, j, k) + tempd10
        wd(i-1, j, k, irho) = wd(i-1, j, k, irho) - rlv(i-1, j, k)*tempd10&
        &       /temp6
        temp5 = w(i, j, k, irho)
        rlvd(i, j, k) = rlvd(i, j, k) + nud/temp5
        wd(i, j, k, irho) = wd(i, j, k, irho) - rlv(i, j, k)*nud/temp5**2
        wd(i+1, j, k, itu1) = wd(i+1, j, k, itu1) + half*nutpd
        wd(i, j, k, itu1) = wd(i, j, k, itu1) + half*nutpd
        wd(i-1, j, k, itu1) = wd(i-1, j, k, itu1) + half*nutmd
        wd(i, j, k, itu1) = wd(i, j, k, itu1) + half*nutmd
        cnudd = ttm*camd + ttp*capd
        wd(i, j, k, itu1) = wd(i, j, k, itu1) - rsacb2*cb3inv*cnudd
    end do
    do ii=0,nx*ny*nz-1
        i = mod(ii, nx) + 2
        j = mod(ii/nx, ny) + 2
        k = ii/(nx*ny) + 2
        ! compute the metrics in eta-direction, i.e. along the
        ! line j = constant.
        voli = one/vol(i, j, k)
        volmi = two/(vol(i, j, k)+vol(i, j-1, k))
        volpi = two/(vol(i, j, k)+vol(i, j+1, k))
        xm = sj(i, j-1, k, 1)*volmi
        ym = sj(i, j-1, k, 2)*volmi
        zm = sj(i, j-1, k, 3)*volmi
        xp = sj(i, j, k, 1)*volpi
        yp = sj(i, j, k, 2)*volpi
        zp = sj(i, j, k, 3)*volpi
        xa = half*(sj(i, j, k, 1)+sj(i, j-1, k, 1))*voli
        ya = half*(sj(i, j, k, 2)+sj(i, j-1, k, 2))*voli
        za = half*(sj(i, j, k, 3)+sj(i, j-1, k, 3))*voli
        ttm = xm*xa + ym*ya + zm*za
        ttp = xp*xa + yp*ya + zp*za
        ! computation of the viscous terms in eta-direction; note
        ! that cross-derivatives are neglected, i.e. the mesh is
        ! assumed to be orthogonal.
        ! furthermore, the grad(nu)**2 has been rewritten as
        ! div(nu grad(nu)) - nu div(grad nu) to enhance stability.
        ! the second derivative in eta-direction is constructed as
        ! the central difference of the first order derivatives, i.e.
        ! d^2/deta^2 = d/deta (d/deta j+1/2 - d/deta j-1/2).
        ! in this way the metric can be taken into account.
        ! compute the diffusion coefficients multiplying the nodes
        ! j+1, j and j-1 in the second derivative. make sure that
        ! these coefficients are nonnegative.
        cnud = -(rsacb2*w(i, j, k, itu1)*cb3inv)
        cam = ttm*cnud
        cap = ttp*cnud
        nutm = half*(w(i, j-1, k, itu1)+w(i, j, k, itu1))
        nutp = half*(w(i, j+1, k, itu1)+w(i, j, k, itu1))
        nu = rlv(i, j, k)/w(i, j, k, irho)
        num = half*(rlv(i, j-1, k)/w(i, j-1, k, irho)+nu)
        nup = half*(rlv(i, j+1, k)/w(i, j+1, k, irho)+nu)
        cdm = (num+(one+rsacb2)*nutm)*ttm*cb3inv
        cdp = (nup+(one+rsacb2)*nutp)*ttp*cb3inv
        if (cdm + cam .lt. zero) then
            c1m = zero
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 0
        else
            c1m = cdm + cam
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 1
        end if
        if (cdp + cap .lt. zero) then
            c1p = zero
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 0
        else
            c1p = cdp + cap
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 1
        end if
        c10 = c1m + c1p
        ! update the residual for this cell and store the possible
        ! coefficients for the matrix in b1, c1 and d1.
        c1md = w(i, j-1, k, itu1)*scratchd(i, j, k, idvt)
        wd(i, j-1, k, itu1) = wd(i, j-1, k, itu1) + c1m*scratchd(i, j, k, &
        &       idvt)
        c1pd = w(i, j+1, k, itu1)*scratchd(i, j, k, idvt)
        wd(i, j+1, k, itu1) = wd(i, j+1, k, itu1) + c1p*scratchd(i, j, k, &
        &       idvt)
        c10d = -(w(i, j, k, itu1)*scratchd(i, j, k, idvt))
        wd(i, j, k, itu1) = wd(i, j, k, itu1) - c10*scratchd(i, j, k, idvt&
        &       )
        c1md = c1md + c10d
        c1pd = c1pd + c10d
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) then
            cdpd = 0.0_8
            capd = 0.0_8
        else
            cdpd = c1pd
            capd = c1pd
        end if
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) then
            cdmd = 0.0_8
            camd = 0.0_8
        else
            cdmd = c1md
            camd = c1md
        end if
        tempd3 = ttp*cb3inv*cdpd
        nupd = tempd3
        nutpd = (one+rsacb2)*tempd3
        tempd4 = ttm*cb3inv*cdmd
        numd = tempd4
        nutmd = (one+rsacb2)*tempd4
        temp4 = w(i, j+1, k, irho)
        tempd5 = half*nupd/temp4
        rlvd(i, j+1, k) = rlvd(i, j+1, k) + tempd5
        wd(i, j+1, k, irho) = wd(i, j+1, k, irho) - rlv(i, j+1, k)*tempd5/&
        &       temp4
        nud = half*numd + half*nupd
        temp3 = w(i, j-1, k, irho)
        tempd6 = half*numd/temp3
        rlvd(i, j-1, k) = rlvd(i, j-1, k) + tempd6
        wd(i, j-1, k, irho) = wd(i, j-1, k, irho) - rlv(i, j-1, k)*tempd6/&
        &       temp3
        temp2 = w(i, j, k, irho)
        rlvd(i, j, k) = rlvd(i, j, k) + nud/temp2
        wd(i, j, k, irho) = wd(i, j, k, irho) - rlv(i, j, k)*nud/temp2**2
        wd(i, j+1, k, itu1) = wd(i, j+1, k, itu1) + half*nutpd
        wd(i, j, k, itu1) = wd(i, j, k, itu1) + half*nutpd
        wd(i, j-1, k, itu1) = wd(i, j-1, k, itu1) + half*nutmd
        wd(i, j, k, itu1) = wd(i, j, k, itu1) + half*nutmd
        cnudd = ttm*camd + ttp*capd
        wd(i, j, k, itu1) = wd(i, j, k, itu1) - rsacb2*cb3inv*cnudd
    end do
    do ii=0,nx*ny*nz-1
        i = mod(ii, nx) + 2
        j = mod(ii/nx, ny) + 2
        k = ii/(nx*ny) + 2
        ! compute the metrics in zeta-direction, i.e. along the
        ! line k = constant.
        voli = one/vol(i, j, k)
        volmi = two/(vol(i, j, k)+vol(i, j, k-1))
        volpi = two/(vol(i, j, k)+vol(i, j, k+1))
        xm = sk(i, j, k-1, 1)*volmi
        ym = sk(i, j, k-1, 2)*volmi
        zm = sk(i, j, k-1, 3)*volmi
        xp = sk(i, j, k, 1)*volpi
        yp = sk(i, j, k, 2)*volpi
        zp = sk(i, j, k, 3)*volpi
        xa = half*(sk(i, j, k, 1)+sk(i, j, k-1, 1))*voli
        ya = half*(sk(i, j, k, 2)+sk(i, j, k-1, 2))*voli
        za = half*(sk(i, j, k, 3)+sk(i, j, k-1, 3))*voli
        ttm = xm*xa + ym*ya + zm*za
        ttp = xp*xa + yp*ya + zp*za
        ! ttm and ttp ~ 1/deltax^2
        ! computation of the viscous terms in zeta-direction; note
        ! that cross-derivatives are neglected, i.e. the mesh is
        ! assumed to be orthogonal.
        ! furthermore, the grad(nu)**2 has been rewritten as
        ! div(nu grad(nu)) - nu div(grad nu) to enhance stability.
        ! the second derivative in zeta-direction is constructed as
        ! the central difference of the first order derivatives, i.e.
        ! d^2/dzeta^2 = d/dzeta (d/dzeta k+1/2 - d/dzeta k-1/2).
        ! in this way the metric can be taken into account.
        ! compute the diffusion coefficients multiplying the nodes
        ! k+1, k and k-1 in the second derivative. make sure that
        ! these coefficients are nonnegative.
        cnud = -(rsacb2*w(i, j, k, itu1)*cb3inv)
        cam = ttm*cnud
        cap = ttp*cnud
        ! compute nutilde at the faces
        nutm = half*(w(i, j, k-1, itu1)+w(i, j, k, itu1))
        nutp = half*(w(i, j, k+1, itu1)+w(i, j, k, itu1))
        ! compute nu at the faces
        nu = rlv(i, j, k)/w(i, j, k, irho)
        num = half*(rlv(i, j, k-1)/w(i, j, k-1, irho)+nu)
        nup = half*(rlv(i, j, k+1)/w(i, j, k+1, irho)+nu)
        cdm = (num+(one+rsacb2)*nutm)*ttm*cb3inv
        cdp = (nup+(one+rsacb2)*nutp)*ttp*cb3inv
        if (cdm + cam .lt. zero) then
            c1m = zero
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 0
        else
            c1m = cdm + cam
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 1
        end if
        if (cdp + cap .lt. zero) then
            c1p = zero
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 0
        else
            c1p = cdp + cap
            myIntPtr = myIntPtr + 1
            myIntStack(myIntPtr) = 1
        end if
        c10 = c1m + c1p
        ! update the residual for this cell and store the possible
        ! coefficients for the matrix in b1, c1 and d1.
        c1md = w(i, j, k-1, itu1)*scratchd(i, j, k, idvt)
        wd(i, j, k-1, itu1) = wd(i, j, k-1, itu1) + c1m*scratchd(i, j, k, &
        &       idvt)
        c1pd = w(i, j, k+1, itu1)*scratchd(i, j, k, idvt)
        wd(i, j, k+1, itu1) = wd(i, j, k+1, itu1) + c1p*scratchd(i, j, k, &
        &       idvt)
        c10d = -(w(i, j, k, itu1)*scratchd(i, j, k, idvt))
        wd(i, j, k, itu1) = wd(i, j, k, itu1) - c10*scratchd(i, j, k, idvt&
        &       )
        c1md = c1md + c10d
        c1pd = c1pd + c10d
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) then
            cdpd = 0.0_8
            capd = 0.0_8
        else
            cdpd = c1pd
            capd = c1pd
        end if
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) then
            cdmd = 0.0_8
            camd = 0.0_8
        else
            cdmd = c1md
            camd = c1md
        end if
        tempd = ttp*cb3inv*cdpd
        nupd = tempd
        nutpd = (one+rsacb2)*tempd
        tempd0 = ttm*cb3inv*cdmd
        numd = tempd0
        nutmd = (one+rsacb2)*tempd0
        temp1 = w(i, j, k+1, irho)
        tempd1 = half*nupd/temp1
        rlvd(i, j, k+1) = rlvd(i, j, k+1) + tempd1
        wd(i, j, k+1, irho) = wd(i, j, k+1, irho) - rlv(i, j, k+1)*tempd1/&
        &       temp1
        nud = half*numd + half*nupd
        temp0 = w(i, j, k-1, irho)
        tempd2 = half*numd/temp0
        rlvd(i, j, k-1) = rlvd(i, j, k-1) + tempd2
        wd(i, j, k-1, irho) = wd(i, j, k-1, irho) - rlv(i, j, k-1)*tempd2/&
        &       temp0
        temp = w(i, j, k, irho)
        rlvd(i, j, k) = rlvd(i, j, k) + nud/temp
        wd(i, j, k, irho) = wd(i, j, k, irho) - rlv(i, j, k)*nud/temp**2
        wd(i, j, k+1, itu1) = wd(i, j, k+1, itu1) + half*nutpd
        wd(i, j, k, itu1) = wd(i, j, k, itu1) + half*nutpd
        wd(i, j, k-1, itu1) = wd(i, j, k-1, itu1) + half*nutmd
        wd(i, j, k, itu1) = wd(i, j, k, itu1) + half*nutmd
        cnudd = ttm*camd + ttp*capd
        wd(i, j, k, itu1) = wd(i, j, k, itu1) - rsacb2*cb3inv*cnudd
    end do
end subroutine saviscous_fast_b

!  differentiation of turbadvection in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: *w *scratch
!   with respect to varying inputs: *w *scratch
!   rw status of diff variables: *w:incr *scratch:in-out
!   plus diff mem management of: w:in scratch:in
  subroutine turbadvection_fast_b(madv, nadv, offset, qq)
!
!       turbadvection discretizes the advection part of the turbulent
!       transport equations. as the advection part is the same for all
!       models, this generic routine can be used. both the
!       discretization and the central jacobian are computed in this
!       subroutine. the former can either be 1st or 2nd order
!       accurate; the latter is always based on the 1st order upwind
!       discretization. when the discretization must be second order
!       accurate, the fully upwind (kappa = -1) scheme in combination
!       with the minmod limiter is used.
!       only nadv equations are treated, while the actual system has
!       size madv. the reason is that some equations for some
!       turbulence equations do not have an advection part, e.g. the
!       f equation in the v2-f model. the argument offset indicates
!       the offset in the w vector where this subsystem starts. as a
!       consequence it is assumed that the indices of the current
!       subsystem are contiguous, e.g. if a 2*2 system is solved the
!       last index in w is offset+1 and offset+2 respectively.
!
    use constants
!     use blockpointers, only : nx, ny, nz, il, jl, kl, vol, sfacei, &
! &   sfacej, sfacek, w, wd, si, sj, sk, addgridvelocities, bmti1, bmti2, &
! &   bmtj1, bmtj2, bmtk1, bmtk2, scratch, scratchd
    use blockpointers, only : addgridvelocities
    use turbmod, only : secondord
    implicit none
!
!      subroutine arguments.
!
    integer(kind=inttype), intent(in) :: nadv, madv, offset
    real(kind=realtype), dimension(2:il, 2:jl, 2:kl, madv, madv), &
&   intent(inout) :: qq
!
!      local variables.
!
    integer(kind=inttype) :: i, j, k, ii, jj, kk, iii
    real(kind=realtype) :: qs, voli, xa, ya, za
    real(kind=realtype) :: uu, dwt, dwtm1, dwtp1, dwti, dwtj, dwtk
    real(kind=realtype) :: uud, dwtd, dwtm1d, dwtp1d, dwtid, dwtjd, &
&   dwtkd
    real(kind=realtype), dimension(madv) :: impl
    intrinsic mod
    intrinsic abs
    integer :: branch
    real(kind=realtype) :: abs23
    real(kind=realtype) :: abs22
    real(kind=realtype) :: abs21
    real(kind=realtype) :: abs20
    real(kind=realtype) :: abs19
    real(kind=realtype) :: abs18
    real(kind=realtype) :: abs17
    real(kind=realtype) :: abs16
    real(kind=realtype) :: abs15
    real(kind=realtype) :: abs14
    real(kind=realtype) :: abs13
    real(kind=realtype) :: abs12
    real(kind=realtype) :: abs11
    real(kind=realtype) :: abs10
    real(kind=realtype) :: abs9
    real(kind=realtype) :: abs8
    real(kind=realtype) :: abs7
    real(kind=realtype) :: abs6
    real(kind=realtype) :: abs5
    real(kind=realtype) :: abs4
    real(kind=realtype) :: abs3
    real(kind=realtype) :: abs2
    real(kind=realtype) :: abs1
    real(kind=realtype) :: abs0
    qs = zero
    do iii=0,nx*ny*nz-1
      i = mod(iii, nx) + 2
      j = mod(iii/nx, ny) + 2
      k = iii/(nx*ny) + 2
! compute the grid velocity if present.
! it is taken as the average of i and i-1,
      voli = half/vol(i, j, k)
      if (addgridvelocities) qs = (sfacei(i, j, k)+sfacei(i-1, j, k))*&
&         voli
! compute the normal velocity, where the normal direction
! is taken as the average of faces i and i-1.
      xa = (si(i, j, k, 1)+si(i-1, j, k, 1))*voli
      ya = (si(i, j, k, 2)+si(i-1, j, k, 2))*voli
      za = (si(i, j, k, 3)+si(i-1, j, k, 3))*voli
      uu = xa*w(i, j, k, ivx) + ya*w(i, j, k, ivy) + za*w(i, j, k, ivz) &
&       - qs
! determine the situation we are having here, i.e. positive
! or negative normal velocity.
      if (uu .gt. zero) then
        uud = 0.0_8
        do 100 ii=1,nadv
! set the value of jj such that it corresponds to the
! turbulent entry in w.
          jj = ii + offset
! check whether a first or a second order discretization
! must be used.
          if (secondord) then
! second order; store the three differences for the
! discretization of the derivative in i-direction.
            dwtm1 = w(i-1, j, k, jj) - w(i-2, j, k, jj)
            dwt = w(i, j, k, jj) - w(i-1, j, k, jj)
            dwtp1 = w(i+1, j, k, jj) - w(i, j, k, jj)
! construct the derivative in this cell center. this is
! the first order upwind derivative with two nonlinear
! corrections.
            dwti = dwt
            if (dwt*dwtp1 .gt. zero) then
              if (dwt .ge. 0.) then
                abs8 = dwt
              else
                abs8 = -dwt
              end if
              if (dwtp1 .ge. 0.) then
                abs20 = dwtp1
              else
                abs20 = -dwtp1
              end if
              if (abs8 .lt. abs20) then
                dwti = dwti + half*dwt
                call pushcontrol2b(0)
              else
                dwti = dwti + half*dwtp1
                call pushcontrol2b(1)
              end if
            else
              call pushcontrol2b(2)
            end if
            if (dwt*dwtm1 .gt. zero) then
              if (dwt .ge. 0.) then
                abs9 = dwt
              else
                abs9 = -dwt
              end if
              if (dwtm1 .ge. 0.) then
                abs21 = dwtm1
              else
                abs21 = -dwtm1
              end if
              if (abs9 .lt. abs21) then
                dwti = dwti - half*dwt
                call pushcontrol2b(0)
              else
                dwti = dwti - half*dwtm1
                call pushcontrol2b(1)
              end if
            else
              call pushcontrol2b(2)
            end if
          else
! 1st order upwind scheme.
            dwti = w(i, j, k, jj) - w(i-1, j, k, jj)
            call pushcontrol2b(3)
          end if
          uud = uud - dwti*scratchd(i, j, k, idvt+ii-1)
          dwtid = -(uu*scratchd(i, j, k, idvt+ii-1))
          call popcontrol2b(branch)
          if (branch .lt. 2) then
            if (branch .eq. 0) then
              dwtd = -(half*dwtid)
              dwtm1d = 0.0_8
            else
              dwtm1d = -(half*dwtid)
              dwtd = 0.0_8
            end if
          else if (branch .eq. 2) then
            dwtd = 0.0_8
            dwtm1d = 0.0_8
          else
            wd(i, j, k, jj) = wd(i, j, k, jj) + dwtid
            wd(i-1, j, k, jj) = wd(i-1, j, k, jj) - dwtid
            goto 100
          end if
          call popcontrol2b(branch)
          if (branch .eq. 0) then
            dwtd = dwtd + half*dwtid
            dwtp1d = 0.0_8
          else if (branch .eq. 1) then
            dwtp1d = half*dwtid
          else
            dwtp1d = 0.0_8
          end if
          dwtd = dwtd + dwtid
          wd(i+1, j, k, jj) = wd(i+1, j, k, jj) + dwtp1d
          wd(i, j, k, jj) = wd(i, j, k, jj) - dwtp1d
          wd(i, j, k, jj) = wd(i, j, k, jj) + dwtd
          wd(i-1, j, k, jj) = wd(i-1, j, k, jj) - dwtd
          wd(i-1, j, k, jj) = wd(i-1, j, k, jj) + dwtm1d
          wd(i-2, j, k, jj) = wd(i-2, j, k, jj) - dwtm1d
 100    continue
      else
        uud = 0.0_8
        do 110 ii=1,nadv
! set the value of jj such that it corresponds to the
! turbulent entry in w.
          jj = ii + offset
! check whether a first or a second order discretization
! must be used.
          if (secondord) then
! second order; store the three differences for the
! discretization of the derivative in i-direction.
            dwtm1 = w(i, j, k, jj) - w(i-1, j, k, jj)
            dwt = w(i+1, j, k, jj) - w(i, j, k, jj)
            dwtp1 = w(i+2, j, k, jj) - w(i+1, j, k, jj)
! construct the derivative in this cell center. this is
! the first order upwind derivative with two nonlinear
! corrections.
            dwti = dwt
            if (dwt*dwtp1 .gt. zero) then
              if (dwt .ge. 0.) then
                abs10 = dwt
              else
                abs10 = -dwt
              end if
              if (dwtp1 .ge. 0.) then
                abs22 = dwtp1
              else
                abs22 = -dwtp1
              end if
              if (abs10 .lt. abs22) then
                dwti = dwti - half*dwt
                call pushcontrol2b(0)
              else
                dwti = dwti - half*dwtp1
                call pushcontrol2b(1)
              end if
            else
              call pushcontrol2b(2)
            end if
            if (dwt*dwtm1 .gt. zero) then
              if (dwt .ge. 0.) then
                abs11 = dwt
              else
                abs11 = -dwt
              end if
              if (dwtm1 .ge. 0.) then
                abs23 = dwtm1
              else
                abs23 = -dwtm1
              end if
              if (abs11 .lt. abs23) then
                dwti = dwti + half*dwt
                call pushcontrol2b(0)
              else
                dwti = dwti + half*dwtm1
                call pushcontrol2b(1)
              end if
            else
              call pushcontrol2b(2)
            end if
          else
! 1st order upwind scheme.
            dwti = w(i+1, j, k, jj) - w(i, j, k, jj)
            call pushcontrol2b(3)
          end if
          uud = uud - dwti*scratchd(i, j, k, idvt+ii-1)
          dwtid = -(uu*scratchd(i, j, k, idvt+ii-1))
          call popcontrol2b(branch)
          if (branch .lt. 2) then
            if (branch .eq. 0) then
              dwtd = half*dwtid
              dwtm1d = 0.0_8
            else
              dwtm1d = half*dwtid
              dwtd = 0.0_8
            end if
          else if (branch .eq. 2) then
            dwtd = 0.0_8
            dwtm1d = 0.0_8
          else
            wd(i+1, j, k, jj) = wd(i+1, j, k, jj) + dwtid
            wd(i, j, k, jj) = wd(i, j, k, jj) - dwtid
            goto 110
          end if
          call popcontrol2b(branch)
          if (branch .eq. 0) then
            dwtd = dwtd - half*dwtid
            dwtp1d = 0.0_8
          else if (branch .eq. 1) then
            dwtp1d = -(half*dwtid)
          else
            dwtp1d = 0.0_8
          end if
          dwtd = dwtd + dwtid
          wd(i+2, j, k, jj) = wd(i+2, j, k, jj) + dwtp1d
          wd(i+1, j, k, jj) = wd(i+1, j, k, jj) - dwtp1d
          wd(i+1, j, k, jj) = wd(i+1, j, k, jj) + dwtd
          wd(i, j, k, jj) = wd(i, j, k, jj) - dwtd
          wd(i, j, k, jj) = wd(i, j, k, jj) + dwtm1d
          wd(i-1, j, k, jj) = wd(i-1, j, k, jj) - dwtm1d
 110    continue
      end if
      wd(i, j, k, ivx) = wd(i, j, k, ivx) + xa*uud
      wd(i, j, k, ivy) = wd(i, j, k, ivy) + ya*uud
      wd(i, j, k, ivz) = wd(i, j, k, ivz) + za*uud
    end do
    qs = zero
    do iii=0,nx*ny*nz-1
      i = mod(iii, nx) + 2
      j = mod(iii/nx, ny) + 2
      k = iii/(nx*ny) + 2
! compute the grid velocity if present.
! it is taken as the average of j and j-1,
      voli = half/vol(i, j, k)
      if (addgridvelocities) qs = (sfacej(i, j, k)+sfacej(i, j-1, k))*&
&         voli
! compute the normal velocity, where the normal direction
! is taken as the average of faces j and j-1.
      xa = (sj(i, j, k, 1)+sj(i, j-1, k, 1))*voli
      ya = (sj(i, j, k, 2)+sj(i, j-1, k, 2))*voli
      za = (sj(i, j, k, 3)+sj(i, j-1, k, 3))*voli
      uu = xa*w(i, j, k, ivx) + ya*w(i, j, k, ivy) + za*w(i, j, k, ivz) &
&       - qs
! determine the situation we are having here, i.e. positive
! or negative normal velocity.
      if (uu .gt. zero) then
        uud = 0.0_8
        do 120 ii=1,nadv
! set the value of jj such that it corresponds to the
! turbulent entry in w.
          jj = ii + offset
! check whether a first or a second order discretization
! must be used.
          if (secondord) then
! second order; store the three differences for the
! discretization of the derivative in j-direction.
            dwtm1 = w(i, j-1, k, jj) - w(i, j-2, k, jj)
            dwt = w(i, j, k, jj) - w(i, j-1, k, jj)
            dwtp1 = w(i, j+1, k, jj) - w(i, j, k, jj)
! construct the derivative in this cell center. this is
! the first order upwind derivative with two nonlinear
! corrections.
            dwtj = dwt
            if (dwt*dwtp1 .gt. zero) then
              if (dwt .ge. 0.) then
                abs4 = dwt
              else
                abs4 = -dwt
              end if
              if (dwtp1 .ge. 0.) then
                abs16 = dwtp1
              else
                abs16 = -dwtp1
              end if
              if (abs4 .lt. abs16) then
                dwtj = dwtj + half*dwt
                call pushcontrol2b(0)
              else
                dwtj = dwtj + half*dwtp1
                call pushcontrol2b(1)
              end if
            else
              call pushcontrol2b(2)
            end if
            if (dwt*dwtm1 .gt. zero) then
              if (dwt .ge. 0.) then
                abs5 = dwt
              else
                abs5 = -dwt
              end if
              if (dwtm1 .ge. 0.) then
                abs17 = dwtm1
              else
                abs17 = -dwtm1
              end if
              if (abs5 .lt. abs17) then
                dwtj = dwtj - half*dwt
                call pushcontrol2b(0)
              else
                dwtj = dwtj - half*dwtm1
                call pushcontrol2b(1)
              end if
            else
              call pushcontrol2b(2)
            end if
          else
! 1st order upwind scheme.
            dwtj = w(i, j, k, jj) - w(i, j-1, k, jj)
            call pushcontrol2b(3)
          end if
          uud = uud - dwtj*scratchd(i, j, k, idvt+ii-1)
          dwtjd = -(uu*scratchd(i, j, k, idvt+ii-1))
          call popcontrol2b(branch)
          if (branch .lt. 2) then
            if (branch .eq. 0) then
              dwtd = -(half*dwtjd)
              dwtm1d = 0.0_8
            else
              dwtm1d = -(half*dwtjd)
              dwtd = 0.0_8
            end if
          else if (branch .eq. 2) then
            dwtd = 0.0_8
            dwtm1d = 0.0_8
          else
            wd(i, j, k, jj) = wd(i, j, k, jj) + dwtjd
            wd(i, j-1, k, jj) = wd(i, j-1, k, jj) - dwtjd
            goto 120
          end if
          call popcontrol2b(branch)
          if (branch .eq. 0) then
            dwtd = dwtd + half*dwtjd
            dwtp1d = 0.0_8
          else if (branch .eq. 1) then
            dwtp1d = half*dwtjd
          else
            dwtp1d = 0.0_8
          end if
          dwtd = dwtd + dwtjd
          wd(i, j+1, k, jj) = wd(i, j+1, k, jj) + dwtp1d
          wd(i, j, k, jj) = wd(i, j, k, jj) - dwtp1d
          wd(i, j, k, jj) = wd(i, j, k, jj) + dwtd
          wd(i, j-1, k, jj) = wd(i, j-1, k, jj) - dwtd
          wd(i, j-1, k, jj) = wd(i, j-1, k, jj) + dwtm1d
          wd(i, j-2, k, jj) = wd(i, j-2, k, jj) - dwtm1d
 120    continue
      else
        uud = 0.0_8
        do 130 ii=1,nadv
! set the value of jj such that it corresponds to the
! turbulent entry in w.
          jj = ii + offset
! check whether a first or a second order discretization
! must be used.
          if (secondord) then
! store the three differences for the discretization of
! the derivative in j-direction.
            dwtm1 = w(i, j, k, jj) - w(i, j-1, k, jj)
            dwt = w(i, j+1, k, jj) - w(i, j, k, jj)
            dwtp1 = w(i, j+2, k, jj) - w(i, j+1, k, jj)
! construct the derivative in this cell center. this is
! the first order upwind derivative with two nonlinear
! corrections.
            dwtj = dwt
            if (dwt*dwtp1 .gt. zero) then
              if (dwt .ge. 0.) then
                abs6 = dwt
              else
                abs6 = -dwt
              end if
              if (dwtp1 .ge. 0.) then
                abs18 = dwtp1
              else
                abs18 = -dwtp1
              end if
              if (abs6 .lt. abs18) then
                dwtj = dwtj - half*dwt
                call pushcontrol2b(0)
              else
                dwtj = dwtj - half*dwtp1
                call pushcontrol2b(1)
              end if
            else
              call pushcontrol2b(2)
            end if
            if (dwt*dwtm1 .gt. zero) then
              if (dwt .ge. 0.) then
                abs7 = dwt
              else
                abs7 = -dwt
              end if
              if (dwtm1 .ge. 0.) then
                abs19 = dwtm1
              else
                abs19 = -dwtm1
              end if
              if (abs7 .lt. abs19) then
                dwtj = dwtj + half*dwt
                call pushcontrol2b(0)
              else
                dwtj = dwtj + half*dwtm1
                call pushcontrol2b(1)
              end if
            else
              call pushcontrol2b(2)
            end if
          else
! 1st order upwind scheme.
            dwtj = w(i, j+1, k, jj) - w(i, j, k, jj)
            call pushcontrol2b(3)
          end if
          uud = uud - dwtj*scratchd(i, j, k, idvt+ii-1)
          dwtjd = -(uu*scratchd(i, j, k, idvt+ii-1))
          call popcontrol2b(branch)
          if (branch .lt. 2) then
            if (branch .eq. 0) then
              dwtd = half*dwtjd
              dwtm1d = 0.0_8
            else
              dwtm1d = half*dwtjd
              dwtd = 0.0_8
            end if
          else if (branch .eq. 2) then
            dwtd = 0.0_8
            dwtm1d = 0.0_8
          else
            wd(i, j+1, k, jj) = wd(i, j+1, k, jj) + dwtjd
            wd(i, j, k, jj) = wd(i, j, k, jj) - dwtjd
            goto 130
          end if
          call popcontrol2b(branch)
          if (branch .eq. 0) then
            dwtd = dwtd - half*dwtjd
            dwtp1d = 0.0_8
          else if (branch .eq. 1) then
            dwtp1d = -(half*dwtjd)
          else
            dwtp1d = 0.0_8
          end if
          dwtd = dwtd + dwtjd
          wd(i, j+2, k, jj) = wd(i, j+2, k, jj) + dwtp1d
          wd(i, j+1, k, jj) = wd(i, j+1, k, jj) - dwtp1d
          wd(i, j+1, k, jj) = wd(i, j+1, k, jj) + dwtd
          wd(i, j, k, jj) = wd(i, j, k, jj) - dwtd
          wd(i, j, k, jj) = wd(i, j, k, jj) + dwtm1d
          wd(i, j-1, k, jj) = wd(i, j-1, k, jj) - dwtm1d
 130    continue
      end if
      wd(i, j, k, ivx) = wd(i, j, k, ivx) + xa*uud
      wd(i, j, k, ivy) = wd(i, j, k, ivy) + ya*uud
      wd(i, j, k, ivz) = wd(i, j, k, ivz) + za*uud
    end do
! initialize the grid velocity to zero. this value will be used
! if the block is not moving.
    qs = zero
    do iii=0,nx*ny*nz-1
      i = mod(iii, nx) + 2
      j = mod(iii/nx, ny) + 2
      k = iii/(nx*ny) + 2
! compute the grid velocity if present.
! it is taken as the average of k and k-1,
      voli = half/vol(i, j, k)
      if (addgridvelocities) qs = (sfacek(i, j, k)+sfacek(i, j, k-1))*&
&         voli
! compute the normal velocity, where the normal direction
! is taken as the average of faces k and k-1.
      xa = (sk(i, j, k, 1)+sk(i, j, k-1, 1))*voli
      ya = (sk(i, j, k, 2)+sk(i, j, k-1, 2))*voli
      za = (sk(i, j, k, 3)+sk(i, j, k-1, 3))*voli
      uu = xa*w(i, j, k, ivx) + ya*w(i, j, k, ivy) + za*w(i, j, k, ivz) &
&       - qs
! this term has unit: velocity/length
! determine the situation we are having here, i.e. positive
! or negative normal velocity.
      if (uu .gt. zero) then
        uud = 0.0_8
        do 140 ii=1,nadv
! set the value of jj such that it corresponds to the
! turbulent entry in w.
          jj = ii + offset
! check whether a first or a second order discretization
! must be used.
          if (secondord) then
! second order; store the three differences for the
! discretization of the derivative in k-direction.
            dwtm1 = w(i, j, k-1, jj) - w(i, j, k-2, jj)
            dwt = w(i, j, k, jj) - w(i, j, k-1, jj)
            dwtp1 = w(i, j, k+1, jj) - w(i, j, k, jj)
! construct the derivative in this cell center. this
! is the first order upwind derivative with two
! nonlinear corrections.
            dwtk = dwt
            if (dwt*dwtp1 .gt. zero) then
              if (dwt .ge. 0.) then
                abs0 = dwt
              else
                abs0 = -dwt
              end if
              if (dwtp1 .ge. 0.) then
                abs12 = dwtp1
              else
                abs12 = -dwtp1
              end if
              if (abs0 .lt. abs12) then
                dwtk = dwtk + half*dwt
                call pushcontrol2b(0)
              else
                dwtk = dwtk + half*dwtp1
                call pushcontrol2b(1)
              end if
            else
              call pushcontrol2b(2)
            end if
            if (dwt*dwtm1 .gt. zero) then
              if (dwt .ge. 0.) then
                abs1 = dwt
              else
                abs1 = -dwt
              end if
              if (dwtm1 .ge. 0.) then
                abs13 = dwtm1
              else
                abs13 = -dwtm1
              end if
              if (abs1 .lt. abs13) then
                dwtk = dwtk - half*dwt
                call pushcontrol2b(0)
              else
                dwtk = dwtk - half*dwtm1
                call pushcontrol2b(1)
              end if
            else
              call pushcontrol2b(2)
            end if
          else
! 1st order upwind scheme.
            dwtk = w(i, j, k, jj) - w(i, j, k-1, jj)
            call pushcontrol2b(3)
          end if
          uud = uud - dwtk*scratchd(i, j, k, idvt+ii-1)
          dwtkd = -(uu*scratchd(i, j, k, idvt+ii-1))
          call popcontrol2b(branch)
          if (branch .lt. 2) then
            if (branch .eq. 0) then
              dwtd = -(half*dwtkd)
              dwtm1d = 0.0_8
            else
              dwtm1d = -(half*dwtkd)
              dwtd = 0.0_8
            end if
          else if (branch .eq. 2) then
            dwtd = 0.0_8
            dwtm1d = 0.0_8
          else
            wd(i, j, k, jj) = wd(i, j, k, jj) + dwtkd
            wd(i, j, k-1, jj) = wd(i, j, k-1, jj) - dwtkd
            goto 140
          end if
          call popcontrol2b(branch)
          if (branch .eq. 0) then
            dwtd = dwtd + half*dwtkd
            dwtp1d = 0.0_8
          else if (branch .eq. 1) then
            dwtp1d = half*dwtkd
          else
            dwtp1d = 0.0_8
          end if
          dwtd = dwtd + dwtkd
          wd(i, j, k+1, jj) = wd(i, j, k+1, jj) + dwtp1d
          wd(i, j, k, jj) = wd(i, j, k, jj) - dwtp1d
          wd(i, j, k, jj) = wd(i, j, k, jj) + dwtd
          wd(i, j, k-1, jj) = wd(i, j, k-1, jj) - dwtd
          wd(i, j, k-1, jj) = wd(i, j, k-1, jj) + dwtm1d
          wd(i, j, k-2, jj) = wd(i, j, k-2, jj) - dwtm1d
 140    continue
      else
        uud = 0.0_8
        do 150 ii=1,nadv
! set the value of jj such that it corresponds to the
! turbulent entry in w.
          jj = ii + offset
! check whether a first or a second order discretization
! must be used.
          if (secondord) then
! store the three differences for the discretization of
! the derivative in k-direction.
            dwtm1 = w(i, j, k, jj) - w(i, j, k-1, jj)
            dwt = w(i, j, k+1, jj) - w(i, j, k, jj)
            dwtp1 = w(i, j, k+2, jj) - w(i, j, k+1, jj)
! construct the derivative in this cell center. this is
! the first order upwind derivative with two nonlinear
! corrections.
            dwtk = dwt
            if (dwt*dwtp1 .gt. zero) then
              if (dwt .ge. 0.) then
                abs2 = dwt
              else
                abs2 = -dwt
              end if
              if (dwtp1 .ge. 0.) then
                abs14 = dwtp1
              else
                abs14 = -dwtp1
              end if
              if (abs2 .lt. abs14) then
                dwtk = dwtk - half*dwt
                call pushcontrol2b(0)
              else
                dwtk = dwtk - half*dwtp1
                call pushcontrol2b(1)
              end if
            else
              call pushcontrol2b(2)
            end if
            if (dwt*dwtm1 .gt. zero) then
              if (dwt .ge. 0.) then
                abs3 = dwt
              else
                abs3 = -dwt
              end if
              if (dwtm1 .ge. 0.) then
                abs15 = dwtm1
              else
                abs15 = -dwtm1
              end if
              if (abs3 .lt. abs15) then
                dwtk = dwtk + half*dwt
                call pushcontrol2b(0)
              else
                dwtk = dwtk + half*dwtm1
                call pushcontrol2b(1)
              end if
            else
              call pushcontrol2b(2)
            end if
          else
! 1st order upwind scheme.
            dwtk = w(i, j, k+1, jj) - w(i, j, k, jj)
            call pushcontrol2b(3)
          end if
          uud = uud - dwtk*scratchd(i, j, k, idvt+ii-1)
          dwtkd = -(uu*scratchd(i, j, k, idvt+ii-1))
          call popcontrol2b(branch)
          if (branch .lt. 2) then
            if (branch .eq. 0) then
              dwtd = half*dwtkd
              dwtm1d = 0.0_8
            else
              dwtm1d = half*dwtkd
              dwtd = 0.0_8
            end if
          else if (branch .eq. 2) then
            dwtd = 0.0_8
            dwtm1d = 0.0_8
          else
            wd(i, j, k+1, jj) = wd(i, j, k+1, jj) + dwtkd
            wd(i, j, k, jj) = wd(i, j, k, jj) - dwtkd
            goto 150
          end if
          call popcontrol2b(branch)
          if (branch .eq. 0) then
            dwtd = dwtd - half*dwtkd
            dwtp1d = 0.0_8
          else if (branch .eq. 1) then
            dwtp1d = -(half*dwtkd)
          else
            dwtp1d = 0.0_8
          end if
          dwtd = dwtd + dwtkd
          wd(i, j, k+2, jj) = wd(i, j, k+2, jj) + dwtp1d
          wd(i, j, k+1, jj) = wd(i, j, k+1, jj) - dwtp1d
          wd(i, j, k+1, jj) = wd(i, j, k+1, jj) + dwtd
          wd(i, j, k, jj) = wd(i, j, k, jj) - dwtd
          wd(i, j, k, jj) = wd(i, j, k, jj) + dwtm1d
          wd(i, j, k-1, jj) = wd(i, j, k-1, jj) - dwtm1d
 150    continue
      end if
      wd(i, j, k, ivx) = wd(i, j, k, ivx) + xa*uud
      wd(i, j, k, ivy) = wd(i, j, k, ivy) + ya*uud
      wd(i, j, k, ivz) = wd(i, j, k, ivz) + za*uud
    end do
end subroutine turbadvection_fast_b

!  differentiation of sasource in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: *w *rlv *scratch
!   with respect to varying inputs: *w *rlv *scratch
!   rw status of diff variables: *w:incr *rlv:incr *scratch:in-out
!   plus diff mem management of: w:in rlv:in scratch:in
subroutine sasource_fast_b()
    !
    !  source terms.
    !  determine the source term and its derivative w.r.t. nutilde
    !  for all internal cells of the block.
    !  remember that the sa field variable nutilde = w(i,j,k,itu1)
    ! use blockpointers
    use blockpointers, only : sectionid
    use constants
    use paramturb
    use section
    use inputphysics
    use inputdiscretization, only : approxsa
    use sa_fast_b, only : cv13, kar2Inv, cw36, cb3Inv
    use flowvarrefstate
    implicit none
    ! local parameters
    real(kind=realtype), parameter :: f23=two*third
    ! local variables.
    integer(kind=inttype) :: i, j, k, nn, ii
    real(kind=realtype) :: fv1, fv2, ft2
    real(kind=realtype) :: fv1d, fv2d, ft2d
    real(kind=realtype) :: ss, sst, nu, dist2inv, chi, chi2, chi3
    real(kind=realtype) :: ssd, sstd, nud, chid, chi2d, chi3d
    real(kind=realtype) :: rr, gg, gg6, termfw, fwsa, term1, term2
    real(kind=realtype) :: rrd, ggd, gg6d, termfwd, fwsad, term1d, &
    &   term2d
    real(kind=realtype) :: dfv1, dfv2, dft2, drr, dgg, dfw
    real(kind=realtype) :: uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz
    real(kind=realtype) :: uuxd, uuyd, uuzd, vvxd, vvyd, vvzd, wwxd, &
    &   wwyd, wwzd
    real(kind=realtype) :: div2, fact, sxx, syy, szz, sxy, sxz, syz
    real(kind=realtype) :: div2d, sxxd, syyd, szzd, sxyd, sxzd, syzd
    real(kind=realtype) :: vortx, vorty, vortz
    real(kind=realtype) :: vortxd, vortyd, vortzd
    real(kind=realtype) :: omegax, omegay, omegaz
    real(kind=realtype) :: strainmag2, strainprod, vortprod
    real(kind=realtype) :: strainmag2d, strainprodd, vortprodd
    real(kind=realtype), parameter :: xminn=1.e-10_realtype
    intrinsic mod
    intrinsic sqrt
    intrinsic exp
    intrinsic min
    intrinsic max
    integer :: branch
    real(kind=realtype) :: temp1
    real(kind=realtype) :: temp0
    real(kind=realtype) :: tempd10
    real(kind=realtype) :: min1
    real(kind=realtype) :: min1d
    real(kind=realtype) :: tempd9
    real(kind=realtype) :: tempd
    real(kind=realtype) :: tempd8
    real(kind=realtype) :: tempd7
    real(kind=realtype) :: tempd6
    real(kind=realtype) :: tempd5
    real(kind=realtype) :: tempd4
    real(kind=realtype) :: tempd3
    real(kind=realtype) :: tempd2
    real(kind=realtype) :: tempd1
    real(kind=realtype) :: tempd0
    real(kind=realtype) :: temp
    real(kind=realtype) :: y1
    real(kind=realtype) :: y1d
    ! set model constants
    cv13 = rsacv1**3
    kar2inv = one/rsak**2
    cw36 = rsacw3**6
    ! determine the non-dimensional wheel speed of this block.
    omegax = timeref*sections(sectionid)%rotrate(1)
    omegay = timeref*sections(sectionid)%rotrate(2)
    omegaz = timeref*sections(sectionid)%rotrate(3)
    ! create switches to production term depending on the variable that
    ! should be used
    if (turbprod .eq. katolaunder) then
        stop
    else
        strainmag2d = 0.0_8
        ssd = 0.0_8
        do ii=0,nx*ny*nz-1
            i = mod(ii, nx) + 2
            j = mod(ii/nx, ny) + 2
            k = ii/(nx*ny) + 2
            ! compute the gradient of u in the cell center. use is made
            ! of the fact that the surrounding normals sum up to zero,
            ! such that the cell i,j,k does not give a contribution.
            ! the gradient is scaled by the factor 2*vol.
            uux = w(i+1, j, k, ivx)*si(i, j, k, 1) - w(i-1, j, k, ivx)*si(i-&
            &         1, j, k, 1) + w(i, j+1, k, ivx)*sj(i, j, k, 1) - w(i, j-1, k, &
            &         ivx)*sj(i, j-1, k, 1) + w(i, j, k+1, ivx)*sk(i, j, k, 1) - w(i&
            &         , j, k-1, ivx)*sk(i, j, k-1, 1)
            uuy = w(i+1, j, k, ivx)*si(i, j, k, 2) - w(i-1, j, k, ivx)*si(i-&
            &         1, j, k, 2) + w(i, j+1, k, ivx)*sj(i, j, k, 2) - w(i, j-1, k, &
            &         ivx)*sj(i, j-1, k, 2) + w(i, j, k+1, ivx)*sk(i, j, k, 2) - w(i&
            &         , j, k-1, ivx)*sk(i, j, k-1, 2)
            uuz = w(i+1, j, k, ivx)*si(i, j, k, 3) - w(i-1, j, k, ivx)*si(i-&
            &         1, j, k, 3) + w(i, j+1, k, ivx)*sj(i, j, k, 3) - w(i, j-1, k, &
            &         ivx)*sj(i, j-1, k, 3) + w(i, j, k+1, ivx)*sk(i, j, k, 3) - w(i&
            &         , j, k-1, ivx)*sk(i, j, k-1, 3)
            ! idem for the gradient of v.
            vvx = w(i+1, j, k, ivy)*si(i, j, k, 1) - w(i-1, j, k, ivy)*si(i-&
            &         1, j, k, 1) + w(i, j+1, k, ivy)*sj(i, j, k, 1) - w(i, j-1, k, &
            &         ivy)*sj(i, j-1, k, 1) + w(i, j, k+1, ivy)*sk(i, j, k, 1) - w(i&
            &         , j, k-1, ivy)*sk(i, j, k-1, 1)
            vvy = w(i+1, j, k, ivy)*si(i, j, k, 2) - w(i-1, j, k, ivy)*si(i-&
            &         1, j, k, 2) + w(i, j+1, k, ivy)*sj(i, j, k, 2) - w(i, j-1, k, &
            &         ivy)*sj(i, j-1, k, 2) + w(i, j, k+1, ivy)*sk(i, j, k, 2) - w(i&
            &         , j, k-1, ivy)*sk(i, j, k-1, 2)
            vvz = w(i+1, j, k, ivy)*si(i, j, k, 3) - w(i-1, j, k, ivy)*si(i-&
            &         1, j, k, 3) + w(i, j+1, k, ivy)*sj(i, j, k, 3) - w(i, j-1, k, &
            &         ivy)*sj(i, j-1, k, 3) + w(i, j, k+1, ivy)*sk(i, j, k, 3) - w(i&
            &         , j, k-1, ivy)*sk(i, j, k-1, 3)
            ! and for the gradient of w.
            wwx = w(i+1, j, k, ivz)*si(i, j, k, 1) - w(i-1, j, k, ivz)*si(i-&
            &         1, j, k, 1) + w(i, j+1, k, ivz)*sj(i, j, k, 1) - w(i, j-1, k, &
            &         ivz)*sj(i, j-1, k, 1) + w(i, j, k+1, ivz)*sk(i, j, k, 1) - w(i&
            &         , j, k-1, ivz)*sk(i, j, k-1, 1)
            wwy = w(i+1, j, k, ivz)*si(i, j, k, 2) - w(i-1, j, k, ivz)*si(i-&
            &         1, j, k, 2) + w(i, j+1, k, ivz)*sj(i, j, k, 2) - w(i, j-1, k, &
            &         ivz)*sj(i, j-1, k, 2) + w(i, j, k+1, ivz)*sk(i, j, k, 2) - w(i&
            &         , j, k-1, ivz)*sk(i, j, k-1, 2)
            wwz = w(i+1, j, k, ivz)*si(i, j, k, 3) - w(i-1, j, k, ivz)*si(i-&
            &         1, j, k, 3) + w(i, j+1, k, ivz)*sj(i, j, k, 3) - w(i, j-1, k, &
            &         ivz)*sj(i, j-1, k, 3) + w(i, j, k+1, ivz)*sk(i, j, k, 3) - w(i&
            &         , j, k-1, ivz)*sk(i, j, k-1, 3)
            ! compute the components of the stress tensor.
            ! the combination of the current scaling of the velocity
            ! gradients (2*vol) and the definition of the stress tensor,
            ! leads to the factor 1/(4*vol).
            fact = fourth/vol(i, j, k)
            if (turbprod .eq. strain) then
                sxx = two*fact*uux
                syy = two*fact*vvy
                szz = two*fact*wwz
                sxy = fact*(uuy+vvx)
                sxz = fact*(uuz+wwx)
                syz = fact*(vvz+wwy)
                ! compute 2/3 * divergence of velocity squared
                div2 = f23*(sxx+syy+szz)**2
                ! compute strain production term
                strainmag2 = two*(sxy**2+sxz**2+syz**2) + sxx**2 + syy**2 + &
                &           szz**2
                strainprod = two*strainmag2 - div2
                ss = sqrt(strainprod)
                call pushcontrol2b(0)
            else if (turbprod .eq. vorticity) then
                ! compute the three components of the vorticity vector.
                ! substract the part coming from the rotating frame.
                vortx = two*fact*(wwy-vvz) - two*omegax
                vorty = two*fact*(uuz-wwx) - two*omegay
                vortz = two*fact*(vvx-uuy) - two*omegaz
                ! compute the vorticity production term
                vortprod = vortx**2 + vorty**2 + vortz**2
                ! first take the square root of the production term to
                ! obtain the correct production term for spalart-allmaras.
                ! we do this to avoid if statements.
                ss = sqrt(vortprod)
                call pushcontrol2b(1)
            else
                call pushcontrol2b(2)
            end if
            ! compute the laminar kinematic viscosity, the inverse of
            ! wall distance squared, the ratio chi (ratio of nutilde
            ! and nu) and the functions fv1 and fv2. the latter corrects
            ! the production term near a viscous wall.
            nu = rlv(i, j, k)/w(i, j, k, irho)
            dist2inv = one/d2wall(i, j, k)**2
            chi = w(i, j, k, itu1)/nu
            chi2 = chi*chi
            chi3 = chi*chi2
            fv1 = chi3/(chi3+cv13)
            fv2 = one - chi/(one+chi*fv1)
            ! the function ft2, which is designed to keep a laminar
            ! solution laminar. when running in fully turbulent mode
            ! this function should be set to 0.0.
            if (useft2sa) then
                ft2 = rsact3*exp(-(rsact4*chi2))
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                ft2 = zero
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            ! correct the production term to account for the influence
            ! of the wall.
            sst = ss + w(i, j, k, itu1)*fv2*kar2inv*dist2inv
            ! add rotation term (userotationsa defined in inputparams.f90)
            if (userotationsa) then
                y1 = sqrt(two*strainmag2)
                if (zero .gt. y1) then
                    min1 = y1
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 0
                else
                    min1 = zero
                    myIntPtr = myIntPtr + 1
                    myIntStack(myIntPtr) = 1
                end if
                sst = sst + rsacrot*min1
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            else
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            end if
            if (sst .lt. xminn) then
                sst = xminn
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
                sst = sst
            end if
            ! compute the function fw. the argument rr is cut off at 10
            ! to avoid numerical problems. this is ok, because the
            ! asymptotical value of fw is then already reached.
            rr = w(i, j, k, itu1)*kar2inv*dist2inv/sst
            if (rr .gt. 10.0_realtype) then
                rr = 10.0_realtype
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
            else
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
                rr = rr
            end if
            gg = rr + rsacw2*(rr**6-rr)
            gg6 = gg**6
            termfw = ((one+cw36)/(gg6+cw36))**sixth
            fwsa = gg*termfw
            ! compute the source term; some terms are saved for the
            ! linearization. the source term is stored in dvt.
            if (approxsa) then
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 0
                term1 = zero
            else
                term1 = rsacb1*(one-ft2)*ss
                myIntPtr = myIntPtr + 1
                myIntStack(myIntPtr) = 1
            end if
            term2 = dist2inv*(kar2inv*rsacb1*((one-ft2)*fv2+ft2)-rsacw1*fwsa&
            &         )
            tempd9 = w(i, j, k, itu1)*scratchd(i, j, k, idvt)
            temp1 = w(i, j, k, itu1)
            term1d = tempd9
            term2d = temp1*tempd9
            wd(i, j, k, itu1) = wd(i, j, k, itu1) + (term1+term2*temp1)*&
            &         scratchd(i, j, k, idvt) + term2*tempd9
            scratchd(i, j, k, idvt) = 0.0_8
            tempd10 = dist2inv*kar2inv*rsacb1*term2d
            ft2d = (1.0_8-fv2)*tempd10
            fv2d = (one-ft2)*tempd10
            fwsad = -(dist2inv*rsacw1*term2d)
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .ne. 0) then
                ft2d = ft2d - ss*rsacb1*term1d
                ssd = ssd + rsacb1*(one-ft2)*term1d
            end if
            termfwd = gg*fwsad
            temp0 = (one+cw36)/(cw36+gg6)
            if (temp0 .le. 0.0_8 .and. (sixth .eq. 0.0_8 .or. sixth .ne. int&
            &           (sixth))) then
            gg6d = 0.0
        else
            gg6d = -(sixth*temp0**(sixth-1)*temp0*termfwd/(cw36+gg6))
        end if
        ggd = 6*gg**5*gg6d + termfw*fwsad
        rrd = (rsacw2*6*rr**5-rsacw2+1.0_8)*ggd
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) rrd = 0.0_8
        tempd8 = kar2inv*dist2inv*rrd/sst
        wd(i, j, k, itu1) = wd(i, j, k, itu1) + tempd8
        sstd = -(w(i, j, k, itu1)*tempd8/sst)
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) sstd = 0.0_8
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .ne. 0) then
            min1d = rsacrot*sstd
            branch = myIntStack(myIntPtr)
            myIntPtr = myIntPtr - 1
            if (branch .eq. 0) then
                y1d = min1d
            else
                y1d = 0.0_8
            end if
            if (.not.two*strainmag2 .eq. 0.0_8) strainmag2d = strainmag2d &
            &             + two*y1d/(2.0*sqrt(two*strainmag2))
        end if
        tempd7 = kar2inv*dist2inv*sstd
        ssd = ssd + sstd
        wd(i, j, k, itu1) = wd(i, j, k, itu1) + fv2*tempd7
        fv2d = fv2d + w(i, j, k, itu1)*tempd7
        branch = myIntStack(myIntPtr)
        myIntPtr = myIntPtr - 1
        if (branch .eq. 0) then
            chi2d = -(exp(-(rsact4*chi2))*rsact3*rsact4*ft2d)
        else
            chi2d = 0.0_8
        end if
        tempd4 = -(fv2d/(one+chi*fv1))
        tempd5 = -(chi*tempd4/(one+chi*fv1))
        fv1d = chi*tempd5
        tempd6 = fv1d/(cv13+chi3)
        chi3d = (1.0_8-chi3/(cv13+chi3))*tempd6
        chi2d = chi2d + chi*chi3d
        chid = chi2*chi3d + 2*chi*chi2d + fv1*tempd5 + tempd4
        wd(i, j, k, itu1) = wd(i, j, k, itu1) + chid/nu
        nud = -(w(i, j, k, itu1)*chid/nu**2)
        temp = w(i, j, k, irho)
        rlvd(i, j, k) = rlvd(i, j, k) + nud/temp
        wd(i, j, k, irho) = wd(i, j, k, irho) - rlv(i, j, k)*nud/temp**2
        call popcontrol2b(branch)
        if (branch .eq. 0) then
            if (strainprod .eq. 0.0_8) then
                strainprodd = 0.0
            else
                strainprodd = ssd/(2.0*sqrt(strainprod))
            end if
            strainmag2d = strainmag2d + two*strainprodd
            div2d = -strainprodd
            tempd = two*strainmag2d
            sxyd = 2*sxy*tempd
            sxzd = 2*sxz*tempd
            syzd = 2*syz*tempd
            tempd0 = f23*2*(sxx+syy+szz)*div2d
            sxxd = tempd0 + 2*sxx*strainmag2d
            syyd = tempd0 + 2*syy*strainmag2d
            szzd = tempd0 + 2*szz*strainmag2d
            vvzd = fact*syzd
            wwyd = fact*syzd
            uuzd = fact*sxzd
            wwxd = fact*sxzd
            uuyd = fact*sxyd
            vvxd = fact*sxyd
            wwzd = two*fact*szzd
            vvyd = two*fact*syyd
            uuxd = two*fact*sxxd
            strainmag2d = 0.0_8
            ssd = 0.0_8
        else
            if (branch .eq. 1) then
                if (vortprod .eq. 0.0_8) then
                    vortprodd = 0.0
                else
                    vortprodd = ssd/(2.0*sqrt(vortprod))
                end if
                vortxd = 2*vortx*vortprodd
                vortyd = 2*vorty*vortprodd
                vortzd = 2*vortz*vortprodd
                tempd1 = two*fact*vortzd
                vvxd = tempd1
                uuyd = -tempd1
                tempd2 = two*fact*vortyd
                uuzd = tempd2
                wwxd = -tempd2
                tempd3 = two*fact*vortxd
                wwyd = tempd3
                vvzd = -tempd3
                ssd = 0.0_8
            else
                wwxd = 0.0_8
                wwyd = 0.0_8
                vvxd = 0.0_8
                vvzd = 0.0_8
                uuyd = 0.0_8
                uuzd = 0.0_8
            end if
            wwzd = 0.0_8
            vvyd = 0.0_8
            uuxd = 0.0_8
        end if
        wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) + si(i, j, k, 3)*wwzd
        wd(i-1, j, k, ivz) = wd(i-1, j, k, ivz) - si(i-1, j, k, 3)*wwzd
        wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) + sj(i, j, k, 3)*wwzd
        wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) + sk(i, j, k, 3)*wwzd
        wd(i, j-1, k, ivz) = wd(i, j-1, k, ivz) - sj(i, j-1, k, 3)*wwzd
        wd(i, j, k-1, ivz) = wd(i, j, k-1, ivz) - sk(i, j, k-1, 3)*wwzd
        wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) + si(i, j, k, 2)*wwyd
        wd(i-1, j, k, ivz) = wd(i-1, j, k, ivz) - si(i-1, j, k, 2)*wwyd
        wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) + sj(i, j, k, 2)*wwyd
        wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) + sk(i, j, k, 2)*wwyd
        wd(i, j-1, k, ivz) = wd(i, j-1, k, ivz) - sj(i, j-1, k, 2)*wwyd
        wd(i, j, k-1, ivz) = wd(i, j, k-1, ivz) - sk(i, j, k-1, 2)*wwyd
        wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) + si(i, j, k, 1)*wwxd
        wd(i-1, j, k, ivz) = wd(i-1, j, k, ivz) - si(i-1, j, k, 1)*wwxd
        wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) + sj(i, j, k, 1)*wwxd
        wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) + sk(i, j, k, 1)*wwxd
        wd(i, j-1, k, ivz) = wd(i, j-1, k, ivz) - sj(i, j-1, k, 1)*wwxd
        wd(i, j, k-1, ivz) = wd(i, j, k-1, ivz) - sk(i, j, k-1, 1)*wwxd
        wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) + si(i, j, k, 3)*vvzd
        wd(i-1, j, k, ivy) = wd(i-1, j, k, ivy) - si(i-1, j, k, 3)*vvzd
        wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) + sj(i, j, k, 3)*vvzd
        wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) + sk(i, j, k, 3)*vvzd
        wd(i, j-1, k, ivy) = wd(i, j-1, k, ivy) - sj(i, j-1, k, 3)*vvzd
        wd(i, j, k-1, ivy) = wd(i, j, k-1, ivy) - sk(i, j, k-1, 3)*vvzd
        wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) + si(i, j, k, 2)*vvyd
        wd(i-1, j, k, ivy) = wd(i-1, j, k, ivy) - si(i-1, j, k, 2)*vvyd
        wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) + sj(i, j, k, 2)*vvyd
        wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) + sk(i, j, k, 2)*vvyd
        wd(i, j-1, k, ivy) = wd(i, j-1, k, ivy) - sj(i, j-1, k, 2)*vvyd
        wd(i, j, k-1, ivy) = wd(i, j, k-1, ivy) - sk(i, j, k-1, 2)*vvyd
        wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) + si(i, j, k, 1)*vvxd
        wd(i-1, j, k, ivy) = wd(i-1, j, k, ivy) - si(i-1, j, k, 1)*vvxd
        wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) + sj(i, j, k, 1)*vvxd
        wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) + sk(i, j, k, 1)*vvxd
        wd(i, j-1, k, ivy) = wd(i, j-1, k, ivy) - sj(i, j-1, k, 1)*vvxd
        wd(i, j, k-1, ivy) = wd(i, j, k-1, ivy) - sk(i, j, k-1, 1)*vvxd
        wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) + si(i, j, k, 3)*uuzd
        wd(i-1, j, k, ivx) = wd(i-1, j, k, ivx) - si(i-1, j, k, 3)*uuzd
        wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) + sj(i, j, k, 3)*uuzd
        wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) + sk(i, j, k, 3)*uuzd
        wd(i, j-1, k, ivx) = wd(i, j-1, k, ivx) - sj(i, j-1, k, 3)*uuzd
        wd(i, j, k-1, ivx) = wd(i, j, k-1, ivx) - sk(i, j, k-1, 3)*uuzd
        wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) + si(i, j, k, 2)*uuyd
        wd(i-1, j, k, ivx) = wd(i-1, j, k, ivx) - si(i-1, j, k, 2)*uuyd
        wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) + sj(i, j, k, 2)*uuyd
        wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) + sk(i, j, k, 2)*uuyd
        wd(i, j-1, k, ivx) = wd(i, j-1, k, ivx) - sj(i, j-1, k, 2)*uuyd
        wd(i, j, k-1, ivx) = wd(i, j, k-1, ivx) - sk(i, j, k-1, 2)*uuyd
        wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) + si(i, j, k, 1)*uuxd
        wd(i-1, j, k, ivx) = wd(i-1, j, k, ivx) - si(i-1, j, k, 1)*uuxd
        wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) + sj(i, j, k, 1)*uuxd
        wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) + sk(i, j, k, 1)*uuxd
        wd(i, j-1, k, ivx) = wd(i, j-1, k, ivx) - sj(i, j-1, k, 1)*uuxd
        wd(i, j, k-1, ivx) = wd(i, j, k-1, ivx) - sk(i, j, k-1, 1)*uuxd
    end do
end if
end subroutine sasource_fast_b

!  differentiation of timestep_block in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: *p *w *radi *radj *radk
!   with respect to varying inputs: *p *w *radi *radj *radk
!   rw status of diff variables: *p:incr *w:incr *radi:in-out *radj:in-out
!                *radk:in-out
!   plus diff mem management of: p:in w:in radi:in radj:in radk:in
  subroutine timestep_block_fast_b(onlyradii)
!
!       timestep computes the time step, or more precisely the time
!       step divided by the volume per unit cfl, in the owned cells.
!       however, for the artificial dissipation schemes, the spectral
!       radii in the halo's are needed. therefore the loop is taken
!       over the the first level of halo cells. the spectral radii are
!       stored and possibly modified for high aspect ratio cells.
!
    use constants
!     use blockpointers, only : ie, je, ke, il, jl, kl, w, wd, p, pd, &
! &   rlv, rlvd, rev, revd, radi, radid, radj, radjd, radk, radkd, si, sj,&
! &   sk, sfacei, sfacej, sfacek, dtl, gamma, vol, addgridvelocities, &
! &   sectionid
    use blockpointers, only : addgridvelocities, sectionid
    use flowvarrefstate, only : timeref, eddymodel, gammainf, pinfcorr&
&   , viscous, rhoinf
    use inputdiscretization, only : adis, dirscaling, &
&   radiineededcoarse, radiineededfine, precond
    use inputphysics, only : equationmode
    use iteration, only : groundlevel, currentlevel
    use section, only : sections
    use inputtimespectral, only : ntimeintervalsspectral
    use utils_fast_b, only : terminate
    implicit none
! the rest of this file can be skipped if only the spectral
! radii need to be computed.
!
!      subroutine argument.
!
    logical, intent(in) :: onlyradii
!
!      local parameters.
!
    real(kind=realtype), parameter :: b=2.0_realtype
!
!      local variables.
!
    integer(kind=inttype) :: i, j, k, ii
    real(kind=realtype) :: plim, rlim, clim2
    real(kind=realtype) :: uux, uuy, uuz, cc2, qsi, qsj, qsk, sx, sy, sz&
&   , rmu
    real(kind=realtype) :: uuxd, uuyd, uuzd, cc2d, qsid, qsjd, qskd
    real(kind=realtype) :: ri, rj, rk, rij, rjk, rki
    real(kind=realtype) :: rid, rjd, rkd, rijd, rjkd, rkid
    real(kind=realtype) :: vsi, vsj, vsk, rfl, dpi, dpj, dpk
    real(kind=realtype) :: sface, tmp
    logical :: radiineeded, doscaling
    intrinsic mod
    intrinsic max
    intrinsic abs
    intrinsic sqrt
    integer :: branch
    real(kind=realtype) :: temp2
    real(kind=realtype) :: temp1
    real(kind=realtype) :: temp0
    real(kind=realtype) :: abs1d
    real(kind=realtype) :: abs0d
    real(kind=realtype) :: tempd
    real(kind=realtype) :: tempd2
    real(kind=realtype) :: tempd1
    real(kind=realtype) :: tempd0
    real(kind=realtype) :: abs2
    real(kind=realtype) :: abs2d
    real(kind=realtype) :: abs1
    real(kind=realtype) :: abs0
    real(kind=realtype) :: temp
! determine whether or not the spectral radii are needed for the
! flux computation.
    radiineeded = radiineededcoarse
    if (currentlevel .le. groundlevel) radiineeded = radiineededfine
! return immediately if only the spectral radii must be computed
! and these are not needed for the flux computation.
    if (.not.(onlyradii .and. (.not.radiineeded))) then
! set the value of plim. to be fully consistent this must have
! the dimension of a pressure. therefore a fraction of pinfcorr
! is used. idem for rlim; compute clim2 as well.
      clim2 = 0.000001_realtype*gammainf*pinfcorr/rhoinf
      doscaling = dirscaling .and. currentlevel .le. groundlevel
! initialize sface to zero. this value will be used if the
! block is not moving.
      sface = zero
!
!           inviscid contribution, depending on the preconditioner.
!           compute the cell centered values of the spectral radii.
!
      select case  (precond)
      case (noprecond)
        do ii=0,ie*je*ke-1
          i = mod(ii, ie) + 1
          j = mod(ii/ie, je) + 1
          k = ii/(ie*je) + 1
! compute the velocities and speed of sound squared.
          uux = w(i, j, k, ivx)
          uuy = w(i, j, k, ivy)
          uuz = w(i, j, k, ivz)
          cc2 = gamma(i, j, k)*p(i, j, k)/w(i, j, k, irho)
          if (cc2 .lt. clim2) then
            cc2 = clim2
myIntPtr = myIntPtr + 1
 myIntStack(myIntPtr) = 0
          else
myIntPtr = myIntPtr + 1
 myIntStack(myIntPtr) = 1
            cc2 = cc2
          end if
! set the dot product of the grid velocity and the
! normal in i-direction for a moving face. to avoid
! a number of multiplications by 0.5 simply the sum
! is taken.
          if (addgridvelocities) sface = sfacei(i-1, j, k) + sfacei(i, j&
&             , k)
! spectral radius in i-direction.
          sx = si(i-1, j, k, 1) + si(i, j, k, 1)
          sy = si(i-1, j, k, 2) + si(i, j, k, 2)
          sz = si(i-1, j, k, 3) + si(i, j, k, 3)
          qsi = uux*sx + uuy*sy + uuz*sz - sface
          if (qsi .ge. 0.) then
            abs0 = qsi
myIntPtr = myIntPtr + 1
 myIntStack(myIntPtr) = 0
          else
            abs0 = -qsi
myIntPtr = myIntPtr + 1
 myIntStack(myIntPtr) = 1
          end if
          ri = half*(abs0+sqrt(cc2*(sx**2+sy**2+sz**2)))
! the grid velocity in j-direction.
          if (addgridvelocities) sface = sfacej(i, j-1, k) + sfacej(i, j&
&             , k)
! spectral radius in j-direction.
          sx = sj(i, j-1, k, 1) + sj(i, j, k, 1)
          sy = sj(i, j-1, k, 2) + sj(i, j, k, 2)
          sz = sj(i, j-1, k, 3) + sj(i, j, k, 3)
          qsj = uux*sx + uuy*sy + uuz*sz - sface
          if (qsj .ge. 0.) then
            abs1 = qsj
myIntPtr = myIntPtr + 1
 myIntStack(myIntPtr) = 0
          else
            abs1 = -qsj
myIntPtr = myIntPtr + 1
 myIntStack(myIntPtr) = 1
          end if
          rj = half*(abs1+sqrt(cc2*(sx**2+sy**2+sz**2)))
! the grid velocity in k-direction.
          if (addgridvelocities) sface = sfacek(i, j, k-1) + sfacek(i, j&
&             , k)
! spectral radius in k-direction.
          sx = sk(i, j, k-1, 1) + sk(i, j, k, 1)
          sy = sk(i, j, k-1, 2) + sk(i, j, k, 2)
          sz = sk(i, j, k-1, 3) + sk(i, j, k, 3)
          qsk = uux*sx + uuy*sy + uuz*sz - sface
          if (qsk .ge. 0.) then
            abs2 = qsk
myIntPtr = myIntPtr + 1
 myIntStack(myIntPtr) = 0
          else
            abs2 = -qsk
myIntPtr = myIntPtr + 1
 myIntStack(myIntPtr) = 1
          end if
          rk = half*(abs2+sqrt(cc2*(sx**2+sy**2+sz**2)))
! compute the inviscid contribution to the time step.
!
!           adapt the spectral radii if directional scaling must be
!           applied.
!
          if (doscaling) then
            if (ri .lt. eps) then
              ri = eps
myIntPtr = myIntPtr + 1
 myIntStack(myIntPtr) = 0
            else
myIntPtr = myIntPtr + 1
 myIntStack(myIntPtr) = 1
              ri = ri
            end if
            if (rj .lt. eps) then
              rj = eps
myIntPtr = myIntPtr + 1
 myIntStack(myIntPtr) = 0
            else
myIntPtr = myIntPtr + 1
 myIntStack(myIntPtr) = 1
              rj = rj
            end if
            if (rk .lt. eps) then
              rk = eps
myIntPtr = myIntPtr + 1
 myIntStack(myIntPtr) = 0
            else
myIntPtr = myIntPtr + 1
 myIntStack(myIntPtr) = 1
              rk = rk
            end if
! compute the scaling in the three coordinate
! directions.
            rij = (ri/rj)**adis
            rjk = (rj/rk)**adis
            rki = (rk/ri)**adis
! create the scaled versions of the aspect ratios.
! note that the multiplication is done with radi, radj
! and radk, such that the influence of the clipping
! is negligible.
            rkd = (one+one/rki+rjk)*radkd(i, j, k)
            rkid = -(rk*one*radkd(i, j, k)/rki**2)
            rjkd = rk*radkd(i, j, k)
            radkd(i, j, k) = 0.0_8
            rjd = (one+one/rjk+rij)*radjd(i, j, k)
            rjkd = rjkd - rj*one*radjd(i, j, k)/rjk**2
            rijd = rj*radjd(i, j, k)
            radjd(i, j, k) = 0.0_8
            rijd = rijd - ri*one*radid(i, j, k)/rij**2
            rkid = rkid + ri*radid(i, j, k)
            if (rk/ri .le. 0.0_8 .and. (adis .eq. 0.0_8 .or. adis .ne. &
&               int(adis))) then
              tempd0 = 0.0
            else
              tempd0 = adis*(rk/ri)**(adis-1)*rkid/ri
            end if
            if (rj/rk .le. 0.0_8 .and. (adis .eq. 0.0_8 .or. adis .ne. &
&               int(adis))) then
              tempd2 = 0.0
            else
              tempd2 = adis*(rj/rk)**(adis-1)*rjkd/rk
            end if
            rkd = rkd + tempd0 - rj*tempd2/rk
            if (ri/rj .le. 0.0_8 .and. (adis .eq. 0.0_8 .or. adis .ne. &
&               int(adis))) then
              tempd1 = 0.0
            else
              tempd1 = adis*(ri/rj)**(adis-1)*rijd/rj
            end if
            rid = tempd1 - rk*tempd0/ri + (one+one/rij+rki)*radid(i, j, &
&             k)
            radid(i, j, k) = 0.0_8
            rjd = rjd + tempd2 - ri*tempd1/rj
branch = myIntStack(myIntPtr)
 myIntPtr = myIntPtr - 1
            if (branch .eq. 0) rkd = 0.0_8
branch = myIntStack(myIntPtr)
 myIntPtr = myIntPtr - 1
            if (branch .eq. 0) rjd = 0.0_8
branch = myIntStack(myIntPtr)
 myIntPtr = myIntPtr - 1
            if (branch .eq. 0) rid = 0.0_8
          else
            rkd = radkd(i, j, k)
            radkd(i, j, k) = 0.0_8
            rjd = radjd(i, j, k)
            radjd(i, j, k) = 0.0_8
            rid = radid(i, j, k)
            radid(i, j, k) = 0.0_8
          end if
          temp2 = sx**2 + sy**2 + sz**2
          abs2d = half*rkd
          if (temp2*cc2 .eq. 0.0_8) then
            cc2d = 0.0
          else
            cc2d = half*temp2*rkd/(2.0*sqrt(temp2*cc2))
          end if
branch = myIntStack(myIntPtr)
 myIntPtr = myIntPtr - 1
          if (branch .eq. 0) then
            qskd = abs2d
          else
            qskd = -abs2d
          end if
          uuxd = sx*qskd
          uuyd = sy*qskd
          uuzd = sz*qskd
          sx = sj(i, j-1, k, 1) + sj(i, j, k, 1)
          sy = sj(i, j-1, k, 2) + sj(i, j, k, 2)
          sz = sj(i, j-1, k, 3) + sj(i, j, k, 3)
          temp1 = sx**2 + sy**2 + sz**2
          abs1d = half*rjd
          if (.not.temp1*cc2 .eq. 0.0_8) cc2d = cc2d + half*temp1*rjd/(&
&             2.0*sqrt(temp1*cc2))
branch = myIntStack(myIntPtr)
 myIntPtr = myIntPtr - 1
          if (branch .eq. 0) then
            qsjd = abs1d
          else
            qsjd = -abs1d
          end if
          uuxd = uuxd + sx*qsjd
          uuyd = uuyd + sy*qsjd
          uuzd = uuzd + sz*qsjd
          sx = si(i-1, j, k, 1) + si(i, j, k, 1)
          sy = si(i-1, j, k, 2) + si(i, j, k, 2)
          sz = si(i-1, j, k, 3) + si(i, j, k, 3)
          temp0 = sx**2 + sy**2 + sz**2
          abs0d = half*rid
          if (.not.temp0*cc2 .eq. 0.0_8) cc2d = cc2d + half*temp0*rid/(&
&             2.0*sqrt(temp0*cc2))
branch = myIntStack(myIntPtr)
 myIntPtr = myIntPtr - 1
          if (branch .eq. 0) then
            qsid = abs0d
          else
            qsid = -abs0d
          end if
          uuxd = uuxd + sx*qsid
          uuyd = uuyd + sy*qsid
          uuzd = uuzd + sz*qsid
branch = myIntStack(myIntPtr)
 myIntPtr = myIntPtr - 1
          if (branch .eq. 0) cc2d = 0.0_8
          temp = w(i, j, k, irho)
          tempd = gamma(i, j, k)*cc2d/temp
          pd(i, j, k) = pd(i, j, k) + tempd
          wd(i, j, k, irho) = wd(i, j, k, irho) - p(i, j, k)*tempd/temp
          wd(i, j, k, ivz) = wd(i, j, k, ivz) + uuzd
          wd(i, j, k, ivy) = wd(i, j, k, ivy) + uuyd
          wd(i, j, k, ivx) = wd(i, j, k, ivx) + uuxd
        end do
      end select
    end if
  end subroutine timestep_block_fast_b

end module blockette_state_b
