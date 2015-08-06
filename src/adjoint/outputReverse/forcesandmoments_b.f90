!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
!  differentiation of forcesandmoments in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: gammainf pinf pref *w *x *(*bcdata.f)
!                lengthref machcoef pointref *xx *rev0 *rev1 *rev2
!                *rev3 *pp0 *pp1 *pp2 *pp3 *rlv0 *rlv1 *rlv2 *rlv3
!                *ssi *ww0 *ww1 *ww2 *ww3 cfp cfv cmp cmv cavitation
!                sepsensor
!   with respect to varying inputs: gammainf pinf pref *rev *p
!                *w *rlv *x *si *sj *sk *(*viscsubface.tau) *(*bcdata.f)
!                veldirfreestream lengthref machcoef pointref *xx
!                *rev0 *rev1 *rev2 *rev3 *pp0 *pp1 *pp2 *pp3 *rlv0
!                *rlv1 *rlv2 *rlv3 *ssi *ww0 *ww1 *ww2 *ww3
!   plus diff mem management of: rev:in p:in w:in rlv:in x:in si:in
!                sj:in sk:in viscsubface:in *viscsubface.tau:in
!                bcdata:in *bcdata.f:in *bcdata.dualarea:in xx:in
!                rev0:in rev1:in rev2:in rev3:in pp0:in pp1:in
!                pp2:in pp3:in rlv0:in rlv1:in rlv2:in rlv3:in
!                ssi:in ww0:in ww1:in ww2:in ww3:in
!
!      ******************************************************************
!      *                                                                *
!      * file:          forcesandmoments.f90                            *
!      * author:        edwin van der weide                             *
!      * starting date: 04-01-2003                                      *
!      * last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine forcesandmoments_b(cfp, cfpd, cfv, cfvd, cmp, cmpd, cmv, cmvd&
& , yplusmax, sepsensor, sepsensord, cavitation, cavitationd)
!
!      ******************************************************************
!      *                                                                *
!      * forcesandmoments computes the contribution of the block        *
!      * given by the pointers in blockpointers to the force and        *
!      * moment coefficients of the geometry. a distinction is made     *
!      * between the inviscid and viscous parts. in case the maximum    *
!      * yplus value must be monitored (only possible for rans), this   *
!      * value is also computed. the separation sensor and the cavita-  *
!      * tion sensor is also computed                                   *
!      * here.                                                          *
!      ******************************************************************
!
  use blockpointers
  use bctypes
  use flowvarrefstate
  use inputphysics
  use bcroutines_b
  use costfunctions
  implicit none
!
!      subroutine arguments
!
  real(kind=realtype), dimension(3) :: cfp, cfv
  real(kind=realtype), dimension(3) :: cfpd, cfvd
  real(kind=realtype), dimension(3) :: cmp, cmv
  real(kind=realtype), dimension(3) :: cmpd, cmvd
  real(kind=realtype) :: yplusmax, sepsensor, cavitation
  real(kind=realtype) :: sepsensord, cavitationd
!
!      local variables.
!
  integer(kind=inttype) :: nn, i, j, ii
  real(kind=realtype) :: pm1, fx, fy, fz, fn, sigma
  real(kind=realtype) :: pm1d, fxd, fyd, fzd
  real(kind=realtype) :: xc, yc, zc, qf(3)
  real(kind=realtype) :: xcd, ycd, zcd
  real(kind=realtype) :: fact, rho, mul, yplus, dwall
  real(kind=realtype) :: factd
  real(kind=realtype) :: scaledim, v(3), sensor, sensor1, cp, tmp, &
& plocal
  real(kind=realtype) :: scaledimd, vd(3), sensord, sensor1d, cpd, tmpd&
& , plocald
  real(kind=realtype) :: tauxx, tauyy, tauzz
  real(kind=realtype) :: tauxxd, tauyyd, tauzzd
  real(kind=realtype) :: tauxy, tauxz, tauyz
  real(kind=realtype) :: tauxyd, tauxzd, tauyzd
  real(kind=realtype), dimension(3) :: refpoint
  real(kind=realtype), dimension(3) :: refpointd
  real(kind=realtype) :: mx, my, mz, qa
  real(kind=realtype) :: mxd, myd, mzd, qad
  logical :: viscoussubface
  intrinsic mod
  intrinsic sqrt
  intrinsic exp
  real(kind=realtype), dimension(3) :: tmp0
  integer :: ad_to
  integer :: ad_from
  integer :: ad_to0
  integer :: ad_from0
  integer :: ad_to1
  integer :: branch
  real(kind=realtype) :: temp3
  real(kind=realtype) :: tempd14
  real(kind=realtype) :: temp2
  real(kind=realtype) :: tempd13
  real(kind=realtype) :: temp1
  real(kind=realtype) :: tempd12
  real(kind=realtype) :: temp0
  real(kind=realtype) :: tempd11
  real(kind=realtype) :: tempd10(3)
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
  real(kind=realtype) :: tmpd0(3)
  integer :: ii1
  real(kind=realtype) :: temp
  real(kind=realtype) :: temp9
  real(kind=realtype) :: temp8
  real(kind=realtype) :: temp7
  real(kind=realtype) :: temp6
  real(kind=realtype) :: temp5
  real(kind=realtype) :: tempd16
  real(kind=realtype) :: temp4
  real(kind=realtype) :: tempd15
!
!      ******************************************************************
!      *                                                                *
!      * begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
! set the actual scaling factor such that actual forces are computed
  scaledim = pref/pinf
! determine the reference point for the moment computation in
! meters.
  refpoint(1) = lref*pointref(1)
  refpoint(2) = lref*pointref(2)
  refpoint(3) = lref*pointref(3)
! initialize the force and moment coefficients to 0 as well as
! yplusmax.
  cfp(1) = zero
  cfp(2) = zero
  cfp(3) = zero
  cfv(1) = zero
  cfv(2) = zero
  cfv(3) = zero
  cmp(1) = zero
  cmp(2) = zero
  cmp(3) = zero
  cmv(1) = zero
  cmv(2) = zero
  cmv(3) = zero
! loop over the boundary subfaces of this block.
bocos:do nn=1,nbocos
!
!        ****************************************************************
!        *                                                              *
!        * integrate the inviscid contribution over the solid walls,    *
!        * either inviscid or viscous. the integration is done with     *
!        * cp. for closed contours this is equal to the integration     *
!        * of p; for open contours this is not the case anymore.        *
!        * question is whether a force for an open contour is           *
!        * meaningful anyway.                                           *
!        *                                                              *
!        ****************************************************************
!
    if ((bctype(nn) .eq. eulerwall .or. bctype(nn) .eq. nswalladiabatic)&
&       .or. bctype(nn) .eq. nswallisothermal) then
! subface is a wall. check if it is a viscous wall.
      viscoussubface = .true.
      if (bctype(nn) .eq. eulerwall) viscoussubface = .false.
! set a bunch of pointers depending on the face id to make
! a generic treatment possible. the routine setbcpointers
! is not used, because quite a few other ones are needed.
      call pushreal8array(ssi, size(ssi, 1)*size(ssi, 2)*size(ssi, 3))
      call pushreal8array(pp2, size(pp2, 1)*size(pp2, 2))
      call pushreal8array(pp1, size(pp1, 1)*size(pp1, 2))
      call setbcpointers(nn, .true.)
      select case  (bcfaceid(nn)) 
      case (imin) 
        call pushreal8(fact)
        fact = -one
        call pushcontrol3b(1)
      case (imax) 
        call pushreal8(fact)
        fact = one
        call pushcontrol3b(2)
      case (jmin) 
        call pushreal8(fact)
        fact = -one
        call pushcontrol3b(3)
      case (jmax) 
        call pushreal8(fact)
        fact = one
        call pushcontrol3b(4)
      case (kmin) 
        call pushreal8(fact)
        fact = -one
        call pushcontrol3b(5)
      case (kmax) 
        call pushreal8(fact)
        fact = one
        call pushcontrol3b(6)
      case default
        call pushcontrol3b(0)
      end select
! loop over the quadrilateral faces of the subface. note that
! the nodal range of bcdata must be used and not the cell
! range, because the latter may include the halo's in i and
! j-direction. the offset +1 is there, because inbeg and jnbeg
! refer to nodal ranges and not to cell ranges. the loop
! (without the ad stuff) would look like:
!
! do j=(bcdata(nn)%jnbeg+1),bcdata(nn)%jnend
!    do i=(bcdata(nn)%inbeg+1),bcdata(nn)%inend
      call pushreal8array(bcdata(nn)%dualarea, size(bcdata(nn)%dualarea&
&                   , 1)*size(bcdata(nn)%dualarea, 2))
      bcdata(nn)%dualarea = zero
      call pushreal8array(bcdata(nn)%f, size(bcdata(nn)%f, 1)*size(&
&                   bcdata(nn)%f, 2)*size(bcdata(nn)%f, 3))
      bcdata(nn)%f = zero
      call pushreal8(qa)
      call pushinteger4(i)
      call pushinteger4(j)
      call pushreal8(xc)
      call pushreal8(plocal)
      call pushreal8(yc)
      call pushreal8(tmp)
      call pushreal8(sensor)
      call pushreal8(zc)
      call pushreal8(sensor1)
      call pushreal8array(v, 3)
      do ii1=1,size(bcdata)
        call pushreal8array(bcdata(ii1)%f, size(bcdata(ii1)%f, 1)*size(&
&                     bcdata(ii1)%f, 2)*size(bcdata(ii1)%f, 3))
      end do
      do ii1=1,size(bcdata)
        call pushreal8array(bcdata(ii1)%dualarea, size(bcdata(ii1)%&
&                     dualarea, 1)*size(bcdata(ii1)%dualarea, 2))
      end do
      call pushreal8array(xx, size(xx, 1)*size(xx, 2)*size(xx, 3))
      call pushreal8array(pp1, size(pp1, 1)*size(pp1, 2))
      call pushreal8array(pp2, size(pp2, 1)*size(pp2, 2))
      call pushreal8array(ssi, size(ssi, 1)*size(ssi, 2)*size(ssi, 3))
      call pushreal8array(ww2, size(ww2, 1)*size(ww2, 2)*size(ww2, 3))
      do ii=0,(bcdata(nn)%jnend-bcdata(nn)%jnbeg)*(bcdata(nn)%inend-&
&         bcdata(nn)%inbeg)-1
        i = mod(ii, bcdata(nn)%inend - bcdata(nn)%inbeg) + bcdata(nn)%&
&         inbeg + 1
        j = ii/(bcdata(nn)%inend-bcdata(nn)%inbeg) + bcdata(nn)%jnbeg + &
&         1
! compute the average pressure minus 1 and the coordinates
! of the centroid of the face relative from from the
! moment reference point. due to the usage of pointers for
! the coordinates, whose original array starts at 0, an
! offset of 1 must be used. the pressure is multipled by
! fact to account for the possibility of an inward or
! outward pointing normal.
        pm1 = fact*(half*(pp2(i, j)+pp1(i, j))-pinf)*scaledim
        xc = fourth*(xx(i, j, 1)+xx(i+1, j, 1)+xx(i, j+1, 1)+xx(i+1, j+1&
&         , 1)) - refpoint(1)
        yc = fourth*(xx(i, j, 2)+xx(i+1, j, 2)+xx(i, j+1, 2)+xx(i+1, j+1&
&         , 2)) - refpoint(2)
        zc = fourth*(xx(i, j, 3)+xx(i+1, j, 3)+xx(i, j+1, 3)+xx(i+1, j+1&
&         , 3)) - refpoint(3)
! compute the force components.
        fx = pm1*ssi(i, j, 1)
        fy = pm1*ssi(i, j, 2)
        fz = pm1*ssi(i, j, 3)
! update the inviscid force and moment coefficients.
        cfp(1) = cfp(1) + fx
        cfp(2) = cfp(2) + fy
        cfp(3) = cfp(3) + fz
        mx = yc*fz - zc*fy
        my = zc*fx - xc*fz
        mz = xc*fy - yc*fx
        cmp(1) = cmp(1) + mx
        cmp(2) = cmp(2) + my
        cmp(3) = cmp(3) + mz
! divide by 4 so we can scatter
        fx = fourth*fx
        fy = fourth*fy
        fz = fourth*fz
! scatter 1/4 of the force to each of the nodes:
        bcdata(nn)%f(i-1, j-1, 1) = bcdata(nn)%f(i-1, j-1, 1) + fx
        bcdata(nn)%f(i, j-1, 1) = bcdata(nn)%f(i, j-1, 1) + fx
        bcdata(nn)%f(i-1, j, 1) = bcdata(nn)%f(i-1, j, 1) + fx
        bcdata(nn)%f(i, j, 1) = bcdata(nn)%f(i, j, 1) + fx
        bcdata(nn)%f(i-1, j-1, 2) = bcdata(nn)%f(i-1, j-1, 2) + fy
        bcdata(nn)%f(i, j-1, 2) = bcdata(nn)%f(i, j-1, 2) + fy
        bcdata(nn)%f(i-1, j, 2) = bcdata(nn)%f(i-1, j, 2) + fy
        bcdata(nn)%f(i, j, 2) = bcdata(nn)%f(i, j, 2) + fy
        bcdata(nn)%f(i-1, j-1, 3) = bcdata(nn)%f(i-1, j-1, 3) + fz
        bcdata(nn)%f(i, j-1, 3) = bcdata(nn)%f(i, j-1, 3) + fz
        bcdata(nn)%f(i-1, j, 3) = bcdata(nn)%f(i-1, j, 3) + fz
        bcdata(nn)%f(i, j, 3) = bcdata(nn)%f(i, j, 3) + fz
! scatter a quarter of the area to each node:
        qa = fourth*sqrt(ssi(i, j, 1)**2+ssi(i, j, 2)**2+ssi(i, j, 3)**2&
&         )
        bcdata(nn)%dualarea(i-1, j-1) = bcdata(nn)%dualarea(i-1, j-1) + &
&         qa
        bcdata(nn)%dualarea(i, j-1) = bcdata(nn)%dualarea(i, j-1) + qa
        bcdata(nn)%dualarea(i-1, j) = bcdata(nn)%dualarea(i-1, j) + qa
        bcdata(nn)%dualarea(i, j) = bcdata(nn)%dualarea(i, j) + qa
! get normalized surface velocity:
        v(1) = ww2(i, j, ivx)
        v(2) = ww2(i, j, ivy)
        v(3) = ww2(i, j, ivz)
        v = v/(sqrt(v(1)**2+v(2)**2+v(3)**2)+1e-16)
! dot product with free stream
        sensor = -(v(1)*veldirfreestream(1)+v(2)*veldirfreestream(2)+v(3&
&         )*veldirfreestream(3))
!now run through a smooth heaviside function:
        sensor = one/(one+exp(-(2*sepsensorsharpness*(sensor-&
&         sepsensoroffset))))
! and integrate over the area of this cell and save:
        sensor = sensor*four*qa
        sepsensor = sepsensor + sensor
        plocal = pp2(i, j)
        tmp = two/(gammainf*pinf*machcoef*machcoef)
        cp = tmp*(plocal-pinf)
        sigma = 1.4
        sensor1 = -cp - sigma
        sensor1 = one/(one+exp(-(2*10*sensor1)))
        sensor1 = sensor1*four*qa
        cavitation = cavitation + sensor1
      end do
!
!          **************************************************************
!          *                                                            *
!          * integration of the viscous forces.                         *
!          * only for viscous boundaries.                               *
!          *                                                            *
!          **************************************************************
!
      if (viscoussubface) then
! replace norm with bcdata norm - peter lyu
!norm => bcdata(nn)%norm
! loop over the quadrilateral faces of the subface and
! compute the viscous contribution to the force and
! moment and update the maximum value of y+.
        do ii=0,(bcdata(nn)%jnend-bcdata(nn)%jnbeg)*(bcdata(nn)%inend-&
&           bcdata(nn)%inbeg)-1
          call pushinteger4(i)
          i = mod(ii, bcdata(nn)%inend - bcdata(nn)%inbeg) + bcdata(nn)%&
&           inbeg + 1
          call pushinteger4(j)
          j = ii/(bcdata(nn)%inend-bcdata(nn)%inbeg) + bcdata(nn)%jnbeg &
&           + 1
! store the viscous stress tensor a bit easier.
          tauxx = viscsubface(nn)%tau(i, j, 1)
          tauyy = viscsubface(nn)%tau(i, j, 2)
          tauzz = viscsubface(nn)%tau(i, j, 3)
          tauxy = viscsubface(nn)%tau(i, j, 4)
          tauxz = viscsubface(nn)%tau(i, j, 5)
          tauyz = viscsubface(nn)%tau(i, j, 6)
! compute the viscous force on the face. a minus sign
! is now present, due to the definition of this force.
          fx = -(fact*(tauxx*ssi(i, j, 1)+tauxy*ssi(i, j, 2)+tauxz*ssi(i&
&           , j, 3))*scaledim)
          fy = -(fact*(tauxy*ssi(i, j, 1)+tauyy*ssi(i, j, 2)+tauyz*ssi(i&
&           , j, 3))*scaledim)
          fz = -(fact*(tauxz*ssi(i, j, 1)+tauyz*ssi(i, j, 2)+tauzz*ssi(i&
&           , j, 3))*scaledim)
! compute the coordinates of the centroid of the face
! relative from the moment reference point. due to the
! usage of pointers for xx and offset of 1 is present,
! because x originally starts at 0.
          call pushreal8(xc)
          xc = fourth*(xx(i, j, 1)+xx(i+1, j, 1)+xx(i, j+1, 1)+xx(i+1, j&
&           +1, 1)) - refpoint(1)
          call pushreal8(yc)
          yc = fourth*(xx(i, j, 2)+xx(i+1, j, 2)+xx(i, j+1, 2)+xx(i+1, j&
&           +1, 2)) - refpoint(2)
          call pushreal8(zc)
          zc = fourth*(xx(i, j, 3)+xx(i+1, j, 3)+xx(i, j+1, 3)+xx(i+1, j&
&           +1, 3)) - refpoint(3)
! update the viscous force and moment coefficients.
          cfv(1) = cfv(1) + fx
          cfv(2) = cfv(2) + fy
          cfv(3) = cfv(3) + fz
          mx = yc*fz - zc*fy
          my = zc*fx - xc*fz
          mz = xc*fy - yc*fx
          cmv(1) = cmv(1) + mx
          cmv(2) = cmv(2) + my
          cmv(3) = cmv(3) + mz
! divide by 4 so we can scatter
          call pushreal8(fx)
          fx = fourth*fx
          call pushreal8(fy)
          fy = fourth*fy
          call pushreal8(fz)
          fz = fourth*fz
! scatter 1/4 of the force to each of the nodes:
          call pushreal8(bcdata(nn)%f(i-1, j-1, 1))
          bcdata(nn)%f(i-1, j-1, 1) = bcdata(nn)%f(i-1, j-1, 1) + fx
          call pushreal8(bcdata(nn)%f(i, j-1, 1))
          bcdata(nn)%f(i, j-1, 1) = bcdata(nn)%f(i, j-1, 1) + fx
          call pushreal8(bcdata(nn)%f(i-1, j, 1))
          bcdata(nn)%f(i-1, j, 1) = bcdata(nn)%f(i-1, j, 1) + fx
          call pushreal8(bcdata(nn)%f(i, j, 1))
          bcdata(nn)%f(i, j, 1) = bcdata(nn)%f(i, j, 1) + fx
          call pushreal8(bcdata(nn)%f(i-1, j-1, 2))
          bcdata(nn)%f(i-1, j-1, 2) = bcdata(nn)%f(i-1, j-1, 2) + fy
          call pushreal8(bcdata(nn)%f(i, j-1, 2))
          bcdata(nn)%f(i, j-1, 2) = bcdata(nn)%f(i, j-1, 2) + fy
          call pushreal8(bcdata(nn)%f(i-1, j, 2))
          bcdata(nn)%f(i-1, j, 2) = bcdata(nn)%f(i-1, j, 2) + fy
          call pushreal8(bcdata(nn)%f(i, j, 2))
          bcdata(nn)%f(i, j, 2) = bcdata(nn)%f(i, j, 2) + fy
          call pushreal8(bcdata(nn)%f(i-1, j-1, 3))
          bcdata(nn)%f(i-1, j-1, 3) = bcdata(nn)%f(i-1, j-1, 3) + fz
          call pushreal8(bcdata(nn)%f(i, j-1, 3))
          bcdata(nn)%f(i, j-1, 3) = bcdata(nn)%f(i, j-1, 3) + fz
          call pushreal8(bcdata(nn)%f(i-1, j, 3))
          bcdata(nn)%f(i-1, j, 3) = bcdata(nn)%f(i-1, j, 3) + fz
          call pushreal8(bcdata(nn)%f(i, j, 3))
          bcdata(nn)%f(i, j, 3) = bcdata(nn)%f(i, j, 3) + fz
! compute the tangential component of the stress tensor,
! which is needed to monitor y+. the result is stored
! in fx, fy, fz, although it is not really a force.
! as later on only the magnitude of the tangential
! component is important, there is no need to take the
! sign into account (it should be a minus sign).
        end do
        call pushinteger4(ii - 1)
        call pushcontrol1b(0)
      else
        call pushcontrol1b(1)
      end if
! compute the local value of y+. due to the usage
! of pointers there is on offset of -1 in dd2wall..
! if forces are tractions we have to divide by the dual area:
      if (forcesastractions) then
        ad_from0 = bcdata(nn)%jnbeg
        call pushinteger4(j)
        do j=ad_from0,bcdata(nn)%jnend
          ad_from = bcdata(nn)%inbeg
          call pushinteger4(i)
          do i=ad_from,bcdata(nn)%inend
            call pushreal8array(bcdata(nn)%f(i, j, :), size(bcdata(nn)%f&
&                         , 3))
            bcdata(nn)%f(i, j, :) = bcdata(nn)%f(i, j, :)/bcdata(nn)%&
&             dualarea(i, j)
          end do
          call pushinteger4(i - 1)
          call pushinteger4(ad_from)
        end do
        call pushinteger4(j - 1)
        call pushinteger4(ad_from0)
        call pushcontrol1b(0)
      else
        call pushcontrol1b(1)
      end if
      call pushreal8array(sk, size(sk, 1)*size(sk, 2)*size(sk, 3)*size(&
&                   sk, 4))
      call pushreal8array(sj, size(sj, 1)*size(sj, 2)*size(sj, 3)*size(&
&                   sj, 4))
      call pushreal8array(si, size(si, 1)*size(si, 2)*size(si, 3)*size(&
&                   si, 4))
      call pushreal8array(x, size(x, 1)*size(x, 2)*size(x, 3)*size(x, 4)&
&                  )
      call pushreal8array(rlv, size(rlv, 1)*size(rlv, 2)*size(rlv, 3))
      call pushreal8array(gamma, size(gamma, 1)*size(gamma, 2)*size(&
&                   gamma, 3))
      call pushreal8array(rev, size(rev, 1)*size(rev, 2)*size(rev, 3))
      call resetbcpointers(nn, .true.)
      call pushcontrol1b(1)
    else
      call pushcontrol1b(0)
    end if
  end do bocos
! currently the coefficients only contain the surface integral
! of the pressure tensor. these values must be scaled to
! obtain the correct coefficients.
  call pushreal8(fact)
  fact = two/(gammainf*pinf*machcoef*machcoef*surfaceref*lref*lref*&
&   scaledim)
  call pushreal8(cfp(1))
  cfp(1) = cfp(1)*fact
  call pushreal8(cfp(2))
  cfp(2) = cfp(2)*fact
  call pushreal8(cfv(1))
  cfv(1) = cfv(1)*fact
  call pushreal8(cfv(2))
  cfv(2) = cfv(2)*fact
  call pushreal8(fact)
  fact = fact/(lengthref*lref)
  call pushreal8(cmp(1))
  cmp(1) = cmp(1)*fact
  call pushreal8(cmp(2))
  cmp(2) = cmp(2)*fact
  call pushreal8(cmv(1))
  cmv(1) = cmv(1)*fact
  call pushreal8(cmv(2))
  cmv(2) = cmv(2)*fact
  factd = cmv(3)*cmvd(3)
  cmvd(3) = fact*cmvd(3)
  call popreal8(cmv(2))
  factd = factd + cmv(2)*cmvd(2)
  cmvd(2) = fact*cmvd(2)
  call popreal8(cmv(1))
  factd = factd + cmp(3)*cmpd(3) + cmv(1)*cmvd(1)
  cmvd(1) = fact*cmvd(1)
  cmpd(3) = fact*cmpd(3)
  call popreal8(cmp(2))
  factd = factd + cmp(2)*cmpd(2)
  cmpd(2) = fact*cmpd(2)
  call popreal8(cmp(1))
  factd = factd + cmp(1)*cmpd(1)
  cmpd(1) = fact*cmpd(1)
  call popreal8(fact)
  tempd5 = factd/(lref*lengthref)
  lengthrefd = lengthrefd - fact*tempd5/lengthref
  factd = cfv(3)*cfvd(3) + tempd5
  cfvd(3) = fact*cfvd(3)
  call popreal8(cfv(2))
  factd = factd + cfv(2)*cfvd(2)
  cfvd(2) = fact*cfvd(2)
  call popreal8(cfv(1))
  factd = factd + cfp(3)*cfpd(3) + cfv(1)*cfvd(1)
  cfvd(1) = fact*cfvd(1)
  cfpd(3) = fact*cfpd(3)
  call popreal8(cfp(2))
  factd = factd + cfp(2)*cfpd(2)
  cfpd(2) = fact*cfpd(2)
  call popreal8(cfp(1))
  factd = factd + cfp(1)*cfpd(1)
  cfpd(1) = fact*cfpd(1)
  call popreal8(fact)
  temp2 = machcoef**2*scaledim
  temp1 = surfaceref*lref**2
  temp0 = temp1*gammainf*pinf
  tempd6 = -(two*factd/(temp0**2*temp2**2))
  tempd7 = temp2*temp1*tempd6
  gammainfd = gammainfd + pinf*tempd7
  pinfd = pinfd + gammainf*tempd7
  machcoefd = machcoefd + scaledim*temp0*2*machcoef*tempd6
  scaledimd = temp0*machcoef**2*tempd6
  revd = 0.0_8
  pd = 0.0_8
  rlvd = 0.0_8
  sid = 0.0_8
  sjd = 0.0_8
  skd = 0.0_8
  do ii1=1,size(viscsubfaced)
    viscsubfaced(ii1)%tau = 0.0_8
  end do
  do ii1=1,size(bcdata)
    bcdatad(ii1)%dualarea = 0.0_8
  end do
  veldirfreestreamd = 0.0_8
  vd = 0.0_8
  refpointd = 0.0_8
  do nn=nbocos,1,-1
    call popcontrol1b(branch)
    if (branch .ne. 0) then
      call popreal8array(rev, size(rev, 1)*size(rev, 2)*size(rev, 3))
      call popreal8array(gamma, size(gamma, 1)*size(gamma, 2)*size(gamma&
&                  , 3))
      call popreal8array(rlv, size(rlv, 1)*size(rlv, 2)*size(rlv, 3))
      call popreal8array(x, size(x, 1)*size(x, 2)*size(x, 3)*size(x, 4))
      call popreal8array(si, size(si, 1)*size(si, 2)*size(si, 3)*size(si&
&                  , 4))
      call popreal8array(sj, size(sj, 1)*size(sj, 2)*size(sj, 3)*size(sj&
&                  , 4))
      call popreal8array(sk, size(sk, 1)*size(sk, 2)*size(sk, 3)*size(sk&
&                  , 4))
      call resetbcpointers_b(nn, .true.)
      call popcontrol1b(branch)
      if (branch .eq. 0) then
        call popinteger4(ad_from0)
        call popinteger4(ad_to1)
        do j=ad_to1,ad_from0,-1
          call popinteger4(ad_from)
          call popinteger4(ad_to0)
          do i=ad_to0,ad_from,-1
            call popreal8array(bcdata(nn)%f(i, j, :), size(bcdata(nn)%f&
&                        , 3))
            temp = bcdata(nn)%dualarea(i, j)
            bcdatad(nn)%dualarea(i, j) = bcdatad(nn)%dualarea(i, j) + &
&             sum(-(bcdata(nn)%f(i, j, :)*bcdatad(nn)%f(i, j, :)/temp))/&
&             temp
            bcdatad(nn)%f(i, j, :) = bcdatad(nn)%f(i, j, :)/temp
          end do
          call popinteger4(i)
        end do
        call popinteger4(j)
      end if
      call popcontrol1b(branch)
      if (branch .eq. 0) then
        call popinteger4(ad_to)
        do ii=ad_to,0,-1
          mxd = cmvd(1)
          myd = cmvd(2)
          mzd = cmvd(3)
          j = ii/(bcdata(nn)%inend-bcdata(nn)%inbeg) + bcdata(nn)%jnbeg &
&           + 1
          call popreal8(bcdata(nn)%f(i, j, 3))
          fzd = bcdatad(nn)%f(i-1, j, 3) + bcdatad(nn)%f(i-1, j-1, 3) + &
&           bcdatad(nn)%f(i, j-1, 3) + bcdatad(nn)%f(i, j, 3)
          call popreal8(bcdata(nn)%f(i-1, j, 3))
          call popreal8(bcdata(nn)%f(i, j-1, 3))
          call popreal8(bcdata(nn)%f(i-1, j-1, 3))
          call popreal8(bcdata(nn)%f(i, j, 2))
          fyd = bcdatad(nn)%f(i-1, j, 2) + bcdatad(nn)%f(i-1, j-1, 2) + &
&           bcdatad(nn)%f(i, j-1, 2) + bcdatad(nn)%f(i, j, 2)
          call popreal8(bcdata(nn)%f(i-1, j, 2))
          call popreal8(bcdata(nn)%f(i, j-1, 2))
          call popreal8(bcdata(nn)%f(i-1, j-1, 2))
          call popreal8(bcdata(nn)%f(i, j, 1))
          fxd = bcdatad(nn)%f(i-1, j, 1) + bcdatad(nn)%f(i-1, j-1, 1) + &
&           bcdatad(nn)%f(i, j-1, 1) + bcdatad(nn)%f(i, j, 1)
          call popreal8(bcdata(nn)%f(i-1, j, 1))
          call popreal8(bcdata(nn)%f(i, j-1, 1))
          call popreal8(bcdata(nn)%f(i-1, j-1, 1))
          call popreal8(fz)
          fzd = cfvd(3) - xc*myd + yc*mxd + fourth*fzd
          call popreal8(fy)
          fyd = xc*mzd + cfvd(2) - zc*mxd + fourth*fyd
          call popreal8(fx)
          fxd = cfvd(1) - yc*mzd + zc*myd + fourth*fxd
          xcd = fy*mzd - fz*myd
          ycd = fz*mxd - fx*mzd
          zcd = fx*myd - fy*mxd
          call popreal8(zc)
          tempd = fourth*zcd
          xxd(i, j, 3) = xxd(i, j, 3) + tempd
          xxd(i+1, j, 3) = xxd(i+1, j, 3) + tempd
          xxd(i, j+1, 3) = xxd(i, j+1, 3) + tempd
          xxd(i+1, j+1, 3) = xxd(i+1, j+1, 3) + tempd
          refpointd(3) = refpointd(3) - zcd
          call popreal8(yc)
          tempd0 = fourth*ycd
          xxd(i, j, 2) = xxd(i, j, 2) + tempd0
          xxd(i+1, j, 2) = xxd(i+1, j, 2) + tempd0
          xxd(i, j+1, 2) = xxd(i, j+1, 2) + tempd0
          xxd(i+1, j+1, 2) = xxd(i+1, j+1, 2) + tempd0
          refpointd(2) = refpointd(2) - ycd
          call popreal8(xc)
          tempd1 = fourth*xcd
          xxd(i, j, 1) = xxd(i, j, 1) + tempd1
          xxd(i+1, j, 1) = xxd(i+1, j, 1) + tempd1
          xxd(i, j+1, 1) = xxd(i, j+1, 1) + tempd1
          xxd(i+1, j+1, 1) = xxd(i+1, j+1, 1) + tempd1
          refpointd(1) = refpointd(1) - xcd
          tauzz = viscsubface(nn)%tau(i, j, 3)
          tauxz = viscsubface(nn)%tau(i, j, 5)
          tauyz = viscsubface(nn)%tau(i, j, 6)
          tempd2 = -(fact*scaledim*fzd)
          ssid(i, j, 1) = ssid(i, j, 1) + tauxz*tempd2
          ssid(i, j, 2) = ssid(i, j, 2) + tauyz*tempd2
          tauzzd = ssi(i, j, 3)*tempd2
          ssid(i, j, 3) = ssid(i, j, 3) + tauzz*tempd2
          tauxy = viscsubface(nn)%tau(i, j, 4)
          tauyy = viscsubface(nn)%tau(i, j, 2)
          tempd4 = -(fact*scaledim*fyd)
          tauyzd = ssi(i, j, 3)*tempd4 + ssi(i, j, 2)*tempd2
          ssid(i, j, 1) = ssid(i, j, 1) + tauxy*tempd4
          tauyyd = ssi(i, j, 2)*tempd4
          ssid(i, j, 2) = ssid(i, j, 2) + tauyy*tempd4
          ssid(i, j, 3) = ssid(i, j, 3) + tauyz*tempd4
          tauxx = viscsubface(nn)%tau(i, j, 1)
          scaledimd = scaledimd - fact*(tauxy*ssi(i, j, 1)+tauyy*ssi(i, &
&           j, 2)+tauyz*ssi(i, j, 3))*fyd - fact*(tauxx*ssi(i, j, 1)+&
&           tauxy*ssi(i, j, 2)+tauxz*ssi(i, j, 3))*fxd - fact*(tauxz*ssi&
&           (i, j, 1)+tauyz*ssi(i, j, 2)+tauzz*ssi(i, j, 3))*fzd
          tempd3 = -(fact*scaledim*fxd)
          tauxzd = ssi(i, j, 3)*tempd3 + ssi(i, j, 1)*tempd2
          tauxyd = ssi(i, j, 2)*tempd3 + ssi(i, j, 1)*tempd4
          tauxxd = ssi(i, j, 1)*tempd3
          ssid(i, j, 1) = ssid(i, j, 1) + tauxx*tempd3
          ssid(i, j, 2) = ssid(i, j, 2) + tauxy*tempd3
          ssid(i, j, 3) = ssid(i, j, 3) + tauxz*tempd3
          viscsubfaced(nn)%tau(i, j, 6) = viscsubfaced(nn)%tau(i, j, 6) &
&           + tauyzd
          viscsubfaced(nn)%tau(i, j, 5) = viscsubfaced(nn)%tau(i, j, 5) &
&           + tauxzd
          viscsubfaced(nn)%tau(i, j, 4) = viscsubfaced(nn)%tau(i, j, 4) &
&           + tauxyd
          viscsubfaced(nn)%tau(i, j, 3) = viscsubfaced(nn)%tau(i, j, 3) &
&           + tauzzd
          viscsubfaced(nn)%tau(i, j, 2) = viscsubfaced(nn)%tau(i, j, 2) &
&           + tauyyd
          viscsubfaced(nn)%tau(i, j, 1) = viscsubfaced(nn)%tau(i, j, 1) &
&           + tauxxd
          call popinteger4(j)
          call popinteger4(i)
        end do
      end if
      call popreal8array(ww2, size(ww2, 1)*size(ww2, 2)*size(ww2, 3))
      call popreal8array(ssi, size(ssi, 1)*size(ssi, 2)*size(ssi, 3))
      call popreal8array(pp2, size(pp2, 1)*size(pp2, 2))
      call popreal8array(pp1, size(pp1, 1)*size(pp1, 2))
      call popreal8array(xx, size(xx, 1)*size(xx, 2)*size(xx, 3))
      do ii1=size(bcdata),1,-1
        call popreal8array(bcdata(ii1)%dualarea, size(bcdata(ii1)%&
&                    dualarea, 1)*size(bcdata(ii1)%dualarea, 2))
      end do
      do ii1=size(bcdata),1,-1
        call popreal8array(bcdata(ii1)%f, size(bcdata(ii1)%f, 1)*size(&
&                    bcdata(ii1)%f, 2)*size(bcdata(ii1)%f, 3))
      end do
      call lookreal8array(v, 3)
      do ii=0,(bcdata(nn)%jnend-bcdata(nn)%jnbeg)*(bcdata(nn)%inend-&
&         bcdata(nn)%inbeg)-1
        i = mod(ii, bcdata(nn)%inend - bcdata(nn)%inbeg) + bcdata(nn)%&
&         inbeg + 1
        j = ii/(bcdata(nn)%inend-bcdata(nn)%inbeg) + bcdata(nn)%jnbeg + &
&         1
! compute the average pressure minus 1 and the coordinates
! of the centroid of the face relative from from the
! moment reference point. due to the usage of pointers for
! the coordinates, whose original array starts at 0, an
! offset of 1 must be used. the pressure is multipled by
! fact to account for the possibility of an inward or
! outward pointing normal.
        pm1 = fact*(half*(pp2(i, j)+pp1(i, j))-pinf)*scaledim
        xc = fourth*(xx(i, j, 1)+xx(i+1, j, 1)+xx(i, j+1, 1)+xx(i+1, j+1&
&         , 1)) - refpoint(1)
        yc = fourth*(xx(i, j, 2)+xx(i+1, j, 2)+xx(i, j+1, 2)+xx(i+1, j+1&
&         , 2)) - refpoint(2)
        zc = fourth*(xx(i, j, 3)+xx(i+1, j, 3)+xx(i, j+1, 3)+xx(i+1, j+1&
&         , 3)) - refpoint(3)
! compute the force components.
        fx = pm1*ssi(i, j, 1)
        fy = pm1*ssi(i, j, 2)
        fz = pm1*ssi(i, j, 3)
! update the inviscid force and moment coefficients.
! divide by 4 so we can scatter
! scatter 1/4 of the force to each of the nodes:
! scatter a quarter of the area to each node:
        qa = fourth*sqrt(ssi(i, j, 1)**2+ssi(i, j, 2)**2+ssi(i, j, 3)**2&
&         )
! get normalized surface velocity:
        v(1) = ww2(i, j, ivx)
        v(2) = ww2(i, j, ivy)
        v(3) = ww2(i, j, ivz)
        tmp0 = v/(sqrt(v(1)**2+v(2)**2+v(3)**2)+1e-16)
        call pushreal8array(v, 3)
        v = tmp0
! dot product with free stream
        sensor = -(v(1)*veldirfreestream(1)+v(2)*veldirfreestream(2)+v(3&
&         )*veldirfreestream(3))
!now run through a smooth heaviside function:
        call pushreal8(sensor)
        sensor = one/(one+exp(-(2*sepsensorsharpness*(sensor-&
&         sepsensoroffset))))
! and integrate over the area of this cell and save:
        plocal = pp2(i, j)
        tmp = two/(gammainf*pinf*machcoef*machcoef)
        cp = tmp*(plocal-pinf)
        sigma = 1.4
        sensor1 = -cp - sigma
        call pushreal8(sensor1)
        sensor1 = one/(one+exp(-(2*10*sensor1)))
        mxd = cmpd(1)
        myd = cmpd(2)
        mzd = cmpd(3)
        sensord = sepsensord
        sensor1d = cavitationd
        qad = four*sensor*sensord + four*sensor1*sensor1d
        sensor1d = four*qa*sensor1d
        call popreal8(sensor1)
        temp9 = -(10*2*sensor1)
        temp8 = one + exp(temp9)
        sensor1d = exp(temp9)*one*10*2*sensor1d/temp8**2
        cpd = -sensor1d
        tmpd = (plocal-pinf)*cpd
        plocald = tmp*cpd
        temp7 = gammainf*pinf*machcoef**2
        tempd9 = -(two*tmpd/temp7**2)
        tempd8 = machcoef**2*tempd9
        pinfd = pinfd + gammainf*tempd8 - tmp*cpd
        gammainfd = gammainfd + pinf*tempd8
        machcoefd = machcoefd + gammainf*pinf*2*machcoef*tempd9
        pp2d(i, j) = pp2d(i, j) + plocald
        sensord = four*qa*sensord
        call popreal8(sensor)
        temp6 = -(2*sepsensorsharpness*(sensor-sepsensoroffset))
        temp5 = one + exp(temp6)
        sensord = exp(temp6)*one*sepsensorsharpness*2*sensord/temp5**2
        vd(1) = vd(1) - veldirfreestream(1)*sensord
        veldirfreestreamd(1) = veldirfreestreamd(1) - v(1)*sensord
        vd(2) = vd(2) - veldirfreestream(2)*sensord
        veldirfreestreamd(2) = veldirfreestreamd(2) - v(2)*sensord
        vd(3) = vd(3) - veldirfreestream(3)*sensord
        veldirfreestreamd(3) = veldirfreestreamd(3) - v(3)*sensord
        call popreal8array(v, 3)
        tmpd0 = vd
        temp3 = v(1)**2 + v(2)**2 + v(3)**2
        temp4 = sqrt(temp3)
        tempd10 = tmpd0/(temp4+1e-16)
        vd = tempd10
        if (temp3 .eq. 0.0_8) then
          tempd11 = 0.0
        else
          tempd11 = sum(-(v*tempd10/(temp4+1e-16)))/(2.0*temp4)
        end if
        vd(1) = vd(1) + 2*v(1)*tempd11
        vd(2) = vd(2) + 2*v(2)*tempd11
        vd(3) = vd(3) + 2*v(3)*tempd11
        ww2d(i, j, ivz) = ww2d(i, j, ivz) + vd(3)
        vd(3) = 0.0_8
        ww2d(i, j, ivy) = ww2d(i, j, ivy) + vd(2)
        vd(2) = 0.0_8
        ww2d(i, j, ivx) = ww2d(i, j, ivx) + vd(1)
        vd(1) = 0.0_8
        qad = qad + bcdatad(nn)%dualarea(i-1, j) + bcdatad(nn)%dualarea(&
&         i-1, j-1) + bcdatad(nn)%dualarea(i, j-1) + bcdatad(nn)%&
&         dualarea(i, j)
        if (ssi(i, j, 1)**2 + ssi(i, j, 2)**2 + ssi(i, j, 3)**2 .eq. &
&           0.0_8) then
          tempd12 = 0.0
        else
          tempd12 = fourth*qad/(2.0*sqrt(ssi(i, j, 1)**2+ssi(i, j, 2)**2&
&           +ssi(i, j, 3)**2))
        end if
        ssid(i, j, 1) = ssid(i, j, 1) + 2*ssi(i, j, 1)*tempd12
        ssid(i, j, 2) = ssid(i, j, 2) + 2*ssi(i, j, 2)*tempd12
        ssid(i, j, 3) = ssid(i, j, 3) + 2*ssi(i, j, 3)*tempd12
        fzd = bcdatad(nn)%f(i-1, j, 3) + bcdatad(nn)%f(i-1, j-1, 3) + &
&         bcdatad(nn)%f(i, j-1, 3) + bcdatad(nn)%f(i, j, 3)
        fyd = bcdatad(nn)%f(i-1, j, 2) + bcdatad(nn)%f(i-1, j-1, 2) + &
&         bcdatad(nn)%f(i, j-1, 2) + bcdatad(nn)%f(i, j, 2)
        fxd = bcdatad(nn)%f(i-1, j, 1) + bcdatad(nn)%f(i-1, j-1, 1) + &
&         bcdatad(nn)%f(i, j-1, 1) + bcdatad(nn)%f(i, j, 1)
        fzd = cfpd(3) - xc*myd + yc*mxd + fourth*fzd
        fyd = xc*mzd + cfpd(2) - zc*mxd + fourth*fyd
        fxd = cfpd(1) - yc*mzd + zc*myd + fourth*fxd
        xcd = fy*mzd - fz*myd
        ycd = fz*mxd - fx*mzd
        zcd = fx*myd - fy*mxd
        pm1d = ssi(i, j, 2)*fyd + ssi(i, j, 1)*fxd + ssi(i, j, 3)*fzd
        ssid(i, j, 3) = ssid(i, j, 3) + pm1*fzd
        ssid(i, j, 2) = ssid(i, j, 2) + pm1*fyd
        ssid(i, j, 1) = ssid(i, j, 1) + pm1*fxd
        tempd13 = fourth*zcd
        xxd(i, j, 3) = xxd(i, j, 3) + tempd13
        xxd(i+1, j, 3) = xxd(i+1, j, 3) + tempd13
        xxd(i, j+1, 3) = xxd(i, j+1, 3) + tempd13
        xxd(i+1, j+1, 3) = xxd(i+1, j+1, 3) + tempd13
        refpointd(3) = refpointd(3) - zcd
        tempd14 = fourth*ycd
        xxd(i, j, 2) = xxd(i, j, 2) + tempd14
        xxd(i+1, j, 2) = xxd(i+1, j, 2) + tempd14
        xxd(i, j+1, 2) = xxd(i, j+1, 2) + tempd14
        xxd(i+1, j+1, 2) = xxd(i+1, j+1, 2) + tempd14
        refpointd(2) = refpointd(2) - ycd
        tempd15 = fourth*xcd
        xxd(i, j, 1) = xxd(i, j, 1) + tempd15
        xxd(i+1, j, 1) = xxd(i+1, j, 1) + tempd15
        xxd(i, j+1, 1) = xxd(i, j+1, 1) + tempd15
        xxd(i+1, j+1, 1) = xxd(i+1, j+1, 1) + tempd15
        refpointd(1) = refpointd(1) - xcd
        tempd16 = fact*scaledim*pm1d
        pp2d(i, j) = pp2d(i, j) + half*tempd16
        pp1d(i, j) = pp1d(i, j) + half*tempd16
        pinfd = pinfd - tempd16
        scaledimd = scaledimd + fact*(half*(pp2(i, j)+pp1(i, j))-pinf)*&
&         pm1d
      end do
      call popreal8array(v, 3)
      call popreal8(sensor1)
      call popreal8(zc)
      call popreal8(sensor)
      call popreal8(tmp)
      call popreal8(yc)
      call popreal8(plocal)
      call popreal8(xc)
      call popinteger4(j)
      call popinteger4(i)
      call popreal8(qa)
      call popreal8array(bcdata(nn)%f, size(bcdata(nn)%f, 1)*size(bcdata&
&                  (nn)%f, 2)*size(bcdata(nn)%f, 3))
      bcdatad(nn)%f = 0.0_8
      call popreal8array(bcdata(nn)%dualarea, size(bcdata(nn)%dualarea, &
&                  1)*size(bcdata(nn)%dualarea, 2))
      bcdatad(nn)%dualarea = 0.0_8
      call popcontrol3b(branch)
      if (branch .lt. 3) then
        if (branch .ne. 0) then
          if (branch .eq. 1) then
            call popreal8(fact)
          else
            call popreal8(fact)
          end if
        end if
      else if (branch .lt. 5) then
        if (branch .eq. 3) then
          call popreal8(fact)
        else
          call popreal8(fact)
        end if
      else if (branch .eq. 5) then
        call popreal8(fact)
      else
        call popreal8(fact)
      end if
      call popreal8array(pp1, size(pp1, 1)*size(pp1, 2))
      call popreal8array(pp2, size(pp2, 1)*size(pp2, 2))
      call popreal8array(ssi, size(ssi, 1)*size(ssi, 2)*size(ssi, 3))
      call setbcpointers_b(nn, .true.)
    end if
  end do
  pointrefd(3) = pointrefd(3) + lref*refpointd(3)
  refpointd(3) = 0.0_8
  pointrefd(2) = pointrefd(2) + lref*refpointd(2)
  refpointd(2) = 0.0_8
  pointrefd(1) = pointrefd(1) + lref*refpointd(1)
  prefd = prefd + scaledimd/pinf
  pinfd = pinfd - pref*scaledimd/pinf**2
end subroutine forcesandmoments_b
