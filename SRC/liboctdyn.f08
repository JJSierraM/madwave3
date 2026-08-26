function factor(from, to) result(output)
   integer, intent(in) :: from, to
   integer :: output, iter

   output = 1
   do iter = from, to
      output = output * iter
   end do

end function factor

! ************************  dwigner  **************************************

subroutine dwigner(dj,jmax2,m2,mp2,beta,ndim)
    implicit none

    integer, intent(in) :: jmax2, m2, mp2, ndim
    real*8, intent(in):: beta
    real*8, intent(inout) :: dj(0:ndim)

    integer, external :: factor

    integer :: jmin2, ene, int_val, eme, ind, indp, indm
    real*8 :: y, xj, xjp, xm, xmp, xjm, xsin, xcos, a, b, fact, sign, facnum, facden1, facden2, factot, facj, fac1jm, fac2jm, facden, facjm

! *     ***********************************************
! *     *   Wigner matrizes: d^J_{M,M'} (cos(beta))   *
! *     *       kept in dj(2*J) to consider           *
! *     *         halfinteger values of J             *
! *     *  input:                                     *
! *     *     jmax2 is 2 * jmax, jmax being           *
! *     *                   the highest  j required   *
! *     *     m2 and mp2: are (m * 2) and (mp *2)     *
! *     *                   respectively              *
! *     *     beta: angle defined between 0 and Pi    *
! *     ***********************************************
   if(jmax2 > ndim)then
         print *, ' ** be carefull with dimensions in dwigner **'
         print *, '    2 jmax > ndim',jmax2,ndim
         stop
   endif
   if(jmax2 < m2 .or. jmax2 < mp2)then
         print *, ' ** jmax2 can NOT be less than  m2 or mp2**'
         print *, '    2 jmax < 2 m or ',jmax2,ndim
         stop
   endif

   dj(:) = 0.0

   jmin2 = max(abs(m2),abs(mp2)) 

   y = dcos(beta)

   ! ***> jmin=0 ---> M=M'=0 ---> d^J_{M,Mp} = P_J(cos(beta))
   ! ********************************************************

   if(jmin2 == 0)then
      dj(jmin2) = 1.d0
      dj(jmin2+2) = y

      do ind = jmin2+4,jmax2,2
         xj = dble(ind)*0.5d0
         dj(ind) = y*(2.d0*xj-1.d0)*dj(ind-2)-(xj-1.d0)*dj(ind-4)
         dj(ind) = dj(ind)/xj
      end do
      
   else
      ! ***> d^J_{M,M'} for J=jmin=max0(iabs(M),iabs(M')) > 0
      ! *a)      Normalization factor
      ene = min(abs(m2),abs(mp2))/2
      eme = max(abs(m2),abs(mp2))/2

      fact = sqrt(dble(factor(jmin2/2+ene+1, jmin2/2+eme)/factor(jmin2/2-eme+1,jmin2/2-ene)))

      ! *b)
      sign = 1.d0
      if (m2-mp2 > 0) sign = (-1)**(int(dble(m2-mp2)*0.5d0+0.5d0))

      xsin = dsin(beta*0.5d0)
      xcos = dcos(beta*0.5d0)

      a = xsin**(abs(m2-mp2)*0.5d0)
      b = xcos**(abs(m2+mp2)*0.5d0)

      ! *c)

      dj(jmin2) = fact*sign*a*b

      ! ***> d^J_{M,M'} for J=jmin+1
      ind = jmin2 + 2
      indp = ind - 2
      if (ind <= jmax2) then 
         xjp = dble(ind)*0.5d0
         xj  = xjp-1.d0
         xm  = dble(m2)*0.5d0
         xmp = dble(mp2)*0.5d0

         facnum  = xjp*(2.d0*xj+1.d0)
         facden1 = xjp*xjp-xm*xm
         facden2 = xjp*xjp-xmp*xmp
         factot  = facnum/dsqrt(facden1*facden2)

         facj = y - xm * xmp /(xj*(xj+1.d0))

         dj(ind) = factot*facj*dj(indp)

         ! ***> d^J_{M,M'} for J>jmin+1

         do ind  = jmin2+4,jmax2,2
            indp = ind-2
            indm = ind-4

            xjp = dble(ind)*0.5d0
            xj  = xjp-1.d0
            xjm = xjp-2.d0

            facnum  = xjp*(2.d0*xj+1.d0)
            facden1 = xjp*xjp-xm*xm
            facden2 = xjp*xjp-xmp*xmp
            factot  = facnum/dsqrt(facden1*facden2)

            facj = y - xm*xmp/(xj*(xj+1.d0))

            fac1jm = xj*xj-xm*xm
            fac2jm = xj*xj-xmp*xmp
            facden = xj*(2.d0*xj+1.d0)
            facjm  = dsqrt(fac1jm*fac2jm)/facden

            dj(ind) = factot*( facj*dj(indp) -facjm*dj(indm) )
         enddo
      endif
   endif
end

! *************************** gauleg **************************

subroutine gauleg(x, w, n)
   implicit none
   integer, intent(in) :: n
   real*8 :: x1 = -1.0, x2 = 1.0
   real*8, dimension(n), intent(out) :: x, w
   integer :: i, j, m
   real*8 :: p1, p2, p3, pp, xl, xm, z, z1
   real*8, parameter :: eps=3.d-14

   m = (n+1)/2
   xm = 0.5d0 * (x2+x1)
   xl = 0.5d0 * (x2-x1)
   
   do i=1,m
      z = cos(3.141592654d0 * (i - 0.25d0) / (n + 0.5d0))
      z1 = 0.0
      do while(abs(z-z1) > eps)
         p1 = 1.0d0
         p2 = 0.0d0
         do j=1,n
            p3 = p2
            p2 = p1
            p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
         end do
         pp = n*(z*p1-p2)/(z*z-1.0d0)
         z1 = z
         z = z1 - p1/pp
      end do
      x(i) = xm - xl*z
      x(n+1-i) = xm + xl*z
      w(i) = (2.0d0*xl)/((1.0d0-z*z)*pp*pp)
      w(n+1-i) = w(i)
   end do
end subroutine gauleg

! ********************  FACTORIAL  *********************************
subroutine factorial
   common/fct/fact(0:10000)
   ! IMPLICIT real*8(A-H,O-Z)

   fact(0)=0.d0
   fact(1)=0.d0
   n=10000
   do i=2, n
      fact(i)=fact(i-1)+dlog(dble(i))
   end do 
end

! *************************  tqli  ****************************

subroutine tqli(d, e, n, np)
   ! implicit real*8(a-h,o-z)
   implicit none 
   real*8, intent(inout) ::  d(np), e(np)
   integer, intent(in) :: n, np

   integer :: l=0, m=0, i=0
   real*8 :: g=0, r=0, s=0, c=0, p=0, f=0, b=0
! ****************************************************
! ***    diagonalization of a tridiagonal matrix   ***
! ***         from numerical recipies              ***
! ****************************************************

   if (n > 1) then
! c overide of silly num recipes convention for off-diagonal
! c elements.  It is now assumed off diags are in E(1)..E(N-1)
! c        DO 11 I=2,N
! c          E(I-1)=E(I)
! c11      CONTINUE
      e(n)=0.d0
      do l=1, n
        do while (1 > 0)
            do m = l, n
            if (e(m) == 0) then
                exit
            end if
            end do 

            if(m == l) exit

            ! if (iter == 500) print *, 'too many iterations'
            ! iter = iter + 1

            g = (d(l+1) - d(l)) / (2.d0 * e(l))
            r = sqrt(g**2 + 1.d0)
            g = d(m) - d(l) + e(l) / (g + sign(r,g))

            s = 1.d0
            c = 1.d0
            p = 0.d0

            do i = m-1, l, -1
                f = s * e(i)
                b = c * e(i)
                if (abs(f) >= abs(g)) then
                    c = g/f
                    r = sqrt(c**2 + 1.d0)
                    e(i+1) = f*r
                    s = 1.d0 / r
                    c = c * s
                else
                    s = f/g
                    r = sqrt(s**2 + 1.d0)
                    e(i+1) = g * r
                    c = 1.d0 / r
                    s = s * c
                end if
                g = d(i+1) - p
                r = (d(i)-g) * s + 2.d0*c*b
                p = s*r
                d(i+1) = g+p
                g = c*r-b
            end do
            d(l) = d(l) - p
            e(l) = g
            e(M) = 0.d0
        end do
      end do
   end if
end subroutine

! **************************  SPLIN  **********************************

subroutine splinqq(f, x, iold, nx, r, ndim, spl)
   implicit none
   real*8, intent(in) :: f(ndim, 2), x(ndim), r
   real*8, intent(out) :: spl
   integer, intent(inout) :: iold
   integer, intent(in) :: nx, ndim

   real*8 :: u(4), ro, yo, yp, ain, c8, c6, hi, xr
   integer :: idol

   if ( r >= x(nx) ) then
      ro = x(nx)
      yo = f(nx, 1)
      yp = f(nx,2)
      ain = yo + yp * ro / 6.0
      c8 = - ain * 3.0 * ro ** 8
      c6 = yo * ro ** 6 - c8 / ro / ro
      spl = c6 / r ** 6 + c8 / r ** 8
   else
      do idol = iold, nx
         if ( r < x(idol) ) then
            exit
         end if
      end do
      iold = idol
      hi = x(idol) - x(idol-1)
      xr = (r-x(idol-1)) / hi
      u(1) = xr * xr * (-2.0*xr + 3.0)
      u(2) = 1.0 - u(1)
      u(3) = hi * xr * xr * (xr - 1.0)
      u(4) = hi * xr * ((xr-2.0) * xr + 1.0)
      spl = u(1) * f(idol,1) + u(2) * f(idol-1,1) + u(3) * f(idol,2) + u(4) * f(idol-1,2)
   end if

end subroutine splinqq

subroutine splset (f, x, nx, ndim)
! C***********************************************************************
! c*        This routine sets the spline interpolation on the grid       *
! c*            (x(i),i=1,nx) for the function (f(i),i=1,nx).            *
! C***********************************************************************
   implicit none
   integer, parameter :: nxmax = 60000
   integer, intent(in) :: nx, ndim
   real*8, intent(inout) :: f(ndim,2)
   real*8, intent(in) :: x(ndim)

   real*8 :: hx(nxmax), rlx(nxmax-1), rmux(nxmax-1), xi(nxmax-1), b(nxmax-1)
   real*8 :: ab(4), yz(4), a(4), p0, pn
   integer:: i, nx2, man

   real*8, external :: dlagra

   if ( nx > nxmax ) then
        print "(//,2X,20('*'),' STOP IN SPLSET ',20('*'),/)"
        print "(2X,'DIMENSION PARAMETER NXMAX = ',I5,' TOO SMALL, ',I5,' REQUIRED')", nxmax, nx
   else
        nx2 = nx-2
        do i = 2, nx
            hx(i-1) = x(i) - x(i-1)
        end do
        do i = 1, nx2
            rlx(i)  = hx(i+1) / (hx(i) + hx(i+1))
            rmux(i) = 1.0 - rlx(i)
        end do
        man = nx - 3
        do i = 1, 4
            a(i) = x(man)
            man = man + 1
        end do

        ! spline-fit of p(x)
        do i = 1, 4
            ab(i) = f(i,1)
            yz(i) = f(nx+i-4,1)      
        end do
        p0 = dlagra(x, ab, 4, 1)
        f(1,2) = p0
        pn = dlagra(a, yz, 4, 4)
        f(nx, 2) = pn
        do i = 1, nx2
            b(i) = 3.0 * rlx(i) / hx(i) * (f(i+1,1) - f(i,1)) + 3.0 * rmux(i) / hx(i+1) * (f(i+2,1) - f(i+1,1))
        end do
        b(1) = b(1) - rlx(1) * p0
        b(nx2) = b(nx2) - rmux(nx2) * pn
        call jordan(rmux, rlx, xi, nx2, b)
        do i = 1, nx2
            f(i+1, 2) = xi(i)
        end do
    end if
end subroutine splset

! *********************************  SCHR *************************

subroutine schr(e0,rmin,rmax,n,maxit,eps,e2,kv,v,p,npunt)
   implicit none
   real*8, intent(in)  :: e0,rmin,rmax,eps,v(npunt)
   real*8, intent(inout):: p(npunt), e2
   integer, intent(in) :: n, maxit, npunt
   integer, intent(inout) :: kv

   integer  :: itry, i, it, xit, m, m1, j, msave, i1, nl
   real*8   :: y(3), h, h2, hv, e, test, de, gn, gi, apr, pm, yin, eold, yout, ym, df, f, dold, schrod, sn

   itry = 0
   h = (rmax - rmin) / (dble(n)-1)
   h2 = h**2
   hv = h2 / 12.0
   e = e0
   test = -1.0
   de = 0.0

   do i = 1, n
      p(i) = 0.0
   end do
   do it = 1, maxit
      xit = it
      p(n) = 1.d-30
      gn = v(n) - e
      gi = v(n-1) - e
      if ( gi < 0.0 ) then
         e = v(n-2)
      end if
      apr = (rmax - h) * sqrt(gi) - rmax * sqrt(gn)
      if ( apr > 50.0 ) then
         apr = 50.0
      end if
      p(n-1) = p(n) * exp(-apr)
      y(1) = (1.d0 - hv * gn) * p(n)
      y(2) = (1.d0 - hv * gi) * p(n-1)

      ! integration
      m = n-2
      do while (.TRUE.)
            y(3) = y(2) + ((y(2) - y(1)) + h2 * gi * p(m+1))
            gi = v(m) - e
            p(m) = y(3) / (1.d0 - hv * gi)
            if ( abs(p(m)) > 1.d+34) then
               m1 = m+1
               pm = p(m1)
               print *, "(2X,'PM = ',E16.8/)"
               do j = m1, n
                  p(j) = p(j) / pm
               end do
               y(1) = y(1) / pm
               y(2) = y(2) / pm
               y(3) = y(3) / pm
               gi = v(m1) - e
               cycle
            end if
            if ( abs(p(m)) <= abs(p(m+1)) .OR. m <= 2) then
               exit
            end if
            y(1) = y(2)
            y(2) = y(3)
            m = m-1
      end do
      pm = p(m) 
      msave = m
      yin = y(2)/pm
      do j = m, n
         p(j) = p(j) / pm
      end do
      p(1) = 1.d0
      y(1) = 0.d0
      gi = v(1) - e
      y(2) = (1.d0 - hv * gi) * p(1)
      do i = 2, m
         y(3) = y(2) + ((y(2) - y(1)) + h2 * gi * p(i-1))
         gi = v(i) - e
         p(i) = y(3) / (1.0 - hv * gi)
         if (abs(p(i)) > 1.d+34) then
            i1 = i-1
            pm = p(i1)
            do j = 1, i1
               p(j) = p(j) / pm
            end do
            y(1) = y(1) / pm
            y(2) = y(2) / pm
            y(3) = y(3) / pm
            gi = v(i1) - e
         else
            y(1) = y(2)
            y(2) = y(3)
         end if
      end do

      pm = p(m)
      if ( pm /= 0 ) then
         yout = y(1) / pm
         ym = y(3) / pm
         do j = 1, m
            p(j) = p(j) / pm
         end do
         df = 0.0
         do j = 1, n
            df = df - p(j) ** 2
         end do
         f = (-yout - yin + 2.d0 * ym) / h2 + (v(m) - e)
         dold = de
      end if
      if ( abs(f) > 1.d+37 ) then
         f = 9.99999d+29
         df = -f
         de = abs(0.0001d0*e)
      else
         de = -f / df
      end if
      eold = e
      e = e + de
      test = max((abs(dold)-abs(de)),test)
      if ( test < 0.0 ) then
         cycle
      end if
      if ( abs(e - eold) < eps ) then
         exit
      end if
   end do
   schrod = 0.0
   kv = 0
   nl = n - 2
   do j = 3, nl
      if ( p(j) < 0 ) then
         if (p(j-1) < 0) then
            cycle
         else if ( p(j-1) == 0 ) then
            if ( p(j+1) < 0) then
               if ( p(j-2) > 0 ) then
                  kv = kv + 1
               end if
            end if
         else
            if ( p(j+1) < 0) then
               if ( p(j-2) >= 0 ) then
                  kv = kv + 1
               end if
            end if
         end if
      else
         if ( p(j-1) < 0) then
            if ( p(j+1) >= 0) then
               if ( p(j-2) < 0 ) then
                  kv = kv + 1
               end if
            end if
         end if
      end if
   end do
   e2 = e

   sn = sqrt(-h * df)
   do j = 1, n
      p(j) = p(j) / sn
   end do
end subroutine schr

subroutine besjot ( l, x, f, df, r )
! c   this program generates the standard and modified versions of the spherical
! c   bessel-functions of first(besjot)- and second(bessen)-kind respectively.
! c   for definitions compare: 'nbs-h1ndbook of mathematical functions',
! c   (abramowitz+stegun,eds./n.y.:1964), sections 10.1.1 on page 437 for
! c   standard versions and ss.10.2.2 + 10.2.3 on p.443 for the modified ones.
! c   l=index(natural numbers including zero), x=argument(real,d.p.), f=output.
! c
! c   the sign of the argument is used to determine the versions:
! c   the outcomes f(=first-kind-functions) and g(=second-kind-f.) must be
! c   divided (resp. multiplied) by the l-th power of the reduction logfac
! c   to get the modified versions, use argument with negative sign }
! c
! c   by formulas 10.1.31 on page 439 loc.cit. and 10.2.7 on p.443 ibid.,
! c   solutions have been tested to be correct to twelve places at least in the
! c   range combining x=1...441 and l=0...340 .
! c
! c   besjot is divided into three parts, corresponding to wether x > l, o
! c   while x < l, beeing 0.5*x*x < 2*l or 0.5*x*x > 2*l respectively .
! c
! c
! c modif. pour x plus grand que l dans besjot  : r=1
! c        qq soit x dans bessen : r=1
! c   pour eviter les overflows ou underflows dans le prog. appele pour 50
! c      version jan. 77
! c
! c  j.m.launay, meudon, france
! c  updated to Fortran 2008: J. Sierra, IFF-CSIC, Spain
   implicit none
   integer, intent(in)  :: l
   real*8, intent(in)   :: x
   real*8, intent(inout):: f, df, r
   
   integer :: i, n, j
   real*8 :: sinix, cosix, w, pi, xr, z, a, f0, g0, f1, f2, b0, b1, b2, y, s0, s1, p0, p1, c0, c1, q

   f = 0.0
   r = 1.0
   if ( x == 0 ) then
      if ( l == 0 ) then
         f = 1.0
      end if
      df = 0.0
      return
   else if ( x < 0 ) then
      sinix = sinh(-x)
      cosix = cosh(-x)
      w = -1.0
   else 
      pi = 6.283185307179586
      xr = mod(x, pi)
      sinix = sin(xr)
      cosix = cos(xr)
      w = 1.0
   end if
   z = 1.0 / abs(x)
   a = dble(l)
   r = a * Z

   if ( abs(x) - a <= 0 ) then
      if ( 0.5 * x * x - 2.0 * a > 0 ) then
         n = a + 25.0 + sqrt(a)
         b0 = 0.0
         b1 = 1.0
         do j = 1, n
            b2 = w * (b1 * dble(2 * (n-j) + 3) * z - b0/r) / r
            b0 = b1
            if ( n-l-j == 0) then
               f = b2
            else if ( n-l+1-j == 0 ) then
               df = b2
            end if
            b1 = b2
         end do
         df = w ** (l-1) * (df/b1) * sinix * z
         f = w ** l * (f/b1) * sinix * z
         df = df * r - (l+1) * f * z
         return 
      else
         y = -w * 0.5 * x * x
         s0 = dble(2*l-1)
         s1 = dble(2*l+1)
         p0 = 1.d0
         p1 = 1.d0
         c0 = 1.d0
         c1 = 1.d0
         do i = 1, 15
            s0 = s0 + 2.d0
            s1 = s1 + 2.d0
            p0 = y * p0 / (s0 * dble(i))
            p1 = y * p1 / (s1 * dble(i))
            c0 = c0 + p0
            c1 = c1 + p1
         end do
         q = 1.0
         if ( l /= 1 ) then
            j = l - 1
            do i = 1, j
               q = q * a / dble(2*i+1)
            end do
            f = q * a / dble(2*l+1) * c1
            df = q * c0 * r - dble(l+1) * f * z
            return 
         end if
            f = c1 / 3.d0
            df = c0*r - 2.d0*f*z
            return
      end if
   else
      f0 = sinix * z
      g0 =-cosix * z
      r = 1.0
      if ( l <= 0 ) then
         f = F0
         df = -g0 - f0 * z
         return
      else
         if ( l-1 <= 0 ) then
            f = w * (f0 - cosix) * z
            df = f0 - 2.0 * f * z
            return
         else 
            f1 = w * (f0 - cosix) * z
            if ( l /= 2 ) then
               j = l-2
               do i = 1, j
                  f2 = w * (f1 * dble(2*i+1) - f0)
                  f0 = F1
                  f1 = f2
               end do
            end if
            f = w * (f1 * dble(2*l-1) * z - f0)
            df = f1 - (l+1) * f * z
            return 
         end if
      end if
   end if
end subroutine besjot

subroutine bessen (l, x, g, dg, r)
   implicit none
   integer, intent(in) :: l
   real*8, intent(in) :: x
   real*8, intent(inout) :: g, r, dg

   integer :: i, j
   real*8 :: a, sinix, cosix, w, pi, xr, z, g0, f0, g1, g2

   g = 0.0
   a = dble(l)
   r = 1.0
   if ( x == 0 ) then
      print *, "ARGUMENT OF SPHERICAL BESSEL-FUNCTION OF SECOND KIND SHOULD NOT BE ZERO"
      return
   else if ( x < 0 ) then
      sinix = sinh(-x)
      cosix = cosh(-x)
      w = -1.0
   else
      pi = 6.283185307179586d0
      xr = mod(x,pi)
      sinix = sin(xr)
      cosix = cos(xr)
      w = +1.d0
   end if
   z = 1.0 / abs(x)
   g0 =-w * cosix * z
   f0 = w * sinix * z
   if ( l <= 0 ) then
      g = g0
      dg = w * f0 - g0 * z
      return
   else
      r = 1.0
      if ( l-1 <= 1 ) then
         g = w * (g0 - sinix) * z
         dg = g0 - 2.0 * g * z
         return 
      else
         g1 = w * (g0 - sinix) * z
         if ( l /= 2 ) then
            j = l-2
            do i = 1, j
               g2 = w * (g1 * dble(2 * i + 1) * z - g0  )
               g0 = g1
               g1 = g2
            end do
         end if
         g = w * (g1 * dble(2*l-1) * z - g0  )
         dg = g1 - dble(l+1) * g * z
      end if
   end if
end subroutine bessen

! *******************************  BESSEL  ***********************************

subroutine besph2 (f, df, g, dg, aa, arg, key, ibug)
! C#######################################################################
! C#    CALCULATES LINEARLY INDEPENDANT SOLUTIONS OF THE EQUATION        #
! C#       2    2        2                                               #
! C#    ( D / DX  - A / X  + 1 ) Y(X) = 0                                #
! C#                       -                                             #
! C#    WHERE A = L*(L+1) WITH L INTEGER .GE. 0                          #
! C#    + (-) SIGN CORRESPONDS TO AN OPEN (CLOSED) CHANNEL               #
! C#    THE SOLUTIONS ARE OBTAINED FROM BESSEL FUNCTIONS OBTAINED        #
! C#    IN BESJOT, BESSEN, BESSIK SUBROUTINES                            #
! C#    SEE ABRAMOWITZ AND STEGUN (CHAP. 10)                             #
! C#    ASYMPTOTIC BEHAVIOUR IS :                                        #
! C#    F # SIN (X-L*PI/2) ; G # -COS (X-L*PI/2) FOR OPEN   CHANNELS     #
! C#    F # SINH (X); G # EXP (-X)               FOR CLOSED CHANNELS     #
! C#---------------------------------------------------------------------#
! C#    AA   : L*(L+1)                                                   #
! C#    ARG  : ARGUMENT VALUE                                            #
! C#           IF POSITIVE THEN ARGUMENT Z IS REAL      (Z = ARG)        #
! C#           IF NEGATIVE THEN ARGUMENT Z IS IMAGINARY (Z = -I*ARG)     #
! C#           WHERE I = (-1)**0.5                                       #
! C#    F,G,DF,DG : REGULAR AND IRREGULAR FUNCTIONS AND THEIR DERIVATIVES#
! C#    KEY  : .LT.0 TO SUPPRESS EXPONENTIAL FACTORS IN BESSIK           #
! C#    IBUG :  >  0 TO PRINT THE OUTPUT                                 #
! C#---------------------------------------------------------------------#
! C#    HAVE BEEN TESTED ON WRONSKIAN RELATION W(F,G) = F*DG-DF*G = 1    #
! C#    FOR THE FOLLOWING VALUES OF THE ARGUMENT AND ORDERS :            #
! C#    0 =< L  < 100  AND 0 =< Z  < 100   (ERROR IN W LESS THAN 10**-12)#
! C#    0 =< L  < 30   AND 0 =< Z  < 100*I (ERROR IN W LESS THAN 10**-10)#
! C#    OUTSIDE THIS RANGE CHECK FOR UNDERFLOWS, OVERFLOWS, DIVIDE CHECKS#
! C#    AND DEXP CAPACITY.                                               #
! C#---------------------------------------------------------------------#
! C#    J.M.L. 08/1981 ; UPDATE : 12/1981                                #
! C#    J.M.LAUNAY, MEUDON, FRANCE                                       #
! C#    J.J.SIERRA, MADRID, SPAIN, 2026 : updated to Fortran 2008        #
! C#######################################################################

   implicit none
   real*8, intent(in) :: aa, arg
   real*8, intent(inout) :: f, df, g, dg
   integer, intent(in):: key, ibug

   real*8 :: ff(2), dff(2), pi, x, sk2, twopim, bi, bk, bim, bkm, bip, bkp
   real*8 :: bj, dbj, r, bn, dbn, rg, vv, wm1
   integer :: l, lp, lm, i

   pi = acos(-1.0)

   if ( aa > -0.25 ) then
      l = -0.5 + sqrt(aa + 0.25)
   else
      print "(' ****** BESPH2 ROUTINE; AA = ',F12.2,' IS NOT L*(L+1) WITH L INTEGER .GE. 0; REQUEST ABORTED ******')", aa
   end if
   
   x = abs(arg)
   sk2 = -1.0
   if ( arg <= 0.0 ) then
      twopim = 2.0/pi
      lp = l+1
      lm = l-1
      call bessik(l,  x, bi,  bk,  key)
      call bessik(lm, x, bim, bkm, key)
      call bessik(lp, x, bip, bkp, key)
      f  = -bi * x
      df = -((l * bim + lp * bip) / (l + lp) * x + bi)
      g  = bk * x * twopim
      dg = (-(l * bkm + lp * bkp) / (l + lp) * x + bk) * twopim
   else
      call besjot(l, x, bj, dbj, r)
      call bessen(l, x, bn, dbn, rg)
      sk2 = 1.0
      if ( l /= 0 ) then
         do i = 1, l
            bj  = bj / r
            dbj = dbj / r
            bn  = bn * rg
            dbn = dbn * rg
         end do
      end if
      f = x * bj
      g = x * bn
      df = x * dbj + bj
      dg = x * dbn + bn
   end if

   ff(1) = f
   ff(2) = g
   dff(1) = df
   dff(2) = dg
   vv = aa / (x*x) - sk2

   if ( ibug > 0 ) then
      wm1 = f * dg - df * g - 1.0
      print "(' BESPH2 ',1F12.4,1F12.4,1P,4D20.12,1D9.1)", aa, arg, f, df, g, dg, wm1
   end if
end subroutine besph2

subroutine bessik (l, x, bi, bk, key)
! C#######################################################################
! C#    CALCULATION OF MODIFIED SPHERICAL BESSEL FUNCTIONS OF THE THIRD  #
! C#    KIND  :                                                          #
! C#    (PI/(2*X))**1/2 * I     (X)  BY FORMULA 10.2.5                   #
! C#                       L+1/2     FROM THE HANDBOOK (P.443)           #
! C#    (PI/(2*X))**1/2 * K     (X)  BY FORMULA 10.2.15                  #
! C#                       L+1/2     FROM THE HANDBOOK (P.444)           #
! C#---------------------------------------------------------------------#
! C#    L     : ORDER (INTEGER)                                          #
! C#    X     : ARGUMENT (POSITIVE REAL NUMBER)                          #
! C#    BI,BK : BESSEL FUNCTIONS  ;                                      #
! C#    KEY  =0: NORMAL, >0: MULTIPLIES BI BY EXP(X), BK BY EXP(-X)      #
! C#                     <0:     ''     BI BY EXP(-X), BK BY EXP(X)      #
! C#######################################################################
   implicit none
   real*8 fct
   common /lgfac/fct(50000)

   real*8, intent(in) :: x
   real*8, intent(inout):: bi, bk 
   integer, intent(in):: l, key

   real*8 :: pi, factor, s2, sx, ci, shx2, alfa, beta, sum, hox, shox, shox1, st, ssum, term
   integer :: kmin, kmax, k, lk

   pi = acos(-1.0)
   factor = log(1.d-30)

   call faclg
   bi = 0.0
   bk = 0.0
   if ( l < 0 ) then
      return
   end if
   s2 = log(2.0)
   sx = log(x)
   hox = 0.5 / x
   ci = l * sx
   shx2 = log(0.5) + 2.0 * sx
   alfa = 0.0
   if ( key < 0 ) then
      alfa = x
   else if ( key > 0) then
      alfa = -x
   end if
   beta = x - alfa

   ! Computation of bk
   sum = hox * exp(-beta)
   if ( l /= 0 ) then
      shox = log(hox)
      shox1 = 0.0
      kmin = 0
      kmax = l
      sum = 0.d0
      do k = kmin, kmax
         shox1 = shox + shox1
         st = fct(l+k+1) - fct(k+1) - fct(l-k+1) + shox1
         if ( st < -180.0 ) then
            exit
         end if
         sum = sum + exp(st - beta)
      end do
   end if
   bk = pi * sum 

   ! Computation of bi
   sum = 0.0
   ssum = 0.0
   k = 0
   lk = l

   do while (ssum + alfa < factor .OR. ssum == 0.0 .OR. term/sum > 0)
      st = k * shx2 - fct(k+1) - (fct(2*lk+2) - lk * s2 - fct(lk+1) )
      k  = k + 1
      lk = l + k
      if ( 2 * lk + 2 > 5000 ) then
         print *,'***> error in bessik: max. value of factorial = ', 2*lk+2
         return
      end if
      term = exp(st - alfa + ci)
      sum = sum + term
      ssum = 0.0
      if ( sum > 0.0 ) then
         ssum = log(sum)
      end if
   end do
   
   bi = sum
end subroutine bessik

! C***********************************************************************
SUBROUTINE JORDAN(MU,LAMBDA,X,N,B)
! C***********************************************************************
   IMPLICIT DOUBLE PRECISION(A-H,L-M,O-Z)
! C
   PARAMETER(NX2MAX=60000)
! C
   DIMENSION MU(N),LAMBDA(N),X(N),B(N)
   DIMENSION PIV(NX2MAX)
! C
      IF(N > NX2MAX) GO TO 999
! C
! C
! C     CALCUL DES PIVOTS
      PIV(1)=2.D0
      DO 10 I=2,N
      PIV(I)=2.D0-LAMBDA(I)*MU(I-1)/PIV(I-1)
   10 B(I)=B(I)-LAMBDA(I)/PIV(I-1)*B(I-1)
! C
      X(N)=B(N)/PIV(N)
      I=N-1
   20 X(I)=(B(I)-X(I+1)*MU(I))/PIV(I)
      I=I-1
      IF(I > 0) GOTO 20
      RETURN
! C
 999  CONTINUE
      PRINT 9000
      PRINT 9999, NX2MAX, N
      STOP
 9000 FORMAT(//,2X,20('*'),' STOP IN JORDAN ',20('*'),/)
 9999 FORMAT(2X,'DIMENSION PARAMETER NX2MAX TOO SMALL , ',I5,'REQUIRED')
END
! C***********************************************************************
DOUBLE PRECISION FUNCTION DLAGRA(X,Y,MIN,IP)
! C***********************************************************************
   IMPLICIT DOUBLE PRECISION(A-H,O-Z)
   DIMENSION X(MIN),Y(MIN)
   DLAGRA=0.D0
   DO 10 I=1,MIN
   IF(I == IP) GOTO 10
   YP=Y(I)
   DO 20 J=1,MIN
      IF(J == IP) GOTO 20
      IF(J == I) GOTO 20
         YP=YP*(X(IP)-X(J))
20 CONTINUE
   DO 30 J=1,MIN
      IF(J == I) GOTO 30
         YP=YP/(X(I)-X(J))
30 CONTINUE
   DLAGRA=DLAGRA+YP
10 CONTINUE
   DO 40 I=1,MIN
      IF(I == IP) GOTO 40
         DLAGRA=DLAGRA+Y(IP)/(X(IP)-X(I))
40 CONTINUE
   RETURN
   END

subroutine faclg
! C#######################################################################
! C#    INITIALISATION OF LOGARITHMS OF FACTORIALS ARRAY                 #
! C#######################################################################
   implicit none
   real*8  :: fct
   integer, save :: ntimes = 0
   integer :: i
   common /lgfac/fct(50000)

   if ( ntimes == 0 ) then
      ntimes = ntimes + 1
      fct(1) = 0.0
      do i = 1, 49999
         fct(i+1) = fct(i) + log(dble(i+1))
      end do      
   end if
end subroutine faclg

!----------------------------------------------------------------
subroutine sinmom(box,npun,npundim,xmred,hbr,pr,p2r)
   implicit real*8(a-h,o-z)

   dimension pr(npundim),p2r(npundim)

   if(npun > npundim)then
      write(6,*)'  npun= ',npun, ' > npundim = ',npundim,' in fftmom'
      stop
   endif   
      
   pi = dacos(-1.d0)
   dpi = 2.d0*pi
   hbrxm = 0.5d0*hbr*hbr/xmred
   ah = (box)/dble(npun-1)
   box = 2.d0*(dble(npun)+1.d0)*ah
   do ir=1,npundim
      pr(ir)=0.d0
      p2r(ir)=0.d0
   enddo

   do ir=1,npun
      iii=ir
      p = dpi*dble(iii)/box
      pr(ir) = p
      p2r(ir) = p*p*hbrxm 
   enddo

   return
end
subroutine noptFFT(nin,nout,ntot)
   integer nin,nout,ntot,i,nmax,imin,idis
   integer nexp2,nexp3,nexp5,nexp7,nexp11
   integer n2,n3,n5,n7,n11
   integer i2,i3,i5,i7,i11

   nmax=10000
   if(nin > nmax.or.ntot > nmax)then
         write(6,*)' nin= ',nin,' , ntot= ',ntot,'  > nmax= ',nmax
         write(6,*)'   in noptFFT of libdyn library'
         stop
   endif

   nout=1
   i=1
   
   nexp2=1
   do n2=0,100
      nexp2=nexp2*2
      if(nexp2 > ntot)go to 2
   enddo
2    continue
   nexp2=n2

   nexp3=1
   do n3=0,100
      nexp3=nexp3*3
      if(nexp3 > ntot)go to 3
   enddo
3    continue
   nexp3=n3

   nexp5=1
   do n5=0,100
      nexp5=nexp5*5
      if(nexp5 > ntot)go to 5
   enddo
5    continue
   nexp5=n5

   nexp7=1
   do n7=0,100
      nexp7=nexp7*7
      if(nexp7 > ntot)go to 7
   enddo
7    continue
   nexp7=n7

   nexp11=1
   do n11=0,100
      nexp11=nexp11*11
      if(nexp11 > ntot)go to 11
   enddo
11   continue
   nexp11=n11

   i=1
   imin=ntot

   do n11=0,nexp11
      i11=11**n11
      do n7=0,nexp7
         i7=7**n7
         do n5=0,nexp5
            i5=5**n5
            do n3=0,nexp3
               i3=3**n3
               do n2=0,nexp2
                  i2=2**n2
                  i=i2*i3*i5*i7*i11
                  idis=i-nin                     
                  if(idis.ge.0.and.idis.lt.imin)then
                        nout=i
                        imin=idis
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo

   if(nout > ntot)nout=ntot

   return
end
!--------------------------------------------------------------------
subroutine fftmom(box,npun,npundim,xmred,hbr,pr,p2r)
   implicit real*8(a-h,o-z)

   dimension pr(npundim),p2r(npundim)

   if(npun > npundim)then
      write(6,*)'  npun= ',npun, ' > npundim = ',npundim,' in fftmom'
      stop
   endif   
      
   pi = dacos(-1.d0)
   dpi = 2.d0*pi
   hbrxm = 0.5d0*hbr*hbr/xmred
   ah = (box)/dble(npun-1)
   box = ah*dble(npun)
   do ir=1,npundim
      pr(ir)=0.d0
      p2r(ir)=0.d0
   enddo

   do ir=1,npun
      if(ir.le.(npun/2))then
         iii=ir-1
      else
         iii=ir-npun-1
      endif
      p = dpi*dble(iii)/box
      pr(ir) = p
      p2r(ir) = p*p*hbrxm 
   enddo

   return
end
! *************************  tresj  **************************

subroutine tresj(j1,j2,j3,m1,m2,m3,coef)
   implicit real*8 (a-h,o-z)

   common/fct/fact(0:10000)

   coef=0.d0

   if(m1+m2+m3 == 0 .and. trian(j1,j2,j3) /= 0.d0)then
         b=fact(j1+m1)+fact(j1-m1)+fact(j2+m2)+fact(j2-m2)+fact(j3+m3)+fact(j3-m3)
         b=0.5d0*b

         k1=j3-j2+m1
         k2=j3-j1-m2
         k3=j1+j2-j3
         k4=j1-m1
         k5=j2+m2

         kmin=max0( -k1 , -k2 , 0 )
         kmax=min0(  k3 ,  k4 , k5 )

         isign=-1
         if(mod(kmin,2) == 0)isign=1

         do k=kmin, kmax 
            a=fact(k)+fact(k3-k)+fact(k4-k)+fact(k5-k)+fact(k1+k)+fact(k2+k)
            coef=coef+isign*dexp(b-a)
            isign=isign*(-1) 
         end do

         isign=-1
         if(mod(j1-j2-m3,2) == 0) isign=1
         coef = coef * trian(j1,j2,j3) * isign
   endif
end
! *******************************  NPLEGM  ********************************

subroutine nplegm(pm,lmax,m,y,ndim)
   implicit real*8(a-h,o-z)
! *     ******************************************
! *     **  normalized                          **
! *     **     Associated Legedre functions     **
! *     **                                      **
! *     ** Recursion formula (Num.Rec. pg.182)  **
! *     **   P_m^m (y) =(2m-1)!! (1-y^2)^(m/2)  **
! *     **         P_{m-1}^m (y)=0              **
! *     **   (l-m)P_l^m = y(2l-1)P_{l-1}^m      **
! *     **              - (l+m-1)P_{l-2}^m      **
! *     **                                      **
! *     **  Input:                              **
! *     **      lmax: maximum l value desired   **
! *     **      m   : m value                   **
! *     **      y   : argument (cos(theta)      **
! *     **  Output:                             **
! *     **      pm  : array with the values..   **
! *     **            pm(0)= Y_0^m              **
! *     **            pm(1)= Y_1^m              **
! *     ****************************************** 

   dimension pm(0:ndim)

   if(lmax > ndim)then
         write(6,*)' ** be carefull with dimensions in nplegm **'
         write(6,*)'     lmax > ndim',lmax,ndim
         stop
   endif
   if(dabs(y) > 1.d0)then
         write(6,*)' ** Bad argument in plegm **'
         stop
   elseif(m.lt.0)then
         write(6,*)' ** m.lt.0 in plegm **'
         stop
   endif 

   do l=0,lmax
      pm(l)=0.d0
   enddo
   faclog=0.d0

! ***>> Compute P_m^m and P_{m+1}^m

   if(m == 0)then
      pm(m)=1.d0
      pm(m+1)=y
   else
      faclog=0.d0
      do i=1,2*m-1,2
         faclog=faclog+dlog(dble(i))
      enddo
         xxx=dsin(dacos(y))
         xxx=xxx**(dble(m))
         xxx=xxx*(-1.d0)**(m)

! c         xxx=(1-y*y)
! c         xxx=xxx**(0.5d0*dble(m))
      pm(m)=xxx
      pm(m+1)=y*dble(2*m+1)*pm(m)
   endif

! ***>> Compute P_l^m  l=m+2,...,lmax

   do l=m+2,lmax
      xlm=dble(l-m)
      pm(l)=y*dble(2*l-1)*pm(l-1)-dble(l+m-1)*pm(l-2)
      pm(l)=pm(l)/xlm
   enddo

! ***>> normalization

   do l=0,lmax
      facden=0.d0
      if(l+m > 0)then
         do i=1,l+m
            facden=facden+dlog(dble(i))
         enddo
      endif
      facnum=0.d0
      if(l-m > 0)then
         do i=1,l-m
            facnum=facnum+dlog(dble(i))
         enddo
      endif
      factor=dexp(0.5d0*(facnum-facden)+faclog)
      if(m == 0)factor=1.d0
      factor=factor*dsqrt(dble(l)+0.5d0)
      pm(l)=pm(l)*factor
   enddo

   return
end

! *************************  coefread  ***************************

subroutine coefread(ifile,nmin,nmax,ntot,nbaslie,cdis,jd,ld,nv1,nv2,eee)
   implicit real*8(a-h,o-z)
   character*40 fcoef

   dimension cdis(nbaslie,nmin:nmax)
   dimension jd(nbaslie),ld(nbaslie),nv1(nbaslie),nv2(nbaslie)

! c      open(ifile,file=fcoef,status='old')
      read(ifile,*)numvec
      if(numvec.lt.nmax)then
         write(6,*)'  ** In coefread: no. of eigenvectors= ',numvec
         write(6,*)'      while desired state are between= ',nmin,nmax
         stop 
      endif

      read(ifile,*)nnn,ntot
      if(ntot > nbaslie)then
         write(6,*)'  ** In coefread: no. of basis fucntions= ',ntot
         write(6,*)'     larger than dimension nbaslie= ',nbaslie
         stop 
      endif
   
! * Reading basis set quantum numbers

      do i=1,ntot
         read(ifile,*)nv1(i),ld(i),nv2(i),iii,jd(i)
      enddo

! * Reading coefficients

      do ibound=1,nmax
         read(ifile,*)kkkk,eee
         write(6,*)'    ',ibound,kkkk,' bound state energy = ',eee
         call flush(6)
         do i=1,ntot
            read(ifile,*)ccc
            if(ibound.ge.nmin)cdis(i,ibound)=ccc
         enddo
         if(ibound.ge.nmin)then
            write(6,"('First and last coef. of wvf =',e15.7,2x,e15.7)")cdis(1,ibound),cdis(ntot,ibound)
         endif

! * Distributions of bound  state

      if(ibound.ge.nmin)then
         write(6,*)'  v distribution '
         m1=nv1(1)
         m2=nv1(ntot)
         do im= m1,m2
            xxx=0.d0
            do i=1,ntot
            if(nv1(i) == im)xxx=xxx+cdis(i,ibound)*cdis(i,ibound)
            enddo
            write(6,*)'      v=',im,' -->',xxx
         enddo
         write(6,*)'  j distribution '
         m1=jd(1)
         m2=jd(ntot)
         do im= m1,m2
            xxx=0.d0
            do i=1,ntot
            if(jd(i) == im)xxx=xxx+cdis(i,ibound)*cdis(i,ibound)
            enddo
            write(6,*)'      j=',im,' -->',xxx
         enddo
         write(6,*)'  n distribution '
         m1=nv2(1)
         m2=nv2(ntot)
         do im= m1,m2
            xxx=0.d0
            do i=1,ntot
            if(nv2(i) == im)xxx=xxx+cdis(i,ibound)*cdis(i,ibound)
            enddo
            write(6,*)'      n=',im,' -->',xxx
         enddo

         endif
      enddo

! c      close(ifile)

      return
end
! ******************************  rfunread  **********************************

subroutine rfunread(ifile,npunr,nvini,nvmax,fread,nuno,fspl,npun,rmis,ah,f,x)
      implicit real*8(a-h,o-z)
      dimension fread(npunr,nvini:nvmax),fspl(npun,nvini:nvmax)
      dimension f(npunr,2),x(npunr)

      CONVL=.52917726D0
      
      write(6,*)ifile,npunr,nvini,nvmax,nuno,npun,rmis,ah
      call flush(6) 
      read(ifile,*)nnn,nuno,nvinir,nvmaxr
      write(6,*)nnn,nuno,nvinir,nvmaxr
      call flush(6)
      if(nuno == 0)then
         if(nnn > npunr)then
            write(6,*)' *** Attention in nfunread ifile= ',ifile,' ***'
            write(6,*)'   number of points in file ',ifile,' = ',nnn
            write(6,*)'   larger than npunlie= ',npunr
            stop
         endif
         if(nvini > nvinir.or.nvmaxr > nvmax)then
            write(6,*)' *** Attention in nfunread ***'
            write(6,*)'   number of functions in file',ifile,'= ',nvinir,nvmaxr
         write(6,*)' are not compatible with dimensions= ',nvini,nvmax
            stop
         endif
      
         do ir=1,nnn
            read(ifile,*)x(ir),(fread(ir,iv),iv=nvinir,nvmaxr)
            x(ir)=x(ir)*convl
         enddo                  
      
         do iv=nvinir,nvmaxr
            xnorm=0.d0
            do ir=1,nnn
               f(ir,1)=fread(ir,iv)/dsqrt(convl)
               f(ir,2)=0.d0
            enddo

            call splset(f,x,nnn,npunr)
            iold=2
            do ir=1,npun
               r=rmis+dble(ir-1)*ah
               if(r.lt.x(1))then
                  fspl(ir,iv)=0.d0
               elseif(r > x(nnn))then
                  fspl(ir,iv)=0.d0
               else
                  call splinqq(f,x,iold,nnn,r,npunr,spl)
! c                  fspl(ir,iv)=splinq(f,x,iold,nnn,r,npunr)
                  fspl(ir,iv)=spl
               endif
               xnorm=xnorm+fspl(ir,iv)*fspl(ir,iv)
            enddo
            write(6,*)iv,xnorm*ah
         
         enddo
   endif

   return
end
! c********************  trian   ****************************

function trian(j1,j2,j3)
   implicit real*8 (a-h,o-z)

   common/fct/fact(0:10000)

   trian=0.d0

   if(j3.ge.iabs(j1-j2).and.j3.le.(j1+j2))then
      cc = fact(j1+j2-j3)+fact(j1-j2+j3)+fact(-j1+j2+j3)-fact(j1+j2+j3+1)
      trian = dexp(cc/2.d0)
   endif

   return
end 

! ************************************************************************
! This function is never called in the code base
double precision function tresjd(j1d,j2d,j3d,m1d,m2d,m3d)
   implicit real*8(a-h,o-z)

   fj1=dble(j1d)*0.5d0
   fj2=dble(j2d)*0.5d0
   fj3=dble(j3d)*0.5d0
   fm1=dble(m1d)*0.5d0
   fm2=dble(m2d)*0.5d0
   fm3=dble(m3d)*0.5d0
   
   tresjd=f3j(fj1,fj2,fj3,fm1,fm2,fm3)

   return
end function tresjd
!---------------------------------------------------------
double precision function seisjd(j1d,j2d,j3d,j4d,j5d,j6d)
   implicit real*8(a-h,o-z)

   fj1=dble(j1d)*0.5d0
   fj2=dble(j2d)*0.5d0
   fj3=dble(j3d)*0.5d0
   fj4=dble(j4d)*0.5d0
   fj5=dble(j5d)*0.5d0
   fj6=dble(j6d)*0.5d0
   
   seisjd=f6j(fj1,fj2,fj3,fj4,fj5,fj6)

   return
end function seisjd

! c-----------------------------------------------------------------------
double precision function f3j (fj1,fj2,fj3, fm1,fm2,fm3)
! c#######################################################################
! c#    calculates 3j coefficients from racah formula                    #
! c#    (messiah: t2, p 910; formula 21) .                               #
! c#    clebsch-gordan coefficients are given by (p. 908, formula 12) :  #
! c#                         j -j +m                |j    j     j|       #
! c#    <j j m m |j m> = (-1) 1  2   (2*j+1)**(0.5) | 1    2     |       #
! c#      1 2 1 2                                   |m    m    -m|       #
! c#                                                | 1    2     |       #
! c#######################################################################
      implicit double precision (a-h,o-z)
      integer t,tmin,tmax
      parameter (nfctmx=5001)
      data tiny,zero,one /0.01d0,0.d0,1.d0/ ,ntimes /1/
      common /cfaclog/ fct(nfctmx)
      if (ntimes  ==  1) call faclog
      ntimes = ntimes+1
      cc = zero
      if (fj3 .gt. (fj1+fj2+tiny))      go to 100
      if (dabs(fj1-fj2)  >  (fj3+tiny)) go to 100
      if (dabs(fm1+fm2+fm3)  >  tiny)   go to 100
      if (dabs(fm1)  >  (fj1+tiny))     go to 100
      if (dabs(fm2)  >  (fj2+tiny))     go to 100
      if (dabs(fm3)  >  (fj3+tiny))     go to 100
      fk1 = fj3-fj2+fm1
      fk2 = fj3-fj1-fm2
      fk3 = fj1-fm1
      fk4 = fj2+fm2
      fk5 = fj1+fj2-fj3
      fk1m = fk1-tiny
      fk2m = fk2-tiny
      fk1p = fk1+tiny
      fk2p = fk2+tiny
      if (fk1m .lt. zero) k1 = fk1m
      if (fk1p  >  zero) k1 = fk1p
      if (fk2m .lt. zero) k2 = fk2m
      if (fk2p  >  zero) k2 = fk2p
      k3 = fk3+tiny
      k4 = fk4+tiny
      k5 = fk5+tiny
      tmin = 0
      if (k1+tmin .lt. 0) tmin = -k1
      if (k2+tmin .lt. 0) tmin = -k2
      tmax = k3
      if (k4-tmax .lt. 0) tmax = k4
      if (k5-tmax .lt. 0) tmax = k5
      n1 = fj1+fj2-fj3+one+tiny
      n2 = fj2+fj3-fj1+one+tiny
      n3 = fj3+fj1-fj2+one+tiny
      n4 = fj1+fm1+one+tiny
      n5 = fj2+fm2+one+tiny
      n6 = fj3+fm3+one+tiny
      n7 = fj1-fm1+one+tiny
      n8 = fj2-fm2+one+tiny
      n9 = fj3-fm3+one+tiny
      n10 = fj1+fj2+fj3+2.d0+tiny
      x = fct(n1)+fct(n2)+fct(n3)+fct(n4)+fct(n5)+fct(n6)+fct(n7)+fct(n8)+fct(n9)-fct(n10)
      x = 0.5d0*x
      do 10  t = tmin,tmax
	 phase = one
	 if (mod(t,2) .ne. 0) phase = -one
	 cc = cc+phase*dexp(-fct(t+1)   -fct(k1+t+1)-fct(k2+t+1)-fct(k3-t+1)-fct(k4-t+1)-fct(k5-t+1)+x)
 10   continue
      fsp = dabs(fj1-fj2-fm3)+tiny
      ns = fsp
      if (mod(ns,2)  >  0) cc = -cc
 100  f3j = cc
      return
      end
! c-----------------------------------------------------------------------
double precision function f6j (fj1,fj2,fj3,fl1,fl2,fl3)
! c#######################################################################
! c#    calculation of 6j-coefficients                                   #
! c#######################################################################
      implicit double precision (a-h,o-z)
      parameter (nfctmx=5001)
      common /cfaclog/ fct(nfctmx)
      data tiny /.01/ ,ntimes /1/
! c
      if (ntimes  ==  1) call faclog
      ntimes = ntimes+1
      d = fdelta (fj1,fj2,fj3)
      d = d*fdelta (fj1,fl2,fl3)
      d = d*fdelta (fl1,fj2,fl3)
      d = d*fdelta (fl1,fl2,fj3)
      f6j = 0.d0
      if (dabs(d)  ==  0.d0) return
! c
      fk1 = fj1+fj2+fj3
      fk2 = fj1+fl2+fl3
      fk3 = fl1+fj2+fl3
      fk4 = fl1+fl2+fj3
      fk5 = fj1+fj2+fl1+fl2
      fk6 = fj2+fj3+fl2+fl3
      fk7 = fj3+fj1+fl3+fl1
      fmin = dmin1 (fk5,fk6,fk7)
      fmax = dmax1 (fk1,fk2,fk3,fk4)
      min = fmin+tiny
      max = fmax+tiny
      k1 = fk1+tiny
      k2 = fk2+tiny
      k3 = fk3+tiny
      k4 = fk4+tiny
      k5 = fk5+tiny
      k6 = fk6+tiny
      k7 = fk7+tiny
      if (min-max) 1000,3,3
 3    if (max) 1000,4,4
 4    if (min) 1000,5,90
 5    k1 = -k1
      k2 = -k2
      k3 = -k3
      k4 = -k4
      bot = fct(k1+1)+fct(k2+1)+fct(k3+1)+fct(k4+1)+fct(k5+1)+fct(k6+1)+fct(k7+1)
      bot = dexp(bot)
      f6j = d/bot
      return
! c
 90   f6j = 0.
      do 100 i = max,min
	 boite = ((-1.)**i)
	 iz = i+1
	 m1 = i-k1
	 m2 = i-k2
	 m3 = i-k3
	 m4 = i-k4
	 m5 = k5-i
	 m6 = k6-i
	 m7 = k7-i
	 dot = fct(iz+1)
	 bot = fct(m1+1)+fct(m2+1)+fct(m3+1)+fct(m4+1)+fct(m5+1)+fct(m6+1)+fct(m7+1)
	 b1 = dot-bot
	 boite = boite*dexp(b1)
	 f6j = f6j+boite
 100  continue
      f6j = f6j*d
 1000 return
! c
end
! c-----------------------------------------------------------------------
double precision function fdelta (fl1,fl2,fl3)
   implicit double precision (a-h,o-z)
   parameter (nfctmx=5001)
   common /cfaclog/ fct(nfctmx)
   data   eps /.01/
! c
   ia=fl1+fl2+fl3 +eps
   a=2.*(fl1+fl2+fl3)+1.
   ib=a +eps
   ib=ib/2
   if(ib-ia)1,6,1
   6 continue
   ik1=fl1+fl2-fl3+eps
   ik2=fl2+fl3-fl1+eps
   ik3=fl3+fl1-fl2+eps
   kk=fl1+fl2+fl3+1+eps
   if(ik1)1,2,2
   2 if(ik2)1,3,3
   3 if(ik3)1,4,4
   4 d1=fct(kk+1)
   d2=fct(ik1+1)+fct(ik2+1)+fct(ik3+1)
   d3 = (d2 - d1) / 2.d0
   fdelta = dexp (d3)
   go to 5
   1 fdelta=0.
   5 return
! c
end
      
subroutine faclog
! c#######################################################################
! c#    initialisation of logarithms of factorials array                 #
! c#######################################################################
   implicit double precision (a-h,o-z)
   parameter (nfctmx=5001)
   common /cfaclog/ fct(nfctmx)
   data ntimes /0/
! c
   ntimes = ntimes+1
   if (ntimes  >  1) return
   fct(1) = 0.d0
   do 10 i = 1,nfctmx-1
   ai = i
   fct(i+1) = fct(i)+dlog(ai)
10   continue
! c
   return
end

! ***********************************************************************

subroutine bndbcele(Eval,fun,potmatrix,xmu,rmis,rfin,nvini,nvmax,j,npun,nelec)
   implicit real*8(a-h,o-z)

! *
! * calculate the bound state of a diatomic system
! * for several coupled electronic states (nelec)
! *   for a particular value of the angular momentum j
! *
! * uses a DVR representation of particles in a box
! *
! *  uses a.u.
! * xmu: reduced mass
! *
   parameter(npunaux=1024,nelecaux=10,ntotaux=nelecaux*npunaux)

   dimension Eval(nvini:nvmax),fun(npun,nelec,nvini:nvmax)      
   dimension potmatrix(npun,nelec,nelec)

   dimension nr(ntotaux),ne(ntotaux),Hmat(ntotaux,ntotaux)
   dimension eigen(ntotaux),T(ntotaux,ntotaux)
   dimension wwork(5*ntotaux)
   dimension ind(ntotaux*5)

! * checking dimensions

   if(npun > npunaux)then
      write(6,*)' npun= ',npun,' > npunaux= ',npunaux
      write(6,*)'  change it in bndbcele '
      stop
   endif

   if(nelec > nelecaux)then
      write(6,*)' nelec= ',nelec,' > nelecaux= ',nelecaux
      write(6,*)'  change it in bndbcele '
      stop
   endif

! * forming basis 

   ii=0
   do ie=1,nelec
   do ir=1,npun
      ii=ii+1
      nr(ii)=ir
      ne(ii)=ie
   enddo
   enddo
   ntot=nelec*npun
   ah=(rfin-rmis)/dble(npun-1)

   xl=0.5d0*dble(j*(j+1))/xmu

   pi=dacos(-1.d0)
   xkinfac=pi*0.5d0/ ( (rfin-rmis)*xmu)
   hbr=1.d0

! * forming H matrix

   ii=0
   do i=1,ntot
   do ip=i,ntot
      ii=ii+1
      Hmat(i,ip)=0.d0
      ie=ne(i)
      iep=ne(ip)
      ir=nr(i)
      irp=nr(ip)

! * rotational term
      if(ir == irp.and.ie == iep)then
         r=rmis+dble(ir-1)*ah

         Hmat(i,ip)= Hmat(i,ip)+ xl/(r*r)
      endif
! * radial kinetic term
      if(ie == iep)then
         cint=0.d0
         menos=(ir-irp)
         mas=ir+irp
         sign=dble( (-1)**(menos))
         if(ir == irp)then
            x1=pi*pi/6.d0
            x2=1.d0/( dble(mas)*dble(mas) )
            cint=x1+x2
         else
            x1=1.d0/( dble(menos)*dble(menos) )
            x2=1.d0/( dble(mas)*dble(mas) )
            cint=x1+x2
         endif
         cint=sign*cint/(ah*ah)

         Tmat=hbr*hbr*cint/xmu
         
         Hmat(i,ip)=Hmat(i,ip)+Tmat
      endif
! * He+Hso potential terms

      if(ir == irp)then
         Hmat(i,ip)=Hmat(i,ip)+potmatrix(ir,ie,iep)
      endif

   enddo
   enddo

! * diagonalization
!  
! c       call dspev('v','l',ntot,Hmat,eigen,T,ntotaux,wwork,inf)

   call diagon(hmat,ntot,ntotaux,T,eigen)

! * keeping desired eigenstates

   iiv=0
   do iv=nvini,nvmax
      iiv=iiv+1
      Eval(iv)=eigen(iiv)
      do ii=1,ntot
         ir=nr(ii)
         ie=ne(ii)
         fun(ir,ie,iv)=T(ii,iiv)/dsqrt(ah)
      enddo
   enddo

   return
end
! ***********************************************************************

subroutine bndbc1ele(Eval,fun,potmatrix,xm, &
                     rmis,rfin,nvini,nvmax,j,npun,nelec, &
                     vnumber,enumber,max_viblevels,ivreal)
   implicit none

   integer,intent(in) :: nelec,j,npun,nvini,nvmax
   integer,intent(in) :: max_viblevels(nelec)
   integer,intent(inout) :: vnumber(nvini:nvmax)
   integer,intent(inout) :: enumber(nvini:nvmax)
   real*8,intent(inout) :: Eval(nvini:nvmax)
   real*8,intent(inout) :: fun(npun,nelec,nvini:nvmax)      
   real*8,intent(in) :: potmatrix(npun,nelec,nelec)
   real*8,intent(in) :: xm,rmis,rfin
   
   real*8 :: alpha(npun),beta(npun),v(npun),p(npun),vv(npun)
   real*8 :: ah,xl,xz,xz1,eps,a,b,r,e0,e2
   integer :: maxit,ivreal,ie,ir,iv,nchan,kv,itry
   integer :: maxvib_real

   Eval(:)=0.d0
   fun(:,:,:)=0.d0
   vnumber(:)=0
   enumber(:)=0
! *
! * calculate the bound state of a diatomic system
! * for several uncoupled electronic states (nelec)
! *   for a particular value of the angular momentum j
! *
! * uses a equispaced radial representation 
! *
! *  uses a.u.
! * xmu: reduced mass
! *
   
   ah=(rfin-rmis)/dble(npun-1)
   xl=0.5d0*dble(j*(j+1))/xm
   xz=0.5d0/xm
   xz1=2.d0*xm
   eps=1.d-10
   maxit=20
!---  > loop in electronic states

   ivreal=nvini-1
   do ie=1,nelec
      V(:)=0.d0
      VV(:)=0.d0
      alpha(:)=0.d0
      beta(:)=0.d0
      p(:)=0.d0
      
! * effective potential
      do ir=1,npun
         r = rmis + dble(ir-1)*ah
         vv(ir) = potmatrix(ir,ie,ie) + xl/(r*r)
         beta(ir)=1.d0
         alpha(ir)=-2.d0-2.d0*xm*vv(ir)*ah*ah
      enddo

! eigenvalues
      call tqli(alpha,beta,npun,npun)

      do iv=1,npun
         alpha(iv)=-xz*alpha(iv)/ah/ah
      enddo

! * orderig of eigenvalues

44      continue
         nchan=0
         do ir=1,npun-1
            if(alpha(ir) > alpha(ir+1))then
               a=alpha(ir)
               b=alpha(ir+1)
               alpha(ir)=b
               alpha(ir+1)=a
               nchan=nchan+1
            endif
         enddo
      if(nchan > 0)go to 44

!     eigenvalues

      if(max_viblevels(ie).lt.0)then
         maxvib_real=nvmax+1
      else
         maxvib_real=max_viblevels(ie)-1
      endif
      do iv=0,min0(npun,maxvib_real)
         if(potmatrix(npun,ie,ie)-alpha(iv+1) > 1.d-2)then
            
            e0=alpha(iv+1)*xz1
            do ir=1,npun
               v(ir)=vv(ir)*xz1
            end do
            p(:)=0.d0
      call schr(e0,rmis,rfin,npun,maxit,eps,e2,kv,itry,v,p,npun)
            e2=e2/xz1
            ivreal=ivreal+1
            if(ivreal.ge.nvini.and.ivreal.le.nvmax)then

               Eval(ivreal)=e2
               do ir=1,npun
                  fun(ir,ie,ivreal)=p(ir)
               end do
!                 write(6,*)' ielec= ',ie,' ivreal=',ivreal
               
               vnumber(ivreal)=iv
               enumber(ivreal)=ie
            else
               write(6,*)' products with ielec,iv= ',ie,iv,' exceed nvmaxprod= ',nvmax
               write(6,*)' increase nvmaxprod or ',' redefine max_viblevels (ielec) in namelist'
            end if ! Energy below ekinmax
         end if              ! eigval < pot(Rmax)
      end do !iv=1,min(npun,max_viblevels(ie))
   enddo                     ! ielec

   if(ivreal.ne.nvmax)then
      write(6,*)' bndbc1ele:  ivreal=',ivreal
      write(6,*)'   while nvmax= ',nvmax
      write(6,*)' end bndbc1ele'
   endif

return
end

! ********************************* l2mat ***************************************

subroutine l2mat(bfl2mat,facmass,iommin,iommax,j,Jtot,iomdim0,iomdim1)
   implicit real*8(a-h,o-z)
   dimension bfl2mat(iomdim0:iomdim1,iomdim0:iomdim1)

! *   l^2 matrix obtained in a body-fixed representation
! *                  for fixed Jtot and j 
! *            and a limited number of Omega projection's
! *                   using parity adapted functions
! *                           facmass = hbr*hbr*0.5d0/xmred
      
   do iom=iommin,iommax
   do jom=iommin,iommax
      bfl2mat(iom,jom)=0.d0
   enddo
   enddo

   xj=dble(j)

   if(iommin.lt.iomdim0.or.iommax > iomdim1)then
      if(idproc == 0)write(6,*)'   problem with Omega:',' j,iommin,iommax= '  ,j,iommin,iommax,'  in l2mat'
      stop
   endif

   do iom=iommin,iommax
   do jom=iommin,iommax
      if(iom == jom)then
         x1=xj*(xj+1.d0)
         x2=dble(Jtot*(Jtot+1)-2*iom*jom)
         xk=(x1+x2)*facmass
         bfl2mat(iom,jom)=xk
      elseif(iabs(iom-jom) == 1)then
         x1=dble(Jtot*(Jtot+1)-iom*jom)
         x2=xj*(xj+1.d0)-dble(iom*jom)
         x12=dsqrt(x1*x2)
         xk=-x12*facmass
         if(iom*jom == 0)xk=xk*dsqrt(2.d0)
         bfl2mat(iom,jom)=xk
      endif
   enddo
   enddo
   
   return
end 

! ********************************* l2mat_Renner ***************************************

subroutine l2mat_Renner(bfl2mat,facmass,iommin,iommax,j,Jtot,lambdaA,isigma,isigmap,iomdim0,iomdim1)
   implicit real*8(a-h,o-z)
   dimension bfl2mat(iomdim0:iomdim1,iomdim0:iomdim1)

! *   l^2 matrix obtained in a body-fixed representation
! *                  for fixed Jtot and j  for electronic \Lambda=1 (Renner -Teller effect)
! *            and a limited number of Omega projection's
! *                   using parity adapted functions
! *                           facmass = hbr*hbr*0.5d0/xmred
      
   do iom=iommin,iommax
   do jom=iommin,iommax
      bfl2mat(iom,jom)=0.d0
   enddo
   enddo

   xj=dble(j)

   if(iommin.lt.iomdim0.or.iommax > iomdim1)then
      if(idproc == 0)write(6,*)'   problem with Omega:',' j,iommin,iommax= '  ,j,iommin,iommax,'  in l2mat'
      stop
   endif

   if(isigma*isigmap == 1)then
      do iom=iommin,iommax
      do jom=iommin,iommax
         if(iom == jom)then
            x1=xj*(xj+1.d0)+dble(iom-lambdaA)**2
            x2=dble(Jtot*(Jtot+1))
            xk=(x1+x2)*facmass
            bfl2mat(iom,jom)=xk
         elseif(iabs(iom-jom) == 1)then
            x1=dble(Jtot*(Jtot+1)-iom*jom)
            x2=xj*(xj+1.d0)-dble((iom-lambdaA)*(jom-lambdaA))
            x12=dsqrt(x1*x2)
            xk=-x12*facmass
            if(iom*jom == 0)xk=xk*dsqrt(2.d0)
            bfl2mat(iom,jom)=xk
         endif
      enddo
      enddo
   elseif(isigma*isigmap == -1)then
      do iom=iommin,iommax
      do jom=iommin,iommax
         if(iom == jom)then
            x1=-dble(2*lambdaA*lambdaA)
               xk=x1*facmass
            bfl2mat(iom,jom)=xk
         endif
      enddo
      enddo
   endif

   return
end 

! ********************** rgauscolini ***********************************
subroutine rgauscolini(rgaus,r,rcol,alpha,xk,factor,l)
   real*8 rgaus,r,rcol,alpha,xk,factor,pepe,f,g,df,dg,arg
   integer l,key

   pepe=dble(l*(l+1))
   arg=(r-rcol)*xk
! c      key=-1
! c      CALL BESPH2(F,DF,G,DG,PEPE,ARG,KEY,0)
! c      g=-g*dsqrt(2.d0)
   g=dcos(arg)*dsqrt(2.d0)

   rgaus=g*factor*dexp(-(r-rcol)*(r-rcol)/alpha/alpha)

   return
end
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
subroutine diagon(Hmat,n,ndim,T,eigen)
!     diagonalize Hmat  --> providing eigen (eigenvalues) and T (eigenvectors)
!            which are ordered in increasing energy
   implicit none
   integer n,ndim,nrot
   real*8 :: Hmat(ndim,ndim),eigen(ndim),T(ndim,ndim)

   call jacobi(Hmat,n,ndim,eigen,T,nrot)

   call eigsrt(eigen,T,n,ndim)

   return
end subroutine diagon
      
!***************************************************
subroutine jacobi(a,n,np,d,v,nrot)
   implicit none
   integer,parameter :: NMAX=50000
   INTEGER :: n,np,nrot
   REAL*8 :: a(np,np),d(np),v(np,np)
   INTEGER :: i,ip,iq,j
   REAL*8 :: c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
   do 12 ip=1,n
      do 11 iq=1,n
         v(ip,iq)=0.d0
11      continue
      v(ip,ip)=1.d0
12    continue
   do 13 ip=1,n
      b(ip)=a(ip,ip)
      d(ip)=b(ip)
      z(ip)=0.d0
13    continue
   nrot=0
   do 24 i=1,50
      sm=0.d0
      do 15 ip=1,n-1
         do 14 iq=ip+1,n
         sm=sm+dabs(a(ip,iq))
14        continue
15      continue
      if(sm == 0.d0)return
      if(i.lt.4)then
         tresh=0.2d0*sm/n**2
      else
         tresh=0.d0
      endif
      do 22 ip=1,n-1
         do 21 iq=ip+1,n
         g=100.d0*dabs(a(ip,iq))
         if((i > 4).and.(dabs(d(ip))+g == dabs(d(ip))).and.(dabs(d(iq))+g == dabs(d(iq))))then
            a(ip,iq)=0.d0
         else if(dabs(a(ip,iq)) > tresh)then
            h=d(iq)-d(ip)
            if(dabs(h)+g == dabs(h))then
               t=a(ip,iq)/h
            else
               theta=0.5d0*h/a(ip,iq)
               t=1.d0/(dabs(theta)+dsqrt(1.d0+theta**2))
               if(theta.lt.0.d0)t=-t
            endif
            c=1.d0/dsqrt(1+t**2)
            s=t*c
            tau=s/(1.d0+c)
            h=t*a(ip,iq)
            z(ip)=z(ip)-h
            z(iq)=z(iq)+h
            d(ip)=d(ip)-h
            d(iq)=d(iq)+h
            a(ip,iq)=0.d0
            do 16 j=1,ip-1
               g=a(j,ip)
               h=a(j,iq)
               a(j,ip)=g-s*(h+g*tau)
               a(j,iq)=h+s*(g-h*tau)
16            continue
            do 17 j=ip+1,iq-1
               g=a(ip,j)
               h=a(j,iq)
               a(ip,j)=g-s*(h+g*tau)
               a(j,iq)=h+s*(g-h*tau)
17            continue
            do 18 j=iq+1,n
               g=a(ip,j)
               h=a(iq,j)
               a(ip,j)=g-s*(h+g*tau)
               a(iq,j)=h+s*(g-h*tau)
18            continue
            do 19 j=1,n
               g=v(j,ip)
               h=v(j,iq)
               v(j,ip)=g-s*(h+g*tau)
               v(j,iq)=h+s*(g-h*tau)
19            continue
            nrot=nrot+1
         endif
21        continue
22      continue
      do 23 ip=1,n
         b(ip)=b(ip)+z(ip)
         d(ip)=b(ip)
         z(ip)=0.d0
23      continue
24    continue
   stop 'too many iterations in jacdiagon'
   return
END subroutine jacobi
! C  (C) Copr. 1986-92 Numerical Recipes Software *1n#!-013.
!*********************************************************************
subroutine eigsrt(d,v,n,np)
!ordering eigenvalues and eigenvectors
   implicit none
   INTEGER:: n,np
   REAL*8 :: d(np),v(np,np)
   INTEGER :: i,j,k
   REAL*8 :: p

   do 13 i=1,n-1
      k=i
      p=d(i)
      do 11 j=i+1,n
         if(d(j).le.p)then
         k=j
         p=d(j)
         endif
11      continue
      if(k.ne.i)then
         d(k)=d(i)
         d(i)=p
         do 12 j=1,n
         p=v(j,i)
         v(j,i)=v(j,k)
         v(j,k)=p
12        continue
      endif
13    continue
   return
! C  (C) Copr. 1986-92 Numerical Recipes Software *1n#!-013.
END subroutine eigsrt


