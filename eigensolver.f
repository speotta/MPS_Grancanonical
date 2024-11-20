! subroutine for the Davidson diagonalization algorithm of a sparse matrix
      subroutine jdqz (alpha, beta, eivec, wanted,
     $     n, shift, eps, kmax, jmax, jmin,
     $     method, m, l,
     $     maxnmv, maxstep, lock, order, testspace, work,lwork,useguess)
c
c     Programmer: Diederik R. Fokkema
c     Modified         : M. van Gijzen
c     Modified 05-24-96: M. Kooper: ra and rb, the Schur matrices of A and B, 
c              added, as well as the vectors sa and sb, which contain the
c              innerproducts of ra with Z and rb with Z. This is added to be
c              enable to compute the eigenvectors in EIVEC
c     Modification 08-27-96: Different computation of eigenvectors, MvG
c

c      USE WFunction


c     .. Parameters ..
c
      implicit none
      integer gmres, cgstab
      parameter ( gmres = 1, cgstab = 2 )
      integer kmax, jmax, jmin, method, m, l, maxnmv, maxstep, order
      integer testspace
      double precision eps, lock
      complex*16 shift
      integer n, lwork
      complex*16 work(n,*), eivec(n,*), alpha(*), beta(*)
      logical wanted,useguess
c
c     .. Local Parameters ..
c
      logical loop, try, found, ok
      integer ldvs, ldzwork, iseed(4)
      parameter (ldvs=50, ldzwork=4*ldvs)
      complex*16 ma(ldvs,ldvs), mb(ldvs,ldvs),
     $     zma(ldvs,ldvs), zmb(ldvs,ldvs),
     $     vsl(ldvs,ldvs), vsr(ldvs,ldvs),
     $     ra(ldvs,ldvs), rb(ldvs, ldvs),
     $     zwork(4*ldvs), aconv(ldvs), bconv(ldvs)

      integer ldqkz
      parameter (ldqkz=ldvs)
      integer ipivqkz(ldqkz)
      complex*16 mqkz(ldqkz,ldqkz), invqkz(ldqkz,ldqkz), f(ldqkz)

      integer i, j, k, info, mxmv, step
      integer d, u, v, w, aux, av, bv, q, z, kz, itmp, tp
      integer solvestep
      real*8 rnrm, deps    !,rwork(3*ldvs)
      real*8 rwork(8*ldvs) !davide: array  dimension changed
      real*8 dtmp
      complex*16 zalpha, zbeta, targeta, targetb, evcond
      complex*16 shifta, shiftb
      complex*16 zero, one
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))


c     new variables for ZGGES and ZGGEV
      integer sdim
      logical bwork (ldvs)
      logical delctg
      external delctg

c
c     .. Functions ..
c
      real*8 dznrm2
      INTRINSIC DCMPLX,DCONJG
c
c     .. Data ..
c
      data iseed /3,3,1966,29/

c     ... Executable Statements
c
c...     Are there errors in the parameters?
      if ( kmax .gt. jmax )
     $   call error('jdqz: kmax greater than jmax')
      if ( jmax .gt. 50 )
     $   call error('jdqz: jmax greater than 50')
      if ( method .ne. 1 .and. method .ne. 2 )
     $   call error('jdqz: illegal choice for solution method')
      if ( order .lt. -2 .or. order .gt. 2 )
     $   call error('illegal value for order, must be between -2 and 2')

c...     d   = rhs, these pointers refer to the columns of the workspace
      d   = 1
c...     Workspace for jdqzmv
      tp  = d+1
c...     u   = pointer to Krylov space GMRES(m) or Bi-CSTAB(l)
      u   = tp+1
c...     v   = pointer to search space JDQZ with max dimension jmax
      if ( method .eq. gmres ) then
	 v   = u+m+1
      else if ( method .eq. cgstab ) then
	 v   = u+2*l+6
      end if
c...     w   = pointer to test subspace JDQZ with max dimension jmax
      w   = v+jmax
c...     av  = pointer to subspace AV with max dimension jmax
      av  = w+jmax
c...     bv  = pointer to subspace BV with max dimension jmax
      bv  = av+jmax
c...     aux =
      aux = bv+jmax
c...     q   = pointer to search Schur basis in JDQZ with max dimension kmax
      q   = aux+jmax
c...     z   = pointer to test Schur basis in JDQZ with max dimension kmax
      z   = q+kmax
c...     kz  = pointer to matrix K^{-1}Z_k
      kz  = z+kmax
      if (kz+kmax-1.gt.lwork) call error ('qz: memory fault')
c
c     --- initialize loop
c

      ok = .true.

      evcond = dcmplx(sqrt(abs(shift)**2+abs(one)**2))
      shifta = shift/evcond
      shiftb = one/evcond

      targeta = shifta
      targetb = shiftb

      zalpha = shifta
      zbeta = shiftb

      step = 0
      deps = dble(one)
      mxmv = 0

      solvestep = 0

      j = 0
      k = 0
		
c
c     --- loop
c
 100  continue
      loop = (k.lt.kmax.and.step.lt.maxstep)
      if (loop) then
	 step = step+1
	 solvestep = solvestep+1
	 
	 if (j.eq.0) then					
           IF (useguess) THEN
              CALL Guess(n,work(1,v+j))
           ELSE
              call zlarnv(2, iseed, n, work(1,v+j))
           ENDIF

	    call zlarnv(2, iseed, n, work(1,w+j))
	    do i=1,n
	       dtmp = dble(work(i,v+j))
	       work(i,v+j) = dcmplx(dtmp,0d0)
	       dtmp = dble(work(i,w+j))
	       work(i,w+j) = dcmplx(dtmp,0d0)
	    enddo	    
	 else
	    mxmv = maxnmv
	    deps = 2d0**(-solvestep)
	    if (j.lt.jmin) then
	       mxmv = 1
	       call zgmres (n, work(1,v+j), work(1,d), m, deps,
     $              mxmv, zalpha, zbeta, k+1,
     $              work(1,kz), work(1,q), invqkz, ldqkz,
     $              ipivqkz, f, work(1,u), work(1,tp) )
	    elseif (method.eq.gmres) then
	       mxmv = m
	       call zgmres (n, work(1,v+j), work(1,d), m, deps,
     $              mxmv, zalpha, zbeta, k+1,
     $              work(1,kz), work(1,q), invqkz, ldqkz,
     $              ipivqkz, f, work(1,u), work(1,tp) )
	    elseif (method.eq.cgstab) then
	       call zcgstabl (n, work(1,v+j), work(1,d), l,
     $              deps, mxmv, zalpha,
     $              zbeta, k+1, work(1,kz), work(1,q), invqkz,
     $              ldqkz, ipivqkz, f, work(1,u), 2*l+6)
	    endif
	 endif
	 j = j+1
		
	 call zmgs (n, j-1, work(1,v), work(1,v+j-1), 1 )
	 call zmgs (n, k, work(1,q), work(1,v+j-1), 1 )

	 if (testspace.eq.1) then
 	    call jdqzmv (n, work(1,v+j-1), work(1,w+j-1), work(1,tp), 
     $           -dconjg(shiftb), dconjg(shifta))
	 elseif (testspace.eq.2) then
 	    call jdqzmv (n, work(1,v+j-1), work(1,w+j-1), work(1,tp),
     $           -dconjg(zbeta), dconjg(zalpha))
	 elseif (testspace.eq.3) then
 	    call jdqzmv (n, work(1,v+j-1), work(1,w+j-1), work(1,tp),
     $           shifta, shiftb)
	 endif

         call zmgs (n, j-1, work(1,w), work(1,w+j-1), 1 )
 	 call zmgs (n, k, work(1,z), work(1,w+j-1), 1 )

         call amul(n, work(1,v+j-1), work(1,av+j-1))
         call bmul(n, work(1,v+j-1), work(1,bv+j-1))
         
         call makemm (n, j, work(1,w), work(1,av), ma, zma, ldvs)
         call makemm (n, j, work(1,w), work(1,bv), mb, zmb, ldvs)

cccccccccccccccccccccccccccccccccccccccccccccc
c      Use ZGGES instead of ZGEGS
cccccccccccccccccccccccccccccccccccccccccccccc
c         call zgegs ('v', 'v', j, zma, ldvs, zmb, ldvs,
c     $        alpha, beta, vsl, ldvs, vsr, ldvs, zwork,
c     $        ldzwork, rwork, info)
         CALL ZGGES('v', 'v', 'n', delctg, j, zma, ldvs, zmb, ldvs,
     *        sdim, alpha, beta, vsl, ldvs, vsr,ldvs, zwork, ldzwork,
     *        rwork, bwork, info)


         try = .true.
 200     continue
         if (try) then
c
c           --- Sort the Petrov pairs ---
c
            call qzsort (targeta, targetb, j, zma, zmb, vsl, vsr,
     $           ldvs, alpha, beta, order)

            zalpha = zma(1,1)
            zbeta = zmb(1,1)

            evcond = dcmplx(sqrt(abs(zalpha)**2+abs(zbeta)**2))
c
c            --- compute new q ---
c
            call zgemv ('n', n, j, one, work(1,v), n, vsr(1,1),
     $           1, zero, work(1,q+k), 1)
            call zmgs (n, k, work(1,q), work(1,q+k), 1 )
c
c           --- compute new z ---
c
            call zgemv ('n', n, j, one, work(1,w), n, vsl(1,1),
     $           1, zero, work(1,z+k), 1)
            call zmgs (n, k, work(1,z), work(1,z+k), 1 )
c
c           --- Make new qkz ---
c
            call zcopy (n, work(1,z+k), 1, work(1,kz+k), 1)
            call precon (n, work(1,kz+k))
            call mkqkz (n, k+1, work(1,q), work(1,kz), mqkz, invqkz,
     $           ldqkz, ipivqkz)
c
c           --- compute new (right) residual= beta Aq - alpha Bq and
c               orthogonalize this vector on Z.
c
            call jdqzmv (n, work(1,q+k), work(1,d), work(1,tp),
     $           zalpha, zbeta)
            call zmgs (n, k, work(1,z), work(1,d), 0 )

            rnrm = dznrm2 (n, work(1,d), 1)/dble(evcond)

            if (rnrm.lt.lock.and.ok) then
               targeta = zalpha
               targetb = zbeta
               ok = .false.
            endif

            found = (rnrm.lt.eps.and.
     $           (j.gt.1.or.k.eq.kmax-1))
            try =  found

            if (found) then

c              --- increase the number of found evs by 1 ---
               k = k+1

c              --- store the eigenvalue
               aconv(k) = zalpha
               bconv(k) = zbeta

               solvestep = 0
               if (k.eq.kmax) goto 100
               call zgemm ('n', 'n', n, j-1, j, one, work(1,v), n,
     $              vsr(1,2), ldvs, zero, work(1,aux), n)
               itmp = v
               v = aux
               aux = itmp
               call zgemm ('n', 'n', n, j-1, j, one, work(1,av), n,
     $              vsr(1,2), ldvs, zero, work(1,aux), n)
               itmp = av
               av = aux
               aux = itmp
               call zgemm ('n', 'n', n, j-1, j, one, work(1,bv), n,
     $              vsr(1,2), ldvs, zero, work(1,aux), n)
               itmp = bv
               bv = aux
               aux = itmp
               call zgemm ('n', 'n', n, j-1, j, one, work(1,w), n,
     $              vsl(1,2), ldvs, zero, work(1,aux), n)
               itmp = w
               w = aux
               aux = itmp
               j = j-1
               call zlacpy ('a', j, j, zma(2,2), ldvs, ma, ldvs)
               call zlacpy ('a', j, j, ma, ldvs, zma, ldvs)
               call zlacpy ('a', j, j, zmb(2,2), ldvs, mb, ldvs)
               call zlacpy ('a', j, j, mb, ldvs, zmb, ldvs)
               call zlaset ('a', j, j, zero, one, vsr, ldvs)
               call zlaset ('a', j, j, zero, one, vsl, ldvs)
               targeta = shifta
               targetb = shiftb
               ok = .true.
               mxmv = 0
               deps = dble(one)
            else if (j.eq.jmax) then
               call zgemm ('n', 'n', n, jmin, j, one, work(1,v), n,
     $              vsr, ldvs, zero, work(1,aux), n)
               itmp = v
               v = aux
               aux = itmp
               call zgemm ('n', 'n', n, jmin, j, one, work(1,av), n,
     $              vsr, ldvs, zero, work(1,aux), n)
               itmp = av
               av = aux
               aux = itmp
               call zgemm ('n', 'n', n, jmin, j, one, work(1,bv), n,
     $              vsr, ldvs, zero, work(1,aux), n)
               itmp = bv
               bv = aux
               aux = itmp
               call zgemm ('n', 'n', n, jmin, j, one, work(1,w), n,
     $              vsl, ldvs, zero, work(1,aux), n)
               itmp = w
               w = aux
               aux = itmp
               j = jmin
               call zlacpy ('a', j, j, zma, ldvs, ma, ldvs)
               call zlacpy ('a', j, j, zmb, ldvs, mb, ldvs)
               call zlaset ('a', j, j, zero, one, vsr, ldvs)
               call zlaset ('a', j, j, zero, one, vsl, ldvs)
            endif
            goto 200
         endif
         goto 100
      endif

c
c...     Did enough eigenpairs converge?
      kmax = k

      if ( wanted ) then
c...        Compute the Schur matrices if the eigenvectors are
c...        wanted, work(1,tp) is used for temporary storage
c...        Compute RA:
         call zlaset ('l', k, k, zero, zero, ra, ldvs)    
         do i = 1, k
            call amul( n, work(1,q+i-1), work(1,tp) )
            call zgemv( 'c', n, i, one, work(1,z), n, work(1,tp), 1, 
     $                  zero, ra(1,i), 1 )
         end do
c...        Compute RB:
         call zlaset ('l', k, k, zero, zero, rb, ldvs)    
         do i = 1, k
            call bmul( n, work(1,q+i-1), work(1,tp) )
            call zgemv( 'c', n, i, one, work(1,z), n, work(1,tp), 1, 
     $                  zero, rb(1,i), 1 )
         end do

c        --- The eigenvectors RA and RB  belonging to the found eigenvalues
c            are computed. The Schur vectors in VR and VS are replaced by the
c            eigenvectors of RA and RB

cccccccccccccccccccccccccccccccccccccccccccccc
c     Use ZGGEV instead of ZGEGV
cccccccccccccccccccccccccccccccccccccccccccccc
c         call zgegv('N','V',k,ra, ldvs, rb, ldvs, alpha, beta,vsl, ldvs,
c     $              vsr,ldvs, zwork,ldzwork,rwork,info)
         CALL ZGGEV('N','V',k,ra, ldvs, rb, ldvs, alpha, beta,vsl, ldvs,
     $        vsr,ldvs, zwork,ldzwork,rwork,info)



c        --- Compute the eigenvectors belonging to the found eigenvalues
c            of A and put them in EIVEC
         call zgemm('n', 'n', n, k, k, one, work(1,q), n,
     $              vsr, ldvs, zero, eivec, n)
      else
c
c...        Store the Schurvectors in eivec:
         call zcopy( k*n, work(1,q), 1, eivec, 1 )
         call zcopy( k, aconv, 1, alpha, 1 )
         call zcopy( k, bconv, 1, beta, 1 )
      end if

      end

      subroutine error (m)
c
c     Coded by Diederik R. Fokkema
c
c     $Id: error.f,v 1.4 1995/07/26 09:26:26 fokkema Exp $
c
c     .. Parameters ..
c
      implicit none
      character m*(*)

ctex@ \begin{manpage}{ERROR} 
ctex@ \subtitle{Name}
ctex@    ERROR --- Type an error message and stop
ctex@
ctex@ \subtitle{Declaration}
ctex@    %declaration
ctex@ \subtitle{Parameters}
ctex@    \variable{m}
ctex@       character string. On entry m must contain the message string.
ctex@
ctex@ \subtitle{Description}
ctex@    This subroutine types an error message and stops.
ctex@
ctex@ \end{manpage}
ctex@ \begin{verbatim}
ctex@    % actual code
ctex@ \end{verbatim}
c
c     .. Local ..
c
c     None.
c
c     .. Called subroutines
c
c     None.
c
c     .. Executable statements
c
      print 10, m
 10   format (/,1x,'Error: ',a,/)
c
c     --- Stop
c
      stop
      end

      subroutine jdqzmv (n, x, y, work, alpha, beta)
c
c     Coded by Diederik Fokkema
c
c     .. Parameters ..
c
      implicit none
      integer n
      complex*16 alpha, beta, x(*), y(*), work(*) 
c
c     .. Local ..
c
      integer i
c     .. Executable statements ..
c
      call amul( n, x, work )
      call bmul( n, x, y )
      do i = 1, n
         y(i) = beta*work(i) - alpha*y(i)
      end do
c
      end

      subroutine makemm (n, k, w, v, m, zm, ldm)
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/08/03 23:33:20 caveman>
c
c     .. Parameters ..
c
      implicit none
      integer n, k, ldm
      complex*16 w(n,*), v(n,*), m(ldm,*), zm(ldm,*)
c
c     .. Local ..
c
      integer i, j
      complex*16 zdotc
c
c     .. Executable statements ..
c
      do i=1,k
         do j=1,k
            if (i.eq.k.or.j.eq.k)
     $           m(i,j) = zdotc (n, w(1,i), 1, v(1,j), 1)
            zm(i,j) = m(i,j)
         enddo
      enddo
c
c     --- Return
c
      end

      subroutine mkqkz (n, k, q, kq, qkz, invqkz, ldqkz, ipiv)
c
c     Coded by Diederik Fokkema
c
c     $Id: mkqkz.f,v 1.1 1995/08/05 09:08:21 caveman Exp $
c
c     Time-stamp: <95/08/03 23:52:52 caveman>
c
c     .. Parameters ..
c
      implicit none
      integer n, k, ldqkz, ipiv(*)
      complex*16 q(n,*), kq(n,*), qkz(ldqkz,*), invqkz(ldqkz,*)
c
c     .. local ..
c
      integer i, j, info
      complex*16 zdotc
c
c     .. Executable statements ..
c
      do i=1,k
         do j=1,k
            if (i.eq.k.or.j.eq.k) qkz(i,j) =
     $           zdotc (n, q(1,i), 1, kq(1,j), 1)
            invqkz(i,j) = qkz(i,j)
         enddo
      enddo
      call zgetrf (k, k, invqkz, ldqkz, ipiv, info)
c
c     --- Return
c
      end


      subroutine myexc (n, s, t, z, q, ldz, ifst, ilst)
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/10/31 23:51:12 caveman>
c
c     .. Parameters ..
c
      implicit none
      integer n, ldz, ifst, ilst
      complex*16 s(ldz,*), t(ldz,*), z(ldz,*), q(ldz,*)
c
c     .. Local ..
c
      logical tlts
      integer k, m1, m2, m3
      complex*16 f, f1, f2, c1, c2, r, sn
      real*8 cs
      intrinsic dconjg
c
c     .. Executable statements ..
c
      if (n.eq.1 .or. ifst.eq.ilst) return
      if (ifst.lt.ilst) then
         m1 = 0
         m2 = -1
         m3 = 1
      else
         m1 = -1
         m2 = 0
         m3 = -1
      endif
      do k = ifst+m1, ilst+m2, m3
         f = max(abs(t(k+1,k+1)),abs(s(k+1,k+1)))
         f1 = t(k+1,k+1)/f
         f2 = s(k+1,k+1)/f
         tlts = .true.
         if (abs(f1).gt.abs(f2)) tlts = .false.
         c1 = f1*s(k,k) - f2*t(k,k)
         c2 = f1*s(k,k+1) - f2*t(k,k+1)
         call zlartg (c2, c1, cs, sn, r)
         call zrot (k+1, s(1,k+1), 1, s(1,k), 1, cs, sn)
         call zrot (k+1, t(1,k+1), 1, t(1,k), 1, cs, sn)
         call zrot (n, q(1,k+1), 1, q(1,k), 1, cs, sn)
         if (tlts) then
            c1 = s(k,k)
            c2 = s(k+1,k) 
         else
            c1 = t(k,k)
            c2 = t(k+1,k) 
         endif
         call zlartg (c1, c2, cs, sn, r)
         call zrot (n-k+1, s(k,k), ldz, s(k+1,k), ldz, cs, sn)
         call zrot (n-k+1, t(k,k), ldz, t(k+1,k), ldz, cs, sn)
         call zrot (n, z(1,k), 1, z(1,k+1), 1, cs, dconjg(sn))
      enddo
c
c     --- Return
c
      end
      subroutine psolve (n, x, nq, q, kz, invqkz, ldqkz, ipiv, f)
c
c     Coded by Diederik Fokkema
c
c     .. Parameters ..
c
      implicit none
      integer n, nq, ldqkz, ipiv(*)
      complex*16 x(*), q(n,*), kz(n,*),
     $     invqkz(ldqkz,*), f(*)
c
c     .. local .. 
c
      integer info
      complex*16 zero, one
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))
c
c     .. Executable Statements ..
c
      call precon( n, x )
      call zgemv ('c', n, nq, one, q, n, x, 1, zero, f, 1)
      call zgetrs('n', nq, 1, invqkz, ldqkz, ipiv, f, ldqkz, info)
      call zgemv ('n', n, nq, -one, kz, n, f, 1, one, x, 1)
c
c     --- Return
c
      end



      subroutine qzsort (ta, tb, k, s, t, z, q, ldz, alpha, beta,
     $     order)
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/08/03 23:34:03 caveman>
c
c     .. Parameters ..
c
      implicit none
      integer k, ldz, order
      complex*16 ta, tb, s(ldz,*), t(ldz,*), z(ldz,*),
     $     q(ldz,*), alpha(*), beta(*)
c
c     .. Local ..
c
      integer i, j, select
c
c     .. Executable statements ..
c
      do i = 1,k
         do j = 1,k
            alpha(j) = s(j,j)
            beta(j) = t(j,j)
         enddo
         j = select (k-i+1, ta, tb, alpha(i), beta(i), order) + i-1
         call myexc (k, s, t, z, q, ldz, j, i)
      enddo
c
c     --- Return
c
      end
      integer function select (n, sa, sb, a, b, order)
c
c     Coded by Diederik Fokkema
c     Modified Martin van Gijzen, test on division by zero
c
c     .. Parameters ..
c
      implicit none
      integer n, order
      complex*16 sa, sb, a(*), b(*)
c
c     .. Local ..
c
      integer i, j
      real*8 dtmp, optval
c
c     .. Executable statements ..
c
      j = 1
      if ( order .le. 0 ) then
         optval =  1.d99
      else
         optval = -1.d99
      end if
      if (order.eq.0) then
c
c...        Nearest to target
c
         do i=1,n
            if ( b(i) .ne. 0.d0 ) then
               dtmp = abs(a(i)/b(i)-sa/sb)
               if (dtmp.lt.optval) then
                  j=i
                  optval = dtmp
               end if
            end if
         enddo
      elseif (order.eq.-1) then
c
c...        Smallest real part
c
         do i=1,n
            if ( b(i) .eq. 0.d0 ) then
               if ( dble(a(i)) .lt. 0.d0 ) then
                  j=i
                  optval = -1.d99
               end if
            else
               dtmp = dble(a(i)/b(i))
               if (dtmp.lt.optval) then
                  j=i
                  optval = dtmp
               end if
            end if
         enddo
      elseif (order.eq.1) then
c
c...        Largest real part
c
         do i=1,n
            if ( b(i) .eq. 0.d0 ) then
               if ( dble(a(i)) .gt. 0.d0 ) then
                  j=i
                  optval = 1.d99
               end if
            else
               dtmp = dble(a(i)/b(i))
               if (dtmp.gt.optval) then
                  j=i
                  optval = dtmp
               end if
            end if
         enddo
      elseif (order.eq.-2) then
c
c...        Smallest imaginari part
c
         do i=1,n
            if ( b(i) .eq. 0.d0 ) then
               if ( dimag(a(i)) .lt. 0.d0 ) then
                  j=i
                  optval = -1.d99
               end if
            else
               dtmp = dimag(a(i)/b(i))
               if (dtmp.lt.optval) then
                  j=i
                  optval = dtmp
               end if
            end if
         enddo
      elseif (order.eq.2) then
c
c...        Largest imaginari part
c
         do i=1,n
            if ( b(i) .eq. 0.d0 ) then
               if ( dimag(a(i)) .gt. 0.d0 ) then
                  j=i
                  optval = 1.d99
               end if
            else
               dtmp = dimag(a(i)/b(i))
               if (dtmp.gt.optval) then
                  j=i
                  optval = dtmp
               end if
            end if
         enddo
      else
         call error ('unknown order in select')
      endif
      select = j
c
c     --- Return
c
      end

      subroutine zcgstabl (n, x, r, l, eps, mxmv,
     $     zalpha, zbeta, nk, kz, qq, invqkz, ldqkz, jpiv, f,
     $     work, lwork)
c
c     Programmer: Diederik R. Fokkema
c
c
      implicit none
c
c     .. Parameters ..
c
      integer l, n, mxmv, nk, lwork, ldqkz, jpiv(*)
      real*8  eps
      complex*16 x(*), r(*), kz(n,*), qq(n,*), work(n,*),
     $     zalpha, zbeta, invqkz(ldqkz,*), f(*)

c
c     .. Local ..
c
      integer mxl
      parameter (mxl = 32)

      integer i, j, k, m, ipiv(mxl), nmv, info

      integer u, q, w, rr, xp, bp
      logical rcomp, xadd

      real*8 maxnorm, delta, bpnorm, tr0norm, trlnorm
      complex*16 varrho, hatgamma

      real*8 rnrm, rnrm0, eps1, dznrm2
      complex*16 alpha, beta, omega, gamma, rho0, rho1,
     $     yr(mxl), ys(mxl), z(mxl,mxl), ztmp

      complex*16 zero, one
      parameter (zero = (0.0d0,0.0d0), one = (1.0d0,0.0d0))

      complex*16 zdotc
      intrinsic dconjg
c
c     .. Executable statements ..
c
      if (mxmv.eq.0) return

      if (l.gt.mxl)
     $     call error ('l exceeds mxl (zcgstabl)')

      u  = 1
      q  = u + (l+1)
      rr = q + (l+1)
      xp = rr + 1
      bp = xp + 1
      w = bp + 1

      if (w.gt.lwork)
     $     call error ('workspace too small (zcgstabl)')

c
c     --- set x to zero and compute first residue
c
      call zzeros (n, x)
      call zscal (n, -one, r, 1)
      call zcopy (n, r, 1, work(1,rr), 1)
      call psolve (n, r, nk, qq, kz, invqkz, ldqkz,
     $     jpiv, f)

c
c     --- initialize loop
c
      nmv = 0

      rnrm0 = dznrm2 (n, r, 1)
      rnrm = rnrm0
      eps1 = eps * rnrm0

      call zcopy (n, x, 1, work(1,xp), 1)
      call zcopy (n, r, 1, work(1,bp), 1)
      maxnorm = 0d0
      bpnorm = rnrm
      rcomp = .false.
      xadd = .false.
      delta = 1d-2

      m = 0
      rho0 = one
      alpha = one
      beta = zero
      omega = one

      call zzeros (n*(l+1), work(1,u))
      call zzeros (n*(l+1), work(1,q))
      call zcopy (n, r, 1, work(1,q), 1)
c
c     --- loop
c
 1000 continue
      m = m + l
c
c     --- BiCG part
c
      rho0 = -omega * rho0
      do k=1,l
         rho1 = zdotc (n, work(1,rr), 1, work(1,q+k-1), 1)
         beta = alpha * (rho1 / rho0)
         rho0 = rho1
         beta = beta
         do j=0,k-1
            call zxpay (n, work(1,q+j), 1, (-beta), work(1,u+j), 1)
         enddo
         call jdqzmv (n, work(1,u+k-1), work(1,u+k), work(1,w),
     $        zalpha, zbeta)
         call psolve (n, work(1,u+k), nk, qq, kz,
     $        invqkz, ldqkz, jpiv, f)
         gamma = zdotc (n, work(1,rr), 1, work(1,u+k), 1)
         alpha = rho0 / gamma
         do j=0,k-1
            call zaxpy (n, (-alpha), work(1,u+j+1), 1, work(1,q+j), 1)
         enddo
         call jdqzmv (n, work(1,q+k-1), work(1,q+k), work(1,w),
     $        zalpha, zbeta)
         call psolve (n, work(1,q+k), nk, qq, kz,
     $        invqkz, ldqkz, jpiv, f)
         call zaxpy (n, alpha, work(1,u), 1, x, 1)
         rnrm = dznrm2 (n, work(1,q), 1)
         maxnorm = max (maxnorm, rnrm)
         nmv = nmv+2
      enddo
c
c     --- MR part + Maintaining the convergence
c
      do i=1,l-1
         do j=1,i
            ztmp = zdotc (n, work(1,q+i), 1, work(1,q+j), 1)
            z(i,j) = ztmp
            z(j,i) = dconjg(ztmp)
         enddo
         yr(i) = zdotc (n, work(1,q+i), 1, work(1,q), 1)
         ys(i) = zdotc (n, work(1,q+i), 1, work(1,q+l), 1)
      enddo
      call zgetrf (l-1, l-1, z, mxl, ipiv, info)
      call zgetrs ('n', l-1, 1, z, mxl, ipiv, yr, mxl, info)
      call zgetrs ('n', l-1, 1, z, mxl, ipiv, ys, mxl, info)
      call zcopy (n, work(1,q), 1, r, 1)
      call zgemv ('n', n, l-1, (-one), work(1,q+1), n, yr, 1, one,
     $     r, 1)
      call zgemv ('n', n, l-1, (-one), work(1,q+1), n, ys, 1, one,
     $     work(1,q+l), 1)
      tr0norm = dznrm2 (n, r, 1)
      trlnorm = dznrm2 (n, work(1,q+l), 1)
      varrho = zdotc (n, work(1,q+l), 1, r, 1) / (tr0norm*trlnorm)
      hatgamma = varrho/abs(varrho) * max(abs(varrho),7d-1)
      hatgamma = (tr0norm/trlnorm)*hatgamma
      yr(l) = zero
      ys(l) = -one
      call zaxpy (l, (-hatgamma), ys, 1, yr, 1)

      omega = yr(l)
      call zgemv ('n', n, l, one, work(1,q), n, yr, 1, one, x, 1)
      call zgemv ('n', n, l, (-one), work(1,u+1), n, yr, 1, one,
     $     work(1,u), 1)
      call zaxpy (n, (-hatgamma), work(1,q+l), 1, r, 1)
      call zcopy (n, r, 1, work(1,q), 1)
c
c     --- reliable update
c
      rnrm = dznrm2 (n, work(1,q), 1)
      maxnorm = max (maxnorm, rnrm)
      xadd = (rnrm.lt.delta*rnrm0.and.rnrm0.lt.maxnorm)
      rcomp = ((rnrm.lt.delta*maxnorm.and.rnrm0.lt.maxnorm).or.xadd)

      if (rcomp) then
         call jdqzmv (n, x, work(1,q), work(1,w),
     $        zalpha, zbeta)
         call psolve (n, work(1,q), nk, qq, kz,
     $        invqkz, ldqkz, jpiv, f)
         call zxpay (n, work(1,bp), 1, -one, work(1,q), 1)
         maxnorm = rnrm
         if (xadd) then
            call zaxpy (n, one, x, 1, work(1,xp), 1)
            call zzeros (n, x)
            call zcopy (n, work(1,q), 1, work(1,bp), 1)
            bpnorm = rnrm
         endif
      endif
         
      if (nmv.lt.mxmv .and. rnrm.gt.eps1) goto 1000

      call zaxpy (n, one, work(1,xp), 1, x, 1)
c
c     --- return
c
      mxmv = nmv
      eps = rnrm/rnrm0

      return
      end

      subroutine zgmres (n, x, r, mxm, eps, mxmv, 
     $     alpha, beta, k, kz, q, invqkz, ldqkz, ipiv, f, v, tp)
c
c     Programmer: Diederik R. Fokkema
c
c
      implicit none
c
c     .. Parameters ..
c
      integer mxm, mxmv, n, k, ldqkz, ipiv(*)
      real*8 eps
      complex*16 x(*), r(*), kz(n,*), q(n,*), v(n,*),
     $     alpha, beta, invqkz(ldqkz,*), f(*), tp(*)

ctex@ \begin{manpage}{ZGMRES} 
ctex@
ctex@ \subtitle{ZGMRES} 
ctex@    ZGMRES -- Generalized Minimal Residual
ctex@    iterative method for solving linear systems $\beta A-\alpha B = -r$.
ctex@    This subroutine in specilized for use in JDQZ.
ctex@ 
ctex@ \subtitle{Declaration}
ctex@ \function{subroutine zgmres (n, x, r, mxm, eps, mxmv, a, ka, b, kb,
ctex@   alpha, beta, k, kz, mqkz, zmqkz, ldvs, q,
ctex@   lu, klu, dlu, v)}
ctex@
ctex@ \subtitle{Parameters}
ctex@    \variable{integer n} 
ctex@       On entry, n specifies the dimension of the matrix A.
ctex@       
ctex@    \variable{x} 
ctex@       complex*16 array of size n. 
ctex@       On exit, x is overwritten by the approximate solution.
ctex@
ctex@    \variable{r}
ctex@       complex*16 array of size n. Before entry, the array r 
ctex@       must contain the righthandside of the linear problem Ax=r. 
ctex@       Changed on exit.
ctex@ 
ctex@    \variable{integer mxm} 
ctex@       On entry, mxm specifies the degree of the Minimal Residual
ctex@       polynomial.
ctex@
ctex@    \variable{{real*8} eps}
ctex@       On entry, eps specifies the stop tolerance. On exit, eps contains
ctex@       the relative norm of the last residual.
ctex@
ctex@    \variable{integer mxmv}
ctex@       On Entry, mxmv specifies the maximum number of matrix 
ctex@       multiplications. On exit, mxmv contains the number of matrix
ctex@       multiplications performed by the method.
ctex@
ctex@    \variable{{complex*16} zalpha}
ctex@       On entry, zalpha specifies $\alpha$. Unchanged on exit.
ctex@ 
ctex@    \variable{{complex*16} zbeta}
ctex@       On entry, zbeta specifies $\beta$. Unchanged on exit.
ctex@ 
ctex@    \variable{integer k}
ctex@       On entry, k specifies the number of columns of the arrays
ctex@       kz and q.
ctex@
ctex@    \variable{z}
ctex@       complex*16 array z, of size (n,k). On entry the array z
ctex@       must contain the preconditioned matrix Z.
ctex@
ctex@    \variable{mqkz}
ctex@       complex*16 array mqkz, of size (ldvs,k). On entry the array 
ctex@       mqkz must contain the matrix Q'*KZ.
ctex@
ctex@    \variable{zmqkz}
ctex@       complex*16 array zmqkz, of size (ldvs,k). Workspace. Used to
ctex@       copy mqkz.
ctex@
ctex@    \variable{q}
ctex@       complex*16 array q, of size (n,k). On entry the array q
ctex@       must contain the preconditioned matrix Q.
ctex@
ctex@    \variable{v}
ctex@       complex*16 array of size (n,mxm+1). Workspace.
ctex@
ctex@ \subtitle{Description}
ctex@    ***
ctex@
ctex@ \subtitle{See Also}
ctex@    ***
ctex@
ctex@ \subtitle{References}
ctex@    ***
ctex@ 
ctex@ \subtitle{Bugs}
ctex@    ***
ctex@
ctex@ \subtitle{Author}
ctex@     Diederik R.\ Fokkema
ctex@
ctex@ \end{manpage}
ctex@ \begin{verbatim}
ctex@    % actual code
ctex@ \end{verbatim}
ctex@
c
c     .. Local ..
c
      logical restrt, loop
      integer maxm
      parameter (maxm = 100)

      integer i, m, m1, nmv
      real*8 rnrm0, rnrm, eps1, c(maxm)
      complex*16 hh(maxm,maxm-1), rs(maxm), s(maxm), y(maxm), rcs
      complex*16 zero, one
      parameter (zero = (0.0d0,0.0d0), one = (1.0d0,0.0d0))

      real*8  dznrm2
      complex*16 ztmp, zdotc
c
c     .. Executable Statements ..
c
      if ( mxm .gt. maxm-1 )
     $     call error ('mxm larger than maxm (zgmres)')
c
c     --- Initialize first residue
c
      call zzeros (n, x)
      call zscal (n, -one, r, 1)
      call psolve (n, r, k, q, kz, invqkz, ldqkz,
     $     ipiv, f)
c
c     --- initialize loop
c
      rnrm0  = dznrm2 (n,r,1)
      rnrm = rnrm0
      eps1  = eps * rnrm

      nmv = 0

      call zcopy (n, r, 1, v(1,1), 1)
c         
c     --- outer loop
c
 1000 restrt = ( nmv .lt. mxmv .and. rnrm .gt. eps1 )
      if ( restrt ) then
         ztmp = one / rnrm
         call zscal (n, ztmp, v(1,1), 1)
         rs(1) = rnrm
c
c     --- inner loop
c
         m = 0
 2000    loop = (nmv.lt.mxmv .and. m.lt.mxm .and. rnrm.gt.eps1)
         if (loop) then
            m  = m + 1
            m1 = m + 1
            call jdqzmv (n, v(1,m), v(1,m1), tp, alpha, beta)
            call psolve (n, v(1,m1), k, q, kz, invqkz,
     $           ldqkz, ipiv, f)
            nmv = nmv+1 
            do i = 1,m
               ztmp = zdotc (n, v(1,i), 1, v(1,m1), 1)
               hh(i,m) = ztmp
               call zaxpy (n, (-ztmp), v(1,i), 1, v(1,m1), 1)
            enddo
            ztmp = dznrm2( n, v(1,m1), 1 )
            hh(m1,m) = ztmp
            call zscal ( n, (one/ztmp), v(1,m1), 1 )
            do i = 1,m-1
               call zrot (1, hh(i,m), 1, hh(i+1,m), 1, c(i), s(i))
            enddo
            call zlartg (hh(m,m), hh(m1,m), c(m), s(m), rcs )
            hh(m,m) = rcs
            hh(m1,m) = zero
            rs(m1) = zero
            call zrot (1, rs(m), 1, rs(m1), 1, c(m), s(m))
            rnrm = abs (rs(m1))
            goto 2000
         endif
c
c     --- compute approximate solution x
c
         call zcopy ( m, rs, 1, y, 1 )
         call ztrsv ( 'u', 'n', 'n', m, hh, maxm, y, 1 )
         call zgemv ( 'n', n, m, one, v, n, y, 1, one, x, 1 )
c
c     --- compute residual for restart
c
         call jdqzmv (n, x, v(1,2), tp, alpha, beta)
         call psolve (n, v(1,2), k, q, kz, invqkz,
     $        ldqkz, ipiv, f)
         call zcopy (n, r, 1, v(1,1), 1)
         call zaxpy (n, -one, v(1,2), 1, v(1,1), 1)
         rnrm = dznrm2 (n, v(1,1), 1)

         goto 1000
      endif
c
c     --- return
c
      eps = rnrm/rnrm0
      mxmv = nmv

      return
      end

      subroutine zmgs (n, k, v, w, job )
c
c     Coded by Diederik Fokkema
c     Modified 05-23-96: M. Kooper: job =1 corrected, array YWORK added
c
c     .. Parameters ..
c
      implicit none
      integer n, k, job
      complex*16 v(n,*), w(*)
c
c     .. Local ..
c
      integer i
      real*8 s0, s1, dznrm2
      complex*16 znrm
c
c     .. Executable statements ..
c
      s1 = dznrm2 (n, w, 1)
      do i=1, k
	 s0 = s1
	 call zortho (n, v(1,i), w, s0, s1, znrm)
      enddo
      if (job.eq.0) then
	 return
      else
         znrm  = 1.d0/s1
         call zscal (n, znrm, w, 1)
      endif
c
c     --- Return
c
      return
      end

      subroutine zortho (n, v, w, s0, s1, znrm)
c
c     Coded by Diederik Fokkema
c
c     .. Parameters ..
c
      implicit none
      integer n
      real*8 s0, s1
      complex*16 v(*), w(*), znrm
c
c     .. Local ..
c
      real*8 kappa, dznrm2
      complex*16 ztmp, zdotc
      parameter (kappa=1d2)
c
c     .. Executable statements ..
c
      znrm = zdotc (n, v, 1, w, 1)
      call zaxpy (n, (-znrm), v, 1, w, 1)
      s1 = dznrm2 (n, w, 1)
      if (s1.gt.s0/kappa) then
	 goto 100
      else
	 s0 = s1
	 ztmp = zdotc (n, v, 1, w, 1)
         znrm = znrm + ztmp
	 call zaxpy (n, (-ztmp), v, 1, w, 1)
	 s1 = dznrm2 (n, w, 1)
	 if (s1.gt.s0/kappa) then
	    goto 100
	 else
	    call error ('zero vector in zmgs')
	 endif
      endif
 100  continue
c
c     --- Return
c
      return
      end

      subroutine zones (n, x)
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/07/30 21:45:46 caveman>
c
c     .. Parameters ..
c
      implicit none
      integer n
      complex*16 x(*)
c
c     .. Local ..
c
      integer i
c
c     .. Executable statements ..
c
      do i=1,n
         x(i) = (1.0d0,0.0d0)
      enddo
c
c     --- Return
c
      end

      subroutine zxpay(n,dx,incx,da,dy,incy)
c
c     modified by:  D.R. Fokkema
c     01/06/94
c
c     a vector plus constant times a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      implicit none
      complex*16 dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = da*dy(iy) + dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = da*dy(i) + dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = da*dy(i) + dx(i)
        dy(i + 1) = da*dy(i + 1) + dx(i + 1)
        dy(i + 2) = da*dy(i + 2) + dx(i + 2)
        dy(i + 3) = da*dy(i + 3) + dx(i + 3)
   50 continue
      return
      end
      subroutine zzeros (n, x)
c
c     Coded by Diederik Fokkema
c
c     $Id$
c
c     Time-stamp: <95/07/30 21:48:00 caveman>
c
c     .. Parameters ..
c
      implicit none
      integer n
      complex*16 x(*)
c
c     .. Local ..
c
      integer i
c
c     .. Executable statements ..
c
      do i=1,n
         x(i) = (0.0d0,0.0d0)
      enddo
c
c     --- Return
c
      end

c     Fictitious function used by ZGEGS
c     In our case (sort='n') this function is not referenced

      logical function delctg (x1,x2)

      complex*16 x1,x2

      DELCTG =.true.

      return
      end
