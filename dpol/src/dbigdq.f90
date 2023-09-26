program dall2020
!! see README.md for old version notes

   use dqmodule
   use wavext

!*================== TYPE INITIALIZATION ========================

   implicit type (dq_real) (a-h, o-z)

   integer*4 jat0, jat, jcs, jatd, jcsd
   integer*8 igs
   integer pam, qam

   real*8 dlog
   real*16 egs16, amm16, dinc16, cp116, cp216, psi216, psi16, cpd1x16, &
   & cpd2x16, cpd116, cpd216, etot16, bscale16

   logical lext, ldump, lstop, lname, ltrid, lreset, lwcof, lpow, lpowf, lgo
   dimension cpd1x16(8), cpd2x16(8), cp116(12150), cp216(12150), cpd116(8), cpd216(8)

   dimension psi2(12150), psi(12150), md(12150), dso(12150), rmass(12), knv(2), &
   & dhm(9, 2), tot(1, 1), ee(12150), v1(1, 1)
   character fnwve*16, cl(8)*1, title*56, an*1, size*5, cz*2, status*8, &
   & line*56, pow*4, date*8, eo*1, fnwv(2)*16, cdump*24, dmp*4, fnin*12, buffer*80
   
   data cl/'S','P','D','F','G','H','I','K'/

!*================== COMMON BLOCKS ========================

   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   & ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), & 
   & fm2(84), f21x(84, 86)

   common/b1b/cpd1x(8), cpd2x(8), nit, next, kst, key1, ibasis, kbasis, linc, &
   & nv1, longnv, ksm(8), mar12, mar1, kono, nbetx, status

   common/d1d/cp1(12150), cp2(12150), cpd1(8), cpd2(8), b(8, 2), der(9, 2), &
   & ndp(12150), ndq(12150), nds(12150), ld1(12150), ld2(12150), pam, qam, maxc, &
   & ia, nblk(9), nb, nblkx(9), nbx, maxbf, nrow(8998), lngth

   common/h1h/ch(25), neig, nh, iz, lrgl, nspn, nv, l1max, nbet, nd, n1, n2, &
   & n1d, n2d, ltrans, iset, namm

   common/r1r/rap, egs, ltrid

   common hb(12150*(12150+1)/2), at(198000), at2(84, 84), atl(84, 84, 4), &
   & atj(84, 84, 4), jat0(84), jat(80, 80), jcs(80, 80), &
   & zsum(84, 86), bin(84, 84), x(84, 9), at0(84, 84), jatd(46, 46), &
   & jcsd(46, 46), atd(38600), ovb(10000)

   common/done/lr12, idone, lext

   common/res/lrglx, fnwve, title

!*================== FORMATS ========================

 1 format(4i3, 1x, a13, a48)
 4 format(1x, 1p, 3d38.31)
33 format(3i3, i5, i3, 1x, a13, a56)
34 format('-----------------------E-N-D----------------------------')
52 format(d12.5, f12.7, i3)
66 format(12h eigenvalues)
71 format(16i5)
72 format(3i5, i2, 3f20.15)
73 format(10(i3, 2i2))

   rmass(1) = 0.0005443205671q0
   rmass(2) = 0.0001370745620q0
   rmass(3) = 0.0000782019500q0
   rmass(4) = 0.00006088119800q0
   rmass(5) = 0.0000498387000q0
   rmass(6) = 0.0000457254400q0
   rmass(7) = 0.0000391848100q0
   rmass(8) = 0.0000343054100q0
   rmass(9) = 0.0000288817300q0
   rmass(10) = 0.0000274462000q0
   rmass(11) = 0.0000238675800q0
   rmass(12) = 0.0000228775500q0

   next = 0
   v1 = dqreal(0.q0)

!*================== EXTRA PARAMETERS ========================
! TODO: consider using command-line arguments to change these params from their defaults.

   !? TRUE to write the large hamiltonian and overlap matrices to disk
   ldump = .true.

   !? TRUE to use Jacobi's ME. FALSE for tridiagonalization
   ltrid = .true.

   !? TRUE to use the power method to calculate derivatives and optimize non-linear parameters
   lpow = .true.

   !? TRUE to exit after a single calculation.
   lstop = .false.

   !? TRUE to reset the input file after each calculation
   lreset = .true.

   !? TRUE to write a dump file contatining the wave function
   lwcof = .false.

   call routedopen(2, file='ermacdq.dat', status='NEW')

   ermac = dqreal(1.q0)
   test = 1.q0
   
 7 continue
 
   ermac = 0.5q0*ermac
   test = test + ermac

   rewind(2)

   call dqwrite(2, test)

   rewind(2)

   read(2,'(a)') buffer

   if (buffer(4:4).ne.'2') go to 7

   close(2, status='DELETE')

   ermac = 2.q0*ermac
   write(*,'(a, D12.4)') 'ERMAC = ', qreal(ermac)

   cdump = 'DUMP'
   ndump = 4
   fnin = 'dall2020.dat'
   
   call routedopen(1, file='notdone', status='UNKNOWN')

   kount = 1
   e = 0.q0

   close(1, status='DELETE')

   call getdate(date)

   call routedopen(4, file='dall2020.out', status='UNKNOWN')
   write(4,*) date

   call routedopen(3, file=fnin, status='OLD')
   
   bscale = 0.q0
   der = dqreal(0.q0)

   read(3,*) nbet, iset, ibasis, amm16, dinc16, nh, bscale16, linc
   nmax = 0

   if (mr12mx.eq.0) mr12mx = 100
   if (mr1mx.eq.0) mr1mx = 100

   klein = .false.
   lgo = .true.
   if (nmax.gt.0) lgo = .false.

   amm=amm16
   dinc=dinc16
   bscale = bscale16

   kmass = 0
   if (amm.lt.0.q0) kmass = -1

   long = 12150*(12150+1)/2
   longnv = long
   nd = 12150
   maxbf = 10000
   n1 = 80
   n2 = 80
   n1d = 36
   n2d = 36
   namm = 0
   
   if (amm.ne.0.q0) namm = 1
   rap = ermac
   ifirst = 0
10 read(3,*, err=99)egs16, nit, neig

   egs=dqreal(egs16*(1.q0 + 1.q-25))
   
   if (.not.lpow) nit = 0
   esav = egs
   if (egs.eq.0.q0) go to 99
   ein = egs

   kontrol = 0
   if (nit.lt.0) then
      !? kontrol = 1: save "sh" and "ovb" matrices
      !? kontrol = 2: generate and save complete eigenvector set
      kontrol = 1
      if (lpow) kontrol = 0
      nit = 0
   endif

   if (egs.eq.dqreal(1.q0)) egs = e

   if (nh.ge.10) then
      nh = nh - 10
   else
      nh = 1
   endif

   !? ltrans = 0:
   !?    kdiv = 1: r1 divided basis set
   !?    kdiv = 2: r12 divided basis set
   !? ltrans = 1 if block 2 basis functions are to be transformed 

   ltrans = 0
   kdiv = 0
   if (kdiv.ge.1) ltrans = 0


   read(3,*) iz, lrgl, nspn, l1max, fnwve, title
   write(*, 1)iz, lrgl, nspn, l1max, fnwve, title
   
   if (kmass.lt.0.and.iz.gt.12) stop 1
   if (kmass.lt.0) amm= rmass(iz)

   cdump = fnwve(1:7)

   if (iz.eq.0) go to 99
   
   line(36:43) = date(1:8)
   status = title(45:53)
   lgo = .true.

   if (nmax.eq.0.and.status(2:3).eq.'OK') lgo = .false.

   if (status(2:4).eq.'NDG') lgo = .true.

   if (status(2:3).eq.'GO') then
      lgo = .true.
      write(*,'(A, I3)') 'Extrapolated alpha''s and beta''s with KOUNT = ', kount
      write(4,'(A, I3)') 'Extrapolated alpha''s and beta''s with KOUNT = ', kount
   endif

   if (iset.lt.2) iset = 1
   if (lrgl.eq.0.and.nspn.eq.0) iset = 0
   if (lrgl.eq.0) ibasis = 0

   !? kbasis introduced to correct basis sets for 1S5D STATE
   !? ibasis = 1 & kbasis = 0 for block 2 with adjusted L's and the same as for block 1. (?)
   !? ibasis = 0 & kbasis = 1 for block 2 with adjusted L's and independent (?)

   ibasis = 0
   kbasis = 1

   if (lrgl.eq.0.and.nspn.eq.0.and.iset.lt.2) iset = 0
   if (lrgl.eq.0) ibasis = 0
   if (lrgl.ge.2) then
      ibasis = 0
      kbasis = 1
   endif
   
   if (iz.le.1) nh = 0
   if (nh.eq.0) ltrans = 0
   

   if (ldump) then
      call routedopen(9, file=cdump(1:ndump)//'2', form='unformatted', status='UNKNOWN')
      call routedopen(10, file=cdump(1:ndump)//'3', form='unformatted', status='UNKNOWN')
   endif
   kst = - nit

25 next = next + 1

   if (next.gt.nit.and.next.gt.1) go to 10

   call basis(md, lgo, ifirst, kdiv, knv, bscale)

   if (.not.lpow.and.nv.gt.(nd+1)/2) then
      write(*,*) 'NVmax = ',(nd+1)/2,' FOR COMPLETE DIAGONALIZATION.'
      stop
   endif

   if (ifirst.eq.1) egs = e
   if (next.eq.1) e = egs
   ifirst = 1
   
   egs16=egs
   !EGS = DQREAL(NINT(2**30*EGS16)/2**30)   !!!!!WE COMMENT IT OUT JUST FOR NOW
   z = dqreal(iz)
   cp11 = 1.q0/(1 + linc)
   cp21 = (z-1.q0)/dqreal(lrgl+neig)/z

   if (nh.ne.0) cp21 = cp2(1)
   if (nh.ne.0) cp11 = cp1(1)

   esh = z*z*(cp11**2 + cp21**2)/2.q0
   e = esh
   lngth = long

   do i=1, nbet + nh
      do k=1, 2
         dhm(i, k) = der(i, k)
      end do
   end do
      
   
   kdiv1 = kdiv + 1
   if (kdiv1.gt.2) kdiv1 = 2

   do kkdiv=1, kdiv1
      nv = knv(kkdiv)
      size = '00000'
      write(size(1:5),'(I5)') nv
      
      write(cz,'(I2)') iz
      n = lrglx+neig+linc
      an = char(n+48)
      if (n.ge.10) an = char(n+55)

      dmp = '.DMP'
      pow = '.POW'
      if (namm.ne.0) pow = '.POL'
      if (namm.ne.0) dmp = '.DML'

      if (kkdiv.eq.2) pow(2:2) = 'V'
      
      is = -dlog10(1.d0*nv) + 6
      is4 = 5
      
      fnwv(kkdiv) = cz//an//char(2*nspn+49)//cl(lrgl+1)//size(is:is4)//pow

      write(*,*) fnwv(kkdiv), kkdiv

      if (iz.lt.10) fnwv(kkdiv) = fnwv(kkdiv)(2:14)

      !? calculate total parity
      ipar = (2*lrgl+1 - 4*((lrgl+1)/2))*(2*linc+1 - 4*((linc+1)/2))

      eo = 'e'
      if (ipar.eq.-1) eo = 'o'

      line = ' '
      line(2:11) = 'Z=' // cz // ' ' // an // ' ' // char(2*nspn+49) // cl(lrgl+1) // eo

      if (.not.lgo) then
         e = ein
         kst = 0
         key1 = 1
         if (lreset) call reset(dhm, ein, kontrol, bscale)
         next = 0
         go to 10
      endif

      if (kkdiv.eq.1) fnwv(2) = fnwv(1)
   end do

   fnwve = fnwv(2)

   if (lwcof) then
      write(*,*) 'creating ', fnwve(1:7)//dmp, nh, iz

      inquire(file=fnwve(1:7)//dmp, exist=lname)
      inquire(file=fnwve(1:7)//pow, exist=lpowf)
      
      if (lname.and.lpowf) then
         write(*,*) fnwve(1:7)//dmp,'  already exists.'
         write(*,*) 'Delete if a new version is needed.'

         der(1, 1) = 0.q0
         kontrol = 0

         write(4, 33) iz, lrgl+linc, nspn, nv, neig, fnwve, ' Exit called, .DMP file already exists.'

         if (lreset) call reset(der, esav, kontrol, bscale)
         stop
      else
         call routedopen(8, file=fnwve(1:7)//dmp, form='UNFORMATTED', status='UNKNOWN')
      endif
   endif

   call routedopen(1, file='notdone', status='UNKNOWN')
   write(1,'(4A2, 1X, A3, 1X, I5, 1X, A14, 1X, A, A)') cz, an, char(2*nspn+49), cl(lrgl+1), pow(2:4), nv, fnwve,'   '
   close(1, status='KEEP')

   call routedopen(1, file='current.dat', status='UNKNOWN')
   write(1,'(4A2, 1X, A3, 1X, I5, 1X, A14, 1X, A/A)') cz, an, char(2*nspn+49), cl(lrgl+1), pow(2:4), nv, fnwve,' running', fnin
   close(1)
   
   call ham(psi, psi2, egs, amm, next, md, .false., ldump, kdiv, bscale)
   
   nvv1 = nv1

   call final(hb, ovb, v1, tot, psi, psi2, dso, e, ee, next, nvv1, md, kontrol, ldump, lwcof, lpow)
   
   if (nit.eq.0) status = ', NOD, dq'
   
   line(36:43) = date(1:8)
   line(44:55) = status
   
   title = line
   if (kst.ne.0)call ham(psi, psi2, e, amm, next, md,.true., ldump, kdiv, bscale)
   
   if (nit.eq.0) status = ', NOD, dq'
   if (nit.eq.0) status = ', NOD, dq'
   line(36:43) = date(1:8)
   line(44:55) = status
   title = line

   !? check for NaN's

   if (e.ne.e) then
      nanc = nanc + 1
      e = ein
      status = 'NaN, dq'
      if (nanc.gt.2) then
         write(*,*) 'TOO MANY NaN''s'
         next = nit
         oksig = 'NaN'
      endif
   endif

   brat = cpd2(1)/beta1
   
   if (nit.eq.0.and.lreset) call reset(der, e, kontrol, bscale)
   
   if (.not.lpow) kontrol = 1
   
   etot = e*z*z/2.q0 - esh

   if (nh.eq.0) etot = e*z*z/2.q0 - esh

   lrgl = lrgl - linc

   write(4, 4) qreal(etot)
   write(*,*) qreal(etot)
   call dqwrite(6, etot)

   if (kst.ne.0.and.next.lt.nit) go to 25
   
   if (ldump) close(10, status='DELETE')

   if (ksm(2).eq.0) then
      do i=2, nb-nh+1
         ksm(i) = ksm(i+1)
      end do
   endif

   write(line(12:33),'(11I2)') (ksm(i), i=1, nb-nh), mar12, mar1, kono

   line(36:43) = date(1:8)
   line(44:55) = status

   write(*, 33) iz, lrgl+linc, nspn, nv, neig, fnwve, line
   write(4, 33) iz, lrgl+linc, nspn, nv, neig, fnwve, line

   call routedopen(7, file=fnwve, status='unknown')

   write(7, 33) iz, lrgl+linc, nspn, nv, neig, fnwve, line
   write(7, 4) qreal(e), qreal(etot), qreal(amm)
   write(7, 71) nb,(nblk(i+1), i=1, nb),(ksm(i), i=1, nb-nh), mar12, mar1, kono, 1-nh, kdiv

   do i=1, nbet
      ii = nbet - i + 1
      cpd1(ii+ibasis) = cpd1(ii)
      cpd2(ii+ibasis) = cpd2(ii)
   end do

   do nbl1=1, nb
      nb1 = nblk(nbl1) + 1
      nb2 = nblk(nbl1+1)

      if (bscale.gt.0.q0.and.nbl1.eq.1+nh) then
         write(7, 72) nb1, nb2, ld1(nb1), ld2(nb1), qreal(cp1(nb1)), qreal(cp2(nb1)), bscale
      else
         write(7, 72) nb1, nb2, ld1(nb1), ld2(nb1), qreal(cp1(nb1)), qreal(cp2(nb1))
      endif
   end do

   write(7, 73) (ndp(i), ndq(i), nds(i), i=1, nv)
   write(7, 4) (qreal(psi(k)), k=1, nv)

   if (lwcof) then
      write(*,*) 'writing eigenvalues'
      write(7, 66)

      write(7, 4) (qreal(ee(i)), i=1, nv)
      write(8) (qreal(ee(i)), i=1, nv)
   endif

   next = 0
   close(7, status='KEEP')

   if (kontrol.ne.0) then
      call routedopen(7, file='spinv.dat', status='UNKNOWN')
      rewind 7

      write(7,*) 'SECOND-ORDER ENERGY CALCULATION'
      title(1:16) = ''''//fnwve(1:2)//fnwve(4:4)//'spn2.mat'''

      if (iz.gt.10) title(1:16) = ''''//fnwve(1:3)//fnwve(5:5)//'spn2.mat'''

      write(7,'(A16, I4)') title(1:16), kontrol

      if (fnwv(1).eq.fnwv(2)) then
         write(7,'(A)') '''A:'//fnwve//''',''SAME'''
      else
         write(7,'(A)') '''A:'//fnwv(1)//''',''A:'//fnwv(2)//''''
      endif
      
      write(7,'(A)') '''EXIT''/'
      close(7, status='KEEP')

      call dumph(hb, ovb, dso, psi, e, nv1, md, ldump, cdump, ndump)

      call routedopen(1, file='notdone', status='UNKNOWN')
      close(1, status='KEEP')
   else
      if (ldump) close(8, status='DELETE')
      if (ldump) close(9, status='DELETE')

      call routedopen(1, file='notdone', status='UNKNOWN')
      close(1, status='DELETE')
   endif

   nblk(2) = 0

   write(6, 71) (nblk(i+1)-nblk(i), i=2, nb), mar12, mar1, kono
   write(4, 71) (nblk(i+1)-nblk(i), i=2, nb), mar12, mar1, kono
   write(6, 34)
   write(4, 34)

   
   do i=1, nbet
      cpd1(i) = cpd1(i+ibasis)
      cpd2(i) = cpd2(i+ibasis) 
   end do

   if (lstop) go to 99
   if (kontrol.eq.0) go to 10

99 continue
   close(3, status='KEEP')
   stop
end program


subroutine hdiagl(h, n, iegen, u, md, nr, n7d)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)
   type (dq_real) h, u
   type (dq_real) yi, yj
   
   common/r1r/rap, egs, ltrid

   logical ltrid, hk
   dimension h(1), u(n, n), x(1800), yi(1800), yj(1800), iq(1800), md(1)
   real*16 ck, dk, ek, gk, ys, rk, sk, tk, uk

   sq2 = sqrt(dqreal(0.5q0))

   if (iegen .eq. 0) then
      do i=1, n
         do j=1, n
            if (i-j .eq. 0) then
               u(i, j)=1.q0
            else
               u(i, j)=0.q0
            end if
         end do
      end do
   end if

   nr = 0
   if (n-1 .le. 0) go to 10

   !? Scan for largest off diagonal element in each row
   !? x(i) contains largest element in i'th row
   !? iq(i) holds second subscript defining position of element

   nmi1 = n-1

   do i=1, nmi1
      x(i) = 0.q0
      ipl1=i+1
      do j=ipl1, n
         ! converting abs from dq to real
         ck = x(i)-abs( h(i+md(j)))
         if (ck .gt. 0) cycle
         x(i) = abs(h(i+md(j)))
         iq(i)=j
      end do
   end do

   hdtest=1.0q38

   !? find maximum in x for pivot element and test for end of problem

20 continue
   
   ipiv = 1

   do i=2, nmi1
      if (x(i).gt.x(ipiv)) ipiv = i
   end do

   xmax = x(ipiv)
   jpiv = iq(ipiv)
   mdj = md(jpiv)
   mdi = md(ipiv)
   
   dk = xmax - hdtest
   if (dk .le. 0) then
      hdimin = abs(h(1))

      do i=2, n
         ek=hdimin-abs(h(i+md(i)))

         if (ek .le. 0) cycle
         
         hdimin = abs(h(i+md(i)))
         imin = i
      end do

      hdtest=hdimin*rap
      gk = hdtest - xmax

      if (gk .ge. 0) go to 10
   end if

   nr = nr+1

   if (mod(nr, 10000).eq.0) write(*,'(A1,$)') '.'

   !? compute tan, sin, and cos, h(i, i), h(j, j)

   y = h(ipiv+mdi) - h(jpiv+mdj) + rap

   xdeb=h(ipiv+mdj)/y
   xd2 = xdeb*xdeb

   tang = 2.q0*xdeb/(1.q0 + sqrt(dqreal(1.q0+4.q0*xd2)))
   cosine = 1.q0/sqrt(dqreal(1.q0 + tang*tang))
   sine = tang*cosine

   delta = -sine*sine*y + 2.q0*sine*cosine*h(ipiv+mdj)
   h(ipiv+mdi)=h(ipiv+mdi)+delta
   h(jpiv+mdj)=h(jpiv+mdj)-delta
   h(ipiv+mdj)=0.q0

   !? pseudo rank the eigen values, adjust sin and cos for computation of h(ik) and u(ik)

   ys = y

   if (ys .gt. 0) then
      htemp = h(ipiv+mdi)
      h(ipiv+mdi) = h(jpiv+mdj)
      h(jpiv+mdj) = htemp

      !? recompute sin and cos

      htemp = abs(cosine)

      if (sine.gt.0.q0) htemp = -htemp

      cosine = abs(sine)
      sine = htemp
   end if

   !? change the other elements of h

   if (ipiv.ne.1) then

      im1 = ipiv - 1
      imdi = mdi
      imdj = mdj

      do i=1, im1
         imdi = imdi + 1
         imdj = imdj + 1

         htemp = h(imdi)
         h(imdi) = cosine*htemp + sine*h(imdj)
         h(imdj) = -sine*htemp + cosine*h(imdj)

         if ((iq(i)-jpiv .ne. 0) .and. (iq(i).ne.ipiv)) cycle

         xi = abs(h(imdi))
         xj = abs(h(imdj))

         if (xi.gt.xj) then 
            x(i) = xi
            iq(i) = ipiv

            cycle
         end if

         x(i) = xj
         iq(i) = jpiv
      end do
   
   end if

   jn = ipiv + 1

   if (ipiv.ne.jpiv-1) then

      ip1 = jn
      jm1 = jpiv - 1

      imdi = ipiv + md(ip1)
      imdj = ipiv + mdj

      do i=ip1, jm1
         imdj = imdj + 1

         htemp = h(imdi)
         h(imdi) = cosine*htemp + sine*h(imdj)
         h(imdj) = -sine*htemp + cosine*h(imdj)

         if (abs(h(imdi)).gt.abs(h(ipiv+md(jn)))) jn = i

         if (iq(i).eq.jpiv) then

            ipl1 = i + 1
            jm = ipl1

            hmax = 0.q0

            ij = i + md(ipl1)

            do j=ipl1, n
               if (abs(h(ij)).gt.hmax) then
                  jm = j
                  hmax = abs(h(ij))
               endif

               ij = ij + j
            end do
            
            x(i) = hmax
            iq(i) = jm

         end if
         
      end do

      imdi = imdi + i   
   
   end if

   if (jpiv.ne.n) then

      jp1 = jpiv + 1
      jq = jp1
      jmdi = jpiv + md(jp1)

      jmdi = jpiv + md(jp1)
      imdi = ipiv + md(jp1)

      do i=jp1, n
         htemp = h(imdi)
         h(imdi) = cosine*htemp + sine*h(jmdi)
         h(jmdi) = -sine*htemp + cosine*h(jmdi)
         if (abs(h(imdi)).gt.abs(h(ipiv+md(jn)))) jn = i
         if (abs(h(jmdi)).gt.abs(h(jpiv+md(jq)))) jq = i

         jmdi = jmdi + i
         imdi = imdi + i
      end do

      x(jpiv) = abs(h(jpiv+md(jq)))
      iq(jpiv) = jq

   end if

   x(ipiv) = abs(h(ipiv+md(jn)))
   iq(ipiv) = jn
   
   !? test for computation of eigenvectors

   if (iegen .ne. 0 ) go to 20

   do i=1, n
      htemp=u(i, ipiv)
      u(i, ipiv)=cosine*htemp+sine*u(i, jpiv)
      u(i, jpiv)=-sine*htemp+cosine*u(i, jpiv)
   end do

   go to 20
   
   !? do final check for large o.d. matrix elements

10 continue

   do i=1, nmi1
      x(i) = 0.0q0
      ipl1=i+1

      do j=ipl1, n
         uk = x(i)-abs( h(i+md(j)))
         if (uk .le. 0) then
            x(i) = abs(h(i+md(j)))
            iq(i)=j
         end if
      end do

      if (x(i).gt.xmax) xmax = x(i)
   end do

   if (xmax.gt.hdtest) go to 20

   write(*,*)

end subroutine

!* calculate factorials, binomial coefficients, and other needed constants
subroutine calc(y1, y2, isum, mp1, lgo)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)

   logical lr12, lgo
   integer*4 jat0, jat, jcs, jatd, jcsd
   type (dq_real) facx(88), sumlj(6), xx(2)

   common/done/lr12, idone, lext

   common hb(12150*(12150+1)/2), at(198000), at2(84, 84), atl(84, 84, 4), &
   &   atj(84, 84, 4), jat0(84), jat(80, 80), jcs(80, 80), &
   &   zsum(84, 86), bin(84, 84), x(84, 9), at0(84, 84), jatd(46, 46),  &
   &   jcsd(46, 46), atd(38600), ovb(12150*(12150+1)/2)

   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), &
   &  fm2(84), f21x(84, 86)

   common/x1x/zsub(84, 86), at3(84, 84)

01 format(' EXIT CALLED, ITOT =', i6)
   
   pi = acos(-1.q0)

   alf = y1
   bet = y2

   xam = (y1-y2)/y1
   xbm = (y2-y1)/y2
   xap = (y1+y2)/y1
   xbp = (y2+y1)/y2

   do i=1, 86
      do j=1, 2
         do k=1, 6
            suml(i, j, k) = 1.q180*(i+j+k)
         end do
      end do
   end do

   do i=1, 84
      do j=1, 86
         f21x(i, j) = 0.q0
         zsub(i, j) = 1.q70
         zsum(i, j) = 1.q80
      end do
   end do

   r(1) = -1.q0
   r(2) = 0.q0
   
   do i=3, 1250
      r(i) = r(i-1) + dqreal(1.q0)/dqreal(i-1)
   end do

   fac(1) = 1.q0
   fac(2) = 1.q0
   xf = 1.q0

   do i=3, 86
      xf = xf*dqreal(i-1)
      fac(i)=xf
   end do
   
   !? note that fac(n) = factorial(n-1)

   if (idone.ne.1.) then

      bin(1, 1) = 1.q0

      do i=2, 84
         im1 = i - 1
         bin(i, i) = 1.q0

         do j=1, im1
            xf = 1.q0

            do k=j, im1
               xf = xf*dqreal(k)
            end do

            bin(j, i) = xf
         end do
      end do

      do i=2, 84
         im1 = (i+1)/2

         do j=1, im1
            k = i - j + 1
            rf = bin(k, i)/fac(j)
            bin(i, j) = rf
            bin(i, k) = bin(i, j)
         end do
      end do

      do i=1, 84
         x(i, 1) = 1.q0
         xff = 1.q0
         do j=2, 9
            xff = xff/dqreal(i+2*j-4)
            x(i, j) = xff
         end do
      end do
      
      idone=1
   
   end if

   !? note that bin(i, j) = binomial coefficient(i-1, j-1) for i > j and
   !?           bin(i, j) = i * (i + 1) * ... * (j-1)      for i < j     

   itot = isum + 4
   itot2 = itot + 2
   itot1 = itot + 1

   if (isum.eq.0) return
   if (itot.gt.82) then
      write(*, 01) itot
      stop
   end if

   y12 = y1 + y2

   call star2(y1, itot2, y1x)
   call star2(y2, itot2, y2x)
   call star2(y12, itot2, yp)

   sum = 0.q0

   do j1=1, itot
      fj1 = dqreal(j1)
      sum = sum + yp(j1)/y1x(j1)
      z1x(j1) = sum*fac(j1)*y1x(j1+1)/y12
   end do
   
   do jp=1, itot1
      jtot = itot1 - jp + 1
      if (jp.eq.1) jtot = itot
      sum = 0.q0

      do j1=1, jtot
         sum = sum + bin(j1, j1+jp-1)*yp(j1)/y1x(j1)
         zsum(j1, jp) = sum*fac(j1)*y1x(j1+1)*yp(jp+1)
      end do

      sum = 0.q0
      sumb = 0.q0

      do j1=1, jtot
         dz = bin(j1, j1+jp-1)*yp(j1)/y2x(j1)
         sum = sum + dz

         if (j1.gt.2) sumb = sumb + dz

         zsub(j1, jp) = sumb*fac(j1)*y2x(j1+1)*yp(jp+1)
         zsum(85-j1, 87-jp) = sum*fac(j1)*y2x(j1+1)*yp(jp+1)
      end do
   end do
   
   lr12 = .false.
   if (.not.lr12) return

   do j1=1, itot
      ia = j1 - 2
      maxbc = isum - ia + 3
      x1 = 1.q0
      sum = 0.q0
      yr = y2/y12

      do j=maxbc, 20000
         sum = sum + x1
         sf = dqreal(ia+j+3)/dqreal(j+1)
         x1 = x1*yr*sf
         if (x1/sum.lt.ermac) go to 22
      end do

      !! ask about this
      read(*,*)

      22 continue

      sum = sum*bin(maxbc+1, ia+maxbc+3)*yp(maxbc+1)/y2x(maxbc+1)
      at0(maxbc-1, j1) = sum
      maxj = maxbc - 2

      do j=1, maxj
         jj = maxbc - j
         sum = sum + bin(jj+1, ia+jj+3)*yp(jj+1)/y2x(jj+1)
         at0(jj-1, j1) = sum
      end do

      maxj = maxj + 1

      do j=1, maxj
         at0(j, j1) = -at0(j, j1)*2.q0*fac(j+1)*y2x(j+2)*yp(ia+4)
      end do
   end do

   !? at0(ib+ic+1, ia+2) is the infinite sum contribution to the remainder
   !? i(r)(ia, ib, ic) after subtraction of the dominant term.
   
   !? log terms for tab.

   xx(1) = y1/(y1+y2)
   xx(2) = y2/(y1+y2)
   is7 = isum + 7
   is6 = isum + 6
   kmax = mp1-1
   
   if (kmax.gt.6) then
      write(*,*) 'KMAX =', kmax
      stop 40
   endif
   
   do i=1, 2

      do k=1, kmax
         sumlj(k) = 0.q0
      end do

      j = is6
      xj = xx(i)**j

      31 continue

      j = j + 1
      xj = xj*xx(i)
      xjj = xj

      do k=1, kmax
         xjj = xjj/dqreal((j+k-1))
         sumlj(k) = sumlj(k) + xjj
      end do

      if (xj/(dqreal(j)*sumlj(1)).gt.ermac) go to 31

      do k=1, kmax
         suml(is7, i, k) = sumlj(k)
      end do
      
      do j=1, is6
         jj = is6 - j + 1
         xj = xx(i)**jj
         do k=1, kmax
            xj = xj/(dqreal(jj+k-1))
            suml(jj, i, k) = suml(jj+1, i, k) + xj
         end do
      end do
   end do

   do i=1, 84
      fp2(i) = 1.q60
      fm2(i) = 1.q60
   end do

   if (lgo) then
      rba = 2.q0/xbp
      do 45 i=1, itot1
      xj = fac(i+2)*rba*rba/(fac(i)*2.q0)
      fp2(i) = xj
      j = 2
      sum = 0.q0
      44 j = j + 1
      xj = xj*rba*dqreal(i+j-1)/dqreal(j)
      sum = sum + xj
      if (abs(xj/sum).gt.ermac) go to 44
      fm2(i) = -sum
      45 fp2(i) = (fp2(i) + sum)
   end if

   rba = min(alf/bet, bet/alf)
   rba2 = rba*rba

   if (abs(xam).ge.1.q-12) then
      rba3 = rba2*rba
      ts1 = rba3/9.q0
      ts2 = rba2/4.q0
      rx = rba2

      do j=4, 1000000, 2
         rx = rx*rba2
         fj = dqreal(j)
         ts2 = ts2 + rx/(fj*fj)
         fj = j + 1.q0
         dts1 = rx*rba/(fj*fj)
         ts1 = ts1 + dts1
         if (dts1/ts1.lt.ermac) exit
      end do
      
      go to 200
   end if

   ts1 = pi*pi/8.q0 - 1.q0
   ts2 = pi*pi/24.q0

   200 continue

   dl1 = 0.q0
   dl2 = 0.q0

   if (xam.ne.0.q0) then
      dl1 = log(dqreal(abs(xap/xam)))
      dl2 = log(dqreal(1.q0-rba2))
   end if

   do k=1, 2
      a = alf
      b = bet

      if (k.ne.1) then
         a = bet
         b = alf
      end if

      abl = log(dqreal(a/b))

      signx = 1.q0
      if (a.lt.b) signx = -1.q0

      ts(1, 1, k) = abl*(0.5q0*(a*a+b*b)*dl1 + a*b*(dl2-1.q0)) + b**3/a + signx*((a*a+b*b)*ts1 - 2.q0*a*b*ts2)
      ts(1, 2, k) = -abl*(b*dl2 + a*dl1 - 2.q0*b)+ 2.q0*signx*(b*ts2 - a*ts1)
      ts(1, 3, k) = abl*dl1 - 2.q0*(b*abl/a - signx*ts1)
      ts(2, 1, k) = -abl*(b*dl1 + a*dl2) - 2.q0*b*b/a - 2.q0*signx*(b*ts1 - a*ts2)
      ts(2, 2, k) = abl*dl2 - 2.q0*signx*ts2

      if (a.lt.b) then
      
         ts(1, 1, k) = ts(1, 1, k) - b**3/a - a**3/b + 0.25q0*pi*pi*(b-a)**2 - a*b*abl*abl - a*b*(2.q0-pi*pi/3.q0)
         ts(1, 2, k) = ts(1, 2, k) + b*abl*abl + (b-a)*pi*pi/2.q0 + b*(2.q0-pi*pi/3.q0) + 2.q0*a*a/b
         ts(1, 3, k) = ts(1, 3, k) + 2.q0*(-b/a - a/b + pi*pi/4.q0)
         ts(2, 1, k) = ts(2, 1, k) + 2.q0*b*b/a - 0.5q0*pi*pi*(b-a) + a*abl*abl + a*(2.q0 - pi*pi/3.q0)
         ts(2, 2, k) = ts(2, 2, k) - abl*abl - pi*pi/6.q0
      
      end if
      
      ts(2, 3, k) = 0.q0
      ts(3, 1, k) = ts(1, 3, k) + 2.q0*b*(abl + 1.q0)/a

   end do
   
end subroutine


subroutine star2(a, n, an)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)

   dimension an(84)
   
   do i=1, n
      an(i) = dqreal(1.q0)/a**(i-1)
   end do

end subroutine

!* calculate the i(0) value needed in the recursion relation
!* if ic = 0, the integral of 1/r12 - 1/r2 is calculated
function t(ia, ib, ic)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)

   integer*4 jat0, jat, jcs, jatd, jcsd
   integer q1, p1, q1p, p1p

   common hb(12150*(12150+1)/2), at(198000), at2(84, 84), atl(84, 84, 4),           &
   &   atj(84, 84, 4), jat0(84), jat(80, 80), jcs(80, 80)                    &
   &   , zsum(84, 86), bin(84, 84), x(84, 9), at0(84, 84), jatd(46, 46),        &
   &   jcsd(46, 46), atd(38600), ovb(12150*(12150+1)/2)

   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84)   &
   &  , ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84)   &
   &  , fm2(84), f21x(84, 86)

   if (ic.le.-2) then
      t = 1.q100
      return
   endif

   if (ic.eq.-1) then
      t = tlog(ia-1, ib-1, ic-1, sumj)
      t = t  + sumj
      return
   endif

   if (ia.lt.0.or.ib.lt.0) then
      t = tab(ia, ib, ic)
      return
   endif

   xt = 0.q0
   icg1 = ishft(ic,-1) + 1
   ic2 = ic + 2
   ibc3 = ib + ic + 3
   iac3 = ia + ic + 3

   do i=1, icg1

      i2 = 2*i
      q1 = ibc3 - i2
      p1 = ia + i2
      q1p = iac3 - i2
      p1p = ib + i2
      xt=xt+bin(ic2, 2*i)*(zsum(q1p, p1p) + zsum(85-q1, 87-p1))

   end do

   t=xt*dqreal(2.q0)/dqreal((ic+1))

end function


function tab(ia, ib, ic)
   use dqmodule

   implicit type (dq_real) (a-h, o-z)

   integer*4 jat0, jat, jcs, jatd, jcsd
   integer p1, q1, p1p, q1p

   common hb(12150*(12150+1)/2), at(198000), at2(84, 84), atl(84, 84, 4), &
   &   atj(84, 84, 4), jat0(84), jat(80, 80), jcs(80, 80), &
   &   zsum(84, 86), bin(84, 84), x(84, 9), at0(84, 84), jatd(46, 46), &
   &   jcsd(46, 46), atd(38600), ovb(12150*(12150+1)/2)

   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), &
   &  fm2(84), f21x(84, 86)
   
   tabx = 0.0q0
   ic1 = ic + 1
   ic2 = ic + 2
   ia1 = ia + 1
   ib1 = ib + 1
   iac3 = ia + ic + 3
   ibc3 = ib + ic + 3
   icg1 = ic/2 + 1

   if (2*icg1-ic1 .gt. 0) then
      iabc3 = ia + ibc3

      if (iabc3.le.0) then
         tab = 1.q120*(ia-ib+10)
         return
      end if

      if (ib1.le.0) then
         do k=1,-ib, 2
            tabx = tabx + bin(ic2, ic2+1-k)*fac(iac3-k)*suml(iabc3, 1, 1-k-ib)*y1x(iabc3+1)
         end do
      end if
      
      if (ia1.le.0) then
         do k=1,-ia, 2
            tabx = tabx + bin(ic2, ic2+1-k)*fac(ibc3-k)*suml(iabc3, 2, 1-k-ia)*y2x(iabc3+1)
         end do
      end if
   end if
   
   tabx1 = 0.q0

   do i=1, icg1
      i2 = 2*i
      q1 = ibc3 - i2
      p1 = ia + i2
      q1p = iac3 - i2
      p1p = ib + i2
      sum1 = 0.q0
      sum2 = 0.q0

      if (q1.gt.0.and.p1.gt.0) sum1 = zsum(85-q1, 87-p1)
      if (q1p.gt.0.and.p1p.gt.0) sum2 = zsum(q1p, p1p)

      tabx1 = tabx1 + bin(ic2, i2)*(sum1+sum2)
   end do

   tabx = tabx + tabx1
   tab = (tabx/dqreal(ic1))*2.q0
   return

end function


function tlog(iax, ibx, ic, sumj)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)
   integer*4 jat0, jat, jcs, jatd, jcsd
   real*16 zzam

   common hb(12150*(12150+1)/2), at(198000), at2(84, 84), atl(84, 84, 4), &
   &   atj(84, 84, 4), jat0(84), jat(80, 80), jcs(80, 80), &
   &   zsum(84, 86), bin(84, 84), xxx(84, 9), at0(84, 84), jatd(46, 46), &
   &   jcsd(46, 46), atd(38600), ovb(12150*(12150+1)/2)
   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84),  &
   &  fm2(84), f21x(84, 86)

   gam = 0.q0
   klog = -3
   sumj = 0.q0

   if (ibx.gt.-2 .or. iax.le.-2) then
      if (ibx+ic.le.-3.and.iax.le.-3) then
         tlog = 0.q0
         return
      endif

      ia = iax
      ib = ibx
      a = alf
      b = bet
      kts = 1
   end if 

   if (ibx.le.-2 .and. iax.gt.-2) then
      ia = ibx
      ib = iax
      a = bet
      b = alf
      kts = 2
   end if

   ab = a/b
   ba = b/a
   apb = a + b
   zap = apb/a
   zbp = apb/b
   zam = (a-b)/a

   abl = log(dqreal(a/b))

   if (abs(abl).lt.1.q-10) abl = 0.q0

   ic2 = ic + 2
   ic3 = ic + 3
   ic4 = ic + 4
   icg = ic/2 + 1

   ia1 = ia + 1
   ia2 = ia + 2
   ia3 = ia + 3
   ia4 = ia + 4
   ia5 = ia + 5

   ib1 = ib + 1
   ib2 = ib + 2

   ibc3 = ib + ic + 3
   ibc5 = ib + ic + 5
   
   rba = min(ba, ab)
   rba2 = rba*rba

   !? compute the j=0 contribution

   sumj = 0.q0

   if (ia.gt.klog.and.ic2.ne.0) sumj = 2.q0*ic2*fac(ibc3)*fac(ia3)*(r(ibc3) - gam - log(dqreal(b/bet))/(a**ia3*b**(ibc3)))
   
   sum1 = 0.q0
   sum1g = 0.q072

   fc = fac(ic2+1)
   x0 = 2.q0*fc/(a**ia3*b**(ibc3))

   if (icg.ge.2) then
      do j1=2, icg
         j2 = 2*(j1-1)

         if (ibc3-j2.le.0) cycle

         x = x0*(fac(ibc3-j2)/fac(ic2-j2))*(fac(j2+ia3)/fac(j2+2))*ba**j2
         sum1g = sum1g + x*(r(ic+2) - r(ic2-j2) + r(ibc3-j2) - gam - log(dqreal(b/bet)))
         sum1 = sum1 + x*(r(ic+2) + r(j2+ia3) - r(j2+2) - gam -log(dqreal(a/bet)))
      end do
   end if

   sum2 = 0.q0

   do m1=1, ib2
      sa = 0.q0
      sb = 0.q0
      mm3 = m1 - 4
      iam = ia + m1
      iam1 = iam + 1
      x0 = 1.q0/(apb**iam1)

      if (m1.gt.3.or.ia4.ne.0) then
         if (mm3.ge.0) then
            if (mm3.gt.0) then
               do l=1, mm3
                  sb = sb + fac(iam1-l)*zbp**l/(fac(mm3-l+1)*dqreal(l))
               end do

               sb = sb*fac(mm3+1)*x0
            end if

            if (ia4.ge.1) then
               do l=1, ia4
                  sa = sa + fac(iam1-l)*zap**l/(fac(ia5-l)*dqreal(l))
               end do

               sa = sa*fac(ia5)*x0
            end if 

            sa = (sa-sb) - fac(iam1)*abl*x0 - fac(ia5)*fac(mm3+1)*f21(ia4, iam1, a, b)/(dqreal(iam1)*a**ia5*b**mm3)
            cycle
         end if

         if (iam.ge.0) then
            if (iam.gt.0) then
               do l=1, iam
                  sa = sa + zap**l/dqreal(l)
               end do
            end if

            sa = (sa-abl)*x0 - (-1.q0)**m1*f21(iam, iam1, a, b)/(dqreal(iam1)*a**iam1)

            if (m1.eq.1) sa = sa - 2.q0*ia2*b*(abl-r(ia3))/(a**ia3)
            if (m1.eq.2) sa = sa + 2.q0*(abl-r(ia3)-1.q0)/(a**ia3)

            sa = sa*fac(iam1)
            cycle
         end if
      end if
      
      sa = ts(m1, ia4+1, kts)
      sum2 = sum2 + sa*bin(ib2, m1)*fac(ibc5-m1)/(b**(ibc5-m1))
   end do

   if (ic2.le.0.and.ia3.gt.0) then
      !? j=0 contribution for ic = -2
      x = 2.q0*fac(ia3)/(a**ia3*b**ib1)
      if (ib1.eq.0) sumj =  x*(abl - r(ia3))
      if (ib1.gt.0) sumj =  x*fac(ib1)
   end if 

   tlog = sum1 + sum2
   
end function


function f21(m1, n1, a, b)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)
   real*16 d

   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), &
   &  fm2(84), f21x(84, 86)

25 format(2i5,' OUT OF RANGE IN F21')

   if (a.eq.alf) then
      ms = m1 + 1
      ns = n1 - m1
   else
      ms = 82 - m1
      ns = 85 - n1 + m1
   end if

   if (ms.lt.0 .or. ns.lt.1) then
      write(*, 25) m1, n1
      stop
   end if

   if (f21x(ms, ns).ne.0.q0) then
      f21 = f21x(ms, ns)
      return 
   end if

   n = n1
   d = a - b

   if (d.lt.0) then
      m = n1 - m1 - 1
      z = - d/b
      c = a/b
   else if (d.gt.0) then
      m = m1
      z = d/a
      c = 1.q0
   else
      f21 = 1.q0
      f21x(ms, ns) = f21
      return
   end if

   sum = 1.q0
   fx = 1.q0

   do j=1, 2000
      fx = fx*z*dqreal(j+m)/dqreal(j+n)
      sum = sum + fx
      if (abs(fx/sum).lt.ermac) exit
   end do

   f21 = sum*c
   f21x(ms, ns) = f21

   return

end function

function tc3(ia, ib, ic)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)

   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), &
   &  fm2(84), f21x(84, 86)

   integer p, p1, p2

01 format(' TC3 ERROR', 3i4)

   p = ia + ib + 2
   p1 = p + 1
   p2 = p + 2

   ia1 = ia + 1
   ia2 = ia1 + 1
   ib1 = ib + 1
   ib2 = ib1 + 1

   sumap = 0.q0
   sumbp = 0.q0
   sumam = 0.q0
   sumbm = 0.q0

   !if (ib1)100, 20, 21

   if (ia1.lt.0 .or. ib1.lt.0) then
      write(*, 01) ia, ib, ic
      stop
   end if

   if (ib1.gt.0) then
      21 do j=1, ib1
         sumbp = sumbp + fac(p1-j)*y2x(j+1)*yp(p2-j)/(dqreal(j)*fac(ib2-j))
      end do
   end if

   if (ia1.gt.0) then
      do j=1, ia1
         sumap = sumap + fac(p1-j)*y1x(j+1)*yp(p2-j)/(dqreal(j)*fac(ia2-j))
      end do
   end if
   
   if (alf.ne.bet) then

      if (ib1.gt.0) then
         do j=1, ib1
            sumbm = sumbm + fac(p1-j)*xbm**j/(fac(ib2-j)*dqreal(j))
         end do
         
         sumbm = sumbm*(-1.q0)**ia1/(bet-alf)**p1*fac(ib2)
      end if
      do j=1, ia1
         sumam = sumam + fac(p1-j)*xam**j/(fac(ia2-j)*dqreal(j))
      end do

      sumam = sumam*(-1.q0)**ib1/(alf-bet)**p1*fac(ia2)
      slog = (-1.q0)**ib*fac(p1)*log(alf/bet)/(alf-bet)**p1
      summ = sumam + sumbm + slog
   end if

   if (abs(summ/max(abs(sumam), abs(sumbm), abs(slog))).le.0.5q0) then
      summ = f21(ia1, p1, alf, bet)/alf
      summ = -summ*fac(ia2)*y1x(ia2)*fac(ib2)*y2x(ib2)/dqreal(p1)
   end if
   
   tc3 = sumbp*fac(ib2) + sumap*fac(ia2) + summ

end function


function tc3x(ia, ib, ic, sumj)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)

   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), &
   &  fm2(84), f21x(84, 86)

01 format(' TC3X ERROR', 3i4)

   ia1 = ia + 1
   ia2 = ia1 + 1
   ib1 = ib + 1
   ib2 = ib1 + 1

   suma = 0.q0
   sumb = 0.q0
   sumjx= 0.q0

   rba = 1.q0/xbp

   if (ib.ne.0 .or. ia1.lt.0) then
      write(*, 1) ia, ib, ic
      stop
   end if

   is = 1
   do j=1, ib1
      is = -is
      n = j + ia1

      c = fac(n)*y2x(ib2-j+1)*yp(n+1)/(fac(j)*dqreal(ib2-j))

      if (is.le.0) then
         x = fm2(n)
         dx = -2.q0*n*(1.q0 + (n+1.q0)*rba)*rba
         if (j.gt.2) then
            sumb = sumb + (x+dx)*c
            cycle
         end if
         sumjx = sumjx + dx*c
         sumb = sumb + x*c
         cycle
      else 
         x = fp2(n)
         dx = 2.q0*(1.q0 + dqreal(n)*rba)

         if ((j-2).le.0) then
            sumjx = sumjx + dx*c
            sumb = sumb + x*c
         else
            sumb = sumb + (x+dx)*c
         end if
      end if
   end do

   if (ia1.gt.0) then
      is = (-1.q0)**ib1

      do j=1, ia1
         n = j + ib1
         if (is.le.0) then
            x = fm2(n) - 2.q0*n*(1.q0 + (n+1.q0)*rba)*rba
            cycle
         end if
         x = fp2(n) + 2.q0*(1.q0 + dqreal(n)*rba)
         suma = suma + x*fac(n)*y1x(ia2-j+1)*yp(n+1)/(fac(j)*dqreal(ia2-j))
      end do
   end if

   slog = (-1.q0)**ib*fac(ia2+ib1)*log(alf/bet)/((alf-bet)**(ia2+ib1))
   tc3x = sumb*fac(ib2) + suma*fac(ia2) + slog
   sumj = sumjx*fac(ib2)

   if (abs((tc3x+sumj)/slog).lt.0.5q0) tc3x = 0.q0
   
end function

!* version of tc3 for ia = -2 or ib = -2.
function tc32(ia, ib, ic)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)
   integer p, p1, p2

   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), &
   &  fm2(84), f21x(84, 86)

01 format(' TC32 ERROR', 3i4)

   ib0 = max0(ia, ib)
   ia0 = min0(ia, ib)

   ilog = 2
   alf0 = alf
   bet0 = bet

   if (ia0.ne.ia) then
      alf0 = bet
      bet0 = alf
      ib0 = ia
      ia0 = ib
      ilog = 1
   endif

   if (ia0.ne.-2 .or. ic.ne.-3 .or. ib0.lt.0) then
      write(*, 01) ia, ib, ic
      stop
   end if

   ib1 = ib0 + 1
   ib2 = ib1 + 1
   sumbp = 0.q0
   sumbm = 0.q0

   tp = (alf0+bet0)/bet0

   if (ib0.gt.0) then
      do j=1, ib0
         sumbp = sumbp + tp**j/dqreal(j)
      end do

      sumbp = sumbp/(alf0+bet0)**ib1
   end if

   sumbp2 = -suml(ib1, ilog, 1)/bet0**ib1
   sumbp = sumbp + 2.q0*sumbp2

   if (alf0.ne.bet0) then
      tm = (bet0-alf0)/bet0

      if (ib0.gt.0)then
         do j=1, ib0
            sumbm = sumbm + tm**j/dqreal(j)
         end do
         sumbm = -sumbm/(bet0-alf0)**ib1*fac(ib1)
      end if

      slog = (-1.q0)**ib0*fac(ib1)*log(alf0/bet0)/(alf0-bet0)**ib1
      summ = sumbm + slog
   end if

   if (abs(summ/max(abs(sumbm), abs(slog))).le.0.5q0) then
      summ = f21(0, ib1, alf0, bet0)/alf0
      summ = summ*fac(ib1)/(dqreal(ib1)*bet0**ib0)
   end if

   tc32 = sumbp*fac(ib1) + summ

end function

!* version of tc3 for ia = -3 or ib = -3
function tc33(ia, ib, ic)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)

   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), &
   &  fm2(84), f21x(84, 86)

   integer p, p1, p2

01 format(' TC33 ERROR', 3i4)

   ib0 = max0(ia, ib)
   ia0 = min0(ia, ib)

   ilog = 2
   isgn = 2
   alf0 = alf
   bet0 = bet

   if (ia0.ne.ia) then
      alf0 = bet
      bet0 = alf
      isgn = -2
      ib0 = ia
      ia0 = ib
      ilog = 1
   endif

   if (ia0.ne.-3.or.ic.ne.-3) then
      write(*, 01) ia, ib, ic
      stop
   end if

   ib1 = ib0 + 1
   ib2 = ib1 + 1

   tc33 = tc3(ia+isgn, ib-isgn, ic) - 2.q0*fac(ib1)*suml(ib0, ilog, 2)/bet0**ib0

end function

!* calculates basis function integrals with cancelling i=1 terms remove
function txx(ia, ib, ic, kk)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)

   integer*4 jat0, jat, jcs, jatd, jcsd
   integer q1, p1, q1p, p1p

   common hb(12150*(12150+1)/2), at(198000), at2(84, 84), atl(84, 84, 4), &
   &   atj(84, 84, 4), jat0(84), jat(80, 80), jcs(80, 80), &
   &   zsum(84, 86), bin(84, 84), x(84, 9), at0(84, 84), jatd(46, 46), &
   &   jcsd(46, 46), atd(38600), ovb(12150*(12150+1)/2)

   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), &
   &  fm2(84), f21x(84, 86)

   xt = 0.q0
   icg1 = ishft(ic,-1) + 1

   ia1 = ia + 1
   ib1 = ib + 1
   ic1 = ic + 1
   ic2 = ic + 2
   ibc3 = ib + ic + 3
   iac3 = ia + ic + 3

   if ((2*icg1-ic1).gt.0) then
      iabc3 = ia + ibc3

      if (iabc3.le.0) then
         txx = 1.12q0
         return
      end if

      if (ib1.le.0) then
         do k=1,-ib, 2
            tabx = tabx + bin(ic2, ic2+1-k)*fac(iac3-k)*suml(iabc3, 1, 1-k-ib)*y1x(iabc3+1)
         end do
      end if

      if (ia1.le.0) then
         do k=1,-ia, 2
            tabx = tabx + bin(ic2, ic2+1-k)*fac(ibc3-k)*suml(iabc3, 2, 1-k-ia)*y2x(iabc3+1)
         end do
      end if
   end if

   do i=1, icg1
      i2 = 2*i
      bn = bin(ic2, i2)
      q1 = ibc3 - i2
      p1 = ia + i2
      q1p = iac3 - i2
      p1p = ib + i2

      if (i.le.1) then
         if (kk.eq.1.and.q1.gt.0.and.p1.gt.0) xt = xt + bn*zsum(85-q1, 87-p1)
         if (kk.eq.2.and.q1p.gt.0.and.p1p.gt.0) xt = xt + bn*zsum(q1p, p1p)
         cycle
      else 
         if (q1p.gt.0.and.p1p.gt.0) xt = xt + bn*zsum(q1p, p1p)
         if (q1.gt.0.and.p1.gt.0) xt = xt + bn*zsum(85-q1, 87-p1)
      end if
   end do

   txx=xt*2.q0/dqreal(ic+1)

end function


subroutine genint(id, nj, ni, kdiv)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)
   integer*4 jat0, jat, jcs, jatd, jcsd
   logical lext, lgo

   common/done/lr12, idone, lext

   common hb(12150*(12150+1)/2), at(198000), at2(84, 84), atl(84, 84, 4), &
   &   atj(84, 84, 4), jat0(84), jat(80, 80), jcs(80, 80), &
   &   zsum(84, 86), bin(84, 84), xx(84, 9), at0(84, 84), jatd(46, 46), &
   &   jcsd(46, 46), atd(38600), ovb(12150*(12150+1)/2)

   common/h1h/ch(25), neig, nh, iz, lrgl, nspn, nv, l1max, nbet, nd, n1, n2, &
   &   n1d, n2d, ltrans, iset, namm

   common/f1f/y1g, y2g, mp1, mp2, mq1, mq2, ms1, ms2, msl, isum, lg
   
   common/p1p/cplt(9, 11), csum(11), lmn(11), lmb(11), ict, na, nc, nad, ncd

   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), &
   &  fm2(84), f21x(84, 86)

   common/x1x/zsub(84, 86), at3(84, 84)

02 format(' MIN. ATD DIMENSIONS ARE', 3i4, i6)

   if (nj.eq.ni) write(*,'(''|'',$)')
   if (nj.ne.ni) write(*,'(''.'',$)')

   n1 = na - 1
   n2 = nc - 1

   icnt = 0

   do  i=1, 82
      icnt = icnt + 1
      at(icnt) = 1.q90
      jat0(i) = 1

      do  j=1, 82
         jat(i, j) = 0
      end do
   end do
   
   iatd = 0

   do i=1, 36
      iatd = iatd + 1
      atd(iatd) = 1.q80
      atd(iatd+19800) = 1.q80

      do j=1, 36
         jatd(i, j) = 1
         jcsd(i, j) = 1
      end do
   end do

   ism = isum
   imin = 2*lrgl - 1
   iminc = imin + mp1 + mq1 + ms1

   mpb = mq2
   mpa = mp2

   jtop = 0
   itop = 0

   if (nj.le.nh) then
      ism = isum + neig - 1
      mpb = mq2 + neig - 1
      jtop = mq2
   end if
  
   if (ni.le.nh) then
      ism = ism + neig - 1
      mpb = mpb + (1-id)*(neig-1)
      mpa = mp2 + id*(neig-1)
      jtop = jtop + (1-id)*(neig-1)
      itop = (1-id)*mq2
   end if

   lgo = .false.
      
   if (y2g/y1g.lt.0.35q0) lgo = .true.

   call calc(y1g, y2g, ism, mp1, lgo)

   ismc = ism + ms1 + mp1 + mq1
   
   do ka=1, mpa
      ia = ka - mp1

      do kb=1, mpb
         ib = kb - mq1
         jat(ka, kb) = icnt

         do kc=1, ms2
            ic = kc - ms1
            it = ia + ib + ic
            
            if (it.gt.ism) exit

            icnt = icnt + 1
            at(icnt) = 1.q80

            if (it.lt.imin) cycle

            if (.not.lext) then
               at(icnt) = t(ia+1, ib+1, ic+1)
               cycle
            end if

            if (ic.gt.-1) at(icnt) =txx(ia+1, ib+1, ic+1, 2) + at0(ib+ic+2, ia+2)
            if (ic.eq.-1) at(icnt) = t(ia+1, ib+1, ic+1)
            if (ic.eq.0) at(icnt) = 0.q0
         end do
      end do
   end do

   icnt0 = icnt + 1

   do i=1, 82
      icnt = icnt + 1
      at(icnt) = 1.q90*icnt
      do j=1, 82
         jcs(i, j) = icnt0
      end do
   end do

   ms2c = ms2 - 2
   mq2c = mpb - 1
   mp2c = mpa - 1
   ismc = ism + ms1 + mp1 + mq1

   do ka=2, mp2c
      ia = ka - mp1

      do kb=2, mq2c
         ib = kb - mq1
         jcs(ka, kb) = icnt

         do kc=1, ms2c
            ic = kc - ms1

            if (ka+kb+kc.gt.ismc) exit

            icnt = icnt + 1

            if (ic.gt.0) xxx = 0.5q0*(txx(ia+2, ib, ic+1, 1) + txx(ia, ib+2, ic+1, 2) - txx(ia, ib, ic+3, 3))
            if (ic.lt.0) xxx = 0.5q0*(t(ia+2, ib, ic+1) + t(ia, ib+2, ic+1) - t(ia, ib, ic+3))
            if (ic.eq.0) xxx = dqreal(0.q0)

            at(icnt) = xxx
         end do
      end do
   end do

   icnt = icnt

   if (icnt.gt.198000) then 
      write(*,*) 'ICNT = ', icnt
      stop 42
   end if

   if (jtop.ne.0) then
      if (mpa.ge.nad .or. jtop.ge.nad .or. ms2.ge.ncd) then
         write(*, 02) mpa, jtop, ms2, 1

         stop
      end if

      do ka=1, mpa
         do kb=1, jtop
            !? pointers and counters for atd array. cos integrals are displaced by 19800
            jatd(ka, kb) = iatd
            jcsd(ka, kb) = iatd + 19800

            do kc=1, ms2
               sum = 0.q0
               sumc = 0.q0

               if (ka+kb+kc.gt.ismc-neig+1) cycle

               iatd = iatd + 1
               atd(iatd+19800) = 1.q90
               atd(iatd) = 1.q90

               if (ka+kb+kc.lt.iminc) cycle

               do n=1, neig
                  sumc = sumc + ch(n)*at(jcs(ka, kb+n-1)+kc)
                  sum = sum + ch(n)*t(ka-mp1+1, kb+n-mq1, kc-ms1+1)
               end do

               at(jat(ka, kb)+kc) = sum
               jj = jcs(ka, kb)

               if (jj.ne.0 .and. kc.le.ms2c) at(jj+kc) = sumc

               sumcd = 0.q0
               sumd = 0.q0

               if (namm.gt.0) then
                  do n=1, neig
                     sumcd = sumcd + dqreal(n-1)*ch(n)*at(jcs(ka, kb+n-1)+kc)
                     sumd = sumd + dqreal(n-1)*ch(n)*t(ka-mp1+1, kb+n-mq1, kc-ms1+1)
                  end do
               endif

               atd(iatd+19800) = sumcd
               atd(iatd) = sumd
            end do
         end do
      end do

      if (iatd.ge.19800) then
         write(*,'('' IATD ='', I6)')iatd
         stop
      end if
   end if

   if (itop.ne.0) then
      do kc=1, ms2
         do ka=1, mpa
            do kb=1, itop
               sum = 0.q0
               sumc = 0.q0

               if (ka+kb+kc.gt.ismc-neig+1) exit !!!-NEIG+1 deleted
               if (ka+kb+kc.lt.iminc) cycle

               do n=1, neig
                  sumc = sumc + ch(n)*at(jcs(ka, kb+n-1)+kc)
                  if (jtop.gt.0)sum = sum + ch(n)*at(jat(ka, kb+n-1)+kc)
                  if (jtop.eq.0)sum = sum + ch(n)*t(ka-mp1+1, kb+n-mq1, kc-ms1+1)
               end do

               at(jat(ka, kb)+kc) = sum
               jj = jcs(ka, kb)

               if (jj.eq.0.or.kc.gt.ms2c) cycle
               
               at(jj+kc) = sumc
            end do
         end do
      end do
   end if

   if (ni.le.nh .and. id.ne.0) then
      do kc=1, ms2
         do kb=1, mq2
            do ka=1, mp2
               sum = 0.q0
               sumc = 0.q0

               if (ka+kb+kc.gt.ismc-neig+1) cycle  !!!-NEIG+1 deleted
               if (ka+kb+kc.lt.iminc) cycle

               do n=1, neig
                  sumc = sumc + ch(n)*at(jcs(ka+n-1, kb)+kc)
                  sum = sum + ch(n)*at(jat(ka+n-1, kb)+kc)
               end do
               
               at(jat(ka, kb)+kc) = sum
               jj = jcs(ka, kb)
               if (jj.eq.0.q0.or.ka.gt.mp2c.or.kb.gt.mq2c.or.kc.gt.ms2c) cycle
               at(jj+kc) = sumc
            end do
         end do
      end do
   end if

end subroutine


subroutine genintx(id, nj, ni, kdiv)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)
   integer*4 jat0, jat, jcs, jatd, jcsd
   logical lext, lgo
   real*16 z16, y1g16, y2g16

   common/done/lr12, idone, lext

   common hb(12150*(12150+1)/2), at(198000), at2(84, 84), atl(84, 84, 4), &
   &   atj(84, 84, 4), jat0(84), jat(80, 80), jcs(80, 80), &
   &   zsum(84, 86), bin(84, 84), xxx(84, 9), at0(84, 84), jatd(46, 46), &
   &   jcsd(46, 46), atd(38600), ovb(12150*(12150+1)/2)

   common/h1h/ch(25), neig, nh, iz, lrgl, nspn, nv, l1max, nbet, nd, n1, n2,    &
   &   n1d, n2d, ltrans, iset, namm

   common/f1f/y1g, y2g, mp1, mp2, mq1, mq2, ms1, ms2, msl, isum, lg
   
   common/p1p/cplt(9, 11), csum(11), lmn(11), lmb(11), ict, na, nc, nad, ncd
   
   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), &
   &  fm2(84), f21x(84, 86)
   
   common/x1x/zsub(84, 86), at3(84, 84)

01 format(' MIN. ATD DIMENSIONS ARE', 3i4, i6)

   y1g16 = y1g
   y2g16 = y2g

   atx = dqreal(0.q0)

   if (nj.eq.ni) write(6,'(''|'',$)')
   if (nj.ne.ni) write(6,'(''.'',$)')

   n1 = na - 1
   n2 = nc - 1
   icnt = 0

   do i=1, 82
      icnt = icnt + 1
      at(icnt) = 1.q90
      jat0(i) = 1

      do j=1, 82
         jat(i, j) = 0
      end do
   end do

   iatd = 0

   do i=1, 36
      iatd = iatd + 1
      atd(iatd) = 1.q80
      atd(iatd+19800) = 1.q80
      do j=1, 36
         jatd(i, j) = 1
         jcsd(i, j) = 1
      end do
   end do

   do i=1, 82
      at(i) = 1.q80*i
      do j=1, 82
         x = 1.q90*i*j
         at2(i, j) = 1.q1*x
         at3(i, j) = 1.q2*x
      end do
   end do

   do k=1, 4
      do i=1, 82
         x = i*k*1.q170
         do j=1, 82
            atj(i, j, k) = x
            atl(i, j, k) = x
         end do
      end do
   end do

   ism = isum
   imin = 2*lrgl - 1

   if (kdiv.ge.1) imin = imin - 3
   iminc = imin + mp1 + mq1 + ms1

   mpb = mq2
   mpa = mp2
   mpl2 = mp2
   mql2 = mq2
   mps1 = max0(mp1, ms1)
   mqs1 = max0(mq1, ms1)

   jtop = 0
   itop = 0

   if (nj.le.nh) then
      ism = isum + neig - 1
      mpb = mq2 + neig - 1
      jtop = mq2
   end if

   if (ni.le.nh) then
      ism = ism + neig - 1
      mpb = mpb + (1-id)*(neig-1)
      mpa = mp2 + id*(neig-1)
      jtop = jtop + (1-id)*(neig-1)
      itop = (1-id)*mq2
   end if

   lgo = .false.
   if (y2g/y1g.lt.0.35q0) lgo = .true.

   call calc(y1g, y2g, ism, mp1, lgo)
   
   z16 = t(0, 0, 0)
   ismc = ism + ms1 + mp1 + mq1
   
   !
   !   SUBROUTINE - TAB IS USED FOR POWERS OF R1 OR R2 < -1 AND TC IS USED
   !   POWERS OF R12 < -1.
   !

   do ka=1, mpa
      ia = ka - mp1

      do kb=1, mpb
         ib = kb - mq1
         jat(ka, kb) = icnt

         do kc=1, ms2
            ic = kc - ms1
            it = ia + ib + ic

            if (it.gt.ism) exit

            icnt = icnt + 1
            at(icnt) = 1.q80

            if (it.lt.imin) cycle
            if (it.gt.isum) exit

            if (ic.le.-3) go to 12

            if (ia.gt.-2 .and. ib.gt.02) then
               z = t(ia+1, ib+1, ic+1)
            else 
               if (ic.eq.-3) go to 15
               
               if (ic.ne.-2) then
                  if (ic.ne.-2) then
                     z = tab(ia+1, ib+1, ic+1)
                     go to 13
                  end if
                  
                  z = tlog(ia, ib, ic, sumj)
                  at2(ka, kb) = z
                  z = z + sumj
                  go to 13

               12 continue 
                  if (lgo) go to 14
               end if

            15 continue

               if (ia.eq.-2.or.ib.eq.-2) then
                  z = tc32(ia, ib, ic)
                  else
                  if (ia.ge.-1.and.ib.ge.-1) z = tc3(ia, ib, ic)
                  if (ia.eq.-3.or.ib.eq.-3) z = tc33(ia, ib, ic)
               end if

               at3(ka, kb) = 0.q0
               go to 13
               14 if (ib.lt.1) go to 15
               z = tc3x(ia, ib, ic, sumj)

               if (z.eq.0.q0) go to 15
               at3(ka, kb) = z
               z = z + sumj
            end if

         13 continue
            at(icnt) = z
         end do
      end do
   end do
   
   if (msl.gt.0) then
      if (msl.gt.4) stop 52

      do kc=1, msl
         ic = (kc-1)*2

         do ka=1, mpl2
            ia = ka - mp1

            do kb=1, mql2
               ib = kb - mq1

               if (ia+ib+ic.lt.imin) cycle
               if (ia+ib.lt.-mp1) cycle
               if (ia+ic.lt.-mps1) cycle
               if (ib+ic.lt.-mqs1) cycle
               if (ia.lt.-2.or.ib.lt.-2) cycle

               if (ia+ib+ic.gt.isum) exit

               atj(ka, kb, kc) = tlog(ia, ib, ic, sumj)/dqreal(ic+2)
               atl(ka, kb, kc) = atj(ka, kb, kc) + sumj/dqreal(ic+2)
            end do
         end do
      end do
   end if 

   icnt0 = icnt + 1

   do i=1, 82
      icnt = icnt + 1
      at(icnt) = 1.q90*icnt
      do  j=1, 82
         jcs(i, j) = icnt0
      end do
   end do

   ms2c = ms2 - 2
   mq2c = mpb - 1
   mp2c = mpa - 1

   ismc = ism + ms1 + mp1 + mq1
   key = 1
   rba = min(y1g/y2g, y2g/y1g)

   if (rba.lt.0.35q0.and.y1g.lt.y2g) key = -1

   do ka=2, mp2c
      ia = ka - mp1
      do kb=2, mq2c
         ib = kb - mq1
         jcs(ka, kb) = icnt
         do kc=1, ms2c
            ic = kc - ms1
            if (ka+kb+kc.gt.ismc+1) exit
            icnt = icnt + 1
            
            if (ia.lt.0 .or. ib.lt.1 .or. at3(ka-1, kb+1).eq.0.q0) then
               at(icnt) = 0.5q0*(at(jat(ka+1, kb-1)+kc) + at(jat(ka-1, kb+1)+kc) - at(jat(ka-1, kb-1)+kc+2))
               cycle
            end if

            at(icnt) = 0.5q0*(at(jat(ka+1, kb-1)+kc) + at3(ka-1, kb+1) - zsub(ib+1, ia+2))
            cycle

            if (ia.gt.-2 .or. ib.gt.-2) then
               at(icnt) = 0.5q0*(at(jat(ka+key, kb-key)+kc) + at2(ka-key, kb+key))
               cycle
            end if

            if (ib.lt.0 .and. ia.lt.0) then
               at(icnt) = 0.5q0*(tabb(ia+2, ib, ic+1, 1) + tabb(ia, ib+2, ic+1, 2) - tabb(ia, ib, ic+3, 3))
               cycle
         
            end if

            if (ic.eq.-1) then
               at(icnt) = 0.5q0*(t(ia+2, ib, ic+1) + t(ia, ib+2, ic+1) - t(ia, ib, ic+3))
            else
               at(icnt) = 0.5q0*(txx(ia+2, ib, ic+1, 1) + txx(ia, ib+2, ic+1, 2) - txx(ia, ib, ic+3, 3))
            end if
         end do
      end do
   end do

   icnt = icnt

   if (icnt.gt.198000) then
      write(6,*) 'ICNT = ', icnt
      stop 42
   end if

   if (jtop.ne.0) then
      if (mpa.ge.nad .or. jtop.ge.nad .or. ms2.ge.ncd) then
         write(6, 01) mpa, jtop, ms2, 1

         stop
      end if
      
      do ka=1, mpa
         do kb=1, jtop
            !? pointers and counters for atd array. cos integrals are displaced by 19800
            jatd(ka, kb) = iatd
            jcsd(ka, kb) = iatd + 19800

            do kc=1, ms2
               sum = 0.q0
               sumc = 0.q0

               if (ka+kb+kc.gt.ismc-neig+2) exit
               
               iatd = iatd + 1

               atd(iatd+19800) = 1.q90
               atd(iatd) = 1.q90

               if (ka+kb+kc.lt.iminc) cycle

               do n=1, neig
                  sumc = sumc + ch(n)*at(jcs(ka, kb+n-1)+kc)
                  sum = sum + ch(n)*t(ka-mp1+1, kb+n-mq1, kc-ms1+1)
               end do

               at(jat(ka, kb)+kc) = sum
               jj = jcs(ka, kb)

               if (jj.ne.0 .and. kc.le.ms2c) at(jj+kc) = sumc

               sumcd = 0.q0
               sumd = 0.q0

               if (namm.gt.0) then
                  do n=1, neig
                     sumcd = sumcd + (n-1.q0)*ch(n)*at(jcs(ka, kb+n-1)+kc)
                     sumd = sumd + (n-1.q0)*ch(n)*t(ka-mp1+1, kb+n-mq1, kc-ms1+1)
                  end do
               end if

               atd(iatd+19800) = sumcd
               atd(iatd) = sumd
            end do
         end do
      end do

      if (iatd.ge.19800) then
         write(6,'('' IATD ='', I6)')iatd
         stop
      end if
   end if
   
   if (itop.ne.0) then
      do kc=1, ms2
         do ka=1, mpa
            do kb=1, itop
               sum = 0.q0
               sumc = 0.q0

               if (ka+kb+kc.gt.ismc-neig+2) exit
               if (ka+kb+kc.lt.iminc) cycle

               do n=1, neig
                  sumc = sumc + ch(n)*at(jcs(ka, kb+n-1)+kc)

                  if (jtop.gt.0)sum = sum + ch(n)*at(jat(ka, kb+n-1)+kc)
                  if (jtop.eq.0)sum = sum + ch(n)*t(ka-mp1+1, kb+n-mq1, kc-ms1+1)
               end do

               at(jat(ka, kb)+kc) = sum
               jj = jcs(ka, kb)

               if (jj.eq.0.or.kc.gt.ms2c) cycle
               
               at(jj+kc) = sumc
            end do
         end do
      end do
   end if

   if (ni.le.nh .and. id.ne.0) then
      do kc=1, ms2
         do kb=1, mq2
            do ka=1, mp2
               sum = 0.q0
               sumc = 0.q0

               if (ka+kb+kc.gt.ismc-neig+2) cycle
               if (ka+kb+kc.lt.iminc) cycle

               do n=1, neig
                  sumc = sumc + ch(n)*at(jcs(ka+n-1, kb)+kc)
                  sum = sum + ch(n)*at(jat(ka+n-1, kb)+kc)
               end do

               at(jat(ka, kb)+kc) = sum
               jj = jcs(ka, kb)

               if (jj.eq.0.or.ka.gt.mp2c.or.kb.gt.mq2c.or.kc.gt.ms2c) cycle

               at(jj+kc) = sumc
            end do
         end do
      end do
   end if

end subroutine

!* same as tab for kk=0. cancelling terms are omitted for kk=1, 2, 3
function tabb(ia, ib, ic, kk)
   
   use dqmodule

   implicit type (dq_real) (a-h, o-z)
   integer*4 jat0, jat, jcs, jatd, jcsd
   integer p1, q1, p1p, q1p

   common hb(12150*(12150+1)/2), at(198000), at2(84, 84), atl(84, 84, 4), &
   &   atj(84, 84, 4), jat0(84), jat(80, 80), jcs(80, 80), &
   &   zsum(84, 86), bin(84, 84), x(84, 9), at0(84, 84), jatd(46, 46), &
   &   jcsd(46, 46), atd(38600), ovb(12150*(12150+1)/2)

   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), &
   &  fm2(84), f21x(84, 86)

   tabx = 0.q0
   ic1 = ic + 1

   if (ic1.eq.0) return

   ic2 = ic + 2
   ia1 = ia + 1
   ib1 = ib + 1
   iac3 = ia + ic + 3
   ibc3 = ib + ic + 3
   icg1 = ic/2 + 1

   if ((icg1+igc1-ic1) .gt. 0) then
      iabc3 = ia + ibc3

      if (iabc3.le.0) then
         tabb = 1.q160*(ia-ib+10)
         return
      end if

      if (ib1.le.0) then
         do k=1,-ib, 2
            tabx = tabx + bin(ic2, ic2+1-k)*fac(iac3-k)*suml(iabc3, 1, 1-k-ib) *y1x(iabc3+1)
         end do
      end if

      if (ia1.le.0) then
         do k=1,-ia, 2
            tabx = tabx + bin(ic2, ic2+1-k)*fac(ibc3-k)*suml(iabc3, 2, 1-k-ia)*y2x(iabc3+1)
         end do
      end if
   end if

   do i=1, icg1
      i2 = i + i
      q1 = ibc3 - i2
      p1 = ia + i2
      q1p = iac3 - i2
      p1p = ib + i2
      sum1 = 0.q0
      sum2 = 0.q0

      if (i2.le.2 .and. kk.ne.0) then
         if (q1.gt.0.and.p1.gt.0.and.kk.eq.1) sum1 = zsum(85-q1, 87-p1)
         if (q1p.gt.0.and.p1p.gt.0.and.kk.eq.2) sum2 = zsum(q1p, p1p)
         cycle
      end if

      if (q1.gt.0.and.p1.gt.0) sum1 = zsum(85-q1, 87-p1)
      if (q1p.gt.0.and.p1p.gt.0) sum2 = zsum(q1p, p1p)
      tabx = tabx + bin(ic2, i2)*(sum1+sum2)
   end do

   tabb = (tabx/dqreal(ic1))*2.q0

end function


function fpl(iat, ibt, ik)
   use dqmodule

   implicit type (dq_real) (a-h, o-z)
   integer*4 jat0, jat, jcs, jatd, jcsd
   dimension p(9, 9)

   common hb(12150*(12150+1)/2), at(198000), at2(84, 84), atl(84, 84, 4), &
   &   atj(84, 84, 4), jat0(84), jat(80, 80), jcs(80, 80), &
   &   zsum(84, 86), bin(84, 84), x(84, 9), at0(84, 84), jatd(46, 46), &
   &   jcsd(46, 46), atd(38600), ovb(12150*(12150+1)/2)
   
   common/f1f/y1g, y2g, mp1, mp2, mq1, mq2, ms1, ms2, msl, isum, lg
   
   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), &
   &  fm2(84), f21x(84, 86)
   
   common/p1p/cplt(9, 11), csum(11), lmn(11), lmb(11), ict, na, nc, nad, ncd

   fpl = 0.q0
   lmax = lmb(ik)
   lmin = lmn(ik)
   ic = ict - ms1

   if (lg.gt.0) then
      ic2 = ishft(ic+2,-1)
      if ((iand(ic, 1)).le.0 .and. ic2.lt.lmax) lmax = ic2
   end if

   if (lmin.gt.lmax) return
   
   ial = iat
   ibl = ibt
   icl = ic/2 + 1
   lm2 = lmax - 2
   idz = 2 - ic

   if ((iand(lmin, 1)).le.0) then
      2 continue
      if (lg.ne.0) then
         p(2, 1) = at(jcs(iat, ibt)+ict)
         if ((lm2-1).le.0) then
            do l=lmin, lmax, 2
               fpl = fpl + cplt(l, ik)*p(l, 1)
            end do
            return
         end if

         id = 1
         xx = 1.q0

         do i=2, lm2, 2
            i1 = i - 1
            id = id + 1

            if ((id+id-idz).eq.0) then
               p(1, i) = atl(ial-i1, ibl-i1, icl+i1)*xx
               x1 = 0.5q0*xx
               x2 = 0.25q0*xx
               p(2, i+1) =  x1*csl(ial-i, ibl-i, icl+i) - x2*at(jcs(iat-i, ibt-i)+ict+2*i)
               id = 2

               go to 63
            end if

            xx = x(ict, id)
            p(1, i) = at(jat(iat-i1, ibt-i1)+ict+i1+i1)*xx
            id = id + 1

            if ((id+id-idz).eq.0) then
               62 p(2, i+1) = csl(ial-i, ibl-i, icl+i)*xx
               x1 = xx
               x2 = 0.q0
               id = 0

               go to 63
            end if

            xx = x(ict, id)
            p(2, i+1) = at(jcs(iat-i, ibt-i)+ict+i+i)*xx
            i1 = i + 2

            do j=3, i1
               k = i - j + 3
               p(j, k) = dqreal(j+j-3)*p(j-1, k+1) + p(j-2, k)
            end do
         end do

         do l=lmin, lmax, 2
            fpl = fpl + cplt(l, ik)*p(l, 1)
         end do
         return
      end if

      !!  IF LG.NE.0, COMPUTE INTEGRAL FOR (R1**IA)*(R2**IB)*(R12**IC)*LN(R12)

      p(2, 1) = csl(ial, ibl, icl)

      if (lm2.le.0) then
         do l=lmin, lmax, 2
            fpl = fpl + cplt(l, ik)*p(l, 1)
         end do
         return
      end if


      p(1, 2) = 0.5q0*atl(ial-1, ibl-1, icl+1) - 0.25q0*at(jat(iat-1, ibt-1)+ict+2)

      x1 = 0.125q0
      x2 = 0.09375q0
      
      p(2, 3) = x1*csl(ial-2, ibl-2, icl+2) - x2*at(jcs(iat-2, ibt-2)+ict+4)
      i = 2
      id = 4

      63 continue 

      do ip=i, lm2, 2
         if ((ip-i).eq.0) then
            i1 = ip - 1
            id = id + 2
            x1 = x1/dqreal(id)
            x2 = (x1+x2)/dqreal(id)

            p(1, ip)= x1*atl(ial-i1, ibl-i1, icl+i1) - x2*at(jat(iat-i1, ibt-i1)+ict+2*i1)

            id = id + 2
            x1 = x1/dqreal(id)
            x2 = (x1+x2)/dqreal(id)

            p(2, ip+1)= x1*csl(ial-ip, ibl-ip, icl+ip) - x2*at(jcs(iat-ip, ibt-ip)+ict+2*ip)
         end if

         i1 = ip + 2

         do j=3, i1
            k = ip - j + 3
            p(j, k) = dqreal(2*j-3)*p(j-1, k+1) + p(j-2, k)
         end do
      end do

      do l=lmin, lmax, 2
         fpl = fpl + cplt(l, ik)*p(l, 1)
      end do
      return

   end if

   if (lg.ne.0) then
      if ((lmax-2).le.0) then
         fpl = cplt(1, ik)*at(jat(iat, ibt)+ict)
         return
      end if

      id = 0
      xx = 1.q0

      do i=3, lmax, 2
         i1 = i - 3
         i2 = i - 2
         id = id + 1

         if ((id+id-idz).eq.0) go to 50

         xx = x(ict, id)
         p(1, i2) = at(jat(iat-i1, ibt-i1)+ict+i1+i1)*xx
         id = id + 1

         if ((id+id-idz).eq.0) go to 51

         xx = x(ict, id)
         p(2, i2+1) = at(jcs(iat-i2, ibt-i2)+ict+i2+i2)*xx
         do j=3, i
            k = i - j + 1
            p(j, k) = dqreal(j+j-3)*p(j-1, k+1) + p(j-2, k)
         end do
      end do

      do l=lmin, lmax, 2
         fpl = fpl + cplt(l, ik)*p(l, 1)
      end do
      return
   end if

   if ((lmax-2).le.0) then
      fpl = cplt(1, ik)*atl(ial, ibl, icl)
      return
   end if

   xx = 1.q0
   i1 = 0
   i2 = 1
   i = 3
   
50 continue

   p(1, i2) = atl(ial-i1, ibl-i1, icl+i1)*xx
   x1 = 0.5q0*xx
   x2 = 0.25q0*xx
   p(2, i2+1) =  x1*csl(ial-i2, ibl-i2, icl+i2) - x2*at(jcs(iat-i2, ibt-i2)+ict+2*i2)
   id = 2

   go to 52

51 continue

   p(2, i2+1) = csl(ial-i2, ibl-i2, icl+i2)*xx
   x1 = xx
   x2 = 0.q0
   id = 0

52 continue 
      
   do ip=i, lmax, 2
      if ((ip-i).ne.0) then
         i1 = ip-3
         i2 = ip-2
         id = id + 2

         x1 = x1/dqreal(id)
         x2 = (x1+x2)/dqreal(id)

         p(1, i2) = atl(ial-i1, ibl-i1, icl+i1)*x1 - x2*at(jat(iat-i1, ibt-i1)+ict+2*i1)

         id = id + 2

         x1 = x1/dqreal(id)
         x2 = (x1+x2)/dqreal(id)

         p(2, i2+1) = csl(ial-i2, ibl-i2, icl+i2)*x1 - x2*at(jcs(iat-i2, ibt-i2)+ict+2*i2)
      end if

      do j=3, ip
         k = ip - j + 1
         p(j, k) = dqreal(j+j-3)*p(j-1, k+1) + p(j-2, k)
      end do
   end do

   do l=lmin, lmax, 2
      fpl = fpl + cplt(l, ik)*p(l, 1)
   end do

end function


function csl(i, j, k)
   
   use dqmodule

   implicit type (dq_real) (a-h, o-z)
   integer*4 jat0, jat, jcs, jatd, jcsd

   common hb(12150*(12150+1)/2), at(198000), at2(84, 84), atl(84, 84, 4), &
   &   atj(84, 84, 4), jat0(84), jat(80, 80), jcs(80, 80), &
   &   zsum(84, 86), bin(84, 84), x(84, 9), at0(84, 84), jatd(46, 46), &
   &   jcsd(46, 46), atd(38600), ovb(12150*(12150+1)/2)
   
   common/f1f/y1g, y2g, mp1, mp2, mq1, mq2, ms1, ms2, msl, isum, lg
   
   if (i.lt.2.or.j.lt.2.or.k.lt.1) then
      csl = 1.9q0
      return
   endif

   if (j.gt.mp1-1) then
      csl = 0.5q0*(atj(i-1, j+1, k) + atl(i+1, j-1, k) - atj(i-1, j-1, k+1))
      return
   end if
   
   csl = 0.5q0*(atl(i-1, j+1, k) + atj(i+1, j-1, k) - atj(i-1, j-1, k+1))

end function


function fpld(iat, ibt, ik)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)
   integer*4 jat0, jat, jcs, jatd, jcsd
   dimension p(9, 9)

   common hb(12150*(12150+1)/2), at(198000), at2(84, 84), atl(84, 84, 4), &
   &   atj(84, 84, 4), jat0(84), jat(80, 80), jcs(80, 80), &
   &   zsum(84, 86), bin(84, 84), x(84, 9), at0(84, 84), jatd(46, 46), &
   &   jcsd(46, 46), atd(38600), ovb(12150*(12150+1)/2)

   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), &
   &  fm2(84), f21x(84, 86)

   common/p1p/cplt(9, 11), csum(11), lmn(11), lmb(11), ict, na, nc, nad, ncd

   fpld = 0.q0
   lmax = lmb(ik)

   if (lmax.eq.0) return

   lmin = lmn(ik)
   ic2 = ict/2

   if (((ict+1)/2).le.0) then
      if (ic2.lt.lmax) lmax = ic2
      if (lmin.gt.lmax) return
   end if

   if (((lmin+1)/2 - lmin/2).le.0) then
      p(2, 1) = atd(jcsd(iat, ibt)+ict)
      lmx = lmax - 2

      if ((lmx-1).le.0) then
         fpld = cplt(2, ik)*p(2, 1)
         return
      end if

      id = 1

      do i=2, lmx , 2
         i1 = i - 1
         id = id + 1

         p(1, i) = atd(jatd(iat-i1, ibt-i1)+ict+i1+i1)*x(ict, id)

         id = id + 1

         p(2, i+1) = atd(jcsd(iat-i, ibt-i)+ict+i+i)*x(ict, id)
         
         i1 = i + 2
         i2 = i+3
         
         do j=3, i1
            k = i2 - j
            p(j, k) = dqreal(j+j-3)*p(j-1, k+1) + p(j-2, k)
         end do
      end do

      do l=lmin, lmax, 2
        fpld = fpld + cplt(l, ik)*p(l, 1)
      end do

      return
   end if
 
   if ((lmax - 2).le.0) then
      fpld = cplt(1, ik)*atd(jatd(iat, ibt)+ict)
      return
   end if

   id = 0

   do i=3, lmax, 2
      i1 = i - 3
      i2 = i - 2
      id = id + 1
      
      p(1, i2) = atd(jatd(iat-i1, ibt-i1)+ict+i1+i1)*x(ict, id)
      
      id = id + 1
      
      p(2, i2+1) = atd(jcsd(iat-i2, ibt-i2)+ict+i2+i2)*x(ict, id)
      
      i1 = i+1
      
      do j=3, i
         k = i1 - j
         p(j, k) = dqreal(j+j-3)*p(j-1, k+1) + p(j-2, k)
      end do
   end do

   do l=lmin, lmax, 2
      fpld = fpld + cplt(l, ik)*p(l, 1)
   end do

end function


subroutine basis(md, lgo, ifirst, kdiv, knv, bscale)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)
   integer ss, pam, qam
   logical ltry, lgo, lksm1
   real*16 cpd116, cpd216, cpd1x16, cpd2x16, cp116, cp216, der16, b16, temp1, temp2
   character status*8
   dimension cpd1x16(8), cpd2x16(8), cp116(12150), cp216(12150), cpd116(8), cpd216(8), & 
   &   il2(8), md(12150), knv(2), &
   &   der16(9, 2), b16(8, 2)
   
   common/b1b/cpd1x(8), cpd2x(8), nit, next, kst, key1, ibasis, kbasis, linc, &
   &   nv1, longnv, ksm(8), mar12, mar1, kono, nbetx, status
   
   common/s1s/faca(80), facb(80), ngv(100, 10), ng(100), &
   &   igrp(698), ngp, nngt
   
   common/r1r/rap, egs, ltrid
   
   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), &
   &  fm2(84), f21x(84, 86)
   
   common/h1h/ch(25), neig, nh, iz, lrgl, nspn, nv, l1mx, nbet, nd, n1, n2, &
   &   n1d, n2d, ltrans, iset, namm
   
   common/f1f/y1g, y2g, mp1, mp2, mq1, mq2, ms1, ms2, msl, isum, lg
   
   common/p1p/cpl(9, 11), csum(11), lmin(11), lmax(11), ss, na, nc, nad, ncd
   
   common/d1d/cp1(12150), cp2(12150), cpd1(8), cpd2(8), b(8, 2), der(9, 2), &
   &   ndp(12150), ndq(12150), nds(12150), ld1(12150), ld2(12150), pam, qam, maxc, &
   &   ia, nblk(9), nb, nblkx(9), nbx, maxbf, nrow(8998), lngth

   common/maxpow/iia(9), jja(9), kka(9), ijka(9)

   save l1max, kono1, kono2

01 format(2f14.5, i2, i1, 6i3)
02 format(4h z =, i3,'  L =', i2,'  S =', i2,'  NV =', i5, '  BASIS SET SIZES ARE', 6i5,$)
44 format(6d12.5)
48 format(6f12.5)
68 format(56h exit called due to looping in basis function generation)
70 format(35h array size of hm will be exceeded. )
96 format("28h tables of l1, l2, p, q, s, k1, k2")
97 format(i5, 3x, 5i5, 2f12.5)

   z = dqreal(iz)
   ltry = .true.
   nbetx = nbet + ibasis

   if (next.le.1) then
      na = n1 + 1
      nc = n2 + 1
      nad = n1d + 1
      ncd = n2d + 1

      linc = 0
      l1max = l1mx

      read(3,*) (cpd1x16(i), cpd2x16(i), i=1, nbet)

      do i=1, nbet
         cpd1x(i)=cpd1x16(i)
         cpd2x(i)=cpd2x16(i)
      end do

      if (ibasis.eq.1 .and. kbasis.eq.0 .and. cpd1x(2).eq.1.q0) then
         !? compress alphas and betas if block 2 is present in list
         backspace(3)
         write(*,*) 'Compressing'
         read(3,*) (cpd1x16(i), cpd2x16(i), i=1, nbet+1)

         do i=1, nbet + 1
            cpd1x(i)=cpd1x16(i)
            cpd2x(i)=cpd2x16(i)
         end do

         do i=2, nbet
            cpd1x(i) = cpd1x(i+1)
            cpd2x(i) = cpd2x(i+1)
         end do
      end if

      read(3,*, err=31) (ksm(i), i=1, nbetx), key1, mar12, mar1, kono

      if (ksm(1).eq.0) then
         do while (.true.)
            !? fill in missing alphas and betas if necessary
         31 if (ltry) then
               ltry = .false.
               ktry = 0
               do while (.true.)
                  write(*,*) 'Fill-in'
                  write(*,'(12i4)') (ksm(i), i=1, nbetx), key1, mar12, mar1, kono
                  backspace(3)
                  backspace(3)

                  ktry = ktry + 1

                  if (ktry.gt.5) stop

                  read(3,'(a80)') line
                  write(*,'(a80)') line
                  backspace(3)
                  read(3,*, err=36) (ksm(i), i=1, nbetx), key1, mar12, mar1, kono

                  if (ksm(1).ne.0) exit

               36 backspace(3)
                  write(*,*) 'Reading error with NBETX =', nbetx
               end do

               backspace(3)

               do i=1, nbet-1
                  ii = nbet - i + 1
                  cpd1x(ii) = cpd1x(ii-1)
                  cpd2x(ii) = cpd2x(ii-1)
               end do
                  
               do i=1, nbet
                  cpd1x16(i)=cpd1x(i)
                  cpd2x16(i)=cpd2x(i)
               end do

               write(*,'(2f12.5)') (cpd1x16(i), cpd2x16(i), i=1, nbet)
               write(*,'(12i4)') (ksm(i), i=1, nbetx), key1, mar12, mar1, kono
            else
               write(*,'(2f12.5)') (cpd1x16(i), cpd2x16(i), i=1, nbet)
               write(*,'(12i4)') (ksm(i), i=1, nbetx), key1, mar12, mar1, kono
               write(*,*) 'Input file error.'
               stop 31
            endif
            
            read(3,*, err=31) (ksm(i), i=1, nbetx), key1, mar12, mar1, kono
            if (ksm(1).ne.0) exit
         end do
      end if

      kss =10

      !? if key1 = 2, use nonlinear parameters from previous calculation
      if (key1.ne.2) then
         do i=1, nbet
            cpd1x16(i)=cpd1x(i)
            cpd2x16(i)=cpd2x(i)
            temp1 = nint(16384*cpd1x16(i))
            cpd1(i)= dqreal(temp1/16384)
            temp2 = nint(16384*cpd2x16(i))
            cpd2(i)= dqreal(temp2/16384)
         end do
      end if

      if (kono.eq.0) kono = 2
      kono2 = kono/10
      kono1 = kono - 10*kono2

      der(1, 1) = 0.q0

      if (key1.eq.1) then
         ifirst = 0

         read(3, 44) ((der16(i, k), k=1, 2), i=1, nbet)
         read(3, 48) ((b16(i, k), k=1, 2), i=1, nbet)
         read(3, 48) (cpd116(i), cpd216(i), i=1, nbet)

         do i=1, nbet
            cpd1(i)=cpd116(i)
            cpd2(i)=cpd216(i)
            
            do k=1, 2
               der(i, k)=der16(i, k)
               b(i, k)=b16(i, k)
            end do
         end do
      end if
   end if

   c = dqreal(iz-1)/dqreal((lrgl+neig)*iz)

   if (iset.eq.1) cpd1(1) = dqreal(1.q0)/dqreal((linc+1))
   if (bscale.gt.0.q0) cpd2(1) = cpd2(1)*bscale

   do i=1, nbet
      ii = nbet - i + 1
      
      cpd116(ii) = cpd1(ii)
      cpd216(ii) = cpd2(ii)
            
      temp1 = nint(16384*cpd116(ii))
      cpd1(ii)=dqreal(temp1/16384)
      temp2 = nint(16384*cpd216(ii))
      cpd2(ii) = dqreal(temp2/16384)
      
      cpd1x(ii) = cpd1(ii)
      cpd2x(ii) = cpd2(ii)
      
      cpd1(ii+ibasis) = cpd1(ii)
      cpd2(ii+ibasis) = cpd2(ii)
   end do

   if (cpd2(1).eq.c) cpd2(1)=cpd2(1)+1.q-12

   l1max = l1mx + 1
   kdiv1 = kdiv + 1
   if (kdiv1.gt.2) kdiv1 = 2

   do kkdiv=1, kdiv1
      j = 1
      if (nh.ne.0) j = 2

      m = 0
      lbasis = 0

      if (lrgl-linc.gt.0) lbasis = max0(ibasis, kbasis)

      do i3=1, l1max
         mmax = 1
         if (i3.eq.1) mmax = nbetx - l1mx

         do mm=1, mmax
            m = m + 1

            if (ksm(m).ne.0) then

               ld1(j)=i3-1

               if (m.eq.2) ld1(j) = ld1(j) + lbasis

               ld2(j)=lrgl-ld1(j)
               mj=100

               if ((ld1(j)-ld2(j)).eq.0)  mj = 0
               npm=ld1(j)+1
               nqm=ld2(j)+1

               k11 = 2 - kkdiv
               kr1 = 2 - kdiv

               if (k11.gt.0) k11 = 1
               if (k11.eq.1) kr1 = 0
               if (k11.eq.-1) kr1 = 1
               
               ksum = k11

               do while (.true.)
                  do k1=k11, ksum
                     ks=ksum-k1+1
                     k2max = ks

                     do k2=1, k2max
                        k1m=k2-1
                        ld1(j)=i3-1

                        if (m.eq.2) ld1(j) = ld1(j) + lbasis

                        ld2(j)=lrgl-ld1(j)
                        ndp(j)=npm+k2-1
                        ndq(j)=nqm+ksum-k1-k2+1
                        nds(j)=k1-1

                        if (kr1.eq.1) then
                           if (ld1(j).gt.0) ndp(j) = ndp(j) -1
                           if (ld2(j).gt.0) ndq(j) = ndq(j) -1
                           nds(j) = nds(j) + 1
                        endif

                        cp1(j)=cpd1(m)
                        cp2(j)=cpd2(m)
                        mj2=2*(mj-ndp(j)+ndq(j))

                        !? removes duplicate of screened hydrogenic term
                        if (nh.ne.0 .and. m.eq.1 .and. ndp(j).eq.npm .and. ndq(j).eq.nqm .and. nds(j).eq.0) cycle

                        if (m.ne.1 .and. nds(j).eq.-1) cycle

                        !? removes terms above upper limit for divided basis sets
                        if (ndp(j).gt.npm+ksm(m)) cycle
                        if (ndq(j).gt.nqm+ksm(m)) cycle
                        if (nds(j).gt.ksm(m)) cycle
                        
                        if (m.eq.1 .and. nds(j).gt.mar12) cycle
                        if (m.eq.1 .and. ndp(j).gt.mar1+1) cycle
                        if (m.eq.2 .and. lbasis.ne.0 .and. ndp(j).gt.mar1+ld1(j)+1) cycle
                        if (m.eq.2 .and. lbasis.ne.0 .and. nds(j).gt.mar12) cycle
                        
                        if (m-lbasis.ge.2 .and. &
                        &  nds(j).gt.kono1 .and. &
                        &  ksum-iabs(ld1(j)-ld2(j)) + iabs(ndp(j)-ndq(j)).gt.ksm(m)+1-kono2 .and. &
                        &  kono.ge.0) cycle
                        
                        k4 = 4
                        km2 = 2
                        
                        if (lrgl.eq.0) km2 = 1
                        if (lrgl.eq.0) k4 = 4000
                        if (lrgl.ne.0.or.neig.gt.2) kss = -1
                        
                        if (m.gt.km2 .and. &
                        &  ksum.le.k4 .and. &
                        &  nspn.eq.1 .and. &
                        &  ld1(j).eq.ld2(j) .and. &
                        &  ndp(j).eq.ndq(j) .and. &
                        &  (neig.eq.1.or.lrgl.eq.0)) cycle
                        
                        if (m.eq.1 .and. &
                        &  j.gt.1 .and. &
                        &  ndp(j)+ndq(j).le.kss .and. &
                        &  nspn.eq.1 .and. &
                        &  ld1(j).eq.ld2(j) .and. &
                        &  ndp(j).eq.ndq(j) .and. &
                        &  (neig.eq.1.or.lrgl.eq.0)) cycle
                        
                        if (ndq(j).ge.50) then
                           write(6, 68)
                           stop
                        end if

                        if ((j-1).gt.0 .and. &
                        &  (iabs(ndp(j)-ndp(j-1))+iabs(ndq(j)-ndq(j-1))+iabs(nds(j)-nds(j-1))).le.0) cycle
         
                        if (mj2.lt.0) cycle
                        
                        j=j+1
                     end do
                  end do

                  ksum=ksum+1
                  
                  if (ksum.gt.ksm(m)+1) exit
               end do
            end if

            il2(m) = j-1
         end do
      end do

      knv(kkdiv) = j - 1
      knv(2) = j - 1
   end do

   nv = j-1
   
   il2(6) = il2(6) - il2(5)
   il2(5) = il2(5) - il2(4)
   il2(4) = il2(4) - il2(3)
   il2(3) = il2(3) - il2(2)
   il2(2) = il2(2) - il2(1)
   
   do j=1, nv
      ndp(j) = ndp(j) - 1 + linc
      ndq(j) = ndq(j) - 1 + linc
      ld1(j) = ld1(j) + linc
      ld2(j) = ld2(j) + linc
   end do
   
   lrgl = lrgl + linc
   call calc(dqreal(1.q0), dqreal(1.q0), 0, 0,.false.)
   !? calculation of screened hydrogenic wavefunctions

   if (nh.ne.0) then
      n = lrgl + neig
      l = linc
      xx= 1.q0
      
      do m=1, neig
         ch(m) = xx*dqreal(2**(l+lrgl+2))*(dqreal(iz-1))**(lrgl+1)*sqrt(dqreal(iz-1) &
         &  *dqreal(fac(n+lrgl+1)/(dqreal(iz)*(fac(n-lrgl))*fac(2*l+2)))) &
         &  /(fac(2*lrgl+2)*dqreal(n)*(dqreal(l+1))**(l+2)*(dqreal(n)*z)**(lrgl+1))
         
         xx = -xx*2.q0*dqreal(iz-1)*dqreal(neig-m)/((2.q0*lrgl+m+1.q0) &
         &  *dqreal(m*n*iz))
      end do

      ld1(1) = ld1(2)
      ld2(1) = ld2(2)
      
      cp1(1) = dqreal(1.q0)/dqreal(linc+1)
      cp2(1) = dqreal(iz-1.q0)/dqreal(n*iz)
      
      ndp(1) = linc
      ndq(1) = lrgl
      nds(1) = 0
   end if

   write(4, 02) iz, lrgl, nspn, nv,(il2(i), i=1, nbetx)
   write(6, 02) iz, lrgl, nspn, nv,(il2(i), i=1, nbetx)
   write(4,*) nv, nd, next
   
   if (nv.gt.nd) then 
      write(4, 70)
      stop
   end if

   if (next.le.1) then
      write(4, 96)

      do j=1, 2
         cp116(j)=cp1(j)
         cp216(j)=cp2(j)
      end do

      write(4, 97)(j, ld1(j), ld2(j), ndp(j), ndq(j), nds(j), cp116(j), cp216(j), j=1, 2)
      if (lgo)write(*, 97)(j, ld1(j), ld2(j), ndp(j), ndq(j), nds(j), cp116(j), cp216(j), j=1, 2)

      k1 = 0

      do k=1, nbetx
         k1 = k1 + il2(k)
         k2 = k1 + 1

         if (k2.gt.nv) k2 = nv

         do j=k1, k2
            cp116(j)=cp1(j)
            cp216(j)=cp2(j)
         end do

         write(4, 97)(j, ld1(j), ld2(j), ndp(j), ndq(j), nds(j), cp116(j), cp216(j), j=k1, k2)
         if (lgo)write(*, 97)(j, ld1(j), ld2(j), ndp(j), ndq(j), nds(j), cp116(j), cp216(j), j=k1, k2)
      end do
   end if
   
   pam = ndp(1)
   qam = ndq(1)
   maxc = nds(1)
   ia = ndp(1) + ndq(1) + nds(1)
   nh1 = nh + 1

   do i=nh1, nv
      if (ndp(i).gt.pam) pam = ndp(i)
      if (ndq(i).gt.qam) qam = ndq(i)
      if (nds(i).gt.maxc) maxc = nds(i)
      ia1 = ndp(i) + ndq(i) + nds(i)
      if (ia1.gt.ia) ia = ia1
   end do

   ia = 2*ia
   nb = 1
   nblk(1) = 0
   nbx = 1
   nblkx(1) = 0

   do i=2, nv
      if (ld1(i).eq.ld1(i-1) .and. cp2(i).eq.cp2(i-1) .and. cp1(i).eq.cp1(i-1)) cycle
      nb = nb + 1
      nblk(nb) = i - 1
      if (cp2(i).eq.cp2(i-1) .and. cp1(i).eq.cp1(i-1)) cycle
      nbx = nbx + 1
      nblkx(nbx) = nb - 1
   end do
   
   nblk(nb+1) = nv
   nblkx(nbx+1) = nb

   do i=1, nbet
      cpd1(i) = cpd1(i+ibasis)
      cpd2(i) = cpd2(i+ibasis)
   end do

   do nbl1=1, nb
         nb1 = nblk(nbl1) + 1
         nb2 = nblk(nbl1+1)
         iia(nbl1) = 0
         jja(nbl1) = 0
         kka(nbl1) = 0
         ijka(nbl1) = 0
      do i=nb1, nb2
         if (ndp(i).gt.iia(nbl1)) iia(nbl1) = ndp(i)
         if (ndq(i).gt.jja(nbl1)) jja(nbl1) = ndq(i)
         if (nds(i).gt.kka(nbl1)) kka(nbl1) = nds(i)
         j = ndp(i) + ndq(i) + nds(i)
         if (j.gt.ijka(nbl1)) ijka(nbl1) = j
      end do  
   end do

   !? arrange block 2 basis functions into groups with the same ndp and nd.
   !! ask about below
   !? array ngv(i, j) stores the intex of the j'th element of the i'th gr.
   !? the matrix elements for each group are nearly equal, and the common 
   !? dominant term for each group will be calculated and stored separately.

   nngt = mar12 + 1

   if (nngt.gt.9) ltrans = 0
   if (ltrans.eq.0) then
      if (ifirst.ne.1) write(6,*)' COMB TRANSFORMATION SUPPRESSED.'
      ngp = 0
   else
      ngp = 1
      nb2 = nblk(3)
      k = nblk(2) + 1
      isum = ndp(k) + ndq(k) + nds(k) - 1
   
   21 isum = isum + 1
   
      k = k - 1
      np = ndp(2) - 1
   
   22 np = np + 1
   23 k = k + 1
   
      if (k.lt.nb2) then
         if (ndp(k)+ndq(k)+nds(k).gt.isum) go to 21
         if (ndp(k).ne.np) go to 23
   
         ic = 0
         do i=k, nb2
            if (ndp(i)+ndq(i)+nds(i).gt.isum) exit
            if (ndp(i).ne.np) cycle
   
            ic = ic + 1
            if (ic.gt.9) stop 23
   
            ngv(ngp, ic) = i
            igrp(i-1) = ngp
         end do
   
         if (ic.le.0) go to 22
   
         ng(ngp) = ic
         ngp = ngp + 1
         go to 22
      end if
   
      ngp = ngp - 1
   
      if (ngp.gt.100) then
         write(6,*) ngp, nb2
         stop 200
      endif
   
      if (nb2-1.gt.698) stop 201
   endif

   !? generate starting addresses of each row of triangular matrix
   
   md(1) = 0

   do i=1, nv
      if (i.gt.1) md(i) = md(i-1) + i - 1
   end do

   !? nv1 = total length of linearized triangular array

   nv1 = md(nv) + nv
   ntest = nv1 + 1

   if (ntest.gt.longnv) write(6,*) ntest, longnv
   if (ntest.gt.longnv) stop 797

end subroutine


subroutine ham(psi, psi2, e, amm, loop, md, lder, ldump, kdiv, bscale)
   
   use dqmodule

   implicit type (dq_real) (a-h, o-z)

   logical lgo, lr12, lder, ltxx, lblk2, lext, ldump
   integer pp, qq, ss, pam, qam
   integer*4 jat0, jat, jcs, jatd, jcsd
   real*16 nspm, nspp, cp116, cp216, d, cpl16, deltr16, deltl16, sum16, sumr16, &
   &   sumll16, pp16, qq16, ss16, sumq16, yx16, psi216, tmaj, tmbj, tmai, tmbi, temp1, &
   &   temp2, temp3, bb16, diff16
   character status*4
   dimension d(6)

   common hm(12150*(12150+1)/2), at(198000), at2(84, 84), atl(84, 84, 4), &
   &   atj(84, 84, 4), jat0(84), jat(80, 80), jcs(80, 80), &
   &   zsum(84, 86), bin(84, 84), fill(84, 9), at0(84, 84), jatd(46, 46), &
   &   jcsd(46, 46), atd(38600), ovb(12150*(12150+1)/2)
   
   common/done/lr12, idone, lext
   
   common/s1s/faca(80), facb(80), ngv(100, 10), ng(100), &
   &   igrp(698), ngp, nngt
   
   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &   ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), &
   &   fm2(84), f21x(84, 86)
   
   common/b1b/cpd1x(8), cpd2x(8), nit, next, kst, key1, ibasis, kbasis, linc, &
   &   nv1, longnv, ksm(8), mar12, mar1, kono, nbetx, status
   
   common/h1h/ch(25), neig, nh, iz, lrgl, nspn, nv, l1max, nbet, nd, n1, n2, &
   &   n1d, n2d, ltrans, iset, namm
   
   common/f1f/y1g, y2g, mp1, mp2, mq1, mq2, ms1, ms2, msl, isum, lg
   
   common/p1p/cpl(9, 11), csum(11), lmin(11), lmax(11), ss, na, nc, nad, ncd
   
   common/d1d/cp1(12150), cp2(12150), cpd(8, 2), b(8, 2), der(9, 2), &
   &   ndp(12150), ndq(12150), nds(12150), ld1(12150), ld2(12150), pam, qam, maxc, &
   &   ia, nblk(9), nb, nblkx(9), nbx, maxbf, nrow(8998), lngth
   
   common/maxpow/iia(9), jja(9), kka(9), ijka(9)
   
   common/q1q/cpl16(9, 11)
   
   dimension psi2(1), psi(1), dhm(9, 2), hmx(2), hmxx(2), md(12150), dhmx(9, 2),  &
   &  d2d(2, 4), df1(2, 2), tt(10), jj(13), fdd(6), fddr(3, 2), df1x(2, 2)
   dimension cpd1x16(8), cpd2x16(8), cp116(12150), cp216(12150), cpd116(8), cpd216(8)
   dimension psum(30)

31 format('ANGULAR COEFFICIENTS FOR', 6i3)
32 format(3i3, 5d16.8)
33 format(2i3,' SUMQ =', 3i4, 5d15.7)
52 format('MI, JBAD =', 2i6)
71 format(' Eold =', f30.24/' Enew =', f30.24/' Diff ='f30.24)
72 format(' LOOP =', i3,', DERIVATIVES ARE'/(1p8d12.4))
73 format(8f12.5)


   nrmax = 8998
   fdd= dqreal(6.0q0)
   fddl1 = 2.0q0
   fddl4= 2.0q0

   nbt1 = nbet + 1

   !? ns1 is the principal quantum number of the inner electron and nl2 the principal
   !? quantum number for the outer electron
   
   idone = 0
   
   !? lg = 0 for log(r12) terms
   
   one = dqreal(1.q0)
   lg = 1
   ns1 = ld1(1) + 1
   nl2 = lrgl + neig
   z = dqreal(iz)
   
   esh = dqreal(-0.5q0*(1.q0/ns1**2 + ((z-1.q0)/(z*dqreal(nl2)))**2))
   s = dqreal((-1.q0)**nspn)
   
   amm4 = amm*dqreal(4.q0)
   cp21 = 0.q0
   
   if (nh.ne.0) cp21 = cp2(1)
   
   if (.not. lder) then
      ibad = 0
      jbad = 0

      do i=1, nv
         psi(i) = 0.q0
         mi = md(i)

         do j=1, i
            hm(mi+j) = 0.q0
               if (hm(mi+j).eq.0.q0) then
                  if (.not.ldump) then
                     ovb(mi+j) = 0.q0
                     if (ovb(mi+j).ne.0.q0) go to 37
                  end if
                  cycle
               end if

            37 write(6, 52) mi, j

            stop
         end do
      end do
   else
      do i=1, 9
         do k=1, 2
            dhmx(i, k) = 0.q0
            dhm(i, k) = 0.q0
         end do
      end do
   end if

   mp1 = 2
   ms1 = 2

   if (kdiv.eq.2) ms1 = 4

   idxx = -1
   
   do idex=1, 2
      ide = 3 - idex
      id = ide - 1
      idxx = idxx + 1
      sx = s**id
      irow = idxx*1.q09
      ndmp = 1
      ls1 = -1
      ls3 = -1

      rewind 8
      rewind 9
      rewind 10

      do nhi=1, nbx
         nbx1 = nblkx(nhi) + 1
         nbx2 = nblkx(nhi+1)

         do nhj=1, nhi
            nbx3 = nblkx(nhj) + 1
            nbx4 = nblkx(nhj+1)

            do nbl1=nbx1, nbx2
               if (nhi.eq.nhj) nbx4 = nbl1
               nb1 = nblk(nbl1) + 1
               nb2 = nblk(nbl1+1)

               do nbl2=nbx3, nbx4
                  nb3 = nblk(nbl2) + 1
                  nb4 = nblk(nbl2+1)
                  if (ld2(nb1).ne.ls1 .or. ld2(nb3).ne.ls3 .or. namm.ne.0) then
                     ls1 = ld2(nb1)
                     ls3 = ld2(nb3)
                     lx1=ld1(nb1)
                     lx2=ld2(nb1)

                     if (ide.ne.1) then
                        lx1=ld2(nb1)
                        lx2=ld1(nb1)
                     end if

                  500 if ((nb1.gt.1.and.nb3.gt.1).or.nh.eq.0) call cross(ld1(nb3), ld2(nb3), lx1, lx2, lrgl, id)
                     if ((nb1.eq.1.or.nb3.eq.1).and.nh.ne.0) call cross(lx1, lx2, ld1(nb3), ld2(nb3), lrgl, id)

                     if (loop.le.0 .and. .not.lder) then
                        write(4, 31) ld1(nb3), ld2(nb3), lx1, lx2, lrgl, id
                        do k=1, 11
                           i1 = lmin(k)
                           i2 = lmax(k)
                           if (i2.eq.0) cycle
                           do i=i1, i2, 2
                           cpl16(i1, k)=cpl(i1, k)
                           end do
                           write(4, 32) k, i1, i2,(cpl16(i, k), i=i1, i2, 2)
                        end do
                     end if

                     l1j = ld1(nb3)*(ld1(nb3) + 1)
                     l2j = ld2(nb3)*(ld2(nb3) + 1)
                     l1j4 = 4*l1j
                     l2j4 = 4*l2j
                     l1i = lx1*(lx1 + 1)
                     l2i = lx2*(lx2 + 1)

                     mldd1 = lmin(1)
                     ldd1 = lmax(1)

                     do k=1, 6
                        if (abs(csum(k)-0.5q0).lt.1.q-10) csum(k) = 0.5q0
                        if (abs(csum(k)+0.5q0).lt.1.q-10) csum(k) = -0.5q0
                        do i=1, ldd1+1

                           if (abs(cpl(i, k)-0.5q0).lt.1.q-10) cpl(i, k) = 0.5q0
                           if (abs(cpl(i, k)+0.5q0).lt.1.q-10) cpl(i, k) = -0.5q0
                        end do
                     end do

                     dd1 = cpl(1, 1)
                     dd2 = cpl(2, 1)
                     ldd4 = lmax(4)
                     mldd5 = lmin(5)
                     ldd5 = lmax(5)
                     dd5 = cpl(1, 5)
                     dd52 = cpl(2, 5)

                     ldel1 = l1j - l1i - ldd1*(ldd1-1)
                     ldel2 = l2j - l2i - ldd1*(ldd1-1)
                     ldel3 = l1i - l1j - ldd1*(ldd1-1)
                     ldel4 = l2i - l2j - ldd1*(ldd1-1)
                  end if

                  !? set lr12 = .true. to calculate dominant terms for block 2
                  lr12 = .false.

                  lblk2 = .false.

                  if (id.eq.0.and.nbl1.eq.2.and.nbl2.eq.2.and.ltrans.ne.0) lblk2 = .true.
                  if (lmax(1).eq.1) lr12 = .true.
                  
                  ltxx = .false.

                  !? set lext = .true. for extended precision calculation of basis function integrals
                  lext = ltxx

                  aj = cp1(nb3)
                  bj = cp2(nb3)
                  ai = dqreal(1-id)*cp1(nb1) + dqreal(id)*cp2(nb1)
                  bi = dqreal(1-id)*cp2(nb1) + dqreal(id)*cp1(nb1)
                  
                  !! ask about this
                  !? use fixed point arithmetic to calculate am = aj - ai, bm - bj - bi e
                  tmaj = aj*16384.q0 + 0.1q0
                  maj=tmaj
                  tmbj = dqreal(iz)*bj*16384.q0 + 0.1q0
                  mbj=tmbj
                  tmai = ai*16384.q0 + 0.1q0
                  mai=tmai
                  tmbi = dqreal(iz)*bi*16384.q0 + 0.1q0
                  mbi=tmbi
                  
                  mam = maj - mai
                  nap = maj + mai
                  mbm = mbj - mbi
                  mbp = mbj + mbi
                  
                  am = mam/16384.q0
                  bm = dqreal(mbm/16384.q0)/z
                  ap = aj + ai
                  bp = bj + bi
                  ap2 = (nap - 2*16384)/16384.q0
                  bp2 = ((mbp - 2*16384*(iz-1))/16384.q0)/z
                  apn2 = (ns1*nap -2*16384)/16384.q0
                  bpn2 = ((nl2*mbp - 2*16384*(iz-1))/16384.q0)/z
                  
                  y1g = ap
                  y2g = bp
                  
                  nl22 = nl2*nl2

                  shift = dqreal(nl22)*( - 0.5q0*(bm*bm + am*am &
                  &   + (dqreal(ns1)*ap+2.q0)*apn2/dqreal(ns1*ns1))) &
                  &   - 0.5q0*(bp*dqreal(nl2)+2.q0*dqreal(iz-1)/z)*bpn2

                  e10 = e*4.q0*16384.q0*16384.q0
                  apx = dqreal(nap)
                  ap2x = apx*apx
                  amx = dqreal(mam)
                  bpx =dqreal(nl2*mbp)
                  bp2x = bpx*bpx
                  bmx = dqreal(nl2*mbm)

                  ca1 = 8.q0*16384.q0*apx
                  ca2 = 4.q0*16384.q0*16384.q0/dqreal(ns1*ns1)
                  cb1 = 8.q0*16384.q0*bpx*dqreal(nl2*(iz-1))
                  cb2 = 4.q0*16384.q0*16384.q0*dqreal((iz-1)*(iz-1))
                  cb3 = dqreal(iz*iz*nl2*nl2)

                  dei = -2.q0*(dqreal(nl22)*(dqreal(ns1)*ai+1.q0)*(dqreal(ns1)*ai-1.q0) &
                  &   / dqreal(ns1**2)+(dqreal(nl2*iz)*bi+dqreal(iz-1))*(dqreal(nl2*iz) &
                  &   * bi-dqreal(iz+1))/dqreal(iz**2)) - 2.q0*e*dqreal(nl22)
                  
                  dej = -2.q0*(dqreal(nl22)*(dqreal(ns1)*aj+1.q0)*(dqreal(ns1)*aj-1.q0) &
                  &   / dqreal(ns1**2)+(dqreal(nl2*iz)*bj+dqreal(iz-1))*(dqreal(nl2*iz)*bj &
                  &   - dqreal(iz+1))/dqreal(iz**2)) - 2.q0*e*dqreal(nl22)
                  
                  ld10 = 2*(l1j + l1i)
                  md20 = -4*16384 + 2*nap
                  d30 = shift
                  le10 = 2*(l2j + l2i)
                  me20 = -4*16384*(iz-1) + 2*mbp
                  qd30 = ai*ap + bi*bp
                  itwo = 0
                  
                  if (namm.gt.0) itwo = 2
                  
                  cjjii = am*ap + bm*bp

                  !? set maximum powers needed and generate basis function integrals
                  
                  mm = 0
                  if (lder) mm = 1
                  mml = 0
                  if (iand(ldd1, 1).eq.0) mml = 1
                  if (namm.gt.0) mml = 1
                  
                  mp1 = 1 - ld1(nb3) - lx1 + ldd1
                  mq1 = 1 - ld2(nb3) - lx2 + ldd1
                  npi = ndp(nb1)
                  nqi = ndq(nb1)
                  
                  if (ide.eq.1) then
                     npi=ndq(nb1)
                     nqi=ndp(nb1)
                  
                  endif

                  if (kdiv.eq.1) then
                     if (npi.lt.lx1) mp1 = mp1 + 1
                     if (nqi.lt.lx2) mq1 = mq1 + 1
                     if (ndp(nb3).lt.ld1(nb3)) mp1 = mp1 + 1
                     if (ndq(nb3).lt.ld2(nb3)) mq1 = mq1 + 1
                     if (npi.lt.lx1.or.ndp(nb3).lt.ld1(nb3)) mp1 = mp1 + 1
                     if (nqi.lt.lx2.or.ndq(nb3).lt.ld2(nb3)) mq1 = mq1 + 1
                  endif

                  msl = 0

                  if (kdiv.eq.2) msl = ldd1 + 1
                  if (mp1.lt.2) mp1 = mp1 + 1
                  if (mq1.lt.2) mq1 = mq1 + 1

                  if (nbl2.ne.1 .or. nh.eq.0) then
                     isum = ijka(nbl1) + ijka(nbl2) + mm

                     mp2 = iia(nbl2)+(1-id)*iia(nbl1) + id*jja(nbl1) + mp1 + mm + mml
                     mq2 = jja(nbl2) + (1-id)*jja(nbl1) + id*iia(nbl1)  + mq1 + mm + mml
                     ms2 = 2*ldd1 + kka(nbl1) + kka(nbl2) + itwo + ms1 - 2

                     if (id+namm.eq.0.and.(nbx2.gt.nbx1.or.nbx4.gt.nbx3)) ms2 = ms2 + 2

                     if (mq2.gt.na-1.or.ms2.gt.nc-1) then 
                        write(6, 20) mq2, ms2
                     20 format(' ARRAY SIZE OF AT EXCEEDED - EXIT CALLED.', 2i4)
                        stop
                     end if

                  else if (nbl1.le.1) then
                     isum = 2*lrgl + mm + 2
                     
                     mp2 = mp1 + 2*lrgl + 1 + mm + mml
                     mq2 = mq1 + 2*lrgl + 1 + mm + mml
                     ms2 = 2*ldd1 + ms1 - 2
                     
                     ap = one/dqreal(ns1) + (1.q0-id)/dqreal(ns1) + id*(iz-1.q0)/dqreal((iz*nl2))
                     y1g = ap
                     
                     bp = (iz-1.q0)/dqreal(iz*nl2) + (1-id)*(iz-1.q0)/dqreal(iz*nl2) + id*1.q0/dqreal(ns1)
                     y2g = bp
                  else
                     isum = ia/2 + lrgl + mm + 1

                     mp2 = (1-id)*pam + id*qam + mp1 + mm + mml
                     mq2 = (1-id)*qam + id*pam + mq1 + lrgl + mm + mml
                     ms2 = 2*ldd1 + maxc + itwo + ms1 - 2

                     ap = one/dqreal(ns1) + dqreal(1-id)*cp1(nb1) + dqreal(id)*cp2(nb1)
                     y1g = ap

                     bp = dqreal((iz-1))/dqreal(iz*nl2) + dqreal(1-id)*cp2(nb1) + dqreal(id)*cp1(nb1)
                     y2g = bp
                  end if

                  !? generate faca(i) and facb(i) in "facab"
                  call facab(ia, y1g, y2g)

                  !! ask about this
                  !? generate basis function integrals for non-linear parameters y1g and ...
                  !? changed from nb3, nb1
                  isum = isum + linc
                  mp2 = mp2 + linc

                  if (nbl1.eq.1.and.nbl2.eq.1) isum = isum + 4

                  call genint(id, nbl2, nbl1, kdiv)

                  !! ask about this
                  !? the following are equivalent for the screened hydrogenic term...?
                  !? xint = 0.5q0 * fac(in + jn + 1 + 2*lrgl) * dqreal(neig + lrgl) ** (in + jn + 1 + 2*lrgl)
                  !? tint = t(1, in + jn - 1 + 2*lrgl, 1)

                  do i=nb1, nb2
                     i1 = nb4

                     if (nb3.eq.nb1) i1 = i

                     if (.not.lder) then
                        if (ldump) then
                           nadd = i1 - nb3 + 1

                           if (irow+nadd.ge.maxbf) then
                              if (idxx.eq.0) nrow(ndmp) = irow

                              if (ndmp.gt.nrmax) then
                                 write(*,*) 'NDMP.GT.NRMAX', ndmp, nrmax
                                 stop
                              endif

                              !! ask about this
                              !   IF IDXX = 1, WRITE BLOCK OF DIRECT + EXCHANGE OVERLAP INTEGRALS TO F
                              !   IF IDXX = 0, WRITE A BLOCK OF DIRECT OVERLAP INTEGRALS TO FILE 8.
                              !   IF IDXX = 1, READ A BLOCK OF DIRECT OVERLAP INTEGRALS FROM FILE 8.
                              !? if idxx = 1, write block of direct + exchange overlap integrals to f...
                              !?              also, read a block of direct overlap integrals from file 8
                              !? if idxx = 0, write a block of direct overlap integrals to file 8

                              if (idxx.eq.1.and.ndmp.gt.1) call dump(ovb, nrow(ndmp-1), 2)
                              call dump(ovb, nrow(ndmp), idxx)
                              
                              irow = 0
                              ndmp = ndmp + 1
                           end if
                        end if
                        
                        ir1 = nb3 - irow - 1
                        irow = irow + nadd

                        if (idxx.ne.0) then
                           do j=nb3, i1
                              psi2(j-nb3+1) = hm(md(i)+j)
                           end do
                        end if
                     end if
                     
                     npi = ndp(i)
                     nqi = ndq(i)
                     nsi = nds(i)
                     
                     if (ide.ne.1) then
                        npi=ndq(i)
                        nqi=ndp(i)
                     end if

                     lsum = -1

                     if (namm.gt.0) then
                        do l=mldd5, ldd5, 2
                           cpl(l, 9) = amm4*(dqreal(npi*nqi)*cpl(l, 5) + dqreal(npi)*cpl(l, 8) &
                           &  + dqreal(nqi)*cpl(l, 7)+ cpl(l, 6))

                           cpl(l, 10) = amm4*(dqreal(npi)*cpl(l, 5) + cpl(l, 7))
                           cpl(l, 11) = amm4*(dqreal(nqi)*cpl(l, 5) + cpl(l, 8))
                        end do

                        do k=9, 11
                           lmin(k) = mldd5
                           lmax(k) = 0
                           do l=mldd5, ldd5, 2
                              if (abs(cpl(l, k)).lt.1.q-10) cpl(l, k) = 0.q0
                              if (cpl(l, k).ne.0.q0) lmax(k) = l
                           end do
                        end do

                        ldd9 = lmax(9)
                        ldd10 = lmax(10)
                        ldd11 = lmax(11)

                        dd9 = cpl(1, 9)
                        dd10 = cpl(1, 10)
                        dd11 = cpl(1, 11)

                        abi = ai*bi*amm4
                        abdd5 = ai*bi*dd5*amm4
                        bdd10 = -bi*dd10
                        add11 = -ai*dd11
                     endif
                     
                     do j=nb3, i1
                        ibf = j - ir1
                        j1 = j-nb3+1

                        npj=ndp(j)
                        nqj=ndq(j)
                        nsj=nds(j)

                        ss = nsi + nsj + ms1

                        npp = npj + npi
                        npm = npj - npi
                        nqp = nqj + nqi
                        nqm = nqj - nqi
                        nsp = nsj + nsi
                        nsm = nsj - nsi
                        nspx = nsp

                        if (nsp.eq.0) then
                           nsm = 0
                           nspx = 1
                        endif
                        
                        nsm2 = nsm + nsm

                        if (namm.gt.0) then
                           cmp = 0.q0
                           ssm = dqreal(nsm)

                           if (ss.gt.ms1) cmp = ssm/dqreal(nsp)

                           cmp2 = 2.q0*cmp
                           cip = 0.q0

                           if (nsi.gt.0) cip = dqreal(nsi)*amm4/dqreal(nsp)

                           cips = dqreal(nspx)*cip
                           qd1 = cip*0.5q0*dqreal(npi*(npp+1)*2 + ldel1)
                           qe1 = cip*0.5q0*dqreal(nqi*(nqp+1)*2 + ldel2)
                           qd2 = - cip*(dqreal(npi)*ap + ai*dqreal(npp+2))
                           qe2 = - cip*(dqreal(nqi)*bp + bi*dqreal(nqp+2))
                           qd3 = cip*qd30
                           qd4 = - amm4*dqreal(nsi)*dqreal(nsi+1)
                           t40 = am*bm*dqreal(nsp)
                           t41 = am*bp + ap*bm
                        endif
                        
                        lgo = .false.

                        lxx = npj + nqj
                     
                        if (lsum.eq.lxx.and.kdiv.eq.0) lgo = .true.
                        lsum = lxx

                        if (nsi.ge.0.and.nsj.ge.0.) then

                           d(1) = 0.5q0*dqreal(nspx*(-(npm*npm+npp-ld10) - (npp*npp+npp)) &
                           &    - nsm2*(l1j-l1i-npm*(npp+1)))
                           
                           d(2) = dqreal((nspx*(md20 + mam*npm + nap*npp) &
                           &    - nsm*(npm*nap + mam*(npp+2))))/16384.q0
                           
                           d3r = dqreal(nsm)*cjjii
                           d(3) = dqreal(nspx)*d30 + d3r*dqreal(nl22)
                           d(6) = dqreal(4*nsi*nsj*nspx)

                           d(4) = 0.5q0*dqreal(nspx*(-(nqm*nqm+nqp-le10) - (nqp*nqp+nqp)) &
                           &    - nsm2*(l2j-l2i-nqm*(nqp+1)))
                           
                           d(5) = ((nspx*(me20 + mbm*nqm + mbp*nqp) &
                           &    - nsm*(nqm*mbp + mbm*(nqp+2)))/16384.q0)/iz
                           
                           ee2 = e*dqreal(2*nspx*nl22)

                        else if (nsj.ge.0) then

                           nsjx = nsj
                           if (nsp.eq.0) nsjx = 0
                           
                           d(1) = -2.q0*dqreal(nspx*(npj*(npj+1) - l1j) &
                           &    - nsjx*(2*npj*(npp+1) + ldel3))
                           
                           d(2) = 4.q0*(dqreal(nspx)*(aj*dqreal(npj+1) - dqreal(1.q0)) &
                           &    - dqreal(nsjx)*(aj*dqreal(npp+2) + dqreal(npj)*ap))
                           
                           d(3) = dej*dqreal(nspx) + 4.q0*dqreal(nsjx*nl22)*(ap*aj + bp*bj)
                           d(6) = dqreal(-4*nsj*(nsj+1)*nspx)
                           
                           d(4) = dqreal(-2*(nspx*(nqj*(nqj+1) - l2j) &
                           &    - nsjx*(2*nqj*(nqp+1) + ldel4)))
                           
                           d(5) = 4.q0*(dqreal(nspx)*(dqreal(iz)*bj*dqreal((nqj+1)-iz+1)) &
                           &    - dqreal(nsjx*iz)*(bj*dqreal(nqp+2) + dqreal(nqj)*bp))/dqreal(iz)
                           
                           ee2 = 0.q0
                        
                        endif

                        if (.not. lder) then

                           !? coding for hamiltonian matrix
                           
                           mdij = md(i) + j

                           if (ldump) then
                              mdij = ibf
                              if (idxx.eq.0) ovb(ibf) = 0.q0
                           endif

                           if (idxx.eq.0) psi2(j1) = 0.q0

                           pp = npi + npj + mp1
                           qq = nqi + nqj + mq1
                           sumq = 0.q0
                           
                           if (j.gt.nh) then
                              if (d(1).ne.0) fdd(1) = fpl(pp-2, qq, 1)

                              if (.not.lgo) then
                                 fdd(2) = fpl(pp-1, qq, 1)
                              else
                                 fdd(2) = fdd(5)
                              end if

                              if (d(4) .ne. 0) fdd(4) = fpl(pp, qq-2, 1)

                              fdd(5) = fpl(pp, qq-1, 1)
                              fdd(3) = fpl(pp, qq, 1)
                              ss = ss - 1

                              if (.not.ltxx) zr12 = 4.q0*(fpl(pp, qq, 1) - fdd(5))/z

                              if (ltxx) zr12 = (txx(npp+1, nqp+1, nsp, 2) &
                              &              - txx(npp+1, nqp, nsp+1, 2))*(2.q0/z)

                              ss = ss + 1

                              if (d(6).ne.0) then
                                 ss = ss - 2
                                 fdd(6) = fpl(pp, qq, 1)
                                 ss = ss + 2
                              end if

                              sum = zr12+((d(1)*fdd(1)+d(2)*fdd(2)+(d(3)-ee2)*fdd(3)/dqreal(nl22)) &
                              &   + d(6)*fdd(6)+(d(4)*fdd(4)+d(5)*fdd(5)))/dqreal(nspx)

                              !? include logarithmic contributions

                              if (nsp.eq.0.and.nsj.ne.nsi) then
                                 lg = 0

                                 if (nsj.gt.0) then
                                    d(1) = 2*((2*npj*(npp+1) + ldel3))
                                    d(2) = -4.q0*(aj*dqreal(npp+2) + dqreal(npj)*ap)
                                    d(4) = 2*((2*nqj*(nqp+1) + ldel4))
                                    d(5) = -4.q0*(bj*dqreal(nqp+2) + dqreal(nqj)*bp)
                                    d(3) = 4.q0*(ap*aj + bp*bj)
                                 else
                                    d(1) = 2*((2*npi*(npp+1) + ldel1))
                                    d(2) = -4.q0*(ai*dqreal(npp+2) + dqreal(npi)*ap)
                                    d(4) = 2*((2*nqi*(nqp+1) + ldel2))
                                    d(5) = -4.q0*(bi*dqreal(nqp+2) + dqreal(nqi)*bp)
                                    d(3) = 4.q0*(ap*ai + bp*bi)
                                 endif

                                 if (d(1).ne.0) fddl1 = fpl(pp-2, qq, 1)
                                 fddl2 = fpl(pp-1, qq, 1)

                                 if (d(4).ne.0) fddl4 = fpl(pp, qq-2, 1)
                                 fddl5 = fpl(pp, qq-1, 1)

                                 fddl3 = fpl(pp, qq, 1)

                                 sum = sum + d(1)*fddl1+d(2)*fddl2+d(3)*fddl3+d(4)*fddl4 &
                                 &   + d(5)*fddl5
                                 lg = 1
                              endif

                              !? 2.* hamiltonian matrix element is accumulated in hm and the overlap integral in
                              !? ovb as linearized triangular arrays

                              fdd(1) = fpl(pp-2, qq, 1)
                              fdd(4) = fpl(pp, qq-2, 1)
                              ss = ss - 2

                              fdd(6) = fpl(pp, qq, 1)
                              ss = ss + 2

                              dr1 = -2.q0*(nspx*(npi*(npi+1) - l1i) &
                              &   - nsi*(2*npi*(npp+1) + ldel1))/nspx
                              
                              dr2 = 4.q0*(ai*dqreal(npi+1) - dqreal(1.q0) &
                              &   - dqreal(nsi)*(ai*dqreal(npp+2) &
                              &   + dqreal(npi)*ap)/dqreal(nspx))
                              
                              dr4 = -2.q0*(nspx*(nqi*(nqi+1) - l2i) &
                              &   - nsi*(2*nqi*(nqp+1) + ldel2))/nspx
                              
                              dr5 = 4.q0*(bi*dqreal(nqi+1) - (z-dqreal(1.q0))/z &
                              &   - dqreal(nsi)*(bi*dqreal(nqp+2) &
                              &   + dqreal(nqi)*bp)/dqreal(nspx))
                              
                              dr3 = dei/dqreal(nl22) + 4.q0*dqreal(nsi)*(ap*ai + bp*bi)/dqreal(nspx)
                              dr6 = - dqreal(4*nsi*(nsi+1))
                              
                              dl1 = -2.q0*(nspx*(npj*(npj+1) - l1j) &
                              &   - nsj*(2*npj*(npp+1) + ldel3))/nspx
                              
                              dl2 = 4.q0*(aj*dqreal(npj+1) - dqreal(1.q0) &
                              &   - dqreal(nsj)*(aj*dqreal(npp+2) + dqreal(npj)*ap)/dqreal(nspx))
                              
                              dl4 = -2.q0*(nspx*(nqj*(nqj+1) - l2j) &
                              &   - nsj*(2*nqj*(nqp+1) + ldel4))/nspx
                              
                              dl5 = 4.q0*(bj*dqreal(nqj+1) &
                              &   - (z-dqreal(1.q0))/z - dqreal(nsj)*(bj*dqreal(nqp+2) &
                              &   + dqreal(nqj)*bp)/dqreal(nspx))
                              
                              dl3 = dej/dqreal(nl22) + 4.q0*dqreal(nsj)*(ap*aj + bp*bj)/dqreal(nspx)
                              dl6 = - dqreal(4*nsj*(nsj+1))
                              
                              sumr  = dr1*fdd(1)+dr2*fdd(2)+dr3*fdd(3)+dr4*fdd(4)+dr5*fdd(5)    &
                              &   +dr6*fdd(6) + zr12
                              
                              sumll = dl1*fdd(1)+dl2*fdd(2)+dl3*fdd(3)+dl4*fdd(4)+dl5*fdd(5)    &
                              &   +dl6*fdd(6) + zr12

                              !? include logarithmic contributions

                              if (nsp.eq.0.and.nsj.ne.nsi) then
                                 lg = 0
                                 nsm = nsj - nsi
                                 nsm2 = 2*nsm
                                 
                                 dr1 = 2.q0*(nsi*(2*npi*(npp+1) + ldel1))/nspx
                                 dr2 = -4.q0*(dqreal(nsi)*(ai*dqreal(npp+2) + dqreal(npi)*ap)/dqreal(nspx))
                                 dr4 = dqreal(2*(nsi*(2*nqi*(nqp+1) + ldel2))/nspx)
                                 dr5 = -4.q0*(dqreal(nsi)*(bi*dqreal(nqp+2) + dqreal(nqi)*bp)/dqreal(nspx))
                                 dr3 = 4.q0*dqreal(nsi)*(ap*ai + bp*bi)/dqreal(nspx)
                                 
                                 dl1 = 2.q0*(nsj*(2*npj*(npp+1) + ldel3))/nspx
                                 dl2 = -4.q0*(dqreal(nsj)*(aj*dqreal(npp+2) + dqreal(npj)*ap)/dqreal(nspx))
                                 dl4 = 2.q0*(nsj*(2*nqj*(nqp+1) + ldel4))/nspx
                                 dl5 = -4.q0*(dqreal(nsj)*(bj*dqreal(nqp+2) + dqreal(nqj)*bp)/dqreal(nspx))
                                 dl3 = 4.q0*dqreal(nsj)*(ap*aj + bp*bj)/dqreal(nspx)
                                 
                                 fddl1 = fpl(pp-2, qq, 1)
                                 fddl2 = fpl(pp-1, qq, 1)
                                 fddl4 = fpl(pp, qq-2, 1)
                                 fddl5 = fpl(pp, qq-1, 1)
                                 fddl3 = fpl(pp, qq, 1)
                                 
                                 fsum = dl1*fddl1 + dl2*fddl2 + dl3*fddl3/dqreal(nl22) + dl4*fddl4 &
                                 &    + dl5*fddl5

                                 sumr  = sumr + dr1*fddl1+dr2*fddl2+dr3*fddl3+dr4*fddl4+dr5*fddl5
                                 sumll = sumll + dl1*fddl1+dl2*fddl2+dl3*fddl3+dl4*fddl4+dl5*fddl5

                                 lg = 1
                              endif

                              jjj = pp + qq - mp1 - mq1 + 3
                              delt = fac(jjj)/(y1g+y2g)**jjj

                              deltr = (sumr-sum)/sum
                              deltl = (sumll-sum)/sum
                              deltrl = ((sumr-sumll)/delt)**2

                              psi2(j1) = psi2(j1) + sx*sum
                              yx = 2.q0*fdd(3)
                              ovb(mdij) = ovb(mdij) + sx*yx

                              deltl16=deltl
                              deltrl16=deltrl

                              !! THE NEXT 100 OR SO LINES OF CODE ARE UNTESTED AS OF WRITING. 
                              !! IF THERE IS A PROBLEM WITH THE PROGRAM AND EVERYTHING ELSE SEEMS 
                              !! TO BE IN ORDER, THIS IS THE LIKELY CAUSE. CONTACT ERIC FOR DETAILS.

                              if (namm.ne.0 .and. ss.eq.ms1 .and. mldd5.eq.1 .and. ldd1.ne.1) then
                                 !? mass polarization contribution for uncorrelated terms
                                 sumq = abdd5*at(jat(pp, qq)+ss) &
                                 &    + dd9*at(jat(pp-1, qq-1)+ss) &
                                 &    + bdd10*at(jat(pp-1, qq)+ss) &
                                 &    + add11*at(jat(pp, qq-1)+ss)
                              else if (namm.ne.0 .and. mldd5.eq.1) then
                                 if (ldd1.ne.1) then
                                    !? general form of mass polarization contribution

                                    sumq = abi*fpl(pp, qq, 5) + fpl(pp-1, qq-1, 9) - bi*fpl(pp-1, qq, 10)   &
                                    &    - ai*fpl(pp, qq-1, 11)
                                    
                                    if (ldd4.gt.0)sumq = sumq - cip*(fpl(pp-2, qq, 4)+fpl(pp, qq-2, 4))
                                    
                                    sumq = sumq + qd1*fpl(pp-2, qq, 1) + qd2*fpl(pp-1, qq, 1)             &
                                    &    + qe1*fpl(pp, qq-2, 1) + qe2*fpl(pp, qq-1, 1) + qd3*fpl(pp, qq, 1)
                                    
                                    if (qd4.ne.0.q0) then
                                       ss = ss - 2
                                       sumq = sumq + qd4*fpl(pp, qq, 1)
                                       ss = ss + 2
                                    end if
                                 else
                                    !? special form of mass polarization contribution for the case ldd1 = 1
                                    
                                    t1 = dqreal(nsp*npm*nqm - nsm*(npp*nqm + nqp*npm))
                                    t2 = am*dqreal(nqm*nsp) - dqreal(nsm)*(ap*dqreal(nqm) + am*dqreal(nqp))
                                    t3 = dqreal(nsp*npm)*bm - dqreal(nsm)*(dqreal(npp)*bm + dqreal(npm)*bp)
                                    t4 = t40 - dqreal(nsm)*t41
                                    
                                    sumq = amm*(dd52*(t1*at(jcs(pp-1, qq-1)+ss)  &
                                    &    - t2*at(jcs(pp, qq-1)+ss) &
                                    &    - t3*at(jcs(pp-1, qq)+ss) &
                                    &    + t4*at(jcs(pp, qq)+ss)) &
                                    &    - dqreal(nsp*(nsm*nsm+nsp))*(0.5q0*at(jat(pp, qq)+ss-2) &
                                    &    + faca(npp2)*facb(nqsp)))/dqreal(nspx)
                     
                                    !? the last continuation line 5 corrects for:
                                    !? a) the calculation of integral of 1/r12 - 1/r2 for ss-2 = 1 and
                                    !? b) the calculation of residual parts of integrals for ss-2 > 1
                                 end if
                              end if

                              psi2(j1) = psi2(j1) - sx*sumq
                              
                           end if

                           if (namm.ne.0) then
                              yx = fpl(pp, qq-1, 1)
                              ss = ss - 1
                              
                              if (i.le.1 .and. j.le.1) then
                                 sum = fpl(pp, qq, 1) - yx
                                 if (namm.gt.0) then
                                    bb = (z-1.q0)/(dqreal(nl2)*z)
                                    sumq = -amm4*dd2*128.q0*bb**5*(1.q0-bb) &
                                    &    **(2*nl2-4)*dqreal(nl2*nl2-1)/(3.q0*(1.q0+bb)**(2*nl2+4))
                                 endif
                              else
                                 if (namm.gt.0) then
                                    ss = ss + 1
                                    sumq = amm4*(aj*bj*fpl(pp, qq, 5) - aj*fpld(pp, qq-1, 5))
                                    
                                    do l=mldd5, ldd5, 2
                                       cpl(l, 9) = amm4*(dqreal(npj*nqj)*cpl(l, 5) &
                                       &         + dqreal(npj)*cpl(l, 8) &
                                       &         + dqreal(nqj)*cpl(l, 7)+ cpl(l, 6))

                                       cpl(l, 10) = amm4*(dqreal(npj)*cpl(l, 5) + cpl(l, 7))
                                       cpl(l, 11) = amm4*(dqreal(nqj)*cpl(l, 5) + cpl(l, 8))
                                    end do

                                    do k=9, 11
                                       lmin(k) = mldd5
                                       lmax(k) = 0
                                       do l=mldd5, ldd5, 2
                                          if (abs(cpl(l, k)).lt.1.q-10) cpl(l, k) = 0.q0
                                          if (cpl(l, k).ne.0.q0) lmax(k) = l
                                       end do
                                    end do
                                    
                                    if (ldd9.gt.0) sumq = sumq + fpl(pp-1, qq-1, 9)
                                    if (ldd10.gt.0) sumq = sumq - bj*fpl(pp-1, qq, 10)
                                    if (ldd11.gt.0) sumq = sumq - aj*fpl(pp, qq-1, 11)
                              
                                    sumq16=sumq

                                    if (abs(sumq).gt.1.q99) write(6, 33) i, j, pp, qq, ss, sumq16
                                    ss = ss - 1
                                 endif

                                 sum = fpl(pp, qq, 1) - yx

                              end if


                              ss = ss + 1
                              ss16=ss
                              sumq16=sumq
                              yx = 2.q0*fpl(pp, qq, 1)
                              
                              psi2(j1) = psi2(j1) + sx*(4.q0*sum/z - e*yx - sumq)
                              ovb(mdij) = ovb(mdij) + sx*yx
                           end if

                           if (i.eq.j) psi(i) = ovb(mdij)

                           cycle
                        end if

                        !? coding for derivatives
                        nspm = 2*(nspx - nsm)
                        nspp = 2*(nspx + nsm)

                        d2d(1, 1) = nspm*(npj+1)
                        d2d(2, 1) = nspm*(nqj+1)
                        d2d(1, 3) = -(nspm*maj)/16384.q0
                        d2d(2, 3) = -((nspm*mbj)/16384.q0)/iz
                        d2d(ide , 2) = nspp*(npi+1)
                        d2d(idex, 2) = nspp*(nqi+1)
                        d2d(ide , 4) = -(nspp*mai)/16384.q0
                        d2d(idex, 4) = -((nspp*mbi)/16384.q0)/iz

                        two = 1.q0
                        if (i.ne.j) two = 2.q0

                        ps2 = psi(i)*psi(j)*sx*two/dqreal(nspx)

                        if (lr12) ps2x = ps2

                        do ipow=1, 2
                           nqz = 2 - ipow
                           npz = ipow - 1
                           pp = npi + npj + mp1 + npz
                           qq = nqi + nqj + mq1 + nqz
                           ipp = pp - mp1

                           jj(8) = pp - 2
                           ipp1 = ipp + 1
                           ipp2 = ipp + 2
                           iqp = qq - mq1

                           jj(9) = qq - 2
                           nqsm = nqm + nsm
                           iqsp = iqp + nsp
                           iqsp1 = iqsp + 1
                           iqsp2 = iqsp + 2
                           sum = 0.q0

                           if (j.gt.nh) then
                              if (lr12) then !! maybe flip this

                                 !? the following evaluates contributions to the derivatives with dominant
                                 !? leading terms subtracted. only "direct" terms are considered.

                                 sumr = 0.q0

                                 if (nsp.le.2) d(6) = 0.q0
                                 if (.not.lblk2) then
                                    tt(5) = (-ap2x*dqreal(ipp2-npm*npm+l1j4) &
                                    &     + dqreal(ipp1)*(amx*(-apx*dqreal(2*npm)+amx*dqreal(ipp2))+ca1-ca2*dqreal(ipp2)))
                                    
                                    tt(6) = (-bp2x*dqreal(iqsp2-nqsm*nqsm+l2j4 &
                                    &     + 4*nsi*nsj) + dqreal(iqsp1)*(bmx*(-bpx*dqreal(2*nqsm)+bmx*dqreal(iqsp2)) + cb1 &
                                    &     - cb2*dqreal(iqsp2)))

                                    tt(3) = dqreal(ipp2*ipp1)
                                    tt(4) = cb3*dqreal(iqsp1*iqsp2)

                                    fdd30 = faca(ipp2)*facb(iqsp2)*dqreal(nspx)/dqreal(16384.q0*16384.q0)

                                    if (ipow.eq.2) then
                                       tt(1) = 2.q0*apx*(apx*dqreal(npm) - amx*dqreal(ipp1))
                                       tt(2) = tt(5) + ap2x
                                    else
                                       tt(1) = 2.q0*bpx*(bpx*dqreal(nqsm) - bmx*dqreal(iqsp1))
                                       tt(2) = tt(6) + bp2x
                                    end if
                                 end if

                                 if (ss.eq.2) then
                                    if (ipow.ne.2) then
                                       df1x(1, 1) = 0.q0
                                       df1x(2, 1) = 0.q0
                                       df1x(1, 2) = 0.q0
                                       df1x(2, 2) = 0.q0
                                    end if
                                 else
                                    if (lgo .or. ipow.ne.1) then

                                       nsp1 = nsp + 1
   
                                       do if=1, 3
                                          fddr(if, 1) = (txx(ipp+if-2, iqp+1, nsp1, 2) &
                                          &           + at0(iqsp2, ipp+if-1))
                                       end do
   
                                       fddr(1, 2) = 0.q0
   
                                       do if=1, 2
                                          fddr(if, 2) = txx(ipp1, iqp+if-2, nsp1, 2) &
                                          &           + at0(iqsp+if-1, ipp+2)
                                       end do
   
                                       fddr(3, 2) = 0.q0
   
                                       if (ss.gt.4) then
                                          fddr(3, 2) = (txx(ipp1, iqp+1, nsp-1, 2) + at0(iqsp, ipp2))
                                       end if
                                    end if
   
                                    if (ipow.ne.1) then
                                       df1x(1, 1) =   df1x(1, 1)   + d2d(1, 1)  *fddr(1, 1)
                                       df1x(ide, 2) = df1x(ide, 2) + d2d(ide, 2)*fddr(1, 1)
                                    else
                                       df1x(1, 1) =                              d2d(1, 3)   *fddr(2, 2)
                                       df1x(2, 1) =    d2d(2, 1)   *fddr(1, 2) + d2d(2, 3)   *fddr(2, 2)
                                       df1x(ide, 2) =                            d2d(ide, 4) *fddr(2, 2)
                                       df1x(idex, 2) = d2d(idex, 2)*fddr(1, 2) + d2d(idex, 4)*fddr(2, 2)
                                    end if
   
                                    sumr = d(1)*fddr(1, 1)+d(2)*fddr(2, 1)+(d(3)-ee2)*fddr(3, 1)/dqreal(nl22) &
                                    &    + d(4)*fddr(1, 2)+d(5)*fddr(2, 2)+d(6)*fddr(3, 2)
                                 end if
            
                                 zr12 = (at(jat(pp, qq)+ss-1) &
                                 &    - at(jat(pp, qq-1)+ss))*(4.q0/z)

                                 hmxx(ipow) = sumr + zr12*dqreal(nspx)

                                 if (lblk2) then
                                    if (ipow.ne.2) then
                                       if (ipow.ne.2) then
                                          df1x(2, 1) = df1x(2, 1) &
                                          &          + fdd30*(tt(5)/tt(3) + (tt(2) + tt(1))/tt(4)   + e10)

                                          df1x(idex, 2) = df1x(idex, 2) + fdd30*(tt(5)/tt(3) &
                                          &             + (tt(2) - tt(1))/tt(4) + e10)


                                       else 
                                          df1x(1, 1) = df1x(1, 1) + fdd30*((tt(2) + tt(1))/tt(3) + tt(6)/tt(4) &
                                          &          + e10)

                                          df1x(ide, 2) = df1x(ide, 2) + fdd30*((tt(2) - tt(1))/tt(3) &
                                          &            + tt(6)/tt(4) + e10)

                                       end if
                                    end if
                                 else
                                    dhmx(nhi, 2) = dhmx(nhi, 2) + 0.5q0*(-hmxx(ide) + df1x(2, 2))*ps2
                                    dhmx(nhj, 2) = dhmx(nhj, 2) + 0.5q0*(-hmxx(1) + df1x(2, 1))*ps2
                                    dhmx(nhi, 1) = dhmx(nhi, 1) + 0.5q0*(-hmxx(idex) + df1x(1, 2))*ps2
                                    dhmx(nhj, 1) = dhmx(nhj, 1) + 0.5q0*(-hmxx(2) + df1x(1, 1))*ps2
                                 end if

                                 if (namm.eq.0) then
                                    if (ipow.eq.2) exit
                                    cycle
                                 end if
                           !
                                 if (ipow.ne.1) then
                                    fqq2 = amm*dd52*at(jcs(pp-1, qq)+ss)
                                    fqq3 = amm*dd52*at(jcs(pp-2, qq)+ss)
                                    fqq4 = amm*dd52*at(jcs(pp-1, qq-1)+ss)
                                    df1(1, 1) = (dqreal(nspx*nqm)-dqreal(nsm)*dqreal(nqp+nqm))*fqq4 &
                                    &         - (dqreal(nspx)*bm-dqreal(nsm)*(bp+bm))*fqq2
                                    df1(2, 1) = (dqreal(nspx*npm)-dqreal(nsm)*dqreal(npp+npm))*fqq3 &
                                    &         - (dqreal(nspx)*am-dqreal(nsm)*(ap+am))*fqq2
                                    if (id.ne.1) then
                                       df1(1, 2) = dqreal(-nspx*nqm+nsm*(nqp-nqm))*fqq4 &
                                       &         - (dqreal(-nspx)*bm+dqreal(nsm)*(bp-bm))*fqq2
                                       df1(2, 2) = dqreal(-nspx*npm+nsm*(npp-npm))*fqq3 &
                                       &         - (dqreal(-nspx)*am+dqreal(nsm)*(ap-am))*fqq2
                                    else
                                       df1(2, 2) = dqreal(-nspx*nqm+nsm*(nqp-nqm))*fqq4 &
                                       &         - (dqreal(-nspx)*bm+dqreal(nsm)*(bp-bm))*fqq2
                                       df1(1, 2) = dqreal(-nspx*npm+nsm*(npp-npm))*fqq3 &
                                       &         - (dqreal(-nspx)*am+dqreal(nsm)*(ap-am))*fqq2
                                    end if
                                 end if

                                 t1 = dqreal(nspx*npm*nqm - nsm*(npp*nqm + nqp*npm))
                                 t2 = am*dqreal(nspx*nqm) - dqreal(nsm)*(ap*dqreal(nqm) + am*dqreal(nqp))
                                 t3 = dqreal(nspx*npm)*bm - dqreal(nsm)*(dqreal(npp)*bm + dqreal(npm)*bp)
                                 t4 = dqreal(nspx)*am*bm - dqreal(nsm)*(am*bp + ap*bm)

                                 sumq = amm*(dd52*(t1*at(jcs(pp-1, qq-1)+ss) &
                                 &    - t2*at(jcs(pp, qq-1)+ss) &
                                 &    - t3*at(jcs(pp-1, qq)+ss) &
                                 &    + t4*at(jcs(pp, qq)+ss)) &
                                 &    - dqreal((nsm*nsm+nsp)*nspx)*(0.5q0*at(jat(pp, qq)+ss-2) &
                                 &    + faca(ipp2)*facb(iqsp)))

                                 !? the last continuation line 5 corrects for
                                 !? a) the calculation of integgral of 1/r12 - 1/42 for ss-2 = 1 and
                                 !? b) the calculation of residual parts of integrals for ss-2 > 1

                                 hmx(ipow) = - sumq
                                 cycle
                              end if

                              if (.not.lgo .or. ipow.ne.1) then
                                 fdd(3) = fpl(pp, qq, 1)

                                 if (pp.gt.2) fdd(1) = fpl(pp-2, qq, 1)
                                 if (qq.gt.2) fdd(4) = fpl(pp, qq-2, 1)

                                 if (ss-2 .ne. 0)  then
                                    ss = ss - 2
                                    fdd(6) = fpl(pp, qq, 1)
                                    ss = ss + 2
                                 end if

                                 if (ipow.ne.1) then
                                    fdd(2) = fdd(5)
                                    fdd(5) = fpl(pp, qq-1, 1)
                                    df1(1, 1) = df1(1, 1) + d2d(1, 1)*fdd(1)

                                    if (id.eq.1) then
                                       df1(2, 2) = df1(2, 2) + d2d(2, 2)*fdd(1)
                                    else
                                       df1(1, 2) = df1(1, 2) + d2d(1, 2)*fdd(1)
                                    end if
                                 end if
                              end if

            
                              if (ipow.eq.1) then
                                 fdd(2) = fpl(pp-1, qq, 1)
                                 fdd(5) = fpl(pp, qq-1, 1)

                                 df1(1, 1) = d2d(1, 3)*fdd(5)
                                 df1(2, 1) = d2d(2, 1)*fdd(4) + d2d(2, 3)*fdd(5)

                                 if (id.ne.1) then
                                    df1(1, 2) =              d2d(1, 4)*fdd(5)
                                    df1(2, 2) = d2d(2, 2)*fdd(4) + d2d(2, 4)*fdd(5)
                                 else 
                                    df1(1, 2) = d2d(1, 2)*fdd(4) + d2d(1, 4)*fdd(5)
                                    df1(2, 2) = d2d(2, 4)*fdd(5)
                                 end if
                              end if

                              if (ipow.ne.1 .or. .not.lgo) then
                                 ss = ss - 1

                                 if (.not.ltxx) zr12 = 4.q0*(fpl(pp, qq, 1) - fdd(5))/z
                                 if (ltxx) zr12 = (txx(ipp1, iqp+1, nsp, 2) &
                                 &              - txx(ipp1, iqp, nsp+1, 2))*(dqreal(2.q0)/z)

                                 !? the above comes from 2*(1/r12 - 1/r2)

                                 ss = ss + 1
                              end if

                              sum = zr12*dqreal(nspx) + (d(1)*fdd(1) + d(2)*fdd(2) &
                              &   + (d(4)*fdd(4)+d(5)*fdd(5)) + (d(3)-ee2)*fdd(3)/dqreal(nl22) & 
                              &   + d(6)*fdd(6))

                              !? the above ariase form terms independent of correlation.
                              
                              !? include logarithmic contributions
                              if (nsp.eq.0.and.nsj.ne.nsi) then
                                 lg = 0
                                 nsm = nsj - nsi

                                 dl1 = -0.5q0*dqreal(nsm2*(l1j-l1i-npm*(npp+1)))
                                 dl2 = dqreal(-(nsm*(npm*map + mam*(npp+2))))/dqreal(16384.q0)
                                 dl3r = dqreal(nsm)*cjjii
                                 dl3 = dl3r*dqreal(nl22)
                                 dl4 = -0.5q0*(nsm2*(l2j-l2i-nqm*(nqp+1)))
                                 dl5 = -dqreal((nsm*(nqm*mbp + mbm*(nqp+2)))/16384.q0)/dqreal(iz)

                                 if (dl1.ne.0q0) fddl1 = fpl(pp-2, qq, 1)

                                 if (lgo) then
                                    fddl2 = fddl5
                                 else
                                    fddl2 = fpl(pp-1, qq, 1)
                                 endif

                                 if (dl4.ne.0q0) fddl4 = fpl(pp, qq-2, 1)

                                 fddl5 = fpl(pp, qq-1, 1)
                                 fddl3 = fpl(pp, qq, 1)

                                 sum = sum + dl1*fddl1 + dl2*fddl2 + dl3*fddl3/dqreal(nl22) &
                                 &   + dl4*fddl4 + dl5*fddl5
                                 
                                 lg = 1
                              endif

                              !? general form of mass polarization contributions

                              if (namm.ne.0) then 
                                 if (ipow.ne.1) then
                                    fqq2 = amm4*dqreal(nspx)*fpl(pp-1, qq, 5)
                                    fqq3 = 0.q0
                                    fqq4 = 0.q0

                                    if (ldd10.gt.0) fqq3 = dqreal(nspx)*fpl(pp-2, qq, 10)
                                    if (ldd11.gt.0) fqq4 = dqreal(nspx)*fpl(pp-1, qq-1, 11)

                                    df1(1, 1) = df1(1, 1) + cips*dqreal(npi)*fdd(1)

                                    if (id.ne.1) then
                                       df1(1, 2)= df1(1, 2) + fqq4 - bi*fqq2 + cips*dqreal(npp+npi+2)*fdd(1)
                                       df1(2, 2) = df1(2, 2) + fqq3 - ai*fqq2
                                    else
                                       df1(2, 2)= df1(2, 2) + fqq4 - bi*fqq2 + cips*dqreal(npp+npi+2)*fdd(1)
                                       df1(1, 2) = df1(1, 2) + fqq3 - ai*fqq2
                                    end if
                                 else
                                    df1(1, 1) = df1(1, 1) - ai*cips*fdd(5)
                                    df1(2, 1) = df1(2, 1) - cips*(bi*fdd(5) - dqreal(nqi)*fdd(4))

                                    if (id.ne.1) then
                                       df1(1, 2) = df1(1, 2) - (ai+ap)*cips*fdd(5)
                                       df1(2, 2) = df1(2, 2) - cips*((bi+bp)*fdd(5) - dqreal(nqp+nqi+2)*fdd(4))
                                    else
                                       df1(2, 2) = df1(2, 2) - (ai+ap)*cips*fdd(5)
                                       df1(1, 2) = df1(1, 2) - cips*((bi+bp)*fdd(5) - dqreal(nqp+nqi+2)*fdd(4))
                                    end if
                                 end if

                                 sumq = abi*fpl(pp, qq, 5) + fpl(pp-1, qq-1, 9) - bi*fpl(pp-1, qq, 10) &
                                 &    - ai*fpl(pp, qq-1, 11)
                                 
                                 if (ldd4.gt.0)sumq = sumq - cip*(fpl(pp-2, qq, 4)+fpl(pp, qq-2, 4))
                                 
                                 sumq = sumq + qd1*fdd(1) + qd2*fdd(2) + qe1*fdd(4) & 
                                 &    + qe2*fdd(5) + qd3*fdd(3) + qd4*fdd(6)
                              end if

                              hmx(ipow) = sum - sumq*dqreal(nspx)

                              cycle
                              
                              yx = fpl(pp, qq-1, 1)
                           end if

                           ss = ss - 1
                           if (i.le.1 .and. j.le.1) then
                              sum = fpl(pp, qq, 2)
                              bb = (z-dqreal(1.q0))/(dqreal(nl2)*z)

                              sumq = -amm4*dd2*128.q0*bb**5*(1.q0-bb)**(2*nl2-4)*dqreal(nl2*nl2-1) &
                              &    / (3.q0*(1.q0+bb)**(2*nl2+4))
                           else
                              ss = ss + 1
                              sumq = amm4*(aj*bj*fpl(pp, qq, 5) - aj*fpld(pp, qq-1, 5))
                        
                              do l=mldd5, ldd5, 2
                                 cpl(l, 9) = amm4*(dqreal(npj*nqj)*cpl(l, 5) + dqreal(npj)*cpl(l, 8) & 
                                 &         + dqreal(nqj)*cpl(l, 7) + cpl(l, 6))

                                 cpl(l, 10) = amm4*(dqreal(npj)*cpl(l, 5) + cpl(l, 7))
                                 cpl(l, 11) = amm4*(dqreal(nqj)*cpl(l, 5) + cpl(l, 8))
                              end do

                              do k=9, 11
                                 lmin(k) = mldd5
                                 lmax(k) = 0
                                 do l=mldd5, ldd5, 2
                                    if (abs(cpl(l, k)).lt.1.q-10) cpl(l, k) = 0.q0
                                    if (cpl(l, k).ne.0.q0) lmax(k) = l
                                 end do
                              end do
                        
                              if (ldd9.gt.0) sumq = sumq + fpl(pp-1, qq-1, 9)
                              if (ldd10.gt.0) sumq = sumq - bj*fpl(pp-1, qq, 10)
                              if (ldd11.gt.0) sumq = sumq - aj*fpl(pp, qq-1, 11)

                              if (abs(sumq).gt.1.q90) write(6, 33) i, j, pp, qq, ss, sumq

                              ss = ss - 1

                              if (ldd1.eq.1) then
                                 ipp = pp - mp1
                                 iqq = qq - mq1
                                 sum = 0.q0

                                 do n=1, neig
                                    sum = sum + 2.q0*dd1*ch(n)*((txx(ipp+1, iqq+n, ss-1, 2) &
                                    &   - txx(ipp+1, iqq+n-1, ss, 2))*2.q0/z - e*t(ipp+1, iqq+n, ss))
                                 end do

                                 hmx(ipow) = (sum - sumq)*dqreal(nspx)
                                 ss = ss + 1
                                 go to 11
                              end if
                              sum = fpl(pp, qq, 1) - yx
                           end if
                     

                           ss = ss + 1
                           yx = 2.q0*fpl(pp, qq, 1)
                           hmx(ipow) = (4.q0*sum/z - e*yx - sumq)*dqreal(nspx)

                        11 continue

                           do if=1, 2
                              do jf=1, 2
                                 df1(if, jf) = 0.q0
                              end do
                           end do
                           
                        end do

                        !? accumulation of derivatives
                        hmi = - hmx(ide) + df1(2, 2)
                        hmj = - hmx(1) + df1(2, 1)

                        dhm(nhi, 2) = dhm(nhi, 2) + hmi*ps2
                        dhm(nhj, 2) = dhm(nhj, 2) + hmj*ps2

                        hmi = - hmx(idex) + df1(1, 2)
                        hmj = - hmx(2) + df1(1, 1)

                        dhm(nhi, 1) = dhm(nhi, 1) + hmi*ps2
                        dhm(nhj, 1) = dhm(nhj, 1) + hmj*ps2
                     end do

                     if (lder) cycle

                     do j=nb3, i1
                        hm(md(i)+j) = psi2(j-nb3+1)
                     end do
                  end do
               end do
            end do
         end do
      end do

      !? empty remaining overlap integraps in ovb
      if (.not.lder) then
         if (ldump) then
            if (idxx.eq.0) nrow(ndmp) = irow
            if (idxx.eq.0) call dump(ovb, nrow(ndmp), 0)
            if (idxx.eq.1) call dump(ovb, nrow(ndmp-1), 2)
         endif
      end if
      
   end do

   write(6,*)

   if (.not.lder) return

   enew = e

   call adder(e, dhmx, psi2)

   if (der(1, 1).eq.0.q0) then
      do i=1, nbet
         do k=1, 2
            dhm(i+nh, k) = dhm(i+nh, k) + dhmx(i+nh, k)
            dhmk = dhm(i+nh, k)

            b(i, k) = cpd(i, k)

            der(i, k) = dhmk
            sc = .01q0

            if (i.eq.1.and.j.eq.1) sc = .001q0

            cpd(i, k) = cpd(i, k) - sign(sc*abs(cpd(i, k)), dhmk)
         end do
      end do

      kst = 1
      kstder = 1
   else
      kstder = 0

      do i=1, nbet
         do k=1, 2
            dhm(i+nh, k) = dhm(i+nh, k) + dhmx(i+nh, k)
            dhmk = dhm(i+nh, k)

            bb = cpd(i, k)
            dold = 10.q0*abs(b(i, k)-bb)

            !? rfac = largest allowed fractional change
            rfac = 0.090q0

            if (dold.gt.rfac*bb) dold = rfac*bb

            diff =  - (b(i, k)-bb)*dhmk/(der(i, k)-dhmk)

            if (diff*dhmk.gt.0.q0) diff = -dqreal(1.1q0)*diff
            if (abs(diff).gt.dold) diff = sign(dold, diff)

            sc = 2.q0/16384.q0

            if (bb.lt.0.1q0) sc = 1.q0/16384.q0
            
            if (abs(diff).lt.sc.or.abs(dhmk).lt.1.q-26/bb) diff = sign(sc, diff)
            if (i.eq.1.and.k.eq.1.and.iset.eq.1) diff = 0.q0
            
            if (abs(diff).gt.8.q0*sc) kstder = 1
            if (diff*dhmk.gt.0.q0) diff = -1.1q0*diff
            if (abs(diff).gt.0.035q0*bb/10.q0) kstder = 10

            diff16=diff
            diff = nint(16384*diff16)/16384.q0
            enew = enew + dhmk*diff
            bb16 = bb
            diff16=diff
            cpd(i, k) = nint(16384*(bb16 + diff16))/16384.q0
            b(i, k) = bb

            der(i, k) = dhmk
         end do
      end do

      if (abs((enew-e)/e).lt.0.5q-16) kst = 0
      if (kstder.eq.0) kst = 0
      if (kstder.eq.10) kst = 1
   end if
 
   write(*, 71) qreal(e), qreal(enew), qreal(e-enew)
   write(4, 71) qreal(e), qreal(enew), qreal(e-enew)

   call reset(dhm, e, 0, bscale)

   write(*, 72) loop,((qreal(dhm(i+nh, k)), k=1, 2), i=1, nbet)
   write(4, 72) loop,((qreal(dhm(i+nh, k)), k=1, 2), i=1, nbet)

   write(*, 73) ((qreal(b(i, k)), k=1, 2), i=1, nbet)
   write(*, 73) ((qreal(cpd(i, k)), k=1, 2), i=1, nbet)
   write(4, 73) ((qreal(b(i, k)), k=1, 2), i=1, nbet)
   write(4, 73) ((qreal(cpd(i, k)), k=1, 2), i=1, nbet)

   cpd(i, k) = dqreal(j)/sc

   do i=1, 10
      write(4,'(''-------------------------------------------------'')')
   end do

end subroutine

!* cross calculates the coefficients of (p(l)(costheta12)f(a, b, c))
!* radial integrals arising from the angular reduction of:
!* 1  :  1
!* 2  :  1 with cpl(1, 1) subtracted from csum
!* 3  :  1 with cpl(1, 1) set equal to zero
!* 4  :  [l*(l+1) - Lmax(Lmax+1)]*cpl(l, 1)
!* 5  :  cos(theta)
!* 6  :  del1.del2
!* 7  :  r2.del1(y)
!* 8  :  r1.del2(y)
!*
!* lmax is the largest contributing value of l+1, and lmin the smallest
!* kde = 1 for exchange terms
subroutine cross(l1px, l2px, l1x, l2x, l, kde)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)

   logical cpl_ik_gt_1em10

   common/p1p/cpl(9, 11), csum(11), lmin(11), lmax(11), ict, na, nc, nad, ncd

70 format(' ERROR IN CROSS, LMAX = ', i4)

   cpl_ik_gt_1em10 = .false.
   
   do k=1, 11
      lmin(k) = 0
      lmax(k) = 0
      csum(k) = 0.q0
      do i=1, 9
         cpl(i, k) = 0.q0
      end do   
   end do

   l1 = l1x
   l2 = l2x
   l1p = l1px
   l2p = l2px
   
   pow = 1.q0
   if (kde.eq.1) pow = (-1.q0)**iabs(l1p+l2p-l)
   
   tt0 = 0.5q0*sqrt((dqreal(2*l1+1)) &
   &   * dqreal(2*l1p+1)*dqreal(2*l2+1)*dqreal(2*l2p+1))*pow
   
   do i2=1, 1
      if (i2.ne.1) then
         l1 = l1px
         l2 = l2px
         l1p = l1x
         l2p = l2x
      end if

      lam1a = l1
      lam1b = l1 + 2
      if (l1.eq.0) lam1a = 2

      lam2a = l2
      lam2b = l2 + 2
      if (l2.eq.0) lam2a = 2

      bl1 = dqreal(l1 + 1)
      
      do i=lam1a, lam1b, 2
         lm1 = i - 1
         if (i.eq.lam1b) bl1 = dqreal(-l1)
         bl2 = dqreal(l2 + 1)

         do j=lam2a, lam2b, 2
            lm2 = j - 1

            if (j.eq.lam2b) bl2 = dqreal(-l2)
            if (lm1+lm2.lt.l.or.iabs(lm1-lm2).gt.l) cycle

            tt = dqreal(-(2*lm1+1)*(2*lm2+1))*tt0*f3j0(l1, 1, lm1) &
            &  * f3j0(l2, 1, lm2)*f6j(l, lm1, lm2, 1, l2, l1)

            ttx = tt*bl1*bl2
            tt1 = tt*bl1
            tt2 = tt*bl2

            lcapa = max0(iabs(l1p-lm1), iabs(l2p-lm2)) + 1
            lcapb = min0(l1p+lm1, l2p+lm2) + 1

            if (lcapa.gt.lcapb) cycle

            do k=lcapa, lcapb, 2
               lc = k - 1
               xx = (-1.q0)**lc*(2*lc+1)*f3j0(l1p, lm1, lc) &
               &  * f3j0(l2p, lm2, lc)*f6j(l, lm1, lm2, lc, l2p, l1p)

               cpl(k, 5) = cpl(k, 5) + xx*tt
               cpl(k, 6) = cpl(k, 6) + xx*ttx
               cpl(k, 7) = cpl(k, 7) + xx*tt1
               cpl(k, 8) = cpl(k, 8) + xx*tt2
            end do
         end do
      end do
   end do
   
   lmin1 = max0(iabs(l1p-l1), iabs(l2p-l2)) + 1
   lmax1 = min0(l1p+l1, l2p+l2) + 1

   if (lmin1.le.lmax1) then
      do i=lmin1, lmax1, 2
         lc = i - 1
         cpl(i, 1) = dqreal((-1.q0)**(l+lc))*f3j0(l1p, l1, lc)*f3j0(l2p, l2, lc)          &
         &   *f6j(l, l1, l2, lc, l2p, l1p)*tt0*dqreal(2*lc+1)
         cpl(i, 2) = cpl(i, 1)
         if (i.gt.1) cpl(i, 3) = cpl(i, 1)
         cpl(i, 4) = dqreal(i*(i-1) - lmax1*(lmax1-1))*0.5q0*cpl(i, 1)
      end do
   end if
   
   do k=1, 4
      do i=1, lmax1
         if (cpl(i, k).ne.0.q0) then
            lmin(k) = i
            lmax(k) = lmax1
            exit
         end if
      end do
   end do
   
   lmin1 = lmin(1)
   lmax(4) = lmax1 - 2

   cpl(lmax1, 4) = 0.q0

   if (lmax(4).lt.lmin(4).or.lmin(4).eq.0) lmax(4) = 0
   lmin5 = lmin1 + 1

   if (lmin5.gt.2) lmin5 = lmin5 - 2
   lmax5 = lmax(1) + 1

   do k=5, 8
      do i=1, lmax5
         j = lmax5 - i + 1
         if (abs(cpl(j, k)).gt.1.q-10) exit
         cpl(j, k) = 0.q0
      end do
      
      lmax(k) = j

      do i=1, lmax5
         if (abs(cpl(i, k)).gt.1.q-10) then
            cpl_ik_gt_1em10 = .true.
            exit
         end if
         cpl(i, k) = 0.q0
      end do

      if (.not.cpl_ik_gt_1em10) then 
         lmax(k) = 0
         i = 0
      end if

      lmin(k) = i
      if (lmax(k).gt.9) then
         write(6, 70) lmax(k)
         stop
      end if
   end do

   do k=1, 4
      do i=lmin1, lmax1, 2
         if (k.eq.2.and.i.eq.1) cycle
         csum(k) = csum(k) + cpl(i, k)
      end do

      if (abs(csum(k)).lt.1.q-10) csum(k) = 0.q0
   end do

   do k=5, 8
      do i=lmin5, lmax5, 2
         csum(k) = csum(k) + cpl(i, k)
      end do

      if (abs(csum(k)).lt.1.q-10) csum(k) = 0.q0
   end do

end subroutine


subroutine dump(ovb, nrow, id)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)

   dimension ovb(nrow)

   if (id.eq.0) write(8) (ovb(i), i=1, nrow)
   if (id.eq.1) read(8) (ovb(i), i=1, nrow)
   if (id.eq.2) write(9) (ovb(i), i=1, nrow)
   if (id.eq.3) read(9) (ovb(i), i=1, nrow)

end subroutine

subroutine facab(ia, a, b)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)

   common/s1s/faca(80), facb(80), ngv(100, 10), ng(100), &
   &   igrp(698), ngp, nngt

   dimension an(84)

   itot = ia + 3
   x = 1.q0

   call star2(a, itot+2, an)
   do i=1, itot
      x = x*dqreal(i)
      faca(i) = x*an(i+2)
   end do

   x = 1.q0

   call star2(b, itot+3, an)
   do i=1, itot
      x = x*dqreal(i)
      facb(i) = x*an(i+2)
   end do

end subroutine

!* calculates the 3j symbol (j1 j2 j3)
!*                          ( 0  0  0)
function f3j0(j1, j2, j3)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)

   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84),  &
   &  fm2(84), f21x(84, 86)

70 format(' J1+J2+J3 NOT EVEN', 3i5)

   jj = j1+j2+j3
   jj2 = jj/2 + 1

   if (2*jj2-2.ne.jj) then
      f3j0 = 0.q0
      write(6, 70)j1, j2, j3
   end if

   f3j0 = dqreal(-(-1.q0)**jj2)*sqrt(del2vw(j1, j2, j3, jsig))*fac(jj2) &
   &    / (fac(jj2-j1)*fac(jj2-j2)*fac(jj2-j3))
   return

end function

!* calculates the 6j symbol (a b c)
!*                          (d e f)
function f6j(ia, ib, ic, id, ie, if)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)

   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), &
   &  fm2(84), f21x(84, 86)

   dimension ipos(4), ineg(3)

70 format("36h ***** fault in call of wbar01, isig=, i2")

   outfac = del2vw(ia, ib, ic, isig)

   if (isig.ne.0) go to 10
   outfac = outfac*del2vw(id, ie, ic, isig)

   if (isig.ne.0) go to 10
   outfac = outfac*del2vw(id, ib, if, isig)

   if (isig.ne.0) go to 10
   outfac = outfac*del2vw(ia, ie, if, isig)

   if (isig.ne.0) go to 10
   outfac = sqrt(outfac)

   ipos(1) = ia+ib+ic
   ipos(2) = id+ie+ic
   ipos(3) = id+ib+if
   ipos(4) = ia+ie+if
   ineg(1) = ia+ib+id+ie
   ineg(2) = ia+ic+id+if
   ineg(3) = ib+ic+ie+if

   kmin = max0(ipos(1), ipos(2), ipos(3), ipos(4))
   kmax = min0(ineg(1), ineg(2), ineg(3))

   sum = 0.q0
   ign = (1-2*(kmin-2*(kmin/2)))

   do k=kmin, kmax
      tot = 1.q0
      do kk=1, 4
         kkk = k-ipos(kk)
         tot = tot*fac(kkk+1)
      end do

      do kk=1, 3
         kkk = ineg(kk)-k
         tot = tot*fac(kkk+1)
      end do

      sum = sum + dqreal(ign)*fac(k+2)/tot
      ign = -ign
   end do
   
   f6j = outfac*sum
   return

10 continue
   if (isig.lt.0) then
      f6j = 0.0q0
      return
   else
      write(4, 70) isig
      f6j = 0.0q0
      return
   end if

end function

!* isig is a signal as follows:
!*    isig = -1 : triangle inequality not satisfied (no write)
!*         = +1 : an angular momentum argument negative (write)
!*         = +2 : sum of arguments non-integral (write)
!*         = 0 : otherwise and del2vw is evaluated
function del2vw(ia, ib, ic, isig)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)

   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), &
   &  fm2(84), f21x(84, 86)

70 format(43h ***** negative angular momentum in del2vw.)
   
   if (ia.lt.0 .or. ib.lt.0 .or. ic.lt.0) then
      write(4, 70)

      isig = +1
      return
   end if

   ism = ia+ib+ic

   if (ia.lt.iabs(ib-ic) .or. ia.gt.ib+ic) then
      isig = -1
      return
   end if

   k1 = ia+ib-ic
   k2 = ib+ic-ia
   k3 = ic+ia-ib

   del2vw = fac(k1+1)*fac(k2+1)*fac(k3+1)/fac(ism+2)
   isig = 0

end function


subroutine eigall(nrow, norder, matrix, eigval, option, eigvec, dpwork,  ierr)

   use dqmodule   

   implicit type (dq_real) (a-h, o-z)
   
   type (dq_real) matrix
   
   dimension matrix(1), eigval(1), eigvec(nrow, 1), dpwork(1)
   
   integer nrow, norder, option, ierr, j, imax, i, i1, n

70 format(/' FATAL ERROR IN SUBROUTINE EIGALL.  SUBROUTINE TERMINATED&
   & DUE TO IMPROPER INPUT:'//' NORDER .GT. NROW.'/)
71 format(/' FATAL ERROR IN SUBROUTINE EIGALL.  SUBROUTINE TERMINATED&
   & DUE TO EISPACK ERROR:'//' IERR SET EQUAL TO ', i10,' BY IMTQL1.'/)
72 format(/' FATAL ERROR IN SUBROUTINE EIGALL.  SUBROUTINE TERMINATED&
   & DUE TO EISPACK ERROR:'//' IERR SET EQUAL TO ', i10,' BY IMTQL2.'/)
   
   !? first executable statement - check input
   if (option .eq. 0) then

   !? find eigenvalues only --------------------------------------------------------
      !?    copy matrix into the first norder*(norder + 1)/2 storage locations
      !?    in dpwork
      j = 1
      imax = norder*(norder+1)/2
      do i = 1, imax
         dpwork(j) = matrix(j)
         j = j+1
      end do

      !?    calculate locations for temporary storage in dpwork
      i1 = 1+norder*(norder+1)/2

      !?    transform the real symmetric matrix whose eigenvalues are wanted to
      !?    tridiagonal form via orthogonal similarity transformations
      nv = norder*(norder+1)/2
      call tred3(norder, nv, dpwork(1), eigval, dpwork(i1), dpwork(i1))

      !?    determine the eigenvalues of the symmetric tridagonal matrix using
      !?    the implicit ql method
      call imtql1(norder, eigval, dpwork(i1), ierr)

      if (ierr .ne. 0) write(6, 71) ierr

      return
   end if

   if (norder .gt. nrow) then
      write(6, 70)
      return
   end if
   
   !?    copy matrix into eigenvector
   k = 1
   do i = 1, norder
      do j = 1, i
         temp = matrix(k)
         k = k+1
         eigvec(i, j) = temp
         eigvec(j, i) = temp
      end do
   end do

   !?    transform the real symmetric matrix whose eigenvalues are wanted to
   !?    tridiagonal form via orthogonal similarity transformations
   call tred2(nrow, norder, eigvec, eigval, dpwork, eigvec)

   !?    determine the eigenvalues and eigenvectors of the symmetric tridagonal 
   !?    matrix using the implicit ql method, and backtransform the eigenvectors
   !?    to the original basis
   call imtql2(nrow, norder, eigval, dpwork, eigvec, ierr)

   if (ierr .ne. 0) write(6, 72) ierr

end subroutine

!* this subroutine is a translation of the algol procedure imtql1,
!* Num. Nath. 12, 377-383 (1968) by Martin and Wilkinson, as modified in
!* Num. Math. 15, 450 (1970) by Dubrulle. Handbook for Auto. Comp.,
!* vol.II-Linear Algebra, 241-248 (1971)
!*
!* this subroutine finds the eigenvalues of a symmetric tridiagonal
!* matrix by the implicit ql method.
!*
!* on input:
!*    n is the order of the matrix
!*    d contains the diagonal elements of the input matrix 
!*    e contains the subdiagonal elements of the input matrix in its last
!*       n-1 positions. e(1) is arbitrary
!*
!* on output:
!*    d contains the eigenvalues in ascending order. if an error exit is
!*       made, the eigenvalues are correct and ordered for indeces 1, 2,
!*       ..., ierr-1, but may not be the smallest eigenvalues
!*    e has been destroyed
!*    ierr is set to:
!*       0  :  normal return
!*       j  :  if the j-th eigenvalue has not been determined after 30
!*             iterations
!*
!* questions and comments should be directed to B.S. Garbow, applied
!* mathematics division, Argonne National Laboratory
subroutine imtql1(n, d, e, ierr)

   use dqmodule
   
   implicit type (dq_real) (a-h, o-z)   
   
   real*16 machep
   integer i, j, l, m, n, ii, mml, ierr
   logical p_ge_d_im1

   dimension d(n), e(n)

   p_ge_d_im1 = .false.
  
   !? machep is a machine dependent parameter specifying the relative precision of
   !? floating point arithmetic. = to the smallest number such that machep+1 != 1
   !? neps is an integer such that machep = (1/2)**neps
   call geteps(machep, neps)
   
   ierr = 0

   if (n .eq. 1) return
   
   do i = 2, n
      e(i-1) = e(i)
   end do
   
   e(n) = 0.0q0
   
   do l = 1, n
      j = 0
   !? look for small sub-diagonal element
      do while (.true.)
         do mm = l, n
            m = mm
            if (m .eq. n) exit
            if (abs(e(m)) .le. machep * (abs(d(m)) + abs(d(m+1)))) exit
         end do
      
         p = d(l)

         if (m .eq. l) exit

         if (j .eq. 180) then
            ierr = l
            return
         end if

         j = j + 1
         !? form shift
         g = (d(l+1) - p) / (2.0q0 * e(l))
         r = sqrt(g*g+1.0q0)
         g = d(m) - p + e(l) / (g + signx(r, g))
         s = 1.0q0
         c = 1.0q0
         p = 0.0q0
         mml = m - l

         do ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            if (abs(f) .ge. abs(g)) then
               c = g / f
               r = sqrt(c*c+1.0q0)
               e(i+1) = f * r
               s = 1.0q0 / r
               c = c * s
            else
               s = f / g
               r = sqrt(s*s+1.0q0)
               e(i+1) = g * r
               c = 1.0q0 / r
               s = s * c
            end if
            
            g = d(i+1) - p
            r = (d(i) - g) * s + 2.0q0 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
         end do
      
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0q0
      end do
   
      !? order eigenvalues
      if (l .ne. 1) then
         do ii = 2, l
            i = l + 2 - ii

            if (p .ge. d(i-1)) then
               p_ge_d_im1 = .true.
               exit
            end if

            d(i) = d(i-1)
         end do
      end if
      
      if (.not.p_ge_d_im1) i = 1
      d(i) = p
   end do

end subroutine

!* this subroutine is a translation of the algol procedure imtql2,
!* Num. Nath. 12, 377-383 (1968) by Martin and Wilkinson, as modified in
!* Num. Math. 15, 450 (1970) by Dubrulle. Handbook for Auto. Comp.,
!* vol.II-Linear Algebra, 241-248 (1971)
!*
!* this subroutine finds the eigenvalues and eigenvectors of a symmetric 
!* tridiagonal matrix by the implicit ql method. the eigenvectors of a
!* full symmetric matrix can also be found if tred2 has been used to 
!* reduce this full matrix to tridiagonal form
!*
!* on input:
!*    nm must be set to the row dimension of two-dimensional array parameters
!8       as declared in the calling program dimension statement
!*    n  is the order of the matrix
!*    d  contains the diagonal elements of the input matrix 
!*    e  contains the subdiagonal elements of the input matrix in its last
!*       n-1 positions. e(1) is arbitrary
!*    z  contains the transformation matrix produced in the reduction by
!*       tred2, if performed. if the eigenvectors of the tridiagonal matirx
!*       are desired, z must contain the identity matrix.
!*
!* on output:
!*    d     contains the eigenvalues in ascending order. if an error exit is
!*          made, the eigenvalues are correct and ordered for indeces 1, 2,
!*          ..., ierr-1, but may not be the smallest eigenvalues
!*    e     has been destroyed
!*    z     contains othonormal eigenvectors of the symmetric tridiagonal (or full)
!*          matrix. if an aeeror exit is made, z contains the eigenvectors
!*          associated with the stored eigenvalues         
!*    ierr  is set to:
!*       0  :  normal return
!*       j  :  if the j-th eigenvalue has not been determined after 30
!*             iterations
!*
!* questions and comments should be directed to B.S. Garbow, applied
!* mathematics division, Argonne National Laboratory
subroutine imtql2(nm, n, d, e, z, ierr)
   
   use dqmodule
   
   implicit type (dq_real) (a-h, o-z)
   
   real*16 machep
   integer i, j, k, l, m, n, ii, nm, mml, ierr
   
   dimension d(n), e(n), z(nm, n)

   !? machep is a machine dependent parameter specifying the relative precision of
   !? floating point arithmetic. = to the smallest number such that machep+1 != 1
   !? neps is an integer such that machep = (1/2)**neps
   call geteps(machep, neps)
   
   ierr = 0
   if (n .eq. 1) then 
      write(*,'(A1)') '|'
      return
   end if
   
   do i = 2, n
      e(i-1) = e(i)
   end do
   
   e(n) = 0.0q0
   
   do l = 1, n
      if (mod(l, 10).eq.0) write(*,'(A1,$)') '.'
      j = 0
      
      !? lok for small sub-diagonal element
      do while (.true.)
         do mm = l, n
            m = mm
            if (m .eq. n) exit
            if (abs(e(m)) .le. machep * (abs(d(m)) + abs(d(m+1)))) exit
         end do
      
         p = d(l)
         if (m .eq. l) cycle

         if (j .eq. 180) then
            ierr = l
            write(*,'(A1)') '|'
            return
         end if

         j = j + 1

         !? form shift
         g = (d(l+1) - p) / (2.0q0 * e(l))
         r = sqrt(g*g+1.0q0)
         g = d(m) - p + e(l) / (g + sign(r, g))
         s = 1.0q0
         c = 1.0q0
         p = 0.0q0
         mml = m - l

         do ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            if (abs(f) .ge. abs(g)) then
               c = g / f
               r = sqrt(c*c+1.0q0)
               e(i+1) = f * r
               s = 1.0q0 / r
               c = c * s
            else
               s = f / g
               r = sqrt(s*s+1.0q0)
               e(i+1) = g * r
               c = 1.0q0 / r
               s = s * c
            end if

            g = d(i+1) - p
            r = (d(i) - g) * s + 2.0q0 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b

            !? form vector
            do k = 1, n
               f = z(k, i+1)
               z(k, i+1) = s * z(k, i) + c * f
               z(k, i) = c * z(k, i) - s * f
            end do
         end do
      
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0q0
      end do
   end do

   !! the following code is unreachable and has been commented out.
   !! it's been cleaned up regardless, in case it ends up being needed
   
   !? order eigenvalues and eigenvectors
   ! do ii = 2, n
   !    i = ii - 1
   !    k = i
   !    p = d(i)
   
   !    do j = ii, n
   !       if (d(j) .ge. p) cycle
   !       k = j
   !       p = d(j)
   !    end do
   
   !    if (k .eq. i) cycle
   !    d(k) = d(i)
   !    d(i) = p
   
   !    do j = 1, n
   !       p = z(j, i)
   !       z(j, i) = z(j, k)
   !       z(j, k) = p
   !    end do
   ! end do

   ! write(*,'(A1)') '|'
   
end subroutine


subroutine geteps(machep, neps)
   
   use dqmodule
   
   implicit type (dq_real) (a-h, o-z)
   
   common/a1a/fac(91), r(1250), y(84), yp(84), y1x(84), y2x(84), z1x(84), &
   &  ts(3, 3, 2), ermac, alf, bet, xap, xbp, xam, xbm, suml(86, 2, 6), fp2(84), &
   &  fm2(84), f21x(84, 86)
   
   integer neps
   
   real*16 machep
   
   neps = 162
   machep = ermac

end

!* this subroutine is a translation of the algol procedre tred2, Num. Math. 11, 
!* 181-195 (1968) by Martin, Reinsch, and Wilkinson
!* Handbook for Auto. Comp., Vol. II - Linear Algebra, 212-226 (1971)
!*
!* this subroutine reduces a real symmetric matrix to a symmetric tridiagonal
!* matrix using and accumulating orthogonal similarity transformations
!*
!* on input :
!*    nm must be set to the row dimension of two-dimenstional array parameters as
!*       declared in the calling program dimenstion statement
!*    n  is the order of the matrix
!*    a  contains the real symmetric input matrix. only the lower triangle of
!*       the matrix needs to be supplied
!*
!* on output :
!*    d contains the diagonal elements of the tridiagonal matrix
!*    e contains the subdiagonal elements of the tridagonal matrix in its last n-1
!*      positions. e(1) is set to zero
!*    z contains the orthogonal transformation matrix produced in the reduciton
!*
!*    a and z may coincide. if distinct, a is unaltered
!*
!* questions and comments should be directed to B.S. Garbow, applied
!* mathematics division, Argonne National Laboratory
subroutine tred2(nm, n, a, d, e, z)

   use dqmodule
   
   implicit type (dq_real) (a-h, o-z)
   
   integer i, j, k, l, n, ii, nm, jp1
   
   dimension a(nm, n), d(n), e(n), z(nm, n)

   logical scale_ne_zero

   scale_ne_zero = .false.

   do i = 1, n
      do j = 1, i
         z(i, j) = a(i, j)
      end do
   end do
   
   if (n .ne. 1) then
      do ii = 2, n
         if (mod(ii, 10).eq.0) write(*,'(A1,$)') '.'

         i = n + 2 - ii
         l = i - 1
         h = 0.0q0
         scale = 0.0q0

         if (l .ge. 2) then
            do k = 1, l
               scale = scale + abs(z(i, k))
            end do
         
            if (scale .ne. 0.0q0) scale_ne_zero = .true.
         end if

         if (.not.scale_ne_zero) then
            e(i) = z(i, l)
         else
            do k = 1, l
               z(i, k) = z(i, k) / scale
               h = h + z(i, k) * z(i, k)
            end do
      
            f = z(i, l)
            g = -signx(sqrt(h), f)
            e(i) = scale * g
            h = h - f * g
            z(i, l) = f - g
            f = 0.0q0
      
            do j = 1, l
               z(j, i) = z(i, j) / h
               g = 0.0q0

               do k = 1, j
                  g = g + z(j, k) * z(i, k)
               end do
         
               jp1 = j + 1

               if (l .ge. jp1) then
                  do k = jp1, l
                     g = g + z(k, j) * z(i, k)
                  end do
               end if

               e(j) = g / h
               f = f + e(j) * z(i, j)
            end do
      
            hh = f / (h + h)

            do j = 1, l
               f = z(i, j)
               g = e(j) - hh * f
               e(j) = g
         
               do k = 1, j
                  z(j, k) = z(j, k) - f * e(k) - g * z(i, k)
               end do
            end do
         end if

         d(i) = h
      end do
   end if

   d(1) = 0.0q0
   e(1) = 0.0q0

   write(*,'(A1,$)') '|'

   do i = 1, n
      if (mod(i, 10).eq.0) write(*,'(A1,$)') '.'

      l = i - 1
      
      if (d(i) .ne. 0.0q0) then
         do j = 1, l
            g = 0.0q0
      
            do k = 1, l
               g = g + z(i, k) * z(k, j)
            end do
      
            do k = 1, l
               z(k, j) = z(k, j) - g * z(k, i)
            end do
         end do
      end if
   
      d(i) = z(i, i)
      z(i, i) = 1.0q0
      if (l .lt. 1) cycle
   
      do j = 1, l
         z(i, j) = 0.0q0
         z(j, i) = 0.0q0
      end do
   
   end do
   
   write(*,'(A1)') '|'
   return

end subroutine

!* this subroutine is a translation of the algol procedure tred3, Num. Math 11,
!* 181-195 (1968) by Matrin, Reinsch, and Wilkinson. Handbook for Auto. Comp.,
!* Vol. II - Linear Algebra, 212-226 (1971)
!* 
!* This subroutine reduces a real symmetric matrix, stored as a one-dimensional
!* array, to a symmetric tridiagonal matrix using orthogonal similarity
!* transformations.
!* 
!* on input : 
!*    n  is the order of the matrix
!*    nv must be set to the dimension of the array parameter "a" as declared
!*       in the calling program dimension statement
!*    a  contains the lower triangle of the quadruple precision symmetric input
!*       matrix, stored row-wise as a one dimensional array in its first n*(n+1)/2
!*       positions
!* on output :
!*    a  contains information about the orthogonal transformations used in the
!*       reduction
!*    d  contains the diagonal elements of the tridiagonal matrix
!*    e  contains the subdiagonal elements of the tridiagonal matrix in its last
!*       n-1 positions. e(1) is set to zero
!*    e2 contains the squares of the corresponding elements of e.
!*       may coincide with e if the squares are not needed
!*
!* questions and comments should be directed to B.S. Garbow, applied
!* mathematics division, Argonne National Laboratory
subroutine tred3(n, nv, a, d, e, e2)
  
   use dqmodule
   
   implicit type (dq_real) (a-h, o-z)
   
   real *16 sum16
   
   integer i, j, k, l, n, ii, iz, jk, nv
   
   dimension a(nv), d(n), e(n), e2(n)

   do ii = 1, n
      i = n + 1 - ii
      l = i - 1
      iz = (i * l) / 2
      
      h = 0.0q0
      scale = 0.0q0
      
      if (l .ge. 1) then
         do k = 1, l
            iz = iz + 1
            d(k) = a(iz)
            scale = scale + abs(d(k))
         end do
      
         if (scale .eq. 0.0q0) then
            e(i) = 0.0q0
            e2(i) = 0.0q0
            
            d(i) = a(iz+1)
            a(iz+1) = scale * sqrt(h)

            cycle
         end if
      end if


      do k = 1, l
         d(k) = d(k) / scale
         h = h + d(k) * d(k)
      end do
   
      e2(i) = scale * scale * h
      f = d(l)
      g = -signx(sqrt(h), f)
      
      e(i) = scale * g
      h = h - f * g
      
      d(l) = f - g
      a(iz) = scale * d(l)
      
      if (l .eq. 1) then
         d(i) = a(iz+1)
         a(iz+1) = scale * sqrt(h)

         cycle
      end if
      
      f = 0.0q0
   
      do j = 1, l
         g = 0.0q0
         jk = (j * (j-1)) / 2
         
         do k = 1, l
            jk = jk + 1
            if (k .gt. j) jk = jk + k - 2
            g = g + a(jk) * d(k)
         end do
         
         e(j) = g / h
         f = f + e(j) * d(j)
      end do
   
      hh = f / (h + h)
      jk = 0

      do j = 1, l
         f = d(j)
         g = e(j) - hh * f
         e(j) = g

         do k = 1, j
            jk = jk + 1
            a(jk) = a(jk) - f * e(k) - g * d(k)
         end do
      end do

      d(i) = a(iz+1)
      a(iz+1) = scale * sqrt(h)
   end do

end subroutine


subroutine dumph(sh, ovb, dso, psi, e, nv1, md, ldump, cdump, ndump)
   
   use dqmodule
   
   implicit type (dq_real) (a-h, o-z)
   
   real *16 sh16
   
   logical ldump, ltrid
   
   common/r1r/rap, egs, ltrid
   
   common/h1h/ch(25), neig, nh, iz, lrgl, nspn, nv, l1max, nbet, nd, n1, n2,    &
   &   n1d, n2d, ltrans, iset, namm
   
   dimension dso(1), sh(nv1), ovb(1), md(12150), psi(1)
   
   character cdump*20
   
   if (ltrans.ne.0.and..not.ldump) call comb(ovb, psi, nv, md,-1, 0)
   if (ltrans.ne.0.and.ldump) call comb(sh, psi, nv, md,-1, 0)
   
   if (ldump) then
      rewind(9)
      do i=1, nv
         mdi = md(i)
         write(9) (sh(mdi+j), j=1, i)
      end do
      rewind(9)
      rewind(8)
      read(8) sh
   endif

   if (ltrans.ne.0) call comb(sh, psi, nv, md,-1, 1)

   de = e - egs
   sum = 0.q0

   do i=1, nv
      mi = md(i)
      mix = mi
      
      if (ldump) then
      
         mix = 0
      read(9) (ovb(j), j=1, i)
      endif
      
      do j=1, i
         mi = mi + 1
         mix = mix + 1
         
         ann = dso(i)*dso(j)
         ovb(mix) = ovb(mix)*ann
         sh(mi) = (sh(mi)*ann - e*ovb(mix))
         
         ps2 = psi(i)*psi(j)
         if (i.ne.j) ps2 = ps2*2.q0

         sum = sum + ps2*sh(mi)
      end do
   end do
   
   sum = sum + egs
   
   write(*,'(A, D12.4)') 'E CHECKSUM ', qreal(sum)

end subroutine


subroutine comb(h, psi, nv, md, itt, khm)

   use dqmodule

   implicit type (dq_real) (a-h, o-z)
   
   common/s1s/faca(80), facb(80), ngv(100, 10), ng(100), &
   &  igrp(698), ngp, nngt

   dimension th(14, 14), tth(14, 14), u(14, 14), psi(1), h(1), md(12150)
   nng = nngt
   
   do j=2, nng
      uj = -1./sqrt((j-1)**2+j-1.q0)
      xx = 0.q0
      nnx = 9*(j-1)
      do i=1, nng
         nnx = nnx + 1
         !? generate inverse transformation if itt < 0
         if (itt.ge.0) then
            if (j.gt.i) u(i, j) = uj
            if (j.lt.i) u(i, j) = 0.q0
         else
            if (j.gt.i) u(j, i) = uj
            if (j.lt.i) u(j, i) = 0.q0
         endif
      end do
      
      u(j, j) = -dqreal(j-1)*uj
   end do
   
   if (itt.ne.1) then
      do i=1, ngp
         nng = ng(i)
         
         do j=1, nng
            jv = ngv(i, j)
            nnj = 10*(j-1)
            nnjj = 10*(j-1)-8
            do k=j, nng
               nnj = nnj + 1
               nnjj = nnjj + 9
               kv = ngv(i, k)
               th(j, k) = h(md(kv)+jv)
               th(k, j) = th(j, k)
            end do
         end do

         uj = 1.q0/sqrt(dqreal(1*nng))

         do j=1, nng
            !? generate inverse transformation if itt < 0
            if (itt.ge.0) then
               u(j, 1) = uj
            else
               u(1, j) = uj
            endif
         end do

         do j=1, nng
            do k=1, nng
               sumh = 0.q0
               nngl = nng
               
               if (k.gt.1.and.itt.ge.0) nngl = k
               
               nnj = j-9
               nnk = 9*(k-1)
               
               do l=1, nngl
                  sumh = sumh + th(j, l)*u(l, k)
               end do
               
               tth(j, k) = sumh
            end do
         end do

         do j=1, nng
            nngl = nng
            if (j.gt.1.and.itt.ge.0) nngl = j
            
            do k=j, nng
               sumh = 0.q0
               nnj = 9*(j-1)
               nnk = 9*(k-1)
               
               do l=1, nngl
                  sumh = sumh + u(l, j)*tth(l, k)
               end do
               
               jv = ngv(i, j)
               kv = ngv(i, k)
               h(md(kv)+jv) = sumh
            end do
         end do
         
         ic = 1
         
         do j=1, nv
            if (j.eq.ngv(i, ic)) then
               ic = ic + 1
               if (ic.gt.nng) ngv(i, ic) = 0
               cycle
            end if
            
            do l=1, nng
               lv = ngv(i, l)
               if (j.lt.lv) then
                  th(l, 1) = h(md(lv)+j)
               else
                  th(l, 1) = h(md(j)+lv)
               end if
            end do
            
            do k=1, nng
               kv = ngv(i, k)
               sumh = 0.q0
               nngl = nng
               
               if (k.gt.1.and.itt.ge.0) nngl = k
               nnk = 9*(k-1)
               nnl = 0
               
               do l=1, nngl
                  sumh = sumh + u(l, k)*th(l, 1)
               end do

               if (j.lt.kv) then
                  h(md(kv)+j) = sumh
               else
                  h(md(j)+kv) = sumh
               end if
            end do
         end do
      end do

      return
   end if

   !? transform the eigenvector
   do i=1, ngp
      nng = ng(i)
      uj = 1.q0/sqrt(dqreal(1*nng))
      
      do j=1, nng
         u(j, 1) = uj
      end do
      
      do j=1, nng
         jv = ngv(i, j)
         th(j, 1) = psi(jv)
      end do
      
      do j=1, nng
         jv = ngv(i, j)
         sum = 0.q0
         nnj = j-9
         
         do k=1, nng
            nnj = nnj + 9
            sum = sum + u(j, k)*th(k, 1)
         end do
         
         psi(jv) = sum
      end do
   end do
   
end subroutine

!* updates the input file for program restart
subroutine reset(dhm, e, kontrol, bscale)
   
   use dqmodule

   implicit type (dq_real) (a-h, o-z)
   
   integer pam, qam

   real *16 cp116, cp216, e16, b16, cpd116, cpd216

   dimension cpd1x16(8), cpd2x16(8), cp116(12150), cp216(12150), cpd116(8), cpd216(8), b16(8, 2)

   common/b1b/cpd1x(8), cpd2x(8), nit, next, kst, key1, ibasis, kbasis, linc, &
   &  nv1, longnv, ksm(8), mar12, mar1, kono, nbetx, status


   common/d1d/cp1(12150), cp2(12150), cpd1(8), cpd2(8), b(8, 2), der(9, 2), &
   &  ndp(12150), ndq(12150), nds(12150), ld1(12150), ld2(12150), pam, qam, maxc, &
   &  ia, nblk(9), nb, nblkx(9), nbx, maxbf, nrow(8998), lngth

   common/h1h/ch(25), neig, nh, iz, lrgl, nspn, nv, l1max, nbet, nd, n1, n2, &
   &  n1d, n2d, ltrans, iset, namm

   common/res/lrglx, fnwve, title

   dimension line(2000), dhm(9, 2)

   character line*82, end*2, fnwve*16, title*56, status*8

70 format(f28.25,',', i2,',', i2)
71 format(i2,',', 3(i1,','),'''', a12,''',''', a52,'''')
72 format(f8.5, 9(',', f8.5))
73 format(10(i2,','))
74 format(1p6d12.4)
75 format(6f12.5)

   iend = 0
   iend2 = 0

   !? check for premature convergence after second iteration
   if (iset.gt.nbet) cpd2(iset-nbet) = cpd2x(iset-nbet)
   if (iset.gt.1.and.iset.le.nbet) cpd1(iset) = cpd1x(iset)

   if (key1.eq.0.and.nit.gt.2.and.next.le.2) kst = - nit
   
   status = ', NDG, dq'
   if (kst.eq.0) status = ', OK , dq'
   if (next.eq.nit) kst = 0
   if (iset.eq.1) cpd1(1) = 1.q0/(linc+1)

   do i=4, 2000
      ii = i - 1
      
      read(3,'(A)', end=01) line(i)
      
      end = line(i)(1:2)
      
      if (end.ne.'0/') cycle
      
      if (iend.ne.0) then
         if (iend2.eq.0) iend2 = ii
      else
         iend = i
      end if
   end do
   
01 continue 

   iend = max0(iend, iend2)

   if (iend.eq.0) then
      ii = ii + 1
      line(ii) = '0/'
   end if

   rewind(3)

   do i=1, 3
      read(3,'(A)') line(i)
   end do

   rewind(3)

   write(3,'(A)') line(1)

   if (kst.eq.0) then
      do i=4, iend
         write(3,'(A)') line(i)
      end do
   end if

   nloop = max0(nit-next, 1)
   
   if (kst.eq.0) nloop = nit
   if (kontrol.ne.0) nloop = -kontrol
   
   e16=e
   
   write(3, 70) e16, nloop, neig
   write(3, 71) iz, lrglx, nspn, l1max, fnwve, title(1:44)//status

   key = 1

   do i=1, nbet
      b16(i, 1)=b(i, 1)
      b16(i, 2)=b(i, 2)
      cpd116(i)=cpd1(i)
      cpd216(i)=cpd2(i)
   end do

   if (der(1, 1).eq.0.q0) key = 0
   if (bscale.gt.0.q0) cpd2(1) = cpd2(1)/bscale
   if (key.eq.1) write(3, 72) (qreal(b(i, 1)), qreal(b(i, 2)), i=1, nbet)
   if (key.eq.0) write(3, 72) (qreal(cpd1(i)), qreal(cpd2(i)), i=1, nbet)
   
   write(3, 73) (ksm(i), i=1, nbetx), key, mar12, mar1, kono

   if (key.ne.0) then
      nh0 = nh
      if (dhm(nbet+nh, 1).eq.0.q0) nh0 = 0
      write(3, 74)((qreal(dhm(i+nh0, k)), k=1, 2), i=1, nbet)

      do i=1, nbet
         do k=1, 2
            b16(i, k)=b(i, k)
         end do
         cpd116(i)=cpd1(i)
         cpd216(i)=cpd2(i)
      end do

      write(3, 75)((qreal(b(i, k)), k=1, 2), i=1, nbet)
      write(3, 75)(qreal(cpd1(i)), qreal(cpd2(i)), i=1, nbet)
   end if

   if (kst.ne.0) iend = 3
   if (iend2.eq.0) write(3,'(''0/'')')

   do i=iend+1, ii
      write(3,'(A)') line(i)
   end do

   endfile(3)
   rewind(3)

   read(3,'(A)')

   if (kst.eq.0) return

   read(3,'(//////)')
   if (nbet.gt.3) read(3,'(//)')

end subroutine

!* subroutine to read the incorrectly ordered overlap integrals form fi(...) and arrange
!* them in the correct order. the loop structure must be i(...) to that in 'ham'
subroutine sort(sh, ovb, dso, psi2, nv, md, key)

   use dqmodule
   
   implicit type (dq_real) (a-h, o-z)
   
   integer pam, qam
   
   dimension sh(1), ovb(1), md(12150), dso(1), psi2(1)
   
   common/d1d/cp1(12150), cp2(12150), cpd(8, 2), b(8, 2), &
   &  der(9, 2), ndp(12150), ndq(12150), nds(12150), ld1(12150), ld2(12150), pam, qam, maxc, &
   &  ia, nblk(9), nb, nblkx(9), nbx, maxbf, nrow(8998), lngth
   
   rewind(9)
   
   loc1 = 131072
   irow = 1.q09
  
   ndmp = 1
   do nhi=1, nbx
      nbx1 = nblkx(nhi) + 1
      nbx2 = nblkx(nhi+1)

      do nhj=1, nhi
         nbx3 = nblkx(nhj) + 1
         nbx4 = nblkx(nhj+1)

         do nbl1=nbx1, nbx2
            if (nhi.eq.nhj) nbx4 = nbl1
            nb1 = nblk(nbl1) + 1
            nb2 = nblk(nbl1+1)

            do nbl2=nbx3, nbx4
               nb3 = nblk(nbl2) + 1
               nb4 = nblk(nbl2+1)

               do i=nb1, nb2
                  i1 = nb4
                  if (nb3.eq.nb1) i1 = i

                  if (key.ne.0) then
                     if (md(i)+i.le.lngth) cycle
                     stop 22
                  end if

                  nadd = i1 - nb3 + 1

                  if (irow+nadd.ge.maxbf) then
                     !? read a block of direct + exchange
                     call dump(ovb, nrow(ndmp), 3)
                     irow = 0
                     ndmp = ndmp + 1
                  end if

                  ir1 = nb3 - irow - 1
                  irow = irow + nadd
                  
                  do j=nb3, i1
                     ibf = j - ir1
                     sh(md(i)+j) = ovb(ibf)/(dso(i)*dso(j))
                  end do
               end do
            end do
         end do
      end do
   end do
end subroutine


subroutine final(sh, ovb, v1, tot, psi, psi2, dso, e, ee, next, nv1, md, kontrol, &
               & ldump, lwcof, lpow)

   use dqmodule
   
   implicit type (dq_real) (a-h, o-z)
   
   logical lr12, lext, ldump, ltrid, lwcof, lpow
   
   integer pam, qam
   
   dimension cpd1x16(8), cpd2x16(8), cp116(12150), cp216(12150), cpd116(8), cpd216(8)
   
   common/r1r/rap, egs, ltrid
   
   common/done/lr12, idone, lext
   
   common/s1s/faca(80), facb(80), ngv(100, 10), ng(100), &
   &  igrp(698), ngp, nngt
   
   common/d1d/cp1(12150), cp2(12150), cpd1(8), cpd2(8), b(8, 2), der(9, 2), &
   &  ndp(12150), ndq(12150), nds(12150), ld1(12150), ld2(12150), pam, qam, maxc, &
   &  ia, nblk(9), nb, nblkx(9), nbx, maxbf, nrow(8998), lngth
   
   common/h1h/ch(25), neig, nh, iz, lrgl, nspn, nv, l1max, nbet, nd, n1, n2, &
   &  n1d, n2d, ltrans, iset, namm
   
   dimension psi2(1), psi(1), dso(nv), sh(nv1), ovb(nv1), md(12150), psid(7), &
   &  tot(nv, nv), v1(nv, nv), ee(nv)
   
   data ni/12150/, klim/6/

70 format('E(0) =', f20.16,' ALPHA =', 5f10.6)
71 format('O(1) =', f20.16,'  BETA =', 5f10.6)
72 format(57x,'  LOOP', i3)
73 format(' EIGENVECTOR IS'/1x, 1p, 7d10.3)
74 format(1h0, 31hhamiltonian matrix elements are )
75 format(2(2i5, d15.7))
76 format(1p5e16.9)
77 format(1x, 1p3d26.19)
78 format(12h eigenvalues )

   esh = e
   kout = 0
   
   write(*, 70) -qreal(esh),(qreal(cpd1(i)), i=1, nbet)
   write(*, 71) qreal(psi(1)),(qreal(cpd2(i)), i=1, nbet)

   call dqwrite(6, psi(1))
   
   do i=1, nv
      dso(i) = sqrt(psi(i))
      if (inan(dso(i), 21).eq.1) write(*,*) i, md(i), qreal(ovb(md(i)+i))
   end do

   nb1 = nblk(2) + 1
   nb2 = nblk(3)

   !? renormalization quantities to make dominant part of overlap integral
   !? unity for block 2 basis functions
   !? set lr12 = .false. if they have not been calculated
   lr12 = .true.
   call dnorm(dso, psi)

   do i=1, nv
      mi = md(i)
      do j=1, i
         mi = mi + 1
         ann = dso(i)*dso(j)
         if (.not.ldump) ovb(mi) = ovb(mi)/ann
         sh(mi) = sh(mi)/ann
         if (inan(sh(mi), 2).eq.1) write(*,*) i, j
      end do
   end do
   
   if (ltrans.ne.0) call comb(sh, psi, nv, md, 0, 1)
   call addom(egs, sh, psi, md,-1)
   if (ldump) then
      !? the following was already commented out prior to cleaning up. 
      !? it's been cleaned up and will remain commented out unless needed

      !? store the hamiltonian matrix elements (int ry) in file 8 as a lineari(...)
      !? triangular array of the form sh(md(i)+j) with i >= j

      ! write(4, 74)
      ! write(4, 75) ((i, j, fetch(sh, md(i)+j), j=1, i), i=1, nv)

      rewind(8)
      write(8) sh

      !? recover the unsorted overlap integraps from file 9 and sort them int(...)
      !? normalized form sh(md(i)+j)/qsqrt(ovb(i)*ovb(j))
      call sort(sh, ovb, dso, psi2, nv, md, 0)
      if (ltrans.ne.0) call comb(sh, psi, nv, md, 0, 0)
      call addom(egs, sh, psi, md, 0)
      rewind(9)
      write(9) sh
   else
      if (ltrans.ne.0) call comb(ovb, psi, nv, md, 0, 0)
      call addom(egs, ovb, psi, md, 0)
   end if
   
   if (.not.lpow) then
         call hall(e, iz, sh, ovb, v1, tot, psi, ee, nv, md)
   else
         write(*,*) 'calling POWER'
         call power(sh, ovb, psi, e, ni, nv, nv1, md, klim, kcyc, ldump)
   endif
   
   !? save transferred wave function coefs. in psi2 for use in adder
   nbb2 = nb2
   do i=1, nbb2
      psi2(i) = psi(i)
   end do
         
   if (ltrans.ne.0) call comb(sh, psi, nv, md, 1, 0)
   
   write(*, 72) next
   write(4, 72) next

   if (.not.lpow) write(*,*) 'E =', qreal(e)

   do i=1, 7
      psid(i) = psi(i)
   end do

   write(*, 73) (qreal(psid(i)), i=1, 7)
   write(4, 73) (qreal(psid(i)), i=1, 7)

   do i=1, nv
      psi(i) = psi(i)/dso(i)
   end do

   if (lpow) return
   
   !? tot(i, j) = i'th component of j'th vector
   
   neigx = neig
   
   !? correct for the 2 3s state where neig = 2
   if (lrgl.eq.0.and.nspn.eq.1) neigx = neig - 1

   do i=1, nv
      do j=1, nv
         tot(i, j) = tot(i, j)/dso(i)
         !? scan for NaN's
         if (tot(i, j).ne.tot(i, j)) then
            write(*,*) 'NaN FOR I, J =', i, j, qreal(tot(i, j)), qreal(dso(i))
            read(*,*)
         end if
      end do
   end do

   do i=1, nv
      ee(i) = ee(i)*dqreal(iz**2)/2.q0 - esh
      psi(i) = tot(i, neigx)
   end do

   e = ee(neigx) + esh
   e = 2.q0*e/dqreal(iz**2)

   !? write the wave function coefficients column-wise if lwcof=.true.
   if (lwcof) then
      rewind(8)
      do j=1, nv
         write(8) (qreal(tot(i, j)), i=1, nv)
      end do
   endif

   write(6, 78)
   write(*, 77) (qreal(ee(i)), i=1, 10)

end subroutine


subroutine addom(e, h, psi, md, khm)

   use dqmodule
   
   implicit type (dq_real) (a-h, o-z)
   
   real *16 alf, bet
   
   integer pam, qam
   
   common/s1s/faca(80), facb(80), ngv(100, 10), ng(100), &
   &  igrp(698), ngp, nngt
   
   common/d1d/cp1(12150), cp2(12150), cpd1(8), cpd2(8), b(8, 2), der(9, 2), &
   &  ndp(12150), ndq(12150), nds(12150), ld1(12150), ld2(12150), pam, qam, maxc, &
   &  ia, nblk(9), nb, nblkx(9), nbx, maxbf, nrow(8998), lngth
   
   common/h1h/ch(25), neig, nh, iz, lrgl, nspn, nv, l1max, nbet, nd, n1, n2, &
   &  n1d, n2d, ltrans, iset, namm
   
   common/done/lr12, idone, lext
   
   dimension psi(1), h(1), md(12150), jj(4), tt(6), ttp(6)
   
   equivalence (jj(1), ja),(jj(2), jap),(jj(3), jb),(jj(4), jbp)
   
   equivalence (tt(1), t1),(tt(2), t2),(tt(3), t3),(tt(4), t4),(tt(5), shift), &
   &  (tt(6), scale)
   
   equivalence (ttp(3), c2)
   
   logical lext
   
   ttp(1)=0.q0
   ttp(2)=1.q0
   ttp(3)=0.q0
   ttp(4)=1.q0
   ttp(5)=0.q0
   ttp(6)=0.q0

   if (ltrans.eq.0) return

   itot = ia + 2
   ns1 = ld1(1) + 1
   nl2 = lrgl + neig
   alf = cp1(2)
   bet = cp2(2)
   l2j = 4*ld2(2)*(ld2(2)+1) + 2
   ma = alf*16384 + 0.1
   mb = iz*bet*16384 + 0.1
   bm = dqreal(mb)
   
   apn2 = dqreal((ns1*ma-16384)/16384.q0)
   shift = - e - (dqreal(ns1*ma+16384)*apn2/16384.q0)/dqreal(ns1*ns1)
   c1 = dqreal((iz-1)*nl2*nl2)*bm*4.q0*16384.q0
   c2 = dqreal(nl2*nl2)*bm*bm
   c3 = dqreal((iz-1)*(iz-1)*16384*16384)
   c4 = dqreal(nl2*nl2*iz*iz*16384*16384)
   nb1 = nblk(2) + 1
   nb2 = nblk(3)

   do i=1, ngp
      iv = ngv(i, 1)
      npi = ndp(iv)
      nqi = ndq(iv)
      nsi = nds(iv)
      nqsi = nqi + nsi
      mi = md(iv)

      do j=i, ngp
         jv = ngv(j, 1)
         npj = ndp(jv)
         nqj = ndq(jv)
         nsj = nds(jv)
         nqsj = nqj + nsj
         npp = npi + npj
         nqsp = nqsi + nqsj
         nqsm = nqsi - nqsj

         ja = min0(npi, npj)
         jap = max0(npi, npj)
         jb = min0(nqsi, nqsj)
         jbp = max0(nqsi, nqsj)
         ngi = ng(i)
         ngj = ng(j)

         if (iv-jv.gt.0) then
            ih = mi + jv
         else
            ih = md(jv) + iv
         end if

         if (khm.ne.0) then
            t1 = dqreal(4.q0*npi*npj)*alf*alf + 4.q0*alf*apn2*dqreal(npp+1)
            t2 = dqreal((npp+1)*(npp+2))
            
            t3 = c2*dqreal(-nqsm*nqsm+nqsp+l2j+(ngi-1)*(ngj-1)) - c1*dqreal(nqsp+1) &
            &  + c3*dqreal((nqsp+2)*(nqsp+1))
            
            t4 = c4*dqreal((nqsp+2)*(nqsp+1))
            ttp(4) = t4
         end if

         scale = dqreal(ngi*ngj)

         x1 = 1.q0

         do k=2*ja+3, ja+jap+2
            x1 = x1*dqreal(k/2)
         end do

         x2 = 1.q0

         do k=ja+jap+3, 2*jap+2
            x2 = x2*dqreal(k/2)
         end do
         
         x3 = 1.q0
         
         do k=2*jb+3, jb+jbp+2
            x3 = x3*dqreal(k/2)
         end do
         
         x4 = 1.q0
         
         do k=jb+jbp+3, 2*jbp+2
            x4 = x4*dqreal(k/2)
         end do
         
         ov = sqrt(scale*(x1*x3)/(x2*x4))
         
         if (khm.ne.0) ov = ov*((t1*t4+t2*t3)/(t2*t4) + shift)

         h(ih) = h(ih) + ov

         if (khm.eq.0) cycle
         
         do k=1, ngi
            isk = k*(k-1)
            
            if (k.eq.1) isk = (ngi-1)*(ngi-1)*ngi
            
            kv = ngv(i, k)
            mk = md(kv)
            l1 = 1
            
            if (k.eq.1) l1 = 2
            
            l2 = ngj
            
            if (i.eq.j) l2 = k
            if (l1.gt.l2) cycle
            
            do l=l1, l2
               isl = l*(l-1)
               
               if (l.eq.1) isl = (ngj-1)*(ngj-1)*ngj
               
               lv = ngv(j, l)
               
               if (kv-lv.gt.0) then
                  ih = mk + lv
               else
                  ih = md(lv) + kv
               end if

               h(ih) = h(ih) + sqrt(dqreal(isk*isl)*x1*x3/(x2*x4))*c2/t4
            end do
         end do
      end do
   end do

end subroutine

!* calculates sqrt(overlap integrals) for diagonal dominant terms. these are used
!* to normalize the basis functions
subroutine dnorm(dso, psi)
   
   use dqmodule
   
   implicit type (dq_real) (a-h, o-z)
   
   integer pam, qam
   
   common/s1s/faca(80), facb(80), ngv(100, 10), ng(100), &
   &  igrp(698), ngp, nngt
   
   common/d1d/cp1(12150), cp2(12150), cpd1(8), cpd2(8), b(8, 2), der(9, 2), &
   &  ndp(12150), ndq(12150), nds(12150), ld1(12150), ld2(12150), pam, qam, maxc, &
   &  ia, nblk(9), nb, nblkx(9), nbx, maxbf, nrow(8998), lngth
   
   dimension psi(1), dso(1)

   alf = 2.q0*cp1(2)
   bet = 2.q0*cp2(2)

   call facab(ia, alf, bet)

   if (ngp.eq.0) return
   
   do i=1, ngp
      iv = ngv(i, 1)
      npi = ndp(iv)
      nqsi = ndq(iv) + nds(iv)
      ov = sqrt(dqreal(2.q0)*faca(2*npi+2)*facb(2*nqsi+2))
      nng = ng(i)
      
      do j=1, nng
         iv = ngv(i, j)
         dso(iv) = ov
      end do
   end do

end subroutine


subroutine hall(e, iz, sh, ovb, v1, tot, x, x2, nv, md)
   
   use dqmodule
   
   implicit type (dq_real) (a-h, o-z)
   
   common/r1r/rap, egs, ltrid
   
   logical ltrid
   
   dimension sh(1), ovb(nv), v1(nv, nv), x(nv), x2(nv), tot(nv, nv), &
   &  md(12150)
   
   real*16 sh16, ovb16
   
70 format(1x, 1p5e15.8)
71 format(' ORTHONORMALIZATION COEFS. ARE', i10)
72 format(' BASIS FUNCTION COEFS. ARE', i10)
   
   if (ltrid) then
      call eigall(nv, nv, ovb, x, 1, v1, x2, ierr)
      do j=1, nv
         ovb(md(j)+j) = x(j)
      end do
   else
      iegen=0
      call hdiagl(ovb, nv, iegen, v1, md, nr, nv)
   endif

   write(*, 71) nr

   do i=1, nv
      x1 = 1.q0/sqrt(ovb(i+md(i)))
      do j=1, nv
         v1(j, i) = v1(j, i)*x1
      end do
   end do

   rewind(8)

   do i=1, nv
      write(8) (v1(j, i), j=1, nv)
   end do

   do j=1, nv
      do k=1, nv
         x(k) = v1(k, j)
      end do

      do i=1, nv
         sum = 0.q0
         do k=1, nv
            if (i.le.k) then
               sum = sum + sh(i+md(k))*x(k)
            else
               sum = sum + sh(md(i)+k)*x(k)
            endif
         end do
         v1(i, j) = sum
      end do
   end do
   
   rewind(8)

   do i=1, nv
      !? retrieve i'th column of original v1(k, i)
      read(8) x
      
      do j=i, nv
         sum = 0.q0
         do k=1, nv
            sum = sum + x(k)*v1(k, j)
         end do
         ovb(i+md(j)) = sum
      end do
   end do
   
   if (ltrid) then
      call eigall(nv, nv, ovb, x2, 1, v1, x, ierr)
      do j=1, nv
         ovb(md(j)+j) = x2(j)
      end do
   else
      iegen=0
      call hdiagl(ovb, nv, iegen, v1, md, nr, nd)
   endif

   write(*, 72) nr

   do i=1, nv
      do j=1, nv
         tot(i, j) = 0.q0
      end do
      x2(i) = ovb(i+md(i)) + egs
   end do

   e = x2(1)

   rewind(8)

   do k=1, nv
      !? retrieve k'th column of original v1(i, k)
      read(8) x
      
      do j=1, nv
         vv = v1(k, j)
         do i=1, nv
            !? tot(i, j) = i'th component of j'th vector
            tot(i, j) = tot(i, j) + vv*x(i)
         end do
      end do
   end do

end subroutine 


subroutine diskw(nv, ndim, sh)
   
   use dqmodule
   
   implicit type (dq_real) (a-h, o-z)
   
   real *16 ssh16
   
   dimension sh(ndim), ssh16(ndim)

   ssh16(ndim) = sh(ndim)

   rewind 8
   write(8) ssh16
   
end


subroutine adder(e, dhm, psi2)

   use dqmodule
   
   implicit type (dq_real) (a-h, o-z)
   
   real *16 alf, bet
   integer pam, qam
   logical lext
   
   common/s1s/faca(80), facb(80), ngv(100, 10), ng(100),                 &
   &   igrp(698), ngp, nngt
   
   common/d1d/cp1(12150), cp2(12150), cpd1(8), cpd2(8), b(8, 2), der(9, 2),  &
   &   ndp(12150), ndq(12150), nds(12150), ld1(12150), ld2(12150), pam, qam, maxc,&
   &   ia, nblk(9), nb, nblkx(9), nbx, maxbf, nrow(8998), lngth
   
   common/h1h/ch(25), neig, nh, iz, lrgl, nspn, nv, l1max, nbet, nd, n1, n2,    &
   &   n1d, n2d, ltrans, iset, namm
   
   common/done/lr12, idone, lext
   
   dimension psi2(1), dhm(9, 2), sum(4), jj(6), tt(10), ttp(10), sumb0(4)

   equivalence (jj(1), ja),(jj(2), jap),(jj(3), jb),(jj(4), jbp)
   
   equivalence (tt(1), t1),(tt(2), t2),(tt(3), t3),(tt(4), t4),(tt(5), shift), &
   &  (tt(6), scale),(tt(7), t1a),(tt(8), t2a),(tt(9), t3b), (tt(10), t4b)

   ttp(1)=0.q0
   ttp(2)=1.q0
   ttp(3)=1.q0
   ttp(4)=1.q0
   ttp(5)=0.q0
   ttp(6)=0.q0
   ttp(7)=0.q0
   ttp(8)=1.q0
   ttp(9)=0.q0
   ttp(10)=1.q0

   if (ltrans.eq.0) return

   ns1 = ld1(1) + 1
   nl2 = lrgl + neig
   alf = cp1(2)
   bet = cp2(2)
   l2j = 4*ld2(2)*(ld2(2)+1) + 2
   ma = alf*16384 + 0.1
   mb = iz*bet*16384 + 0.1
   bm = dqreal(mb)

   apn2 = dqreal(ns1*ma-16384)/16384.q0
   shift = - e - (dqreal(ns1*ma+16384)*apn2/16384.q0/dqreal(ns1*ns1))
   c1 = dqreal((iz-1)*nl2*nl2)*bm*4.q0*16384.q0
   c2 = dqreal(nl2*nl2)*bm*bm
   c3 = dqreal((iz-1)*(iz-1))*16384.q0*16384.q0
   c4 = dqreal(nl2*nl2*iz*iz)*16384.q0*16384.q0
   nb1 = nblk(2) + 1
   nb2 = nblk(3)

   sum(1) = 0.q0
   sum(2) = 0.q0

   do i=1, ngp
      iv = ngv(i, 1)
      npi = ndp(iv)
      nqsi = ndq(iv) + nds(iv)
      nsi = nds(iv)

      do j=i, ngp
         jv = ngv(j, 1)
         npj = ndp(jv)
         nqsj = ndq(jv) + nds(jv)
         npp = npi + npj
         nqsp = nqsi + nqsj
         nqsm = nqsi - nqsj
         nsj = nds(jv)
         
         ja = min0(npi, npj)
         jap = max0(npi, npj)
         jb = min0(nqsi, nqsj)
         jbp = max0(nqsi, nqsj)
         ngi = ng(i)
         ngj = ng(j)

         t1 = (4*npi*npj)*alf*alf + 4.q0*alf*apn2*dqreal(npp+1)
         t1a = (4*npi*npj+2*npp)*alf*alf + 4.q0*alf*apn2*dqreal(npp+2)
         t2 = dqreal((npp+1)*(npp+2))
         t2a = dqreal((npp+2)*(npp+3))
         t3 =  c2*dqreal(-nqsm*nqsm+nqsp+l2j+(ngi-1)*(ngj-1)) - c1*dqreal(nqsp+1)  &
         &   + c3*dqreal(nqsp+2)*dqreal(nqsp+1)
         t3b =  c2*dqreal(-nqsm*nqsm+nqsp+l2j+(ngi-1)*(ngj-1)) - c1*dqreal(nqsp+2)  &
         &   + c3*dqreal(nqsp+3)*dqreal(nqsp+2)
         t4 = c4*dqreal(nqsp+2)*dqreal(nqsp+1)
         t4b = c4*dqreal(nqsp+3)*dqreal(nqsp+2)

         scale = dqreal(ngi*ngj)
         
         x1 = 1.q0
         do k=2*ja+3, ja+jap+2
            x1 = x1*dqreal(k/2)
         end do
         
         x2 = 1.q0   
         do k=ja+jap+3, 2*jap+2
            x2 = x2*dqreal(k/2)
         end do
         
         x3 = 1.q0
         do k=2*jb+3, jb+jbp+2
            x3 = x3*dqreal(k/2)
         end do
         
         x4 = 1.q0        
         do k=jb+jbp+3, 2*jbp+2
            x4 = x4*dqreal(k/2)
         end do
         
         ov = sqrt(scale*(x1*x3)/(x2*x4))
         pss = psi2(iv)*psi2(jv)

         if (iv.ne.jv) pss = 2.q0*pss

         sum(1) = sum(1) + ov*pss*dqreal(npp+3)*(t1a/t2a + t3/t4 + shift)
         sum(2) = sum(2) + ov*pss*dqreal(nqsp+3)*(t1/t2 + t3b/t4b + shift)

         sumb0(1) = 0.q0

         jj(5) = 1

         do k=1, ngi
            isk = k*(k-1)
            
            if (k.eq.1) isk = (ngi-1)*(ngi-1)*ngi
            
            kv = ngv(i, k)
            l1 = 1
            
            if (k.eq.1) l1 = 2
            
            l2 = ngj
            
            if (i.eq.j) l2 = k
            if (l1.gt.l2) cycle
            
            do l=l1, l2
               isl = l*(l-1)
               
               if (l.eq.1) isl = (ngj-1)*(ngj-1)*ngj
               
               lv = ngv(j, l)
               pss = psi2(kv)*psi2(lv)
               
               if (kv.ne.lv) pss = 2.q0*pss
               
               sumb0(1) = sumb0(1) + pss*sqrt(dqreal(isk*isl)*x1*x3/(x2*x4))
            end do
         end do

         sum(1) = sum(1) + sumb0(1)*dqreal(npp+3)*c2/t4
         sum(2) = sum(2) + sumb0(1)*dqreal(nqsp+3)*c2/t4b
      end do
   end do

   dhm(2, 1) = dhm(2, 1) - sum(1)/alf
   dhm(2, 2) = dhm(2, 2) - sum(2)/bet

end subroutine

!* finds the closest eigenvalue to the gauss egs using the power method
subroutine power(sh, ovb, x, e, ni, n, nv1, md, ncyc, kcyc, ldump)

   use dqmodule
   
   implicit type (dq_real) (a-h, o-z)
   
   logical ldump, ltrid, lext
   real *16 ssh16, dxm16, et16, chk16, sumh16, e16
   
   common/r1r/rap, egs, ltrid
   
   common/h1h/ch(25), neig, nh, iz, lrgl, nspn, nv, l1max, nbet, nd, n1, n2, &
   &  n1d, n2d, ltrans, iset, namm
   
   common/done/lr12, idone, lext
   
   dimension x(1), y(12150), xx(12150), sh(nv1), ovb(nv1), md(12150), ee(2)
   dimension ssh16(nv1)


70 format(' ITER =', i2,' DXM =', d10.3,' FOR', i4,' E =', d32.25, d12.3)
71 format(1h0,'SYSTEM HAS NOT CONVERGED AFTER ', i3, 7h cycles)
72 format(' DX =', d10.3, i4,'  DE =', d12.5, d32.25, d12.3)
   
   write(4,*) 'Starting POWER, LDUMP = ', ldump

   do i=1, n
      x(i)=0.q0
   end do

   x(1) = 1.q0
   e=0.q0
   dxm = 1.q60

   if (.not.ldump) then
      rewind(8)
      write(8) (sh(i), i=1, nv1)
   endif

   do k=1, ncyc
      if (k.ne.1) then
         if (ldump) then
            rewind 9
            read(9) sh
         endif
         write(4,*) 'OK1 READ(9) SH, K =', k
      end if

      kcyc = k

      do i=1, n
         y(i) = x(i)/x(1)
      end do

      if (ldump) call yax(xx, sh, y, md, n)
      if (.not.ldump) call yax(xx, ovb, y, md, n)

      write(4,*) 'OK2, DONE YAX'

      if (ldump) then
         if (k.ne.1) then
            rewind(10)
            read(10) sh
            write(4,*) 'OK3 READ(10) SH'
         else         
            rewind(8)
            read(8) sh
            write(4,*) 'OK3 READ(8) SH'
         end if
      endif

      write(4,*) 'OK4 calling solve', kcyc
      call solve(sh, xx, x, chk, n, md, nd, kcyc)
      write(4,*) 'done solve', kcyc

      if (k.le.1 .and. ldump) write(10) sh
      
      dxmp = dxm

      if (.not.dxm.eq.dxm)then
            print*, "NAN DETECTED"
      end if

      es=e

      sum = 0.q0
      sum2 = 0.q0
      dxmp = dxm
      dxm = 0.q0

      do i=1, n
         xx1 = x(i)/x(1)
         dx = abs((xx1-y(i))/xx1)

         if (dxm.lt.dx) then
            im = i
            dxm = dx
         end if
         
         sum  = sum  + x(i)*y(i)
         sum2 = sum2 + x(i)*x(i)
      end do

      e = sum/sum2
      et = e + egs

      dxm16=dxm
      et16=et
      chk16=chk
      
      write(*, 70) kcyc, dxm16, im, et16, chk16
      write(4, 70) kcyc, dxm16, im, et16, chk16

      if (kcyc.lt.2) cycle

      if (dxm.ge.dxmp .or. dxm.lt.1.q-26) go to 01
   end do

   write(*, 71) ncyc
   write(4, 71) ncyc

01 continue

   rewind(8)
   if (ldump) rewind(9)
   write(4,*) 'OK5 reading 8', nv1
   read(8) (sh(i), i=1, nv1)
   write(4,*) 'done 8'

   write(4,*) 'OK6 calling yax', n
   call yax(xx, sh, x, md, n)
   write(4,*) 'done yax'

   write(4,*) 'OK7 reading 9'
   if (ldump) read(9) sh
   write(4,*) 'done 9'

   write(4,*) 'OK8 calling yax2', ldump
   if (ldump) call yax(y, sh, x, md, n)
   if (.not.ldump) call yax(y, ovb, x, md, n)
   write(4,*) 'done yax2'

   dxm = 0.q0
   sumh = 0.q0
   sumo = 0.q0

   write(4,*) 'starting loop 32'

   do i=1, n
      dx = x(i)*(xx(i) - e*y(i))
      sumh = sumh + dx
      dy = x(i)*y(i)
      sumo = sumo + dy
      dx = abs(dx/dy)
      
      if (dx.le.dxm) cycle
      
      dxm = dx
      im = i
   end do

   sumh = sumh/sumo
   sumo = sqrt(sumo)

   do i=1, n
      x(i) = x(i)/sumo
   end do

   e = e + egs
   e16=e
   chk16=chk
   sumh16=sumh
   dxm16=dxm

   write(*, 72) dxm16, im, sumh16, e16, chk16
   write(4, 72) dxm16, im, sumh16, e16, chk16
   
end subroutine


subroutine yax(xx, ovb, y, md, n)

   use dqmodule
   
   implicit type (dq_real) (a-h, o-z)
   
   dimension xx(1), ovb(1), y(1), md(12150)
   
   do i=1, n
      sum = 0.q0
      mi = md(i)
      do j=1, n
         if (i.ge.j)sum  = sum  + y(j)*ovb(mi+j)
         if (i.lt.j)sum  = sum  + y(j)*ovb(md(j)+i)
      end do
      xx(i) = sum
   end do
   
end subroutine

!* solves the equation ax = f by the square root method
subroutine solve(a, f, x, chk, n, md, nd, kcyc)

   use dqmodule
   
   implicit type (dq_real) (a-h, o-z)
   
   common/done/lr12, idone, lext
   
   save isgn, ss, a1

   logical lext, ldot
   
   dimension f(1), x(1), ak(12150), ss(12150), isgn(12150), &
   &  a1(12150), md(12150), a(1)

   ldot = .false.

   if (kcyc.le.1) then
      do i=1, n
         a1(i) = a(md(i)+1)
         isgn(i) = 1
      end do

      ss(1) = dqreal(1.q0)/sqrt(abs(a(1)))

      if (a(1).lt.0.q0) isgn(1) = -1
      
      do j=2, n
         a(md(j)+1) = a(md(j)+1)*ss(1)
      end do

      do i=2, n
         sum = 0.q0
         im1 = i - 1
         ip1 = i + 1
         m = md(i)
         
         do l=1, im1
            ak(l) = a(m+l)*dqreal(isgn(l))
            sum = sum + a(m+l)*ak(l)
         end do

         ax = a(m+i) - sum
         ss(i) = dqreal(1.q0)/sqrt(abs(ax))

         if (ax.lt.0.q0) isgn(i) = -1
         if (ip1.gt.n) exit
         
         do j=ip1, n
            sum = 0.q0
            m = md(j)

            do l=1, im1
               sum = sum + ak(l)*a(m+l)
            end do

            a(m+i) = (a(m+i) - sum)*ss(i)
         end do

         if (inan(sum, 4).eq.1) write(*,*) l, m

         if (mod(i, 100).eq.0) then
            if (mod(i, 1000).eq.0) then
               write(*,'(A1,$)') 'x'
            elseif (mod(i, 500).eq.0) then
               write(*,'(A1,$)') '|'
            else
               write(*,'(A1,$)') '.'
            endif

            ldot= .true.
         endif
      end do
   end if

   ak(1) = f(1)*ss(1)

   do i=2, n
      im1 = i - 1
      sum = 0.q0
      m = md(i)

      do l=1, im1
         sum = sum + a(m+l)*ak(l)*dqreal(isgn(l))
      end do
      
      ak(i) = (f(i) - sum)*ss(i)
   end do
   
   x(n) = ak(n)*ss(n)*dqreal(isgn(n))

   do ip=2, n
      i = n - ip + 1
      ip1 = i + 1
      sum = 0.q0

      do l=ip1, n
         sum = sum + a(md(l)+i)*x(l)
      end do

      x(i) = (ak(i) - sum)*ss(i)*dqreal(isgn(i))
   end do

   sum = 0.q0

   do i=1, n
      sum = sum + x(i)*a1(i)
   end do

   chk = sum/f(1) - 1.q0
   
   if (ldot) write(*,*)
   
end subroutine


function signx(a, b)
   
   use dqmodule
   
   implicit type (dq_real) (a-h, o-z)
   
   signx = abs(a)
   if (b.lt.0.q0) signx = -signx !! this line was previously "sgnx = -sgnx" typo?
   
end function

!* detect and signal NaN's
function inan(x, isig)
   
   use dqmodule
   
   implicit type (dq_real) (a-h, o-z)
   
   if (x.ne.x) then
      write(*,'(a, i2)') 'NaN found, ISIG =', isig
      read(*,*)
      inan = 1
   else
      inan = 0
   end if
   
end function


subroutine getdate(date)

   integer ndmy(3)
   character date*8

70 format(i2.2,'/', i2.2,'/', i2.2)
   
   call idate(ndmy)
   
   nd = ndmy(1)
   nm = ndmy(2)
   ny = ndmy(3)
   
   if (ny.gt.2000) ny = ny - 2000
   
   write(date, 70) nm, nd, ny
   
end subroutine


subroutine dqwritx(n, x)
   use dqmodule
   use wavext

   implicit type (dq_real) (a-h, o-z)
   
   character buffer*80

   if (abs(x).gt.1.q50.or.abs(x).lt.1.q-20) return

   call routedopen(2, file='ermacdq.dat', status='NEW')
   rewind(2)

   call dqwrite(2, x)
   rewind(2)

   read(2,'(a)') buffer
   do i=45, 50
      if (buffer(i:i).ne.buffer(i+1:i+1)) then
         if (n.ne.6) write(*,*) 'n =', n
         call dqwrite(6, x)
         exit
      end if
   end do

   close(2, status='DELETE')
   
end subroutine
