c    This code was written by Karl Taylor
c    lightly edited by Stephen Po-Chedley (increase nmont, add intent to ss, jcnt/icnt)
c    I'm hoping he'll help me add some documentation here. 




      subroutine solvmid(alon, alat, nmon, conv, dt, tmin, tmax,
     &     bbmin, maxiter, a, c, obsmean, ss, icnt, jcnt)
      implicit none
      integer nmont, nmon12
      parameter (nmont=500*12, nmon12=12) ! PJD Oceanonly 1870-2014 - 171025
c      parameter (nmont=149*12, nmon12=12) ! PJD Oceanonly 1870-2014 - 170412
c      parameter (nmont=149*12, nmon12=12) ! PJD Oceanonly 1870-2014 - 160414
c      parameter (nmont=147*12, nmon12=12) ! PJD Oceanonly 1870-2014
      integer nmon, maxiter
      integer, intent(inout) :: icnt, jcnt
      real conv, tmin, tmax, dt, bbmin, alon, alat
      real obsmean(nmon), a(nmon), c(nmon)
      real, intent(out) :: ss(nmon)
      integer i, n, imethod, n1, n2, nn, jj, i1, i2, i3, nnn, jend,
     &      kk, j, k, kkk, nm, np
      integer jbeg(nmont)
      real relax, residmax, resid, dxm, dxp, s1, s2, addmax, addmin
      real r(nmont), avg(nmont), aa(nmont), bb(nmont), cc(nmont),
     &     add(nmont)
      double precision s(nmont), sum
      imethod = 1
c ???  check following value
      relax = 1.0
      if (nmon .gt. nmont) then
        print*, 'error-- nmont not declared large enough in '
        print*, 'subroutine solvmid'
        stop
      endif
c    check for occurance where obs monthly means are consecutively
c      at upper and lower limits. If so, smooth data, being careful
c      to preserve annual mean.
      do 48 n=1,nmon
        add(n) = 0.0
   48 continue
      n2 = nmon
      do 50 n=1,nmon
        n1 = n2
        n2 = n
        if ((obsmean(n2)-obsmean(n1)) .gt. (tmax-tmin-2.*dt)) then
          add(n1) = add(n1) +
     &           (0.5*(obsmean(n2)-obsmean(n1)-tmax+tmin) + dt)/c(n1)
          add(n2) = add(n2) -
     &           (0.5*(obsmean(n2)-obsmean(n1)-tmax+tmin) + dt)/a(n2)
        elseif ((obsmean(n1)-obsmean(n2)) .gt. (tmax-tmin-2.*dt)) then
          add(n1) = add(n1) -
     &           (0.5*(obsmean(n1)-obsmean(n2)-tmax+tmin) + dt)/c(n1)
          add(n2) = add(n2) +
     &           (0.5*(obsmean(n1)-obsmean(n2)-tmax+tmin) + dt)/a(n2)
        endif
   50 continue
      nn = 0
      addmax = 0.0
      addmin = 0.0
      do 51 n=1,nmon
        if (add(n) .ne. 0.0) then
           if (jcnt .lt. 5000) then
              print*, 'jump from extreme to extreme at time ', n
              write(9,*) 'jump from extreme to extreme at time ', n
           endif
           if (jcnt .eq. 5000) then
              print*, ' **********************************'
              print*, ' No more jump messages will be written'
              print*, ' **********************************'
           endif
           if (jcnt .lt. 5) then
          print*,
     &     'monthly means go from one limit to the other in 1 month'
          print*, 'alat = ', alat, ' alon = ', alon, ' n = ', n,
     &         ' add = ', add(n)
             n1 = mod((n+nmon-2), nmon) + 1
             n2 = mod(n, nmon) + 1
          print*, 'observed mean for 3 cells: ', obsmean(n1),
     &          obsmean(n), obsmean(n2)
           endif
          if (jcnt .le. 2)
     &        print*, (obsmean(n1), n1=1,nmon)
          addmax = amax1(addmax, add(n))
          addmin = amin1(addmin, add(n))
          obsmean(n) = obsmean(n) + add(n)
          nn = nn + 1
        endif
   51 continue
      icnt = icnt + nn
      if (nn .gt. 0) then
        jcnt = jcnt + 1
c        if (nmon .eq. nmon12) then
c          print*, 'Climatology: '
c          write(9,*) 'Climatology: '
c        endif
        if (jcnt .le. 1000) then
          print*,  nn,
     &      ' monthly values smoothed at lat,lon', alat, alon
          print*, 'max added = ', addmax,
     &             '  max subtracted = ', addmin
          write(9,*)  nn,
     &      ' monthly values smoothed at lat,lon', alat, alon
          write(9,*) 'max added = ', addmax,
     &             '  max subtracted = ', addmin
        endif
        if (jcnt .eq. 50) then
          print*, ' '
          if (nmon .eq. nmon12) then
            print*,
     &      'No more warnings will be printed concerning smoothing '//
     &       'of climatological data'
          else
            print*,
     &      'No more warnings will be printed concerning smoothing '//
     &       'of monthly data'
          endif
          print*, ' '
          write(9,*) ' '
          if (nmon .eq. nmon12) then
            write(9,*)
     &      'No more warnings will be printed concerning smoothing '//
     &       'of climatological data'
          else
            write(9,*)
     &      'No more warnings will be printed concerning smoothing '//
     &       'of monthly data'
          endif
          write(9,*) ' '
        endif
      endif
c    check if all are le tmin or all are ge tmax
      if (obsmean(1) .le. (tmin+0.01*dt)) then
        do 80 i=2,nmon
          if (obsmean(i) .gt. (tmin+0.01*dt)) go to 99
   80   continue
        do 85 i=1,nmon
          ss(i) = tmin
   85   continue
c        if (nmon .eq. nmon12) print*, 'Climatology: '
c        print*, 'all values were at minimum at this grid cell:'
c        print*, 'latitude = ', alat, ' longitude = ', alon
        return
      elseif (obsmean(1) .ge. (tmax-0.01*dt)) then
        do 90 i=2,nmon
          if (obsmean(i) .lt. (tmax-0.01*dt)) go to 99
   90   continue
        do 95 i=1,nmon
          ss(i) = tmax
   95   continue
c        if (nmon .eq. nmon12) print*, 'Climatology: '
c        print*, 'all values were at maximum at this grid cell:'
c        print*, 'latitude = ', alat, ' longitude = ', alon
        return
      endif
   99 jj = 0
      do 100 i=1,nmon
        i1 = i
        i2 = mod(i,nmon) + 1
        i3 = mod((i+1), nmon) + 1
        if (((obsmean(i1) .le. tmin+0.01*dt) .and.
     &       (obsmean(i2) .le. tmin+0.01*dt) .and.
     &       (obsmean(i3) .gt. tmin+0.01*dt)) .or.
     &      ((obsmean(i1) .ge. tmax-0.01*dt) .and.
     &       (obsmean(i2) .ge. tmax-0.01*dt) .and.
     &       (obsmean(i3) .lt. tmax-0.01*dt))) then
          jj = jj + 1
          jbeg(jj) = i2
        endif
  100 continue
      if (jj .eq. 0) then
c       simple cyclic treatment
c         latest approximation of means (given mid-month values)
         nnn = 0
  105    nnn = nnn + 1
         sum = 0.0
         residmax = 0.0
         do 110 n = 1, nmon
           nm = mod((n+nmon-2), nmon) + 1
           np = mod(n, nmon) + 1
           bb(n) = 0.0
           avg(n) = 0.0
           if (nnn .lt. imethod) then
             call approx(tmin, tmax, a(n), c(n), ss(nm), ss(n),
     &                ss(np), aa(n), bb(n), cc(n), avg(n))
           else
             call numer(conv, tmin, tmax, bbmin, a(n), c(n), ss(nm),
     &                ss(n), ss(np), aa(n), bb(n), cc(n), avg(n))
           endif
           r(n) = obsmean(n) - avg(n)
           sum = sum + r(n)**2
           residmax = amax1(residmax, abs(r(n)))
110      continue
         resid = dsqrt(sum)/nmon
         if (residmax .gt. conv) then
           if (nnn .gt. maxiter*0.5) then
             print*, 'iteration = ', nnn, ' residual = ', resid,
     &          ' maximum residual = ', residmax
           endif
           if (nnn .gt. maxiter*0.9) then
             print*, ' '
             print*, 'latitude = ', alat, ' longitude = ', alon
             do 1234 n=1,nmon
               write(*,'(8(1pe10.2))') obsmean(n), avg(n), r(n),
     &              s(n), ss(n), aa(n), bb(n), cc(n)
 1234        continue
           endif
           if (nnn .gt. maxiter) then
             print*, 'latitude = ', alat, ' longitude = ', alon
             print*, 'does not converge'
             write(9,*) 'latitude = ', alat, ' longitude = ', alon
             write(9,*) 'does not converge'
c             call exit(1)
           else
c            solve for new estimate of mid-month values
             call cyclic(alon, alat, aa, bb, cc, cc(nmon), aa(1), r, s,
     &                nmon)
             do 120 n=1,nmon
               ss(n) = ss(n) + relax*s(n)
  120        continue
c            if ss exceeds tmax or tmin, then it should exceed it no
c                more than absolutely necessary:
             do 130 n=1,nmon
               nm = mod((n+nmon-2), nmon) + 1
               np = mod(n, nmon) + 1
               if (ss(n) .gt. tmax) then
                 if (ss(nm) .le. tmax) then
                   dxm = (ss(n)-tmax)/((ss(n)-ss(nm))*a(n))
                 else
                   dxm = 0.0
                 endif
                 if (ss(np) .le. tmax) then
                   dxp = (ss(n)-tmax)/((ss(n)-ss(np))*c(n))
                 else
                   dxp = 0.0
                 endif
                 if ((dxm .gt. 0.5) .and. (dxp .gt. 0.5)) then
                   s1 = tmax + (tmax-ss(nm))*a(n)/(2.-a(n))
                   s2 = tmax + (tmax-ss(np))*c(n)/(2.-c(n))
                   ss(n) = amin1(s1, s2)
                 elseif ((dxm .eq. 0.0) .and. (dxp .eq. 0.0)) then
                   ss(n) = tmax
                 elseif ((dxp .eq. 0.0) .and. (dxm .gt. 0.5)) then
                   ss(n) = tmax + (tmax-ss(nm))*a(n)/(2.-a(n))
                 elseif ((dxm .eq. 0.0) .and. (dxp .gt. 0.5)) then
                   ss(n) = tmax + (tmax-ss(np))*c(n)/(2.-c(n))
                 endif
               elseif (ss(n) .lt. tmin) then
                 if (ss(nm) .ge. tmin) then
                   dxm = (ss(n)-tmin)/((ss(n)-ss(nm))*a(n))
                 else
                   dxm = 0.0
                 endif
                 if (ss(np) .ge. tmin) then
                   dxp = (ss(n)-tmin)/((ss(n)-ss(np))*c(n))
                 else
                   dxp = 0.0
                 endif
                 if ((dxm .gt. 0.5) .and. (dxp .gt. 0.5)) then
                   s1 = tmin + (tmin-ss(nm))*a(n)/(2.-a(n))
                   s2 = tmin + (tmin-ss(np))*c(n)/(2.-c(n))
                   ss(n) = amin1(s1, s2)
                 elseif ((dxm .eq. 0.0) .and. (dxp .eq. 0.0)) then
                   ss(n) = tmin
                 elseif ((dxp .eq. 0.0) .and. (dxm .gt. 0.5)) then
                   ss(n) = tmin + (tmin-ss(nm))*a(n)/(2.-a(n))
                 elseif ((dxm .eq. 0.0) .and. (dxp .gt. 0.5)) then
                   ss(n) = tmin + (tmin-ss(np))*c(n)/(2.-c(n))
                 endif
               endif
  130        continue
             go to 105
           endif
         endif
      else
c        treat independent segments
         do 300 j=1,jj
           jend = jbeg(j)
  150      jend = jend + 1
           i1 = mod((jend-1), nmon) + 1
           i2 = mod(jend, nmon) + 1
           if (((obsmean(i1) .le. tmin+0.01*dt) .and.
     &          (obsmean(i2) .le. tmin+0.01*dt)) .or.
     &         ((obsmean(i1) .ge. tmax-0.01*dt) .and.
     &          (obsmean(i2) .ge. tmax-0.01*dt))) then
c            calculate values for interval jbeg(j) to jend
c
c            latest approximation of means (given mid-month values)
             nnn = 0
  205        nnn = nnn + 1
             kk = jend - jbeg(j) + 1
             n = jbeg(j)
             avg(1) = obsmean(n)
             r(1) = 0.0
             n = mod((jend-1), nmon) + 1
             avg(kk) = obsmean(n)
             r(kk) = 0.0
             sum = 0.0
             residmax = 0.0
             do 210 k = 2, kk-1
               nm = mod((k+jbeg(j)-3), nmon) + 1
               n  = mod((k+jbeg(j)-2), nmon) + 1
               np = mod((k+jbeg(j)-1), nmon) + 1
               bb(k) = 0.0
               avg(k) = 0.0
               if (nnn .lt. imethod) then
                 call approx(tmin, tmax, a(n), c(n), ss(nm), ss(n),
     &                ss(np), aa(k), bb(k), cc(k), avg(k))
               else
                 call numer(conv, tmin, tmax, bbmin, a(n), c(n), ss(nm),
     &                ss(n), ss(np), aa(k), bb(k), cc(k), avg(k))
               endif
               r(k) = obsmean(n) - avg(k)
               sum = sum + r(k)**2
               residmax = amax1(residmax, abs(r(k)))
210          continue
             resid = dsqrt(sum)/(kk-2)
             if (residmax .gt. conv) then
               if (nnn .gt. maxiter*0.9) then
                  print*, 'iter = ', nnn, ' kk = ', kk, ' residual = ',
     &               resid, ' maximum residual = ', residmax
c                  print*, ss(nm), ss(n), ss(np)
               endif
               if (nnn .gt. maxiter*0.9) then
                 print*, ' '
                 print*, 'latitude = ', alat, ' longitude = ', alon
                 do 2234 k=1,kk
                   n  = mod((k+jbeg(j)-2), nmon) + 1
                   write(*,'(8(1pe10.2))') obsmean(n), avg(k), r(k),
     &                  s(k), ss(n), aa(k), bb(k), cc(k)
 2234            continue
               endif
               if (nnn .gt. maxiter) then
                 print*, 'latitude = ', alat, ' longitude = ', alon
                 print*, 'does not converge'
                 write(9,*) 'latitude = ', alat, ' longitude = ', alon
                 write(9,*) 'does not converge'
c                 call exit(1)
               else
c                solve for new estimate of mid-month values
                 kkk = kk - 2
                 call tridag(alon, alat, aa(2), bb(2), cc(2), r(2),
     &                 s(2), kkk)
                 do 220 k=2,kk-1
                   n  = mod((k+jbeg(j)-2), nmon) + 1
                   ss(n) = ss(n) + relax*s(k)
  220            continue
c               if ss exceeds tmax or tmin, then it should exceed it no
c                    more than absolutely necessary:
                 n  = mod((jbeg(j)-1), nmon) + 1
                 np  = mod(jbeg(j), nmon) + 1
                 if (obsmean(n) .ge. (tmax-0.01*dt)) then
                   ss(n) =
     &                amax1(tmax, (tmax + (tmax-ss(np))*c(n)/(2.-c(n))))
                 else
                   ss(n) =
     &                amin1(tmin, (tmin + (tmin-ss(np))*c(n)/(2.-c(n))))
                 endif
                 nm  = mod((jend+nmon-2), nmon) + 1
                 n  = mod((jend-1), nmon) + 1
                 if (obsmean(n) .ge. (tmax-0.01*dt)) then
                   ss(n) =
     &                amax1(tmax, (tmax + (tmax-ss(nm))*a(n)/(2.-a(n))))
                 else
                   ss(n) =
     &                amin1(tmin, (tmin + (tmin-ss(nm))*a(n)/(2.-a(n))))
                 endif
                 do 230 k=2,kk-1
                   nm = mod((k+jbeg(j)+nmon-3), nmon) + 1
                   n  = mod((k+jbeg(j)-2), nmon) + 1
                   np = mod((k+jbeg(j)-1), nmon) + 1
                   if (ss(n) .gt. tmax) then
                     if (ss(nm) .le. tmax) then
                       dxm = (ss(n)-tmax)/((ss(n)-ss(nm))*a(n))
                     else
                       dxm = 0.0
                     endif
                     if (ss(np) .le. tmax) then
                       dxp = (ss(n)-tmax)/((ss(n)-ss(np))*c(n))
                     else
                       dxp = 0.0
                     endif
                     if ((dxm .gt. 0.5) .and. (dxp .gt. 0.5)) then
                       s1 = tmax + (tmax-ss(nm))*a(n)/(2.-a(n))
                       s2 = tmax + (tmax-ss(np))*c(n)/(2.-c(n))
                       ss(n) = amin1(s1, s2)
                     elseif ((dxm .eq. 0.0) .and. (dxp .eq. 0.0)) then
                       ss(n) = tmax
                     elseif ((dxp .eq. 0.0) .and. (dxm .gt. 0.5)) then
                       ss(n) = tmax + (tmax-ss(nm))*a(n)/(2.-a(n))
                     elseif ((dxm .eq. 0.0) .and. (dxp .gt. 0.5)) then
                       ss(n) = tmax + (tmax-ss(np))*c(n)/(2.-c(n))
                     endif
                   elseif (ss(n) .lt. tmin) then
                     if (ss(nm) .ge. tmin) then
                       dxm = (ss(n)-tmin)/((ss(n)-ss(nm))*a(n))
                     else
                       dxm = 0.0
                     endif
                     if (ss(np) .ge. tmin) then
                       dxp = (ss(n)-tmin)/((ss(n)-ss(np))*c(n))
                     else
                       dxp = 0.0
                     endif
                     if ((dxm .gt. 0.5) .and. (dxp .gt. 0.5)) then
                       s1 = tmin + (tmin-ss(nm))*a(n)/(2.-a(n))
                       s2 = tmin + (tmin-ss(np))*c(n)/(2.-c(n))
                       ss(n) = amax1(s1, s2)
                     elseif ((dxm .eq. 0.0) .and. (dxp .eq. 0.0)) then
                       ss(n) = tmin
                     elseif ((dxp .eq. 0.0) .and. (dxm .gt. 0.5)) then
                       ss(n) = tmin + (tmin-ss(nm))*a(n)/(2.-a(n))
                     elseif ((dxm .eq. 0.0) .and. (dxp .gt. 0.5)) then
                       ss(n) = tmin + (tmin-ss(np))*c(n)/(2.-c(n))
                     endif
                   endif
  230            continue
                 go to 205
               endif
             endif
             go to 300
           else
             go to 150
           endif
           if ((ss(kk) .gt. 700.0) .or. (ss(kk) .lt. -700.0)) then
             print*, 'exceeds 700 at lat = ', alat,
     &              ' alon = ', alon
             print*, 'kk = ', kk, ' obs = ', obsmean(kk-1), obsmean(kk),
     &               obsmean(kk+1)
             print*, 'kk = ', kk, '  ss = ', ss(kk-1), ss(kk),
     &              ss(kk+1)
           endif
  300    continue
c        fill in values where consecutive means are outside limits
         do 250 i=1,nmon
           i1 = mod((i-2+nmon), nmon) + 1
           i2 = mod((i-1), nmon) + 1
           i3 = mod(i, nmon) + 1
           if ((obsmean(i1) .le. (tmin+0.01*dt)) .and.
     &         (obsmean(i2) .le. (tmin+0.01*dt)) .and.
     &         (obsmean(i3) .le. (tmin+0.01*dt))) then
             ss(i2) = tmin
           elseif ((obsmean(i1) .ge. (tmax-0.01*dt)) .and.
     &             (obsmean(i2) .ge. (tmax-0.01*dt)) .and.
     &             (obsmean(i3) .ge. (tmax-0.01*dt))) then
             ss(i2) = tmax
           endif
  250    continue
      endif
           return
      end

      subroutine numer(conv, tmin, tmax, bbmin, a, c, ssm, ss, ssp, aa,
     &                 bb, cc, avg)
c *********************************************************************
      implicit none
      real conv, tmin, tmax, bbmin, a, c, ssm, ss, ssp, aa, bb, cc, avg
      real ssmm, ssmp, sssm, sssp, sspm, sspp, r
      real amean
      avg = amean(tmin,tmax,a,c,ssm,ss,ssp)
      ssmm = ssm - conv
      ssmp = ssm + conv
      sssm = ss  - conv
      sssp = ss  + conv
      sspm = ssp - conv
      sspp = ssp + conv
      aa = (amean(tmin,tmax,a,c,ssmp,ss,ssp) -
     &      amean(tmin,tmax,a,c,ssmm,ss,ssp)) / (2.*conv)
      bb = (amean(tmin,tmax,a,c,ssm,sssp,ssp) -
     &      amean(tmin,tmax,a,c,ssm,sssm,ssp)) / (2.*conv)
      cc = (amean(tmin,tmax,a,c,ssm,ss,sspp) -
     &      amean(tmin,tmax,a,c,ssm,ss,sspm)) / (2.*conv)
      aa = amin1(aa, bb)
      cc = amin1(cc, bb)
C     the following ensure that the diagonal elements of the Jacobian dominate
      if (bb .lt. bbmin) then
        bb=bbmin
        aa = amin1(aa, bb)
        cc = amin1(cc, bb)
      endif
C       if (bb .lt. bbmin) then
C         bb = bbmin
C         r = 0.2*bbmin
C         aa = amax1(r, aa)
C         cc = amax1(r, cc)
C       endif
      return
      end

      subroutine approx(tmin, tmax, a, c, ssm, ss, ssp, aa, bb, cc, avg)
c *********************************************************************
      implicit none
      real tmin, tmax, a, c, ssm, ss, ssp, aa, bb, cc, avg
      real dx, dy
      if (ss .le.  tmin) then
        if (ssm .le. tmin) then
          avg = avg + tmin*0.5
          aa = a/32.
          bb = bb + 0.125 - a/32.
        elseif (ssm .ge. tmax) then
          dx = (ss-tmin)/((ss-ssm)*a)
          dy = (ss-tmax)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + tmin*0.5
            aa = a/32.
            bb = bb + 0.125 - a/32.
          elseif (dy .le. 0.5) then
            avg = avg + tmin*dx + tmax*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
            aa = a/16.
            bb = bb + 0.25 - a/16.
          else
            avg = avg +
     &              tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*a*(ss-ssm))
            aa = a/16.
            bb = bb + 0.25 - a/16.
          endif
        else
          dx = (ss-tmin)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + tmin*0.5
            aa = a/32.
            bb = bb + 0.125 - a/32.
          else
            avg = avg +
     &              tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*a*(ss-ssm))
            aa = a/16.
            bb = bb + 0.25 - a/16.
          endif
        endif
        if (ssp .le. tmin) then
          avg = avg + tmin*0.5
          cc = c/32.
          bb = bb + 0.125 - c/32.
        elseif (ssp .ge. tmax) then
          dx = (ss-tmin)/((ss-ssp)*c)
          dy = (ss-tmax)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + tmin*0.5
            cc = c/32.
            bb = bb + 0.125 - c/32.
          elseif (dy .le. 0.5) then
            avg = avg + tmin*dx + tmax*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
            cc = c/16.
            bb = bb + 0.25 - c/16.
          else
            avg = avg +
     &              tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*c*(ss-ssp))
            cc = c/16.
            bb = bb + 0.25 - c/16.
          endif
        else
          dx = (ss-tmin)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + tmin*0.5
            cc = c/32.
            bb = bb + 0.125 - c/32.
          else
            avg = avg +
     &              tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*c*(ss-ssp))
            cc = c/16.
            bb = bb + 0.25 - c/16.
          endif
        endif
      elseif (ss .ge. tmax) then
        if (ssm .ge. tmax) then
          avg = avg + tmax*0.5
          aa = a/32.
          bb = bb + 0.125 - a/32.
        elseif (ssm .le. tmin) then
          dx = (ss-tmax)/((ss-ssm)*a)
          dy = (ss-tmin)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + tmax*0.5
            aa = a/32.
            bb = bb + 0.125 - a/32.
          elseif (dy .le. 0.5) then
            avg = avg + tmax*dx + tmin*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
            aa = a/16.
            bb = bb + 0.25 - a/16.
          else
            avg = avg +
     &              tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*a*(ss-ssm))
            aa = a/16.
            bb = bb + 0.25 - a/16.
          endif
        else
          dx = (ss-tmax)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + tmax*0.5
            aa = a/32.
            bb = bb + 0.125 - a/32.
          else
            avg = avg +
     &        tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*a*(ss-ssm))
            aa = a/16.
            bb = bb + 0.25 - a/16.
          endif
        endif
        if (ssp .ge. tmax) then
          avg = avg + tmax*0.5
          cc = c/32.
          bb = bb + 0.125 - c/32.
        elseif (ssp .le. tmin) then
          dx = (ss-tmax)/((ss-ssp)*c)
          dy = (ss-tmin)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            cc = c/32.
            bb = bb + 0.125 - c/32.
            avg = avg + tmax*0.5
          elseif (dy .le. 0.5) then
            avg = avg + tmax*dx + tmin*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
            cc = c/16.
            bb = bb + 0.25 - c/16.
          else
            avg = avg +
     &              tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*c*(ss-ssp))
            cc = c/16.
            bb = bb + 0.25 - c/16.
          endif
        else
          dx = (ss-tmax)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + tmax*0.5
            cc = c/32.
            bb = bb + 0.125 - c/32.
          else
            avg = avg +
     &              tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*c*(ss-ssp))
            cc = c/16.
            bb = bb + 0.25 - c/16.
          endif
        endif
      else
        if (ssm .le. tmin) then
          dx = (ss-tmin)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssm)*a)
            aa = a/16.
            bb = bb + 0.25 - a/16.
          else
            avg = avg + tmin*(.5-dx) + dx*0.5*(tmin + ss)
            aa = a/16.
            bb = bb + 0.25 - a/16.
          endif
        elseif (ssm .ge. tmax) then
          dx = (ss-tmax)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssm)*a)
            aa = a/16.
            bb = bb + 0.25 - a/16.
          else
            avg = avg + tmax*(.5-dx) + dx*0.5*(tmax + ss)
            aa = a/16.
            bb = bb + 0.25 - a/16.
          endif
        else
          avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssm)*a)
          aa = a/8.
          bb = bb + 0.5 - a/8.
        endif
        if (ssp .le. tmin) then
          dx = (ss-tmin)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssp)*c)
            cc = c/16.
            bb = bb + 0.25 - c/16.
          else
            avg = avg + tmin*(.5-dx) + dx*0.5*(tmin + ss)
            cc = c/16.
            bb = bb + 0.25 - c/16.
          endif
        elseif (ssp .ge. tmax) then
          dx = (ss-tmax)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssp)*c)
            cc = c/16.
            bb = bb + 0.25 - c/16.
          else
            avg = avg + tmax*(.5-dx) + dx*0.5*(tmax + ss)
            cc = c/16.
            bb = bb + 0.25 - c/16.
          endif
        else
          avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssp)*c)
          cc = c/8.
          bb = bb + 0.5 - c/8.
        endif
      endif
      return
      end
      SUBROUTINE tridag(alon,alat,a,b,c,r,u,n)
      implicit none
      INTEGER n,nmax
      REAL alon,alat,a(n),b(n),c(n),r(n)
      double precision u(n)
      PARAMETER (nmax=500*12) ! PJD Oceanonly 1870-2014 - 171025
c      PARAMETER (nmax=149*12) ! PJD Oceanonly 1870-2014 - 170412
c      PARAMETER (nmax=149*12) ! PJD Oceanonly 1870-2014 - 160414
c      PARAMETER (nmax=147*12) ! PJD Oceanonly 1870-2014
      INTEGER j
      REAL bet, gam(nmax)
      if (nmax .lt. n) then
           print*, 'Error nmax not declared large enough'
           print*, 'in tridag'
           stop
      endif
      if(b(1).eq.0.) then
          print*, 'longitude = ', alon, '  latitude = ', alat
c          pause 'tridag: rewrite equations'
      endif
      bet=b(1)
      u(1)=r(1)/bet
      if (n .gt. 1) then
        do 11 j=2,n
          gam(j)=c(j-1)/bet
          bet=b(j)-a(j)*gam(j)
          if(bet.eq.0.) then
            print*, 'longitude = ', alon, '  latitude = ', alat
c            pause 'tridag failed'
          endif
          u(j)=(r(j)-a(j)*u(j-1))/bet
11      continue
        do 12 j=n-1,1,-1
          u(j)=u(j)-gam(j+1)*u(j+1)
12      continue
      endif
      return
      END

      SUBROUTINE cyclic(alon,alat,a,b,c,alpha,beta,r,x,n)
      implicit none
      INTEGER n,nmax
      real alon,alat,alpha,beta,a(n),b(n),c(n),r(n)
      double precision x(n)
      PARAMETER (nmax=500*12) ! PJD Oceanonly 1870-2014 - 171025
c      PARAMETER (nmax=149*12) ! PJD Oceanonly 1870-2014 - 170412
c      PARAMETER (nmax=149*12) ! PJD Oceanonly 1870-2014 - 160414
c      PARAMETER (nmax=147*12) ! PJD Oceanonly 1870-2014
CU    USES tridag
      INTEGER i
      REAL fact,gamma,bb(nmax),u(nmax)
      double precision z(nmax)
c      if(n.le.2)pause 'n too small in cyclic'
c      if(n.gt.nmax)pause 'nmax too small in cyclic'
      gamma=-b(1)
      bb(1)=b(1)-gamma
      bb(n)=b(n)-alpha*beta/gamma
      do 11 i=2,n-1
        bb(i)=b(i)
11    continue
      call tridag(alon,alat,a,bb,c,r,x,n)
      u(1)=gamma
      u(n)=alpha
      do 12 i=2,n-1
        u(i)=0.
12    continue
      call tridag(alon,alat,a,bb,c,u,z,n)
      fact=(x(1)+beta*x(n)/gamma)/(1.+z(1)+beta*z(n)/gamma)
      do 13 i=1,n
        x(i)=x(i)-fact*z(i)
13    continue
      return
      END


      function amean(tmin, tmax, a, c, ssm, ss, ssp)
c *********************************************************************
      implicit none
      real amean
      real tmin, tmax, a, c, ssm, ss, ssp
      real dx, dy, avg
      avg = 0.0
      if (ss .le.  tmin) then
        if (ssm .le. tmin) then
          avg = avg + tmin*0.5
        elseif (ssm .ge. tmax) then
          dx = (ss-tmin)/((ss-ssm)*a)
          dy = (ss-tmax)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + tmin*0.5
          elseif (dy .le. 0.5) then
            avg = avg + tmin*dx + tmax*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
          else
            avg = avg +
     &              tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*a*(ss-ssm))
          endif
        else
          dx = (ss-tmin)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + tmin*0.5
          else
            avg = avg +
     &              tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*a*(ss-ssm))
          endif
        endif
        if (ssp .le. tmin) then
          avg = avg + tmin*0.5
        elseif (ssp .ge. tmax) then
          dx = (ss-tmin)/((ss-ssp)*c)
          dy = (ss-tmax)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + tmin*0.5
          elseif (dy .le. 0.5) then
            avg = avg + tmin*dx + tmax*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
          else
            avg = avg +
     &              tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*c*(ss-ssp))
          endif
        else
          dx = (ss-tmin)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + tmin*0.5
          else
            avg = avg +
     &              tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*c*(ss-ssp))
          endif
        endif
      elseif (ss .ge. tmax) then
        if (ssm .ge. tmax) then
          avg = avg + tmax*0.5
        elseif (ssm .le. tmin) then
          dx = (ss-tmax)/((ss-ssm)*a)
          dy = (ss-tmin)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + tmax*0.5
          elseif (dy .le. 0.5) then
            avg = avg + tmax*dx + tmin*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
          else
            avg = avg +
     &              tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*a*(ss-ssm))
          endif
        else
          dx = (ss-tmax)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + tmax*0.5
          else
            avg = avg +
     &        tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*a*(ss-ssm))
          endif
        endif
        if (ssp .ge. tmax) then
          avg = avg + tmax*0.5
        elseif (ssp .le. tmin) then
          dx = (ss-tmax)/((ss-ssp)*c)
          dy = (ss-tmin)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + tmax*0.5
          elseif (dy .le. 0.5) then
            avg = avg + tmax*dx + tmin*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
          else
            avg = avg +
     &              tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*c*(ss-ssp))
          endif
        else
          dx = (ss-tmax)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + tmax*0.5
          else
            avg = avg +
     &              tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*c*(ss-ssp))
          endif
        endif
      else
        if (ssm .le. tmin) then
          dx = (ss-tmin)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssm)*a)
          else
            avg = avg + tmin*(.5-dx) + dx*0.5*(tmin + ss)
          endif
        elseif (ssm .ge. tmax) then
          dx = (ss-tmax)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssm)*a)
          else
            avg = avg + tmax*(.5-dx) + dx*0.5*(tmax + ss)
          endif
        else
          avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssm)*a)
        endif
        if (ssp .le. tmin) then
          dx = (ss-tmin)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssp)*c)
          else
            avg = avg + tmin*(.5-dx) + dx*0.5*(tmin + ss)
          endif
        elseif (ssp .ge. tmax) then
          dx = (ss-tmax)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssp)*c)
          else
            avg = avg + tmax*(.5-dx) + dx*0.5*(tmax + ss)
          endif
        else
          avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssp)*c)
        endif
      endif
      amean = avg
      return
      end                  