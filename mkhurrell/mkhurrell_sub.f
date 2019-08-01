c    This code was written by Karl Taylor
c    lightly edited by Stephen Po-Chedley (increase nmont, add intent to ss, jcnt/icnt)
c    I'm hoping he'll help me add some documentation here. 

c      icnt (out) = number of jumps from max to min or min to max
c      niter (out) = number of iterations required to converge (summed over all segments)
c      notconverg (out) = # of segments where convergence cannot be achieved
c      jj = number of isolated segments (separated by times when at least two 
c             consecutive values are very near either the max or min allowed).
c      jumps(nmon) (out) time (index) when obs values jump from one limit to the other
c      jcnt (inout) = accumulator counting number of cells where at least one jump occurs
c                used in this subroutine to determine if diagnostics should be 
c                printed).  Python code should set to a negative integer before
c                first calling solvmid
c      obsmean may be modified by this subroutine (if the time series contains 
c              consecutive values that jump from near minimum to near maximum
c              values allowed, or vice versa)
c       
      


      subroutine solvmid(alon, alat, nmon, conv, dt, tmin, tmax,
     &     bbmin, maxiter, a, c, obsmean, ss, icnt, niter, notconverg, 
     &     jj, resid, residmax, jumps, jcnt)
      implicit none
      integer nmon12
      parameter (nmon12=12)
c      parameter (nmont=500*12, nmon12=12) ! PJD Oceanonly 1870-2014 - 171025
c      parameter (nmont=149*12, nmon12=12) ! PJD Oceanonly 1870-2014 - 170412
c      parameter (nmont=149*12, nmon12=12) ! PJD Oceanonly 1870-2014 - 160414
c      parameter (nmont=147*12, nmon12=12) ! PJD Oceanonly 1870-2014
      integer nmon, maxiter
      integer, intent(inout) :: jcnt
      integer, intent(out) :: icnt, niter, notconverg, jj
      integer, intent(out) :: jumps(nmon)
      real conv, tmin, tmax, dt, bbmin, alon, alat
      real obsmean(nmon), a(nmon), c(nmon)
      real, intent(out) :: ss(nmon)
      real, intent(out) :: resid, residmax
      integer i, n, imethod, n1, n2, i1, i2, i3, nnn, jend,
     &      kk, j, k, kkk, nm, np
      integer jbeg(nmon)
      real relax, residmax1, resid1, dxm, dxp, s1, s2, addmax, addmin
      real r(nmon), avg(nmon), aa(nmon), bb(nmon), cc(nmon),
     &     add(nmon)
      double precision s(nmon), sum
      imethod = 1
c ???  check following value
      relax = 1.0
c      if (nmon .gt. nmont) then
c        print*, 'error-- nmont not declared large enough in '
c        print*, 'subroutine solvmid'
c        stop
c      endif
c
c    check for occurance where obs monthly means are consecutively
c      at upper and lower limits. If so, smooth data, being careful
c      to preserve annual mean.
c    also initialize ss = obsmean
c
      niter = 0
      notconverg = 0
      do 40, n=1,nmon
        jumps(n) = 0
        aa(n) = 1.e20
        bb(n) = 1.e20
        cc(n) = 1.e20
  40  continue
      resid = 0.
      residmax = 0.
c
      if (jcnt .lt. 0) then
        jcnt = 0
        write(9,'("       lat    lon   #jumps   #iter  failed    ",
     &    "segs    resid  residmax     jump times")')
      endif
c
      do 48 n=1,nmon
        ss(n)=obsmean(n)
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
      icnt = 0
      addmax = 0.0
      addmin = 0.0
      do 51 n=1,nmon
        if (add(n) .ne. 0.0) then
            jumps(icnt+1) = n
c           if (jcnt .lt. 5000) then
c              print*, 'jump from extreme to extreme at time ', n
c              write(9,*) 'jump from extreme to extreme at time ', n
c           endif
c           if (jcnt .eq. 5000) then
c              print*, ' **********************************'
c              print*, ' No more jump messages will be written'
c              print*, ' **********************************'
c           endif
c           if (jcnt .lt. 5) then
c          print*,
c     &     'monthly means go from one limit to the other in 1 month'
c          print*, 'alat = ', alat, ' alon = ', alon, ' n = ', n,
c     &         ' add = ', add(n)
c             n1 = mod((n+nmon-2), nmon) + 1
c             n2 = mod(n, nmon) + 1
c          print*, 'observed mean for 3 cells: ', obsmean(n1),
c     &          obsmean(n), obsmean(n2)
c           endif
c          if (jcnt .le. 2)
c     &        print*, (obsmean(n1), n1=1,nmon)
          addmax = amax1(addmax, add(n))
          addmin = amin1(addmin, add(n))
          obsmean(n) = obsmean(n) + add(n)
          icnt = icnt + 1
        endif
   51 continue
c
      if (icnt .gt. 0) then
        jcnt = jcnt + 1
c        if (nmon .eq. nmon12) then
c          print*, 'Climatology: '
c          write(9,*) 'Climatology: '
c        endif
c        if (jcnt .le. 1000) then
c          print*,  icnt,
c     &      ' monthly values smoothed at lat,lon', alat, alon
c          print*, 'max added = ', addmax,
c     &             '  max subtracted = ', addmin
c          write(9,*)  icnt,
c     &      ' monthly values smoothed at lat,lon', alat, alon
c          write(9,*) 'max added = ', addmax,
c     &             '  max subtracted = ', addmin
c        endif
c        if (jcnt .eq. 50) then
c          print*, ' '
c          if (nmon .eq. nmon12) then
c            print*,
c     &      'No more warnings will be printed concerning smoothing '//
c     &       'of climatological data'
c          else
c            print*,
c     &      'No more warnings will be printed concerning smoothing '//
c     &       'of monthly data'
c          endif
c          print*, ' '
c          write(9,*) ' '
c          if (nmon .eq. nmon12) then
c            write(9,*)
c     &      'No more warnings will be printed concerning smoothing '//
c     &       'of climatological data'
c          else
c            write(9,*)
c     &      'No more warnings will be printed concerning smoothing '//
c     &       'of monthly data'
c          endif
c          write(9,*) ' '
c        endif
      endif
c    check if all are le tmin or all are ge tmax
      if (obsmean(1) .le. (tmin+0.01*dt)) then
        do 80 i=2,nmon
          if (obsmean(i) .gt. (tmin+0.01*dt)) go to 99
c      All values are at tmin.          
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
c     Find time intervals when there are two consecutive values 
c            less than minimum followed by a value greater than maximum 
c            OR intervals when there are two consecutive values 
c            greater than maximum followed by a value less than maximum
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
c         jj is counter for number of intervals satisfying criteria
          jj = jj + 1
c         jbeg is is pointer set to time marking beginning of an interval
c              that can be treated independent of some portion of the
c              entire time-series.
          jbeg(jj) = i2
        endif
  100 continue
      if (jj .eq. 0) then
c       simple cyclic treatment
c         latest approximation of means (given mid-month values)
         nnn = 0
  105    nnn = nnn + 1
         niter = niter + 1
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
           sum = sum + abs(r(n))
           residmax = amax1(residmax, abs(r(n)))
110      continue
         resid = sum/nmon
         if (residmax .gt. conv) then
           if (nnn .gt. maxiter*0.95) then
             print*, 'iteration = ', nnn, ' residual = ', resid,
     &          ' maximum residual = ', residmax
           endif
           if (nnn .gt. maxiter*0.99) then
c             print*, ' '
c             print*, 'latitude = ', alat, ' longitude = ', alon
             do 1234 n=1,nmon
               write(*,'(i5, 8(1pe10.2))') n, obsmean(n), avg(n), r(n),
     &              s(n), ss(n), aa(n), bb(n), cc(n)
 1234        continue
           endif
           if (nnn .gt. maxiter) then
             notconverg = notconverg + 1
c             print*, 'latitude = ', alat, ' longitude = ', alon
c             print*, 'does not converge'
c             write(9,*) 'latitude = ', alat, ' longitude = ', alon
c             write(9,*) 'does not converge'
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
c             find end of interval that is independent of part of the
c                full time series.
           if (((obsmean(i1) .le. tmin+0.01*dt) .and.
     &          (obsmean(i2) .le. tmin+0.01*dt)) .or.
     &         ((obsmean(i1) .ge. tmax-0.01*dt) .and.
     &          (obsmean(i2) .ge. tmax-0.01*dt))) goto 204
            goto 150
c            calculate values for interval jbeg(j) to jend
c
c            latest approximation of means (given mid-month values)
  204        nnn = 0
  205        nnn = nnn + 1
             niter = niter + 1
             kk = jend - jbeg(j) + 1
             n = jbeg(j)
             avg(1) = obsmean(n)
             r(1) = 0.0
             n = mod((jend-1), nmon) + 1
             avg(kk) = obsmean(n)
             r(kk) = 0.0
             sum = 0.0
             residmax1 = 0.0
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
               sum = sum + abs(r(k))
               residmax1 = amax1(residmax1, abs(r(k)))
210          continue
             resid1 = sum/(kk-2)
             if (residmax1 .gt. conv) then
c               if (nnn .gt. maxiter*0.9) then
c                  print*, 'iter = ', nnn, ' kk = ', kk, ' residual = ',
c     &               resid1, ' maximum residual = ', residmax1
c                  print*, ss(nm), ss(n), ss(np)
c               endif
               if ((nnn .gt. maxiter*0.99) .and. (residmax1 .gt. 30.))
     &           then
c                 print*, ' '
c                 print*, 'latitude = ', alat, ' longitude = ', alon
                 write(*,'(f8.1, i5, 2f8.1)') residmax1, nnn, alat, alon
                 do 2234 k=1,kk
                   n  = mod((k+jbeg(j)-2), nmon) + 1
                   write(*,'(i5, 8(1pe10.2), i5, i5)') n, obsmean(n), 
     &                  avg(k), r(k), s(k), ss(n), aa(k),
     &                  bb(k), cc(k), j, jj
 2234            continue
               endif
               if (nnn .gt. maxiter) then
                 notconverg = notconverg + 1
c                 print*, 'latitude = ', alat, ' longitude = ', alon
c                 print*, 'does not converge'
c                 write(9,*) 'latitude = ', alat, ' longitude = ', alon
c                 write(9,*) 'does not converge'
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
c
c                treat first sample of segment:
c
                 n  = mod((jbeg(j)-1), nmon) + 1
                 np  = mod(jbeg(j), nmon) + 1
                 if (obsmean(n) .ge. (tmax-0.01*dt)) then
                   ss(n) =
     &                amax1(tmax, (tmax + (tmax-ss(np))*c(n)/(2.-c(n))))
                 else
                   ss(n) =
     &                amin1(tmin, (tmin + (tmin-ss(np))*c(n)/(2.-c(n))))
                 endif
c
c                treat last sample of segment:
c                 
                 nm  = mod((jend+nmon-2), nmon) + 1
                 n  = mod((jend-1), nmon) + 1
                 if (obsmean(n) .ge. (tmax-0.01*dt)) then
                   ss(n) =
     &                amax1(tmax, (tmax + (tmax-ss(nm))*a(n)/(2.-a(n))))
                 else
                   ss(n) =
     &                amin1(tmin, (tmin + (tmin-ss(nm))*a(n)/(2.-a(n))))
                 endif
c
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
c               need another iteration
                 go to 205
               endif
             endif
c              calculation has converged
c
c             go to 300
c           else
c             go to 150
c           endif
c
            resid = resid + sum
            residmax = amax1(residmax, residmax1)
c          finished loop over independent segments.
  300    continue
c
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
c
        resid = resid/nmon
c
c     end of if/else distinguishing between cyclic case and 
c         independent segments case.
      endif
c
      if (notconverg .gt. 0) then
        write(9, '("***", f7.1, f7.1, i8, i8, i8, i8, 1pe10.2, 
     &   1pe10.2, 3i6)')  
     &       alat, alon, icnt, niter, notconverg, jj, resid, residmax,
     &         jumps(1), jumps(2), jumps(3)
      elseif (icnt .gt. 0) then
        write(9, '("   ", f7.1, f7.1, i8, i8, i8, i8, 1pe10.2, 
     &   1pe10.2 )')  
     &       alat, alon, icnt, niter, notconverg, jj, resid, residmax,
     &         jumps(1), jumps(2), jumps(3)
      endif
      if (icnt .gt. 0) then
        print*, 'icnt= ',icnt,'  jumps at times' ,(jumps(n), n=1,icnt)
      endif
c
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
C      if (bb .lt. bbmin) then
C        bb=bbmin
C        aa = amin1(aa, bb)
C        cc = amin1(cc, bb)
C      endif
      if (bb .lt. bbmin) then
        bb = bbmin
        r = 0.2*bbmin
        aa = amax1(r, aa)
        cc = amax1(r, cc)
      endif
c      if (bb .lt. bbmin) bb = bbmin
c      r = 2*bb
c      aa = amin1(aa, r)
c      cc = amin1(cc, r)
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
      INTEGER n
      REAL alon,alat,a(n),b(n),c(n),r(n)
      double precision u(n)
c      PARAMETER (nmax=500*12) ! PJD Oceanonly 1870-2014 - 171025
c      PARAMETER (nmax=149*12) ! PJD Oceanonly 1870-2014 - 170412
c      PARAMETER (nmax=149*12) ! PJD Oceanonly 1870-2014 - 160414
c      PARAMETER (nmax=147*12) ! PJD Oceanonly 1870-2014
      INTEGER j
      REAL bet, gam(n)
c      if (nmax .lt. n) then
c           print*, 'Error nmax not declared large enough'
c           print*, 'in tridag'
c           stop
c      endif
      if(b(1).eq.0.) then
c          print*, 'longitude = ', alon, '  latitude = ', alat
c          pause 'tridag: rewrite equations'
      endif
      bet=b(1)
      u(1)=r(1)/bet
      if (n .gt. 1) then
        do 11 j=2,n
          gam(j)=c(j-1)/bet
          bet=b(j)-a(j)*gam(j)
          if(bet.eq.0.) then
c            print*, 'longitude = ', alon, '  latitude = ', alat
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
      INTEGER n
      real alon,alat,alpha,beta,a(n),b(n),c(n),r(n)
      double precision x(n)
c      PARAMETER (nmax=500*12) ! PJD Oceanonly 1870-2014 - 171025
c      PARAMETER (nmax=149*12) ! PJD Oceanonly 1870-2014 - 170412
c      PARAMETER (nmax=149*12) ! PJD Oceanonly 1870-2014 - 160414
c      PARAMETER (nmax=147*12) ! PJD Oceanonly 1870-2014
CU    USES tridag
      INTEGER i
      REAL fact,gamma,bb(n),u(n)
      double precision z(n)
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
