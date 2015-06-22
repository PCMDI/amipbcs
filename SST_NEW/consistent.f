C WRAPIT -pg consistent_coast.f

C NCLFORTSTART
      subroutine ssticejh (nlat, mlon, ice, sst, zmsg )
      implicit      none
      integer       nlat, mlon
      real          ice(mlon,nlat), sst(mlon,nlat), zmsg
C NCLEND

c  NOTES:  
c   - land values for sst do not necessarily coincide with land values 
c       from ice analysis

      integer    id, jd, krecl
      parameter (id=360,jd=180)
      parameter (krecl=id*jd*4)

      real    sstnew(id,jd),icenew(id,jd)
      real*8  ix 
      real    sstmax(91),icemax(91)
      real    si, sstm

      integer i,j,ic

c *************************************************************
c  This section modifies sea ice data to be consistent with SST: 
c  [Jim Hurrell    sst-ice.consistency.f]
c *************************************************************
c
c  create a maximum SST allowed for a particular sea-ice concentration
c
c  Function created by Jim Hack 2/4/02
c
c  sstmax = 9.328 * (0.729-ice**3) - 1.8
c
c  icemax(1)  = 0.00   sstmax(1) = 5.0
c  icemax(16) = 0.15   sstmax(16) = 4.97    ! si cutoff in orginal data
c  icemax(90) = 0.89   sstmax(90) = -1.57
c  icemax(91) = 0.90   sstmax(91) = -1.8

      ic = 0
      do ix = 0.0, 0.90, 0.01
        ic = ic + 1
        icemax(ic) = ix
        sstmax(ic) = 9.328 * (0.729-ix**3) - 1.8
      enddo

      do j = 1, jd
        do i = 1,id 
           sstnew(i,j) = sst(i,j)
           icenew(i,j) = ice(i,j)
        end do
      end do

      do j = 1, jd
        do i = 1,id 

c (a) first, do not allow sst < -1.8 
c (b) set all values with ice > 90% c to -1.8
c
c THE FOLLOWING IS DONE IN NCL  ... 2006
c
           if (sstnew(i,j).ne.zmsg .and. 
     +         sstnew(i,j).lt.-1.8) sstnew(i,j) = -1.8

           if (ice(i,j).ne.zmsg.and.sstnew(i,j).ne.zmsg) then
             if (ice(i,j).ge.90.0) then
                 sstnew(i,j) = -1.8
             endif
           endif

c adjust the ice concentration data

           if (sstnew(i,j).ne.zmsg) then
               if (icenew(i,j).ne.zmsg.and.icenew(i,j).gt.0.0) then

                if (icenew(i,j).lt.90.0) then   ! Don't adjust values > 90%
                  si   = icenew(i,j) * 0.01     ! Convert to fraction from %
                  sstm = 9.328 * (0.729-si**3) - 1.8
                  if (sstnew(i,j).gt.sstm) then
                    if (sstnew(i,j).gt.sstmax(1)) then
                      icenew(i,j) = 0.0
                    else
                      do ic = 1, 90
                       if (sstnew(i,j).lt.sstmax(ic).and.sstnew(i,j).gt.
     *                     sstmax(ic+1)) then
                            icenew(i,j) = icemax(ic) * 100.
                       endif
                      enddo
                    endif
                  endif
                endif

               endif
              endif

c             force no sea ice < 15%, as in HadISST data  
c             DJS: UNCOMMENTED 16 Aug 2006

              if (icenew(i,j).ne.zmsg) then
                if (icenew(i,j).lt.15.0) icenew(i,j) = 0.0
C DJS           if (icenew(i,j).gt.99.9) then
C DJS             write (*,*) iyr,imn,i,j,ice(i,j)
C DJS           endif
              endif
        enddo      
      enddo         

c DJS: PUT NEW VALUES BACK TO SST/ICE FOR RETURN TO NCL

      DO J = 1, JD
        DO I = 1,ID 
           SST(I,J) = SSTNEW(I,J)
           ICE(I,J) = ICENEW(I,J)
        END DO
      END DO

c DJS: THE FOLLOWING IS SOMETHING LIKE WHAT JIM HAD IN hadissst+oiv2.f

c djs DO J = 1, JD
c djs   DO I = 1, ID
c djs     IF (ICE(I,J).NE.ZMSG) THEN
c djs       IF (ICE(I,J).LT.0.0.OR.ICE(I,J).GT.100.0) THEN
c djs           WRITE (*,*) 'OIv2 ',ICE(I,J),I,J
c djs       ENDIF
c djs     ENDIF
c djs   ENDDO
c djs ENDDO

      return
      end

