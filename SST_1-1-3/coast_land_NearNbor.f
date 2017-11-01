C NCLFORTSTART     
      subroutine coastnn (s,mlon,nlat,ntim,smsg,npassx,npassy)
      implicit none
      integer  mlon, nlat, ntim, npassx, npassy
      real     s(mlon,nlat,ntim), smsg
C NCLEND
      real     stmp(mlon,nlat)
      integer  ml, nl, nt, np
c c c integer  kmsg, new, nwe

c physically the idea is that SSTs don't 
c .   vary much in the east-west direction over small distances
c all this does is set the 1st coastal point to the 
c .   nearest left/right neighbor

      do nt=1,ntim
        do np=1,npassx

          do nl=1,nlat
            do ml=1,mlon
               stmp(ml,nl) = s(ml,nl,nt)
            end do
          end do

          do nl=1,nlat
c                                ! west-to-east [Gulf, Kurishio]
            do ml=mlon,2,-1 
               if (stmp(ml  ,nl).ne.smsg .and.stmp(ml-1,nl).eq.smsg)then
                   s(ml-1,nl,nt) = s(ml,nl,nt) 
               end if
            end do
c                                ! east-to-west
            do ml=1,mlon-1
               if (stmp(ml,nl)  .ne.smsg .and.stmp(ml+1,nl).eq.smsg)then
                   s(ml+1,nl,nt) = s(ml,nl,nt) 
               end if
            end do
          end do

        end do

        do np=1,npassy

          do nl=1,nlat
            do ml=1,mlon
               stmp(ml,nl) = s(ml,nl,nt)
            end do
          end do

          do ml=1,mlon
c                                ! south-to-north
            do nl=1,nlat-1
               if (stmp(ml  ,nl).ne.smsg .and.stmp(ml,nl+1).eq.smsg)then
                   s(ml,nl+1,nt) = s(ml,nl,nt)
               end if   
            end do
c                                ! north-to-south
            do nl=nlat,2,-1
               if (stmp(ml  ,nl).ne.smsg .and.stmp(ml,nl-1).eq.smsg)then
                   s(ml,nl-1,nt) = s(ml,nl,nt)
               end if   
            end do
             
             
          end do
        end do
      end do

      return
      end
