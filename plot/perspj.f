c$Id: perspj.f,v 1.1 2006/11/20 20:33:12 rlt Exp $
      subroutine perspj(xp,ntyp,numnp,errv)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute perspective projection of coordinates

c      Inputs:
c         xp(3,*) - Global coordinates
c         ntyp(*) - Active node indicator
c         numnp   - Total number of nodal points

c      Outputs:
c         xp(3,*) -Perspective projected coordinates
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'ppers.h'
      include  'pview.h'

      logical   errv
      integer   numnp,i,n, ntyp(*)
      real*8    t1(3),xp(3,numnp), alpha

      save

c     Loop over data points and find projection

      do n=1,numnp
        if(ntyp(n).ge.0) then
          lview = .false.
          do i=1,3
            t1(i) = xp(i,n) - e(i)
          end do ! i

          do i=1,3
            xp(i,n) = xlbda(i,1)*t1(1) + xlbda(i,2)*t1(2)
     &              + xlbda(i,3)*t1(3)
          end do ! i

          alpha = -enorm/(t1(1)*q(1,3)+t1(2)*q(2,3)+t1(3)*q(3,3))
          if(alpha.lt.0.d0) then
            errv = .true.
            if(ior.lt.0) then
              write(*,2000)
              return
            else
              write(iow,2000)
              call plstop
            endif
          elseif(xp(3,n).gt.zview) then
            lview = .true.
            errv  = .false.
          else
            errv  = .false.
          endif
          do i=1,3
            xp(i,n) = alpha * xp(i,n)
          end do ! i
        endif
      end do ! n

2000  format(//1x,' Point too close, choose another one!'//)

      end
