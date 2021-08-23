c$Id: prot01.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
      subroutine prot01 (tn ,ta ,tl ,ds ,vn ,v1 ,an ,a1 ,dm ,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Rotational Transformation Updates: (Delta-T)
c         Rotational updates for 5DOF shell using the
c         Hardwired to the Energy-Momentum Method

c      INPUT variables
c         tn = 3x3 Lambda at t_n.
c         vn = 3x1 Spatial rotational velocity at t_n.
c         an = 3x1 Spatial rotational acceleration at t_n.
c         dm = 2x1 Material Delta-T at t_n+1 (BASIC SOLVER VARIABLE)
c         isw=     Task switch := 1 Begining  of the time step
c                              := 2 Iteration within time step
c                              := 3 Back to Beginig  time step

c      OUTPUT variables
c         ta = 3x3 Lambda at t_n+a.
c         tl = 3x3 Lambda at t_n+1.
c         ds = 3x1 Spatial  Delta-t.
c         v1 = 3x1 Spatial director velocity at t_n+1.
c         a1 = 3X1 Spatial director acceleration at t_n+1 (NOT NEEDED)

c      IMPORTANT routine variables
c         du = 3x1 Spatial rotational increment at t_n+a.
c         dl = 3x3 Delta-Lambda.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'ddata.h'
      include  'fdata.h'
      include  'tdata.h'

      integer   isw  , i , j
      real*8    dtr    ,tn(3,3),ta(3,3),tl(3,3),dl(3,3),dm(2),ds(3)
      real*8    du(3)  ,vn(3)  ,v1(3)  ,an(3),a1(3), tt(3,3)

      save

c     INITIALIZE THE ARRAYS:

      if(isw.eq.0) then
c       N.B. Shell uses initial values from directors.
        return
      endif

c     AT THE BEGINING OF A TIME STEP:

      if (fl(9) .and. (isw.eq.1 .or. isw.eq.3)) then

c     Move {v,a}n+1 -> {v,a}n

        if (isw.eq.1) then
          call pmove ( v1 , vn , 3 )
          call pmove ( a1 , an , 3 )
        else
          call pmove ( vn , v1 , 3 )
          call pmove ( an , a1 , 3 )
        endif

c       Zero increment

        call pzero ( du , 3 )
        call pzero ( ds , 3)

c     AT AN ITERATION WITHIN A TIME STEP

      else

c     Calculate Spatial Increments:

        do i = 1 , 3
          ds(i) = tl(i,1)*dm(1) + tl(i,2)*dm(2)
        end do ! i

c     Calculate du:

        du(1) = tl(2,3)*ds(3) - tl(3,3)*ds(2)
        du(2) = tl(3,3)*ds(1) - tl(1,3)*ds(3)
        du(3) = tl(1,3)*ds(2) - tl(2,3)*ds(1)
      endif

c     Compute rotation matrix dl associated with du

      call lamrot (du, dl)

c     Apply Delta.Lambda to Lambda-n+1,k:

      call pmove( tl, tt, 9 )
      do j = 1 , 3
        do i = 1 , 3
          tl(i,j) = dl(i,1)*tt(1,j) + dl(i,2)*tt(2,j)
     &            + dl(i,3)*tt(3,j)
        end do ! i
      end do ! j

c    Compute Lambda-n+alpha

      if (fl(9))then
        call pmove(ta,tt,9)
        do i=1,3
          du(i) = du(i)*theta(3)
        end do ! i
        call lamrot (du, dl)
        do i = 1 , 3
          do j = 1 , 3
            ta(i,j) = dl(i,1)*tt(1,j) + dl(i,2)*tt(2,j)
     &              + dl(i,3)*tt(3,j)
          end do ! j
        end do ! i
      else
        call pmove(tl,ta,9)
      endif

c     DYNAMIC UPDATES: Update Velocities:

      if (fl(9)) then
        dtr  = 2.d0/dt
        do i = 1 , 3
          v1(i) = dtr*(tl(i,3)-tn(i,3))-vn(i)
        end do ! i
      endif

      end
