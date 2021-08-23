c$Id: dyna09.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine dyna09(du,urate,nneq,ndf,ndp,ndo,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Perform Static and 1st order ODE using BDF2 method

c      Inputs:
c         du(*)             Increment to displacement
c         urate(nneq,*)     Rate vectors - fixed by ALGO
c         nneq              numnp * ndf
c         ndf               Number of DOF/node
c         ndp(*)            Partition dof's
c         ndo(*)            Order dof's
c         isw               Control switch
c                            1  STARTING update: begining of time step
c                            2  UPDATE at an iteration within time step
c                            3  BACK solution to begining of time step

c      Outputs:
c         urate(nneq,nn)    Rate vectors fixed by ALGO
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'part0.h'
      include  'tdata.h'

      integer   i, n, nneq,ndf,isw, ndp(*),ndo(*)
      real*8    du(*),urate(nneq,*)

      save

c     Backward Difference Formula (BDF2): Initialize at start of step

      if(isw.eq.1) then

        do i = 1,ndf
          if(ndp(i).eq.npart) then
            do n = i,nneq,ndf
              urate(n,6) = urate(n,1)
              urate(n,7) = urate(n,2)
            end do ! n
            if(ndo(i).ge.2) then
              do n = i,nneq,ndf
                urate(n,5) = urate(n,4)
                urate(n,4) = urate(n,3)
                urate(n,3) = 0.0d0
                urate(n,1) = urate(n,4)*c3
                urate(n,2) = urate(n,1)
              end do ! n
            elseif(ndo(i).eq.1) then
              do n = i,nneq,ndf
c               urate(n,5) = urate(n,4)
c               urate(n,4) = urate(n,3)
c               urate(n,3) = 0.0d0
                urate(n,1) = urate(n,4)*c3
                urate(n,2) = 0.0d0
              end do ! n
            endif
          endif
        end do ! i

c     Backward Difference Formula BDF2: Updates in iterations

      elseif(isw.eq.2) then

        do i = 1,ndf
          if(ndp(i).eq.npart) then
            if(ndo(i).ge.2) then
              do n = i,nneq,ndf
                urate(n,1) = urate(n,1) + c1*du(n)
                urate(n,2) = urate(n,1)
c               urate(n,3) = urate(n,3) + c2*du(n)
              end do ! n
            elseif(ndo(i).eq.1) then
              do n = i,nneq,ndf
                urate(n,1) = urate(n,1) + c1*du(n)
                urate(n,2) = 0.0d0
c               urate(n,3) = urate(n,3) + c2*du(n)
                urate(n,3) = 0.0d0
              end do ! n
            endif
          endif
        end do ! i

c     Backward Difference Formula BDF2: Back up to start of step

      elseif(isw.eq.3) then

        do i = 1,ndf
          if(ndp(i).eq.npart) then
            do n = i,nneq,ndf
              urate(n,1) = urate(n,6)
              urate(n,2) = urate(n,7)
              urate(n,3) = urate(n,4)
              urate(n,4) = urate(n,5)
              urate(n,5) = urate(n,5) ! WIP: Need to figure out
            end do ! n
          endif
        end do ! i

      endif

      end
