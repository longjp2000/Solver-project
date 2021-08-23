c$Id: dyna01.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine dyna01(du,urate,nneq,ndf,ndp,ndo,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Perform static and ODE updates using Newmark formulas.

c      Inputs:
c         du(*)             Increment to displacement
c         urate(nneq,*)     Rate vectors - fixed by ALGO
c         nneq              numnp * ndf
c         ndf               Number of DOF/node
c         ndp(*)            Partition dof's
c         ndo(*)            Order dof's: 0 = static, n = n-th order.
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
      real*8    ur1,c6, du(*),urate(nneq,*)

      save

c     Newmark-Beta updates: urate(n,1) = velocity
c                           urate(n,2) = acceleration

c     (1) Update solution vectors to begin a step

      if(isw.eq.1) then

c       Initialize step

        c6 = dt*c1
        do i = 1,ndf
          if(ndp(i).eq.npart .and. ndo(i).ge.1) then
            do n = i,nneq,ndf
              urate(n,6) =   urate(n,1)
              urate(n,7) =   urate(n,2)
              ur1        = - c6*urate(n,1) + c3*urate(n,2)
              urate(n,1) =   c4*urate(n,1) + c5*urate(n,2)
              urate(n,2) =   ur1
            end do ! n
          endif
        end do ! i

c     (2) Update with in time step

      elseif(isw.eq.2) then

        do i = 1,ndf
          if(ndp(i).eq.npart) then
            if(ndo(i).ge.2) then
              do n = i,nneq,ndf
                urate(n,1) = urate(n,1) + c2*du(n)
                urate(n,2) = urate(n,2) + c1*du(n)
              end do ! n
            endif
            if(ndo(i).eq.1) then
              do n = i,nneq,ndf
                urate(n,1) = urate(n,1) + c2*du(n)
                urate(n,2) = 0.0d0
              end do ! n
            endif
          endif
        end do ! i

c     (3) Backup solution vectors to reinitiate a step

      elseif(isw.eq.3) then

        do i = 1,ndf
          if(ndp(i).eq.npart) then
            do n = i,nneq,ndf
              urate(n,1) = urate(n,6)
              urate(n,2) = urate(n,7)
            end do ! n
          endif
        end do ! i

      endif

      end
