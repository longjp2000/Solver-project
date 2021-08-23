c$Id: dyna08.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine dyna08(du,urate,nneq,ndf,ndp,ndo,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Perform 2nd order ODE update using explicit Central
c               Difference method.

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

      include  'comblk.h'
      include  'part0.h'
      include  'pointer.h'
      include  'tdata.h'

      integer   i, n, nod,nneq,nneq2,ndf,isw, ndp(*),ndo(*)
      real*8    ub, du(*),urate(nneq,*)

      save

c     Update solution vectors at start of step

      if(abs(isw).eq.1) then

        nneq2 = nneq + nneq
        do i = 1,ndf
          if(ndp(i).eq.npart .and. ndo(i).ge.2) then
            nod = 0
            do n = i,nneq,ndf

c             Compute values from forced inputs for fixed dof
c                  RIXT = np(100)
c             if(  mr(np(100)+nod).eq.0  .or. isw.lt.0        ) then
c                  ID = np(31)
                if(mr(np(31)+n-1).gt.0       .or. isw.lt.0) then

c                 Save velocity and acceleration from t_n

                  if(isw.gt.0) then
                    urate(n,6)  = urate(n,1)
                    urate(n,7)  = urate(n,2)
                  endif

c                 Update velocity and set acceleration zero

                  urate(n,1)  = urate(n,1) + c1*urate(n,2)
                  urate(n,2)  = 0.0d0

                  ub          = dt*urate(n,1)
                  du(n)       = du(n)         + ub
                  du(n+nneq)  =                 ub
                  du(n+nneq2) =                 ub

c                 Make u_n+alpha = u_n+1

                  urate(n,3)  = du(n)

c               Compute values from current solution state

                else
c                              FTN = np(30)
                  du(n)       = hr(np(30)+n-1)
                  du(n+nneq)  = du(n) - hr(np(30)+n+nneq-1)
                  du(n+nneq2) = du(n+nneq)

c       WARNING: Boundary velocity and accelerations are NOT computed.
c                Do not use consistent mass or problems with damping
c                with specified boundary motions!

                endif
c             endif
              nod = nod + 1
            end do ! n
          endif
        end do ! i

c     Update solution vectors within step

      elseif(abs(isw).eq.2) then

        do i = 1,ndf
          if(ndp(i).eq.npart .and. ndo(i).ge.2) then
            do n = i,nneq,ndf
              urate(n,2) = urate(n,2) + du(n)
c             du(n)      = 0.0d0
            end do ! n
          endif
        end do ! i

c     Backup solution vectors

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
