c$Id: dyna05.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine dyna05(du,urate,nneq,ndf,ndp,ndo,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Perform Static and ODE updates using conserving form for
c               rate terms.

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

      include  'ddata.h'
      include  'part0.h'
      include  'tdata.h'

      integer   i, n, nneq,ndf,isw, ndp(*),ndo(*)
      real*8    cn3,cn4,cn5,cn6, ur1,ur2, du(*),urate(nneq,*)

      save

c     Update solution vectors to begin a step

c     Three-parameter algorithm in Conservation form

c          u(n,1)     = d(t-n+1)
c          u(n,2)     = d(t-n+a) - d(t-n);
c          u(n,3)     = dd(t_n+1)
c                     = [primary dependent variable for solver]
c          urate(n,1) = v(t-n+1)
c          urate(n,2) = a(t-n+1)
c          urate(n,3) = u(t-n+a)
c          urate(n,4) = v(t-n+a);
c          urate(n,5) = a(t-n+a)

c     Starting condition based on d(t-n+1) = d(t-n)

      if(isw.eq.1) then

        cn3 = 1.d0 - 0.5d0/theta(1)
        cn4 = 1.d0 - theta(2)/theta(1)
        cn5 = dt*(1.d0 - 0.5d0*theta(2)/theta(1))
        cn6 = 1.d0/(theta(1)*dt)

        do i = 1,ndf
          if(ndp(i).eq.npart) then
            if(ndo(i).ge.1) then
              do n = i,nneq,ndf
                urate(n,6) = urate(n,1)
                urate(n,7) = urate(n,2)
                ur1        =  cn4*urate(n,1) + cn5*urate(n,2)
                ur2        = -cn6*urate(n,1) + cn3*urate(n,2)
                urate(n,3) =  du(n)
                urate(n,4) = (1.d0-theta(3))*urate(n,1) + theta(3)*ur1
                urate(n,5) = (ur1-urate(n,1))/dt
                urate(n,1) =  ur1
                urate(n,2) =  ur2
              end do
            else
              do n = i,nneq,ndf
                urate(n,3) =  du(n)
              end do
            endif
          endif
        end do

c     Conserving algorithm solution

      elseif(isw.eq.2) then

        do i = 1,ndf
          if(ndp(i).eq.npart) then
            if(ndo(i).ge.2) then
              do n = i,nneq,ndf
                urate(n,1) = urate(n,1) + c2*du(n)
                urate(n,2) = urate(n,2) + c1*du(n)
                urate(n,3) = urate(n,3) + c3*du(n)
                urate(n,4) = urate(n,4) + c4*du(n)
                urate(n,5) = urate(n,5) + c5*du(n)
              end do
            elseif(ndo(i).eq.1) then
              do n = i,nneq,ndf
                urate(n,1) = urate(n,1) + c2*du(n)
                urate(n,2) = 0.0d0
                urate(n,3) = urate(n,3) + c3*du(n)
                urate(n,4) = urate(n,4) + c4*du(n)
                urate(n,5) = 0.0d0
              end do
            else
              do n = i,nneq,ndf
                urate(n,3) = urate(n,3) + c3*du(n)
              end do
            endif
          endif
        end do

c     Backup solution vectors to reinitiate a step

      elseif(isw.eq.3) then

        do i = 1,ndf
          if(ndp(i).eq.npart) then
            do n = i,nneq,ndf
              urate(n,1) = urate(n,6)
              urate(n,2) = urate(n,7)
            end do
          endif
        end do

      endif

      end
