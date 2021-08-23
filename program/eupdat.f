c$Id: eupdat.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine eupdat(urate,du,eul,nneq,fdyn,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Perform element static and dynamic UPDATE

c      Inputs:
c         urate(nneq,nn)  - Rate vectors fixed by ALGO
c                           nn = 1:  U_n+1
c                           nn = 2: DU_n+a
c                           nn = 3:  V_n+1
c                           nn = 4:  A_n+1
c                           nn = 5:  U_n+a
c                           nn = 6:  V_n+a
c                           nn = 7:  A_n+a
c         du(nneq)        - Displacement increment from SOLVE
c         nneq            - Number of internal equations
c         fdyn            - Flag: true for dynamics
c         isw             - Control switch
c                            1  STARTING update: begining of time step
c                            2  UPDATE at an iteration within time step
c                            3  BACK solution to begining of time step

c      Outputs:
c         eul(nneq,5)     - Solution parameters
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'ddata.h'
      include  'part0.h'
      include  'tdata.h'
      include  'tdatb.h'
      include  'tdato.h'

      logical   fdyn
      integer   n, nneq,isw
      real*8    dtsave, urate(nneq,*),du(*),eul(nneq,5)

      save

c     Update solution vectors to begin a step

      if(isw.eq.1) then

c       Save time step for updates

        dtn = dt

c       Update rate terms

c       Newmark-beta updates

        if(noi.eq.1) then

          call dyna01(urate,urate(1,3),nneq,1,ndfp,ndfo,1)

c       Backward Euler update

        elseif(noi.eq.2) then

          call dyna02(urate,urate(1,3),nneq,1,ndfp,ndfo,1)

c       Conserving HHT update (JCS/Doblare version)

        elseif(noi.eq.3) then

          call dyna03(urate,urate(1,3),nneq,1,ndfp,ndfo,1)

c       Newmark explicit update

        elseif(noi.eq.4) then

          do n = 1,nneq
            du(n)      = dt*urate(n,3) + c1*urate(n,4)
            urate(n,1) = urate(n,1)    + du(n)
            urate(n,2) =                 du(n)
            urate(n,3) = urate(n,3)    + dt*urate(n,4)
          end do ! n

c       Three-parameter algorithm in Conservation form

        elseif(noi.eq.5) then

          call dyna05(urate,urate(1,3),nneq,1,ndfp,ndfo,1)

c       STATICS:  Generalized Mid-point configuration

        elseif(noi.eq.6) then

          call dyna06(urate,urate(1,3),nneq,1,ndfp,ndfo,1)

c       FIRST ORDER :  Generalized Mid-point configuration

        elseif(noi.eq.7) then

          call dyna07(urate,urate(1,3),nneq,1,ndfp,ndfo,1)

        endif

      elseif(isw.eq.2) then

c       STEP 1: Update displacement and its increments within step.

        do n = 1,nneq

          urate(n,1) = urate(n,1) + cc1*du(n)
          urate(n,2) = urate(n,2) + cc2*du(n)

c         Static return values

          eul(n,1)   = urate(n,1)
          eul(n,2)   = urate(n,2)
          eul(n,3)   = du(n)
          eul(n,4)   = 0.0d0
          eul(n,5)   = 0.0d0

        end do ! n

c       STEP 2: For Dynamics, update rate terms [urate-vectors]

        if(fdyn) then

c         Set time step to value for time step

          dtsave = dt
          dt     = dtn

c         Newmark-beta update for velocity and acceleration

          if(noi.eq.1) then

            call dyna01(du,urate(1,3),nneq,1,ndfp,ndfo,2)

c         Backward Euler solution

          elseif(noi.eq.2) then

            call dyna02(du,urate(1,3),nneq,1,ndfp,ndfo,2)

c         Conserving HHT solution

          elseif(noi.eq.3) then

            call dyna03(du,urate(1,3),nneq,1,ndfp,ndfo,2)

c         Newmark explicit Euler solution

          elseif(noi.eq.4) then

            call dyna04(du,urate(1,3),nneq,1,ndfp,ndfo,2)

c         Conserving algorithm solution

          elseif(noi.eq.5) then

            call dyna05(du,urate(1,3),nneq,1,ndfp,ndfo,2)

c         STATIC generalized mid-point rule

          elseif(noi.eq.6) then

            call dyna06(du,urate(1,3),nneq,1,ndfp,ndfo,2)

c         FIRST ORDER generalized mid-point rule

          elseif(noi.eq.7) then

            call dyna07(du,urate(1,3),nneq,1,ndfp,ndfo,2)

          endif

c         Set return values

          if(nrk.gt.0) then
            do n = 1,nneq
              eul(n,1) = urate(n,nrk+2)
            end do ! n
          endif

          if(nrc.gt.0) then
            do n = 1,nneq
              eul(n,4) = urate(n,nrc+2)
            end do ! n
          endif

          if(nrm.gt.0) then
            do n = 1,nneq
              eul(n,5) = urate(n,nrm+2)
            end do ! n
          endif

        endif

c       Restore time step

        dt = dtsave

c     Backup solution vectors to reinitiate a step

      elseif(isw.eq.3) then

        if(fdyn)then

c         Newmark-beta updates

          if(noi.eq.1) then

            call dyna01(urate(1,2),urate(1,3),nneq,1,ndfp,ndfo,3)

c         Backward Euler updates

          elseif(noi.eq.2) then

            call dyna02(urate(1,2),urate(1,3),nneq,1,ndfp,ndfo,3)

          endif
        endif

c       Back up solution vectors: u(*)

        do n = 1,nneq
          urate(n,1) = urate(n,1) - urate(n,2)
          urate(n,2) = 0.d0
        end do ! n
      endif

      end
