c$Id: gamma1.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      double precision function gamma1(id,pu,pr,du,t,s)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute solution energy for step size 's'

c      Inputs:
c         id(*)  - Equation numbers for each dof
c         pu     - Pointer to solution vectors
c         du(*)  - Last increment to solution
c         s      - Step size

c      Outputs:
c         gamma1 - Solution energy for current step

c      Scratch:
c         pr     - Pointer to residual array
c         t(*)   - Temporary storage for solution arrays
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'counts.h'
      include  'ddata.h'
      include  'fdata.h'
      include  'ndata.h'
      include  'pfeapb.h'
      include  'pointer.h'
      include  'prlod.h'
      include  'sdata.h'
      include  'tdata.h'
      include  'comblk.h'

      logical   fa,tr
      integer   n,pr,pu, id(*), gneq,pneq
      real*8    dot, s, du(*),t(nneq,*), gammap,tdata

      save

      data      fa,tr/.false.,.true./

c     Check for parallel solution

      if(pfeap_on) then
        pneq = numpeq
        gneq = numnp*ndf
      else
        pneq = neq
        gneq = neq
      endif

c     Increment counter on RHS forms

      iform = iform + 1

c     Get search displacement

c     Move quantities for saves

      call pmove(hr(pu)  ,t     ,3*nneq)
      call pmove(    du  ,t(1,4),  nneq)
      if(np(42).ne.0) call pmove(hr(np(42)),t(1,5),nneq*nrt)

c     Multiply increment by current step size

      do n = 1,gneq
        du(n) = du(n)*s
      end do ! n

c     Update with step control

      call update(id,hr(np(30)),hr(pu),hr(np(42)),du,fl(9),2)

c     Compute residual

      call pload(id,hr(pu),hr(np(30)),hr(pr),prop,tr,fa)
      call formfe(pu,pr,pr,pr,fa,tr,fa,fa,6,1,numel,1)

c     Restore quantities from saves

      call pmove( t     ,hr(pu),3*nneq)
      call pmove( t(1,4),    du,  nneq)
      if(np(42).ne.0) call pmove(t(1,5),hr(np(42)),nneq*nrt)

c     Compute value of gamma

      gammap = dot (du,hr(pr),pneq)
      if(pfeap_on) then
        call pfeapsr(gammap,tdata,1)
      endif
      gamma1 = gammap

      end
