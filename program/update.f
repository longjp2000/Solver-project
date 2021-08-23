c$Id: update.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine update(id,f,u,urate,du,fdyn,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Static and Dynamic UPDATE algorithms for F E A P

c      Inputs:
c       id (ndf,numnp,*)   - ID-array (vector storage)
c       f  (ndf,numnp,2)   - Nodal force/displ: 1 = current; 2 = previous
c       u (3*nneq)         - Displacement vectors
c       urate(ndf,numnp,*) - Rate vectors fixed by ALGO
c       du(nneq)           - Displacement increment from SOLVER
c       nneq               - numnp * ndf
c       ndf                - Number of DOF/node
c       fdyn               - Flag: true for dynamics
c       isw                - Control switch
c                            1  STARTING update: begining of time step
c                            2  UPDATE at an iteration within time step
c                            3  BACK solution to begining of time step

c      Outputs:
c       u (ndf,numnp,3)    - Displacement vectors
c       urate(ndf,numnp,*) - Rate vectors fixed by ALGO
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'complx.h'
      include   'counts.h'
      include   'crotas.h'
      include   'fdata.h'
      include   'iofile.h'
      include   'part0.h'
      include   'ptest.h'
      include   'sdata.h'
      include   'tdatb.h'
      include   'pointer.h'
      include   'comblk.h'

      logical    fdyn
      integer    i, j, n, isw, id(ndf,numnp,*)
      real*8     ub,cc
      real*8     f(ndf,numnp,*),u(ndf,numnp,*),urate(ndf,numnp,*),du(*)

      save

c     ISW = 1: UPDATE SOLUTION VECTORS TO BEGIN STEP

      if(isw.eq.1) then

c       Compute and output norms

        if(pfr) then
          call udnorm(urate)
        endif

c       Update rate terms

        call dynlib(u,urate, 1)

c     ISW = 2: UPDATES WITHIN A TIME STEP

      elseif(isw.eq.2) then

c       STEP 1: Update displacement and its increments within step.

        call updatl(id,mr(np(100)),du,u,hr(np(41)),f,ndf,numnp)

c       Do complex imaginary updates

        if(cplxfl) then
          call cupdat(id,u(1,1,4),du(neq+1),nneq,ndf)
        endif

c       Update Lagrange multiplier variables

        if(ndl.gt.0) then
          call uplagm(du,hr(np(213)),mr(np(211)),mr(np(32)),mr(np(33)))
        endif

c       Rigid body transformations and update correctors

        call rigidb(3,isw,fdyn)

c       STEP 2: For Dynamics, update rate terms [urate-vectors]

        if(fdyn) then
          call dynlib(u,urate, 2)
        endif

c       Contact update

        call contact (314)

c     ISW = 3: BACKUP SOLUTION VECTORS: Reinitiate a transient step

      elseif(isw.eq.3) then

        if(fdyn)then
          call dynlib(u,urate, 3)
        endif

c       Back up solution vectors: u(*)

        cc = 1.d0/cc2
        do i = 1,ndf
          if(ndfp(i).eq.npart) then
            do n = 1,numnp
              ub = u(i,n,2)*cc
              j  = id(i,n,1)

c             Compute values from current solution

              if (j.gt.0) then
                u(i,n,1) = u(i,n,1) - ub
                u(i,n,2) = 0.0d0
                u(i,n,3) = 0.0d0

c             Compute values from forced inputs

              else
                du(neq-j) = 0.0d0
                u(i,n,3)  = 0.0d0
                u(i,n,2)  = 0.0d0
                f(i,n,1)  = f(i,n,2)
                u(i,n,1)  = f(i,n,2)
              endif
            end do ! n
          endif
        end do ! i
      endif

c     Rigid body transformations and update correctors

      call rigidb(4,isw,fdyn)

c     Update transformation matrices

      if (frotas) then
        call updrot(u(1,1,3),ndf,hr(np(82)),mr(np(81)),numnp,isw)
      endif

      end
