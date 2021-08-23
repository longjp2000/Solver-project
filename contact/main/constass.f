c$Id: constass.f,v 1.1 2006/11/20 20:32:46 rlt Exp $
      subroutine constass (ixl,ida,nnod,ndof,ilm,lnod,nlag,size,
     &                     tanm,resv)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0
c               Robert L Taylor          October 28, 1996            1.1

c      Acronym: CONtact STiffness ASSembling

c      Purpose: Call assembling subroutine in blind way for the user

c      Inputs :
c         ixl(*)  - IX Local vector with involved nodes
c         ida(*)  - ID Active vector with active contact dof
c         nnod    - # of nodes in ixl
c         ndof    - # of dof in ida
c         ilm(*)  - Nodes for Lagrange Multipliers
c         lnod    - Number of multiplier nodes
c         nlag    - Number LAGrange multipliers/node
c         size    - SIZE of tangent matrix (dimension of array)
c         tanm(*) - TANgent Matrix
c         resv(*) - RESidual Vector

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_contac.h'
      include  'c_tanfl.h'
      include  'compas.h'
      include  'corset.h'
      include  'counts.h'
      include  'eqsym.h'
      include  'ndata.h'
      include  'part0.h'
      include  'pointer.h'
      include  'ptdat4.h'
      include  'rigid1.h'
      include  'sdata.h'
      include  'tdatb.h'
      include  'comblk.h'

      logical   rel
      integer   ixl(*),ida(*),ilm(*),nnod,ndof,lnod,nlag,size
      integer   idl(100),nrot(2),i,i1,i2,iu,ix,j,k1,k2, idlag
      real*8    tanm(*),resv(*), umax,du,tol

      save

      call cdebug0 ('      constass ->',-1)

c     Use defined idl vector if big enough

      if (size.le.100) then

c       Get involved DOF

        do i = 1,size
          idl(i) = 0
        end do ! i
        call defdof (ixl,ida,idl,nnod,ndof,mr(np(31)))

c       Add Lagrange multipler equations to idl(*)

        if(lnod.gt.0) then
          i1 = nnod*ndof
          do i2 = 1,lnod
            if(ilm(i2).gt.0) then
              iu = idlag(mr(np(224)),ilm(i2)) - nlag
              do i = 1,nlag
                idl(i1+i) = iu + i
              end do ! i
              i1 = i1 + nlag
            endif
          end do ! i2

        endif

c       Transform stiffness/residual for sloping boundary in 1-2 plane

        if(anglefl) then
          call pangl(ixl,nnod,hr(np(46)),hr(np(45)),nrot(1))
        endif
        if(eulerfl) then
          call peule(ixl,nnod,hr(np(243)),hr(np(242)),nrot(2))
        endif
        if(nrot(1)+nrot(2).gt.0) then
          call ptlocal(du,resv,tanm,.false.,ndof,1,nnod,size,nrot,2)
        endif

c       Modify for specified boundary displacements

        umax = 0.0d0
        do i1 = 1,nnod
          do i2 = 1,ndof
            iu = ndf*(ixl(i1) - 1) + ida(i2) - 1
            umax = max(umax,abs(hr(np(30)+iu)))
          end do
        end do

        k1  = 0
        k2  = 0
        tol = umax*1.d-08
        do i1 = 1,nnod
          do i2 = 1,ndof
            iu = ndf*(ixl(i1) - 1) + ida(i2) - 1
            du = (hr(np(30)+iu) - hr(np(40)+iu))*cc3
            if(idl(k2+i2).le.0 .and. abs(du).gt.tol) then
              do j = 1,size
                resv(j) = resv(j) - tanm(k1+j)*du
              end do
            endif
            k1 = k1 + size
          end do
          k2 = k2 + ndof
        end do

c       Check for rigid body interface nodes

        rel = .false.
        if(rbody .and.nrbprt.eq.npart) then
          iu = 0
          ix = 0
          call pzero(hr(np(44)), ndm*nnod)
          call pzero(hr(np(41)),ndof*nnod)
          do i = 1,nnod
            if(mr(np(100)+ixl(i)-1).ne.0) then
              rel  = .true.
              do j = 1,ndm
                hr(np(44)+ix+j) = hr(np(43)+(ixl(i)-1)*ndm+j)
              end do
              do j = 1,ndof
                hr(np(41)+iu+ida(j)) = hr(np(40)+(ixl(i)-1)*ndf+ida(j))
              end do
            endif
            iu = iu + ndof
            ix = ix + ndm
          end do
        endif

c       Transform and assemble rigid body part

        if(rel) then

          call rasbly(tanm,resv,hr(np(44)),hr(np(41)),idl,
     &                mr(np(20+npart)),ixl,mr(np(100)),mr(np(96)),
     &                mr(np(99)),hr(np(95)),ndm,ndof,nnod,nnod,size,
     &                lafl,uafl,dbfl,ddfl,hr(np(26)),
     &                hr(nal),hr(nau),hr(na))

c       Assemble for real arithmetic

        else

          call dasble (tanm,resv,idl,mr(np(20+npart)),size,neqs,
     &                 uafl,dbfl,hr(np(26)),hr(nal),hr(nau),hr(na))
        endif

c       Store time history contact plot data

        if (ncplts.gt.0) then

          do i1 = 1,nnod
            do j = 1, ncplts
              if (icpl(1,j).eq.ixl(i1)) then
                do i2 = 1,ndof
                  if(icpl(2,j).eq.ida(i2)) then
                    cpl(j) = cpl(j) - resv(ndof*(i1-1) + i2)
                  endif
                end do ! i2
              end if
            end do ! j
          end do ! i1

        end if ! ncplts

c     Allocate scratch idl vector
c     WARNING still to be done

      else
        write (*,*) 'CONSTASS - idl vector too small'
        call plstop()
      endif

      end

      integer function idlag(ida,i2)

      include  'cdata.h'
      integer   ida(numnp,3), i2, ie

      ie    = ida(i2,2)
      if(ie.ne.0) then
        if(ida(ie,1).ne.0) then
          idlag = ida(ie,1) + ida(i2,3)
        else
          idlag = -1
        endif
      else
        idlag = -2
      endif

      end
