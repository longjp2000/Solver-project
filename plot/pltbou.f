c$Id: pltbou.f,v 1.1 2006/11/20 20:33:12 rlt Exp $
      subroutine pltbou(id,x,angl,ip,ib, ndm,ndf,numnp,nbou)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Display boundary conditions on screen for 1-3 dof

c      Inputs:
c         id(ndf,*) - Boundary condition indicator array
c         x(ndm,*)  - Nodal coordinates of mesh
c         angl(*)   - Angle for sloping boundaries
c         ip(*)     - Active node indicators
c         ib(2,*)   - Node rotation dof
c         ndm       - Dimension of x array
c         ndf       - Number dof/node
c         numnp     - Number of nodes in mesh
c         nbou      - Component to display ( 0 = all)

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata1.h'
      include  'pdata4.h'
      include  'pdatay.h'
      include  'pointer.h'
      include  'rigid1.h'
      include  'comblk.h'

      logical   zoom,bc0
      integer   ndm,ndf,numnp,nbou, nsy, n, i
      integer   id(ndf,numnp,*),ip(*),ib(3,*), nn(3),nb(3)
      real*8    dx1, x1,x2,x3, cs,sn,ang, x(ndm,*),angl(*)

      save

c     Plot boundary restraints (lines = fixed)

      bc0 = nbou.eq.0
      dx1 = .006d0/scale
      do nsy = 1,nsym
        lsym = isym(nsy)
        call pltsym(x,ndm,numnp,lsym)
        do n = 1,numnp
          if(ip(n).gt.0) then
            if(angl(n).ne.0.0d0) then
              ang = angl(n)*0.017453292d0
              cs = cos(ang)*dx1
              sn = sin(ang)*dx1
            else
              cs = dx1
              sn = 0.0d0
            endif
            if(zoom(x(1,n),ndm) .and. mr(np(190)+n-1).ge.0 .and.
     &        (.not.rbody .or. (rbody.and.mr(np(100)+n-1).eq.0))) then
              x1 = x(1,n)
              x2 = x(2,n)
              x3 = x(3,n)
              do i = 1,3
                if(ib(i,n).gt.0) then
                  nn(i) = ib(i,n)
                else
                  nn(i) = i
                endif
                nb(i) = id(nn(i),n,2)
              end do ! i
              if (nb(1).gt.0 .and. ( bc0 .or. nbou.eq.nn(1))) then
                call plotl(x1+cs, x2+sn, x3, 3)
                call plotl(x1-cs, x2-sn, x3, 2)
              endif
              if (ndf.ge.2 .and. ndm.ge.2) then
                if(nb(2).gt.0 .and. (bc0 .or. nbou.eq.nn(2))) then
                  call plotl(x1-sn, x2+cs, x3, 3)
                  call plotl(x1+sn, x2-cs, x3, 2)
                endif
              endif
              if (ndf.ge.3 .and. ndm.ge.2)then
                if(nb(3).gt.0 .and. (bc0 .or. nbou.eq.nn(3))) then
                  call plotl(x1,x2, x3+dx1, 3)
                  call plotl(x1,x2, x3-dx1, 2)
                endif
              endif
              if (ndf.ge.4 .and. ndf.ge.nbou .and. ndm.ge.2) then
                i = mod(nbou-1,3) + 1
                if    (id(nbou,n,2) .gt. 0 .and. i.eq.1 ) then
                  call plotl(x1+cs, x2+sn, x3, 3)
                  call plotl(x1-cs, x2-sn, x3, 2)
                elseif(id(nbou,n,2) .gt. 0 .and. i.eq.2 ) then
                  call plotl(x1-sn, x2+cs, x3, 3)
                  call plotl(x1+sn, x2-cs, x3, 2)
                elseif(id(nbou,n,2) .gt. 0 .and. i.eq.3 ) then
                  call plotl(x1,x2, x3+dx1, 3)
                  call plotl(x1,x2, x3-dx1, 2)
                endif
              endif
            endif
          endif

        end do ! n

        call pltsym(x,ndm,numnp,lsym)

      end do ! nsy

      end
