c$Id: ck2dblk.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine ck2dblk(ix,xl,shp,nel,ndm, err)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Check 2d block for correct numbering

c      Inputs:
c         ix(*)     - List of nodes on block
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh
c         nel       - Number of block nodes

c      Outputs:
c         None

c      Scratch:
c         shp(*)    - Storage for shape functions
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      logical   err
      integer   ndm,nel, ineg, i,l, xn(9),yn(9),ic(16),ix(*)
      real*8    xsj, ss(2),shp(*),xl(ndm,*),jac(16)

      save

      data      xn/-1, 1,1,-1, 0,1,0,-1,0/
      data      yn/-1,-1,1, 1,-1,0,1, 0,0/

      if(nel.ge.4) then
        ineg = 0
        err  = .false.
        do l = 1,nel
          ss(1) = xn(l)
          ss(2) = yn(l)
          call  shp2d (ss,xl,shp,xsj,ndm,nel,ix,.false.)
          if(xsj.le.0.0d0) then
            ineg      = ineg + 1
            ic(ineg)  = ix(l)
            jac(ineg) = xsj
          endif
        end do ! l
        if(ineg.gt.0) then
          err = .true.
          write(iow,2001) (ic(i),jac(i),i=1,ineg)
          if(ior.lt.0) then
            write(*,2001) (ic(i),jac(i),i=1,ineg)
          endif
          call iprint(ix,1,4,1,'IXL-nodes')
          call mprint(xl,2,4,ndm,'XL-coord')
        endif
      endif

2001  format(/5x,'Block has zero or negative jacobian':,
     &        ' at nodes:'/
     &      (10x,'Node =',i3,' Jacobian =',1p,1e12.5))

      end
