c$Id: cktris.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine cktris(ix,xl,shp,ndm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Check isoparametric triangles elements for valid data

c      Inputs:
c         ix(*)     - Nodes connected to element
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of problem

c      Outputs:
c         None

c      Scratch:
c         shp(*)    - Storage for shape functions
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'comblk.h'

      integer   ndm, ineg,  i,l, ic(2,7),ix(*)
      real*8    xsj, shp(*),xl(ndm,*),el(3,7)

      save

      data el/0.0d0,0.0d0,1.0d0, 1.0d0,0.0d0,0.0d0, 0.0d0,1.0d0,0.0d0,
     &        0.5d0,0.0d0,0.5d0, 0.5d0,0.5d0,0.0d0, 0.0d0,0.5d0,0.5d0,
     &        0.3333333333333d0, 0.3333333333333d0, 0.3333333333333d0/

c     Check element for input errors

      ineg = 0
      do l = 1,nel
        if(ix(l).gt.0) then
          if(mr(np(190)+ix(l)-1).lt.0) then
            ineg       = ineg + 1
            ic(1,ineg) = l
            ic(2,ineg) = abs(ix(l))
          endif
        endif
      end do ! l

      if(ineg.gt.0) then
        write(iow,2000) n,(ic(1,i),ic(2,i),i=1,ineg)
        if(ior.lt.0) then
          write(*,2000) n,(ic(1,i),ic(2,i),i=1,ineg)
        endif
      else
        do l = 1,nel
          call  trishp(el(1,l),xl,ndm,nel-4,xsj,shp)
          if(xsj.le.0.0d0) then
            ineg       = ineg + 1
            ic(1,ineg) = l
            ic(2,ineg) = abs(ix(l))
          endif
        end do ! l

        if(ineg.gt.0) then
          write(iow,2001) n,(ic(1,i),ic(2,i),i=1,ineg)
          if(ior.lt.0) then
            write(*,2001) n,(ic(1,i),ic(2,i),i=1,ineg)
          endif
        endif

c       Try to fix element by reversing order of nodes on 'ix'

        if(ineg.eq.nel) then
          if(ior.lt.0) write(*,2002) n
          l     = ix(2)
          ix(2) = ix(3)
          ix(3) = l
          if(nel.eq.6 .or. nel.eq.7) then
            l     = ix(6)
            ix(6) = ix(4)
            ix(4) = l
          endif
        endif

      endif

c     Formats

2000  format(' >Element',i9,' coordinates not input for nodes:'/
     &      ('                    Local =',i3,' Global =',i9))

2001  format(' >Element',i9,' has negative jacobian at nodes:'/
     &      ('                    Local =',i3,' Global =',i9))

2002  format(' >Element',i9,' Reverse numbers to fix negative jacobian')
      end
