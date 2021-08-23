c$Id: c2dplot.f,v 1.1 2006/11/20 20:32:46 rlt Exp $
      subroutine c2dplot (x,ix1,ix2,pen1,pen2,ifsurf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact 2d GEOmetry PLot

c      Purpose: Plot geometry of 2D contact pair

c      Inputs:
c         x(*)    - Nodal coordinates in deformed position
c         ix1(*)  - Element nodal connection list for surface 1
c         ix2(*)  - Element nodal connection list for surface 2

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_geom.h'
      include  'cdata.h'
      include  'sdata.h'
      include  'pdatay.h'

      integer   nsy,ke,n1,n2
      integer   ix1(dnope1,*),ix2(dnope2,*),pen1,pen2,ifsurf
      real*8    x(3,*)

      save

      call cdebug0 ('  c2dplot',-1)

c     Reflect coordinates for symmetry plots

      do nsy = 1,nsym
        lsym = isym(nsy)
        call pltsym(x,3,numnp,lsym)

c       Plot slave surfaces
c       WARNING - the current version work for 2-node elements

        if (ifsurf.ne.2) then
          call pppcol(pen1,0)

          if(dnope1.gt.1) then
            do ke = 1, neps1
              n1 = ix1(1,ke)
              n2 = ix1(2,ke)
              call plotl(x(1,n1),x(2,n1),x(3,n1),3)
              call plotl(x(1,n2),x(2,n2),x(3,n2),2)
            end do ! ke
          else
            do ke = 1, neps1
              n1 = ix1(1,ke)
              n2 = 1
              call pltnod(x(1,n1),n2,ndm,1,0,1,1)
            end do ! ke
          endif
        endif

c       Plot master surfaces

        if (ifsurf.ne.1) then
          call pppcol(pen2,0)

          if(dnope2.gt.1) then
            do ke = 1, neps2
              n1 = ix2(1,ke)
              n2 = ix2(2,ke)
              call plotl(x(1,n1),x(2,n1),x(3,n1),3)
              call plotl(x(1,n2),x(2,n2),x(3,n2),2)
            end do ! ke
          else
            do ke = 1, neps1
              n1 = ix2(1,ke)
              n2 = 1
              call pltnod(x(1,n1),n2,ndm,1,0,1,1)
            end do ! ke
          endif
        endif

c       Reset coordinates

        call pltsym(x,3,numnp,lsym)
      end do

      end
