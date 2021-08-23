c$Id: damp1d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine damp1d(d,xl,s,ndf,ndm,nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute mass matrix for plane and axisymmetric problems

c      Inputs:
c         d(*)      - Material set parameters
c         xl(ndm,*) - Nodal coordinates for element
c         ndf       - Number dof/node
c         ndm       - Spatial dimension of mesh
c         nst       - Size of element arrays

c      Outputs:
c         s(nst,*)  - Consistent or interpolated mass
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'pmod2d.h'

      integer   ndf,ndm,nst, ii,i1, jj,j1, l,lint
      real*8    d(*),xl(ndm,*),s(nst,nst), shp(2,4)
      real*8    sg(2,4), xsj,dv, aj1,xx

      save

c     Compute damping matrix

      if(nint(d(182)).eq.1) then
        lint = nel
        call int1dn(lint,sg)
      else
        lint = nint(d(5))
        call int1d(lint,sg)
      endif

      do l = 1,lint

c       Compute shape functions

        call shp1d(sg(1,l),xl,shp,ndm,nel,xsj)
        dv = sg(2,l)*abs(xsj)*d(70)
        if(stype.eq.3) then
          xx = 0.0d0
          do jj = 1,nel
            xx = xx + shp(2,jj)*xl(1,jj)
          end do ! jj
          dv = dv*xx
        end if

c       Compute db = c*shape*dv

        j1 = 1
        do jj = 1,nel

c         Compute damping matrix

          aj1      = shp(2,jj)*dv
          i1  = 1
          do ii = 1,nel
            s(i1,j1) = s(i1,j1) + shp(2,ii)*aj1
            i1       = i1 + ndf
          end do ! ii
          j1 = j1 + ndf
        end do ! jj
      end do ! l

      end
