c$Id: mass1d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine mass1d(d,xl,s,p,ndf,ndm,nst)

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
c         p(nst)    - Diagonal (lumped) mass
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'pmod2d.h'

      integer   ndf,ndm,nst, ii,i1, jj,j1, l,lint
      real*8    d(*),xl(ndm,*),s(nst,nst),p(nst), shp(2,4)
      real*8    sg(2,4), xsj,dv, aj1,xx,lfac,cfac

      save

c     Compute mass matrix

      if(nint(d(182)).eq.1) then
        lint = nel
        call int1dn(lint,sg)
      else
        lint = nint(d(5))
        call int1d(lint,sg)
      endif
      cfac = d(7)
      lfac = 1.d0 - cfac

      do l = 1,lint

c       Compute shape functions

        call shp1d(sg(1,l),xl,shp,ndm,nel,xsj)
        dv = sg(2,l)*abs(xsj)*d(4)
        if(stype.eq.3) then
          xx = 0.0d0
          do jj = 1,nel
            xx = xx + shp(2,jj)*xl(1,jj)
          end do ! jj
          dv = dv*xx
        end if

c       Compute db = rho*shape*dv

        j1 = 1
        do jj = 1,nel

c         Compute lumped mass matrices

          aj1      = shp(2,jj)*dv
          p(j1)    = p(j1)    + aj1
          s(j1,j1) = s(j1,j1) + aj1*lfac

c         Compute consistent mass matrix

          aj1 = aj1*cfac
          i1  = 1
          do ii = 1,nel
            s(i1,j1) = s(i1,j1) + shp(2,ii)*aj1
            i1       = i1 + ndf
          end do ! ii
          j1 = j1 + ndf
        end do ! jj
      end do ! l

      end
