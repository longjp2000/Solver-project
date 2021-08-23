c$Id: slcn1d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine slcn1d(sig,shp,xsj,p,s,se,lint,nel,nes)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Project element variables to nodes

c      Inputs:
c        sig(nes,*) - Stresses at quadrature points
c        shp(2,4,*) - Shape functions at quadrature points
c        xsj(*)     - Volume element at quadrature points
c        lint       - Number of quadrature points
c        nel        - Number nodes on element
c        nes        - Dimension of stress array

c      Outputs:
c        p(nen)   - Weights for 'lumped' projection
c        s(nen,*) - Integral of variables
c        se(nen)  - Error projectors
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'sdata.h'
      include  'prstrs.h'
      include  'strnum.h'
      include  'pointer.h'
      include  'comblk.h'

      integer   nel,nes, i,l,lint
      real*8    p(*),s(nen,*),se(*),xsj(*) ,sig(nes,*)
      real*8    shp(2,4,*), xg

      save

c     Lumped and consistent projection routine

      do l = 1,lint

c       Compute lumped projection and assemble stress integrals

        do i = 1,nel

          xg   = shp(2,i,l)*xsj(l)
          p(i) = p(i) + xg

c         Stress projections

          s(i,1) = s(i,1) + sig(1,l)*xg
          s(i,2) = s(i,2) + sig(2,l)*xg
          s(i,3) = s(i,3) + sig(3,l)*xg

c         Error estimation projection

          se(i)  = se(i)  + erav*xg

        end do ! i
      end do ! l

      iste = 3

      end
