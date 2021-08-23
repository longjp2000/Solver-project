c$Id: stcn2z.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine stcn2z(xl,sig,shp,xsj,lint,ndm,nel,nen)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: ZZ Element Projection

c      Inputs:
c         xl(ndm,*)    - Element node coordinates
c         sig(nen,*)   - Element stresses at quadrature points
c         shp(3,nen,*) - Shape functions at quadrature points
c         xsj(*)       - Jacobian at quadrature point (time weight)
c         lint         - Number quadrature points
c         ndm          - Mesh spatial dimension
c         nel          - Number nodes on element
c         nen          - Dimension for stress and shape functions

c      Outputs:
c         ek(10,10)    - Contribution to projection matrix
c         est(30,10)   - Contribution for each projection component
c                        N.B. Returned in 'zzcom1.h'
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'debugs.h'
      include   'strnum.h'
      include   'zzcom1.h'

      integer    ndm,nel,nen, i,j,l,lint, nps
      real*8     xsj(*),xl(ndm,*),shp(3,nen,*),sig(nen,*),xx(10), cons

c     Compute zz-projection and assemble local stress integrals

      if(debug) then
        call mprint(sig,4,lint,nen,'SIG-ST')
      endif

      nps = 3
      if(nel.eq.6 .or. nel.eq.8 .or. nel.eq.9) then
        nps = 6
      elseif(nel.eq.16) then
        nps = 10
      end if

      xx(1) =  1.d0
      do l = 1,lint
        xx(2) = -xnodz(1)
        xx(3) = -xnodz(2)
        do i = 1,nel
          xx(2) = xx(2) + shp(3,i,l)*xl(1,i)
          xx(3) = xx(3) + shp(3,i,l)*xl(2,i)
        end do ! i

c       Quadratic and cubic projections

        if(nps.ge.6) then
          xx( 4) = xx(2)*xx(2)
          xx( 5) = xx(2)*xx(3)
          xx( 6) = xx(3)*xx(3)
        endif

c       Cubic projection

        if(nps.eq.10) then
          xx( 7) = xx(4)*xx(2)
          xx( 8) = xx(5)*xx(2)
          xx( 9) = xx(5)*xx(3)
          xx(10) = xx(6)*xx(3)
        endif

c       Accumulate matrix and projection contributions

        do i = 1,nps
          cons = xx(i)*xsj(l)
          do j = 1,nps
            ek(j,i)  = ek(j,i) + xx(j)*cons
          end do ! j
          do j = 1,6
            est(j,i) = est(j,i) + cons*sig(j,l)
          end do ! j
        end do ! i

      end do ! l

      iste = 6

      end
