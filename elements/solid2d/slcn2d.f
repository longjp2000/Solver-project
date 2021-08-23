c$Id: slcn2d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine slcn2d(ix,sig,shp,xsj,p,s,se,lint,nel,nes)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Project element variables to nodes

c      Inputs:
c        ix(*)    - Nodal connections for element
c        sig(nes,*) - Stresses at quadrature points
c        shp(3,nes,*) - Shape functions at quadrature points
c        xsj(*)       - Volume element at quadrature points
c        lint         - Number of quadrature points
c        nel          - Number nodes on element
c        nes          - Dimension of stress array

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

      integer   nel,nes, i,l,lint, ix(*),ixl(9)
      real*8    p(*),s(nen,*),se(*),xsj(*) ,sig(nes,*)
      real*8    shp(3,nes,*), shpl(3,9),sg(3,9), el(4,7), xg

      save

      data      ixl /1,2,3,4,5,6,7,8,9/

      if(nel.eq.8) then
        if(lint.eq.9) then
          call int2d(3,l,sg)
        else
          call int2d(2,l,sg)
        endif
        do l = 1,lint
          call meanx(ix,hr(np(44)),ndm)
          call shp2d(sg(1,l),hr(np(44)),shpl,xsj(l),ndm,9,ixl,.false.)
          do i = 1,nel
            shp(3,i,l) = shpl(3,i)
          end do ! i
        end do ! l
      elseif(nel.eq.6 .or. nel.eq.7) then
        call tint2d(7,l,el)
        do l = 1,lint
          shp(3,1,l) =  el(1,l)
          shp(3,2,l) =  el(2,l)
          shp(3,3,l) =  el(3,l)
          shp(3,4,l) = (el(1,l) + el(2,l))*0.5d0
          shp(3,5,l) = (el(2,l) + el(3,l))*0.5d0
          shp(3,6,l) = (el(3,l) + el(1,l))*0.5d0
          xsj(l)     = 1.d0
        end do ! l
      endif

c     Lumped and consistent projection routine

      do l = 1,lint

c       Compute lumped projection and assemble stress integrals

        do i = 1,nel

          xg   = shp(3,i,l)*xsj(l)
          p(i) = p(i) + xg

c         Stress projections

          s(i,1) = s(i,1) + sig(1,l)*xg
          s(i,2) = s(i,2) + sig(2,l)*xg
          s(i,3) = s(i,3) + sig(3,l)*xg
          s(i,4) = s(i,4) + sig(4,l)*xg
          s(i,5) = s(i,5) + sig(5,l)*xg
          s(i,6) = s(i,6) + sig(6,l)*xg

c         Error estimation projection

          se(i)  = se(i)  + erav*xg

        end do ! i
      end do ! l

      iste = 6

      end

      subroutine meanx(ix,xl,ndm)

      implicit   none
      integer    ndm, ix(*)
      real*8     xl(ndm,*)

c     Check midside of edges

      if(ix(5).eq.0) then
        xl(1,5) = 0.5*(xl(1,1)+xl(1,2))
        xl(2,5) = 0.5*(xl(2,1)+xl(2,2))
      endif
      if(ix(6).eq.0) then
        xl(1,6) = 0.5*(xl(1,2)+xl(1,3))
        xl(2,6) = 0.5*(xl(2,2)+xl(2,3))
      endif
      if(ix(7).eq.0) then
        xl(1,7) = 0.5*(xl(1,3)+xl(1,4))
        xl(2,7) = 0.5*(xl(2,3)+xl(2,4))
      endif
      if(ix(8).eq.0) then
        xl(1,8) = 0.5*(xl(1,4)+xl(1,1))
        xl(2,8) = 0.5*(xl(2,4)+xl(2,1))
      endif

c     Compute center node location

      xl(1,9) = 0.50d0*(xl(1,5)+xl(1,6)+xl(1,7)+xl(1,8))
     &        - 0.25d0*(xl(1,1)+xl(1,2)+xl(1,3)+xl(1,4))
      xl(2,9) = 0.50d0*(xl(2,5)+xl(2,6)+xl(2,7)+xl(2,8))
     &        - 0.25d0*(xl(2,1)+xl(2,2)+xl(2,3)+xl(2,4))
      end
