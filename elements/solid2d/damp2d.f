c$Id: damp2d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine damp2d(d,xl,ix,s,ndf,ndm,nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute damping matrix for plane/axisymmetric problems

c      Inputs:
c         d(*)      - Material set parameters
c         xl(ndm,*) - Nodal coordinates for element
c         ix(*)     - Element nodal connections
c         ndf       - Number dof/node
c         ndm       - Spatial dimension of mesh
c         nst       - Size of element arrays

c      Outputs:
c         s(nst,*)  - Damping matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'pmod2d.h'

      logical   quad
      integer   ndf,ndm,nst, i,ii,i1, jj,j1, l,lint, ncp, ix(*)
      real*8    d(*),xl(ndm,*),s(nst,nst), shp(3,9)
      real*8    sg(3,16),el(4,7), xsj,dv, aj1,xx

      save

c     Number of mass components

      if(stype.eq.8) then
        ncp = 3
      else
        ncp = 2
      endif

c     Compute mass matrix

      if(nel.eq.3) then
        if(nint(d(182)).gt.0) then
          call tint2dn(nel,lint,el)
        else
          l =  1
          call tint2d (l,lint,el)
        endif
        quad = .false.
      elseif(nel.eq.6 .or. nel.eq.7) then
        if(nint(d(182)).gt.0) then
          call tint2dn(nel,lint,el)
        else
          l =  7
          call tint2d (l,lint,el)
        endif
        quad = .false.
      else
        quad = .true.
        if(nint(d(182)).gt.0) then
          call int2dn(nel,lint,sg)
        else
          call int2d (l,lint,sg)
        endif
      endif

      do l = 1,lint

c       Compute shape functions

        if(quad) then
          call shp2d(sg(1,l),xl,shp,xsj,ndm,nel,ix,.false.)
          dv = sg(3,l)*abs(xsj)*d(70)
        else
          call trishp(el(1,l),xl,ndm,nel-4,xsj,shp)
          dv = el(4,l)*abs(xsj)*d(70)
        endif
        if(stype.eq.3 .or. stype.eq.8) then
          xx = 0.0d0
          do jj = 1,nel
            xx = xx + shp(3,jj)*xl(1,jj)
          end do ! jj
          dv = dv*xx
        end if

c       Compute db = c*shape*dv

        j1 = 0
        do jj = 1,nel

c         Compute damping matrix

          aj1 = shp(3,jj)*dv
          i1  = 0
          do ii = 1,nel
            do i = 1,ncp
              s(i1+i,j1+i) = s(i1+i,j1+i) + shp(3,ii)*aj1
            end do ! i
            i1 = i1 + ndf
          end do ! ii
          j1 = j1 + ndf
        end do ! jj
      end do ! l

      end
