c$Id: damp3d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine damp3d(d,xl,s,ndf,ndm,nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute damping matrix for 3-d brick elements

c      Inputs:
c         d(*)      - Material set parameters
c         xl(ndm,*) - Nodal coordinates for element
c         ndf       - Number dof/node
c         ndm       - Spatial dimension of mesh
c         nst       - Size of element arrays

c      Outputs:
c         s(nst,*)  - Damping matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'pmod2d.h'

      integer   ndf,ndm,nst, i,ii,i1, jj,j1, l,lint, ord
      real*8    d(*),xl(ndm,*),s(nst,nst)
      real*8    shp(4,27),sg(4,125),sv(5,16), xsj,dv, aj1

      save

c     Compute mass matrix

      if(nel.eq.4) then
        ord = 1
        if(nint(d(182)).gt.0) then
          call tint3dn(nel,lint,sv)
        else
          l =  2
          call tint3d (l,lint,sv)
        endif
      elseif(nel.eq.10) then
        ord =  2
        l   =  3
        call tint3d(l,lint,sv)
      elseif(nel.eq.11) then
        ord =  2
        if(nint(d(182)).gt.0) then
          call tint3dn(nel,lint,sv)
        else
          l =  4
          call tint3d(l,lint,sv)
        endif
      else
        ord = 0
        if(nint(d(182)).gt.0) then
          call int3dn(nel,lint,sg)
        else
          l = nint(d(5))
          call int3d(l,lint,sg)
        endif
      endif

c     Loop over quadrature points

      do l = 1,lint

c       Compute shape functions

        if(ord.eq.0) then
          call shp3d(sg(1,l),xsj,shp,xl,ndm,nel)
          dv = sg(4,l)*xsj*d(70)
        else
          call tetshp(sv(1,l),xl,ndm,ord,xsj,shp)
          dv = sv(5,l)*xsj*d(70)
        endif

c       Compute db = rho*shape*dv

        j1 = 0
        do jj = 1,nel

c         Compute consistent damping matrix

          aj1 = shp(4,jj)*dv
          i1  = 0
          do ii = 1,nel
            do i = 1,3
              s(i1+i,j1+i) = s(i1+i,j1+i) + shp(4,ii)*aj1
            end do ! i
            i1 = i1 + ndf
          end do ! ii
          j1 = j1 + ndf
        end do ! jj
      end do ! l

      end
