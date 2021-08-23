c$Id: mass2d.f,v 1.2 2006/12/12 02:16:29 rlt Exp $
      subroutine mass2d(d,xl,ul,ix,s,p,ndf,ndm,nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1.  Change 'nel.eq.4' to 'nel.eq.3'                 11/12/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute mass matrix for plane and axisymmetric problems

c      Inputs:
c         d(*)      - Material set parameters
c         xl(ndm,*) - Nodal coordinates for element
c         ul(ndm,*) - Nodal displacements for element
c         ix(*)     - Element nodal connections
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

      logical   quad
      integer   ndf,ndm,nst, i,ii,i1, jj,j1, l,lint, ncp, ix(*)
      real*8    d(*),xl(ndm,*),ul(ndf,*),s(nst,nst),p(nst), shp(3,16)
      real*8    sg(3,25),el(4,7), xsj,dv, aj1,xx,uu,lfac,cfac
      real*8    rfac(3)

      save

c     Number of mass components

      if(stype.eq.8) then
        ncp = 3
      else
        ncp = 2
      endif

c     Set default radius type

      do i = 1,ncp
        rfac(i) = 1.d0
      end do ! i

c     Compute mass matrix

      if(nel.eq.3) then
        if(d(182).gt.0.0d0) then
          call tint2dn(nel,lint,el)
        else
          l = 1
          call tint2d (l,lint,el)
        endif
        quad = .false.
      elseif(nel.eq.6 .or. nel.eq.7) then
        if(d(182).gt.0.0d0) then
          call tint2dn(nel,lint,el)
        else
          l = 7
          call tint2d (l,lint,el)
        endif
        quad = .false.
      else
        quad = .true.
        if(nint(d(182)).gt.0) then
          call int2dn(nel,lint,sg)
        else
          l = nint(d(5))
          call int2d (l,lint,sg)
        endif
      endif
      cfac = d(7)
      lfac = 1.d0 - cfac

      do l = 1,lint

c       Compute shape functions

        if(quad) then
          call shp2d(sg(1,l),xl,shp,xsj,ndm,nel,ix,.false.)
          dv = sg(3,l)*abs(xsj)*d(4)
        else
          call trishp(el(1,l),xl,ndm,nel-4,xsj,shp)
          dv = el(4,l)*abs(xsj)*d(4)
        endif
        if(stype.eq.3 .or. stype.eq.8) then
          xx = 0.0d0
          do jj = 1,nel
            xx = xx + shp(3,jj)*xl(1,jj)
          end do ! jj
          dv = dv*xx
        end if

c       Check for finite deformation torsion case

        if(dtype.ne.0 .and. stype.eq.8) then
          uu = 0.0d0
          do jj = 1,nel
            uu = uu + shp(3,jj)*ul(1,jj)
          end do ! jj
          rfac(3) = (xx +uu)**2
        endif

c       Compute db = rho*shape*dv

        j1 = 0
        do jj = 1,nel

c         Compute lumped mass matrices

          aj1 = shp(3,jj)*dv
          do i = 1,ncp
            p(j1+i)      = p(j1+i)      + rfac(i)*aj1
            s(j1+i,j1+i) = s(j1+i,j1+i) + rfac(i)*aj1*lfac
          end do ! i

c         Compute consistent mass matrix

          aj1 = aj1*cfac
          i1  = 0
          do ii = 1,nel
            do i = 1,ncp
              s(i1+i,j1+i) = s(i1+i,j1+i) + rfac(i)*shp(3,ii)*aj1
            end do ! i
            i1 = i1 + ndf
          end do ! ii
          j1 = j1 + ndf
        end do ! jj
      end do ! l

      end
