c$Id: strn3m.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine strn3m(shp,xl,ul,theta,ndm,ndf,nel,nen,eps)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute mixed strain for near incompressible formulation.

c      Inputs:
c        shp(4,*)      = Shape functions
c        xl(ndm,nen)   = Nodal coordinates
c        ul(ndf,nen,*) = Nodal solution parameters
c        theta         = Volume change (mixed form)
c        ndm           = Spatial dimension of mesh
c        ndf           = DOF/node (max)
c        nel           = Number nodes on element (4 or 9)
c        nen           = Max nodes/element (dimension uses only)

c      Outputs:
c        eps(9,*)      = Mixed strain at point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'elcoor.h'
      include  'pconstant.h'
      include  'pmod2d.h'

      integer   ndm,      ndf,         nel,         nen,     k
      real*8    theta,    dtheta
      real*8    shp(4,*), xl(ndm,nen), ul(ndf,nen), eps(9,*)

c     Compute strain tensor for constitutive equations

      do k = 1,6
        eps(k,1) = 0.0d0
      end do ! k
      do k = 1,3
        xref(k) = 0.0d0
        xcur(k) = 0.0d0
      end do ! k

      do k = 1,nel
        eps(1,1)  = eps(1,1) + shp(1,k)*ul(1,k)
        eps(2,1)  = eps(2,1) + shp(2,k)*ul(2,k)
        eps(3,1)  = eps(3,1) + shp(3,k)*ul(3,k)
        eps(4,1)  = eps(4,1) + shp(1,k)*ul(2,k)
     &                       + shp(2,k)*ul(1,k)
        eps(5,1)  = eps(5,1) + shp(2,k)*ul(3,k)
     &                       + shp(3,k)*ul(2,k)
        eps(6,1)  = eps(6,1) + shp(3,k)*ul(1,k)
     &                       + shp(1,k)*ul(3,k)
        xref(1)   = xref(1) + shp(4,k)*xl(1,k)
        xref(2)   = xref(2) + shp(4,k)*xl(2,k)
        xref(3)   = xref(3) + shp(4,k)*xl(3,k)
        xcur(1)   = xcur(1) + shp(4,k)*ul(1,k)
        xcur(2)   = xcur(2) + shp(4,k)*ul(2,k)
        xcur(3)   = xcur(3) + shp(4,k)*ul(3,k)
      end do ! k

c     Set current coords

      do k = 1,3
        xcur(k) = xcur(k) + xref(k)
      end do ! k

c     Correct strains and incremental strains for mixed formulation

      dtheta = one3*(theta - eps(1,1)-eps(2,1)-eps(3,1))
      eps(1,1) = eps(1,1) + dtheta
      eps(2,1) = eps(2,1) + dtheta
      eps(3,1) = eps(3,1) + dtheta

      end
