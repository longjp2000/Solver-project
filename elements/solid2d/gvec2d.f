c$Id: gvec2d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine gvec2d(xl,ul,c,ndm,ndf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute integrals for mixed B-bar 4 node element

c      Inputs:
c         xl(ndm,*)  - Nodal coordinates for element
c         ul(ndf,*)  - Nodal solution for element
c         c          - Plane strain/axisymmetric indicator
c         ndm        - Spatial dimension of mesh
c         ndf        - Number dof/node

c      Outputs:
c         none       - Output is through common /elm2d/
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'elm2d.h'

      integer   i,n,ndm,ndf
      real*8    c,crpsi,creta,crpxe,crsum,reta,rpsi,rpxe,zeta,zpsi,zpxe
      real*8    gpsi1,geta1,gpxe1,gpsi2,geta2,gpxe2,xj0,xj1,xj2,xj0c,vol
      real*8    xl(ndm,*),ul(ndf,*),cr(4)

      save

c     Set up cr array (r for axisymmetric, 1's for plane)

      do i = 1,4
        cr(i) = 1.0d0
        if(c.ne.0.0d0) cr(i) = xl(1,i)
      end do ! i

c     Set up sums

      reta = - xl(1,1) - xl(1,2) + xl(1,3) + xl(1,4)
      rpsi = - xl(1,1) + xl(1,2) + xl(1,3) - xl(1,4)
      rpxe =   xl(1,1) - xl(1,2) + xl(1,3) - xl(1,4)
      zeta = - xl(2,1) - xl(2,2) + xl(2,3) + xl(2,4)
      zpsi = - xl(2,1) + xl(2,2) + xl(2,3) - xl(2,4)
      zpxe =   xl(2,1) - xl(2,2) + xl(2,3) - xl(2,4)
      crpsi = - cr(1) + cr(2) + cr(3) - cr(4)
      creta = - cr(1) - cr(2) + cr(3) + cr(4)
      crpxe = + cr(1) - cr(2) + cr(3) - cr(4)
      crsum = + cr(1) + cr(2) + cr(3) + cr(4)

c     Compute jacobian constants

      xj0 = (rpsi*zeta - reta*zpsi)
      xj1 = (rpsi*zpxe - rpxe*zpsi)
      xj2 = (rpxe*zeta - reta*zpxe)
      vol = xj0*crsum + (xj1*crpsi + xj2*creta)/3.0

c     Modify terms to form volumetric matrix

      xj0c = xj0*c/vol
      crsum = crsum/vol
      crpsi = crpsi/vol
      creta = creta/vol
      crpxe = crpxe/vol

c     Form g-vector constants

      gpsi1 = zeta*crsum + (zpxe*crpsi + xj1*c/vol)/3.0
      geta1 =-zpsi*crsum - (zpxe*creta - xj2*c/vol)/3.0
      gpxe1 = (zeta*creta - zpsi*crpsi)/3.0
      gpsi2 =-reta*crsum - (rpxe*crpsi)/3.0
      geta2 = rpsi*crsum + (rpxe*creta)/3.0
      gpxe2 = (rpsi*crpsi - reta*creta)/3.0

c     Form g-vector for each shape function

      g(1,1) = - gpsi1 - geta1 + gpxe1 + xj0c
      g(1,2) = + gpsi1 - geta1 - gpxe1 + xj0c
      g(1,3) = + gpsi1 + geta1 + gpxe1 + xj0c
      g(1,4) = - gpsi1 + geta1 - gpxe1 + xj0c
      g(2,1) = - gpsi2 - geta2 + gpxe2
      g(2,2) = + gpsi2 - geta2 - gpxe2
      g(2,3) = + gpsi2 + geta2 + gpxe2
      g(2,4) = - gpsi2 + geta2 - gpxe2

c     Compute trace of strain (or its increment)

      trep = 0.0d0
      do n = 1,4
        trep = trep + g(1,n)*ul(1,n) + g(2,n)*ul(2,n)
      end do ! n

      end
