c$Id: geopar2.f,v 1.1 2006/11/20 20:32:46 rlt Exp $
      subroutine geopar2 (ndm,ndf,x,u,nel1,nod1,ix1,ix2,cxs,
     &                    cp0,ch1,ch2,ch3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise          August 27, 1997            1.0

c      Acronym: GEOmetrical PARameters subroutine

c      Purpose: Compute and store all geometrical data

c      Inputs :
c         ndm     - Space dimension of mesh
c         ndf     - Number dof/node
c         x(*)    - nodal coordinates
c         u(*)    - Current nodal solution vectors
c         nel1    - Current contact element of surf. 1
c         nod1    - Current node of contact element nel1
c         ix1(*)  - Element nodal connection list for surface 1
c         ix2(*)  - Element nodal connection list for surface 2
c         cxs(*)  - Coordinates of slavenode S
c         cp0(*)  - Contactpair control data
c         ch1(*)  - Contact history variables (old)
c         ch2(*)  - Contact history variables (current)

c      Outputs:
c         ch2(*)  - Contact history variables (current)
c         ch3(*)  - Contact history variables (static)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_geom.h'
      include  'c_keyh.h'
      include  'c_pair.h'
      include  'c_tole.h'
      include  'counts.h'

      integer   ndm,ndf,nel1,nod1,ix1(dnope1,*),ix2(dnope2,*)
      integer   masts,lnc,istgt,n1,n2,ns,nb,na,n1s,n2s,ifal
      integer   istgno,ke,mastso
      real*8    x(ndm,*),u(ndf,*),cxs(*),cp0(nr0,n0c3:*)
      real*8    ch1(*),ch2(*),ch3(*),cx1(2),cx2(2),cxa(2),cxb(2),cxc(2)
      real*8    d21,s21,c21,gt,csi,gn,area,gnbar,s21c,c21c,csic
      real*8    dtd,d21o,d21oi,csio,gto,cx1s(2),cx2s(2),tdo

      save

      call cdebug0 ('      geopar2',-1)

c     Get main geometrical variables

      masts = ch2(p1(1))
      if(masts.eq.0) then
        s21   = 0.0d0
        c21   = 1.0d0
        d21   = 1.0d0
        csi   = 0.0d0
        gn    = 0.0d0
        gt    = 0.0d0
        area  = 1.0d0
        go to 100
      endif

        lnc = ch2(p1(2))
      istgt = ch2(p1(3))

c     Find geometrical basic parameters

      n1 = ix2(1,masts)
      n2 = ix2(2,masts)
      cx1(1) = x(1,n1) + u(1,n1)
      cx1(2) = x(2,n1) + u(2,n1)
      cx2(1) = x(1,n2) + u(1,n2)
      cx2(2) = x(2,n2) + u(2,n2)

c     Tangent vector t
c     REMARK - s21 and c21 normalized at end to minimize numerical errors

      s21 = (cx2(2) - cx1(2))
      c21 = (cx2(1) - cx1(1))

c     Length of master segment

      d21 = sqrt(c21*c21 + s21*s21)

c     Normal projection of S on to 21
c     REMARK - positive when gap is open

      gn = (cxs(1) - cx1(1)) * s21 - (cxs(2) - cx1(2)) * c21
      gn = gn / d21

c     Tangential projection of S on to 21

      gt = (cxs(1) - cx1(1)) * c21 + (cxs(2) - cx1(2)) * s21
      gt = gt / d21

c     Normalized tangential projection

      csi = gt / d21

c     Now normalize s21 and c21

      s21 = s21 / d21
      c21 = c21 / d21

c     Slavenode in special positions -> update parameters

      if (istgt.gt.1) then

c       Set csi to closest node

        if (lnc.eq.1) then
          csic = 0.d0
          cxc(1) = cx1(1)
          cxc(2) = cx1(2)
        else
          csic = 1.d0
          cxc(1) = cx2(1)
          cxc(2) = cx2(2)
        endif

c       Compute new parameters

        gnbar = sqrt((cxs(1) - cxc(1))*(cxs(1) - cxc(1)) +
     &               (cxs(2) - cxc(2))*(cxs(2) - cxc(2)))

        gn = gnbar * sign(1.d0,gn)

c       If penetration and case 4 or 5 check if inside tolerance

        if ((istgt.eq.4) .and. (gn.lt.0.d0)) then

c         If outside tolerance change sign to consider it open

          if (abs(csi*d21).gt.tlopen) then
            gn = abs(gn)
          endif

        else if ((istgt.eq.5) .and. (gn.lt.0.d0)) then

c         If outside tolerance change sign to consider it open

          if (abs((csi-1.d0)*d21).gt.tlopen) then
            gn = abs(gn)
          endif
        endif

        s21c =  (cxs(1) - cxc(1)) / gn
        c21c = -(cxs(2) - cxc(2)) / gn

c       Store private variables

        ch2(p1(13)) = s21c
        ch2(p1(14)) = c21c
        ch2(p1(15)) = csic

      endif

c     Set tangential displacement if requested

      if (iffric.eq.1) then
        istgno = ch1(p1(4))
        mastso = ch1(p1(1))
        d21o   = ch1(p1(7))
        csio   = ch1(p1(8))
        gto    = ch1(p1(10))
        tdo    = ch1(p1(16))

c       Check tangential solution mode

c       Master segment not changed from previous step

        if (masts.eq.mastso) then
          dtd = (csi-csio)*d21o
          d21oi = d21o

c       Master changed -> compute total sliding,
c       Set old segm. length = to new one for iteration

        elseif (masts.gt.mastso) then
          dtd = (1.d0 - csio)*d21o + csi*d21
c         do ke = mastso+1,mastso-1
          do ke = mastso+1,masts
            n1s = ix2(1,ke)
            n2s = ix2(2,ke)
            cx1s(1) = x(1,n1s) + u(1,n1s)
            cx1s(2) = x(2,n1s) + u(2,n1s)
            cx2s(1) = x(1,n2s) + u(1,n2s)
            cx2s(2) = x(2,n2s) + u(2,n2s)
            dtd = dtd + sqrt((cx2s(1)-cx1s(1))*(cx2s(1)-cx1s(1)) +
     &                       (cx2s(2)-cx1s(2))*(cx2s(2)-cx1s(2)))
            d21oi = d21
          end do
        else
          dtd = -(csio*d21o + (1.d0 - csi)*d21)
          do ke = masts+1,mastso-1
            n1s = ix2(1,ke)
            n2s = ix2(2,ke)
            cx1s(1) = x(1,n1s) + u(1,n1s)
            cx1s(2) = x(2,n1s) + u(2,n1s)
            cx2s(1) = x(1,n2s) + u(1,n2s)
            cx2s(2) = x(2,n2s) + u(2,n2s)
            dtd = dtd - sqrt((cx2s(1)-cx1s(1))*(cx2s(1)-cx1s(1)) +
     &                       (cx2s(2)-cx1s(2))*(cx2s(2)-cx1s(2)))
            d21oi = d21
          end do
        endif

c       Store private parameters

        ch2(p1(12)) = dtd
        ch3(p3(1)) = d21oi

c       Tangential initialization

      elseif (iffric.eq.2) then
        ch2(p1(17)) = csi
      endif

c     Set contact area
c     Set left and right nodes with respect to slave one

      if(dnope1.gt.1) then
        if (nod1.eq.1) then
          na = ix1(dnope1-1,nel1)
          if (na.eq.0) then
            na = nel1
          endif
          na = ix1(1,na)
          nb = ix1(2,nel1)
        else
          na = ix1(1,nel1)
          nb = ix1(dnope1  ,nel1)
          if (nb.eq.0) then
            nb = nel1
            endif
          nb = ix1(2,nb)
        endif

c       If contact area is linearized pick current coord
c       REMARK - still to do --> constant contact area

        ifal = 0

c       Constant contact area (original)

        if ((ifal.eq.0) .or. (ndf.eq.0)) then
          cxa(1) = x(1,na)
          cxa(2) = x(2,na)
          cxb(1) = x(1,nb)
          cxb(2) = x(2,nb)

c         Correct node slave coordinates

          ns = ix1(nod1,nel1)
          cxs(1) = x(1,ns)
          cxs(2) = x(2,ns)

c         Variable contact area (current)

        else
          cxa(1) = x(1,na) + u(1,na)
          cxa(2) = x(2,na) + u(2,na)
          cxb(1) = x(1,nb) + u(1,nb)
          cxb(2) = x(2,nb) + u(2,nb)
        endif

c       Compute area

        area = (sqrt((cxa(1)-cxs(1))*(cxa(1)-cxs(1))   +
     &               (cxa(2)-cxs(2))*(cxa(2)-cxs(2)))  +
     &          sqrt((cxb(1)-cxs(1))*(cxb(1)-cxs(1))   +
     &               (cxb(2)-cxs(2))*(cxb(2)-cxs(2)))) * 0.5d0
      else
        area = 1.d0
      endif

c     Store variables

100   ch2(p1(5))  = s21
      ch2(p1(6))  = c21
      ch2(p1(7))  = d21
      ch2(p1(8))  = csi
      ch2(p1(9))  = gn
      ch2(p1(10)) = gt
      ch2(p1(11)) = area

      end
