c$Id: fbody2d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine fbody2d(d,xl,ul, r,s, ndm,ndf,nst, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute body force loads
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'elbody.h'
      include   'eldata.h'

      integer    ndm,ndf,nst,isw, ii,ll,lint
      real*8     d(*),xl(ndm,*),ul(ndf,*),r(ndf,*),s(nst,nst)
      real*8     shp(2,3), sg(2,3), body(3)
      real*8     le,dx,dy,pn

c     Set body force loading factors

      if(isw.eq.15) then
        body(1) = bodyf(1)
        body(2) = bodyf(2)
        le      = sqrt((xl(1,2)-xl(1,1))**2
     &                 + (xl(2,2)-xl(2,1))**2)*0.5d0
        r(1,1)  = r(1,1) + body(1)*le
        r(2,1)  = r(2,1) + body(2)*le
        r(1,2)  = r(1,2) + body(1)*le
        r(2,2)  = r(2,2) + body(2)*le
      else
        call sbodyf(d, body)

c       2-node element

        if(nel.eq.2) then

          pn     = 0.5d0*dm*d(10)
          dx     = (xl(1,2) - xl(1,1))
     &           + (ul(1,2) - ul(1,1))*d(68)
          dy     = (xl(2,2) - xl(2,1))
     &           + (ul(2,2) - ul(2,1))*d(68)
          le     = sqrt((xl(1,2)-xl(1,1))**2
     &                + (xl(2,2)-xl(2,1))**2)*0.5d0
          r(1,1) = r(1,1)  + body(1)*le - pn*dy
          r(2,1) = r(2,1)  + body(2)*le + pn*dx
          r(1,2) = r(1,2)  + body(1)*le - pn*dy
          r(2,2) = r(2,2)  + body(2)*le + pn*dx

c         Tangent for follower forces (unsymmetric)

          if(d(68).gt.0.0d0) then
            s(1    ,2    ) = -pn
            s(1    ,2+ndf) =  pn
            s(2    ,1    ) =  pn
            s(2    ,1+ndf) = -pn
            s(1+ndf,2    ) = -pn
            s(1+ndf,2+ndf) =  pn
            s(2+ndf,1    ) =  pn
            s(2+ndf,1+ndf) = -pn
          endif

c       3-node element

        else
          lint = nel
          if(nint(d(182)).gt.0) then
            call int1dn(lint, sg)
          else
            call int1d(lint, sg)
          endif
          do ll = 1,lint
            call shp1d(sg(1,ll),xl,shp,ndm,nel,dx)
            dx = sg(2,ll)*dx
            do ii = 1,nel
              r(1,ii)   = r(1,ii) + body(1)*shp(2,ii)*dx
              r(2,ii)   = r(2,ii) + body(2)*shp(2,ii)*dx
            end do ! ii
          end do ! ll
        endif
      endif

      end
