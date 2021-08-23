c$Id: rprint.f,v 1.1 2006/11/20 20:33:12 rlt Exp $
      subroutine rprint(dr,ndf,nfl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Computes range and profile of plot value requested

c      Inputs:
c         dr(ndf,*) - Values for plot (N.B. checks dr(1,i) values)
c         ndf       - Dimension of dr-array

c      Outputs:
c         nfl       - Returns -nfl if all values are < 1.0d-08
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'fdata.h'
      include  'iofile.h'
      include  'pdata2.h'
      include  'pointer.h'
      include  'prmptd.h'
      include  'rpdata.h'
      include  'comblk.h'

      include  'p_int.h'

      integer   ndf,nfl, i,n
      real*8    drv,rs, dr(ndf,*),pr(9)

      save

c     Compute profile of values

      do n = 1,9
        pr(n) = 0.0d0
      end do ! n
      drv   = 100.d0/numnp
      rmx   = dr(1,1)
      rmn   = dr(1,1)
      fp(1) = np(190) - 1
      do n = 1,numnp
        if(mr(fp(1)+n) .ge. 0) then
          rmx = max(rmx,dr(1,n))
          rmn = min(rmn,dr(1,n))
        endif
      end do ! n

c     Check range for contour outputs

      if(abs(rmx-rmn).gt.1.d-5*max(abs(rmx),abs(rmn))) then
        do n = 1,numnp
          if(mr(fp(1)+n) .ge. 0) then
            rs = (dr(1,n) - rmn)/(rmx - rmn)
            do i = 1,9
              if(rs.ge.0.1d0*i) pr(i) = pr(i) + drv
            end do ! i
          endif
        end do ! n
      else
        if(max(abs(rmx),abs(rmn)).gt.1.d-8) then
          rmn = rmn - 0.0999*max(abs(rmx),abs(rmn))
          rmx = rmx + 0.1001*max(abs(rmx),abs(rmn))
        else
          rmn = -1.0d-08 + 1.0d-12
          rmx =  1.0d-08 + 1.0d-12
          nfl = -abs(nfl)
        endif
      endif

c     Output range for min/max

      if(prompt) then
        if(pfr) write(iow,2000) rmn,rmx
        if(ior.lt.0 .and. .not.defalt) then
          write(*,2000) rmn,rmx
        endif
      endif

2000  format('    Minimum is ',1p,e10.2,' Maximum is ',1p,e10.2:/
     &  22x,'10%   20%   30%   40%   50%   60%   70%   80%   90%'/
     &       '    Profile above is:',9f6.1)

      end
