c$Id: etime.f,v 1.1 2006/11/20 20:34:56 rlt Exp $
      real function etime(tt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Timing information for output file

c      Inputs:
c         none

c      Outputs:
c         tt(*)   - CPU and System time (no system time for DOS)
c         etime   - Total CPU + System time
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'etime1.h'

      real      tt(2)
      integer*2 ihr, imin, isec, ihth

      save

      call gettim (ihr,imin,isec,ihth)
      tt(2) = 0.0
      tt(1) = 0.01*ihth + isec + 60.*(imin + 60.*ihr) - tim0
      if(tt(1) .lt. 0.0) then
        tt(1) = tt(1) + 86400.0
        tim0  = tim0  - 86400.0
      endif
      etime = tt(1) + tt(2)

      end
