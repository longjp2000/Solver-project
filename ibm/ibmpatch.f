c$Id: ibmpatch.f,v 1.1 2006/11/20 20:33:52 rlt Exp $

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Patch routine for IBM AIX computers
c-----[--.----+----.----+----.-----------------------------------------]

      function etime(tt)

c     Timing on the IBM RS/6000

      real etime, tt(2)

      integer mclock

      tt(1) = real( mclock() )/100.0
      tt(2) = 0.0d0
      etime = tt(1) + tt(2)

      end

      subroutine flush(io)

      integer io

      end

      subroutine getlog(uname)

      character uname*8
      uname = 'FEAPuser'

      end

      subroutine fdate (cdate)

      character cdate*24
      cdate = 'Undated Output'

      end

      integer function  lnblnk (argv)

      integer   i
      character argv*14

      do i = 14,1,-1
        if(argv(i:i).ne.' ') then
          lnblnk = i
          return
        endif
      end do ! i
      lnblnk =  0

      end
