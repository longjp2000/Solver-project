c$Id: ireshis.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine ireshis(intel, n1,n2)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Initialize t_n+1 history variables for interfaces

c      Inputs:
c         intel(8,*)  - Interface element/history pointer array
c         n1          - Pointer in ix to t_n   data
c         n2          - Pointer in ix to t_n+1 data

c      Outputs:
c         none        - Output is retained in blank common
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'hdata.h'
      include   'pointer.h'
      include   'comblk.h'

      integer    intel(8,*), n1,n2, n, nh

      save

      n = 1
      do while (intel(1,n).ne.0)
        nh1 = np(214) + intel(5+n1,n)
        nh2 = np(214) + intel(5+n2,n)
        if(nh2.ne.nh1) then
          do nh = 1, abs(nh2 - nh1)
            hr(nh+nh2) = hr(nh+nh1)
          end do ! nh
        endif
        n = n + 1
      end do ! while

      end
