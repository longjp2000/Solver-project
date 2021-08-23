c$Id: lnblnk.f,v 1.1 2006/11/20 20:44:57 rlt Exp $
      integer function  lnblnk (argv)

c      * * F E A P * * A Finite Element Analysis Program
c                        -      -       -        -
c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Determine last non-blank character in name

c      Inputs:  argv   = character name (14 characters max)

c      Outputs: lnblnk = length of name.
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit  none

      integer   i
      character argv*14

c     Find first non-blank character in name

      do i = 14,1,-1
        if(argv(i:i).ne.' ') then
          lnblnk = i
          return
        endif
      end do ! i
      lnblnk =  0

      end
