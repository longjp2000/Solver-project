      subroutine pdriver()

c      * * F E A P * * A Finite Element Analysis Program
c                        -      -       -        -
c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Set graphics driver for UNIX/Linux

c      Inputs: none

c      Outputs: none
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'pdata2.h'

c     Set driver for UNIX/Linux

      idev = 1

      end
