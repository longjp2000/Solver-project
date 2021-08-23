c$Id: pinitm.f,v 1.1 2006/11/20 20:33:40 rlt Exp $
      subroutine pinitm()

c      * * F E A P * * A Finite Element Analysis Program
c                        -      -       -        -
c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Initialize memory use for compilers which do not support
c               malloc and free.

c      Inputs:  none

c      Outputs: none
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'allotd.h'
      include   'psize.h'

      integer    mrmax
      parameter (mrmax = 80 000 000 )

      integer    mr
      common     mr( mrmax )

c     Set maximum memory available to solve problems

      mmax      = 1
      maxm      = mrmax
      ipoint(1) = 1

      end
