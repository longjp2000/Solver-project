c$Id: pinitm.f,v 1.1 2006/11/20 20:33:21 rlt Exp $
      subroutine pinitm()

c      * * F E A P * * A Finite Element Analysis Program
c                        -      -       -        -
c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Initialize memory use

c      Inputs:  none

c      Outputs: none
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit none

      include 'comblk.h'

      save

c     Initialize pointers for real and imaginary storage

      call initm(hr(1),mr(1))

      end
