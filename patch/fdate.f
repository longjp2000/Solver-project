c$Id: fdate.f,v 1.1 2006/11/20 20:44:57 rlt Exp $
      subroutine fdate(cdate)

c      * * F E A P * * A Finite Element Analysis Program
c                        -      -       -        -
c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Get 'date' for labels in outputs -- Fortran 90 compatible

c      Inputs: none

c      Outputs: cdate
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      character  cdate*24

      call date_and_time(date = cdate)

      end
