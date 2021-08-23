c$Id: getlog.f,v 1.1 2006/11/20 20:44:57 rlt Exp $
      subroutine getlog(uname)

c      * * F E A P * * A Finite Element Analysis Program
c                        -      -       -        -
c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: User name for header in PostScript files

c      Inputs: none

c      Outputs: User name (dummy -  returned)
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      character  uname*8

      uname = ' - ' !  Insert user name

      end
