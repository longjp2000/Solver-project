c$Id: uglobl.f,v 1.1 2006/11/20 20:33:27 rlt Exp $
      subroutine uglobl(type,td)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Dummy user global parameter routine

c      Inputs:
c         type(2)   - Character array describing user global command
c         td(5)     - Real array of data for global command

c      Outputs:
c         N.B. Users must provide output via common blocks, etc.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      character type(2)*15
      real*8    td(5)

      end
