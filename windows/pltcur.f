c$Id: pltcur.f,v 1.2 2006/12/09 21:17:29 rlt Exp $
      subroutine pltcur()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Replace MSFLIB by DFLIB                          08/12/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Turn cursor on for text inputs

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      use DFLIB

      implicit none

      integer  status

      status = displaycursor($GCURSORON)

      end
