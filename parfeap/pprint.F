c$Id: pprint.F,v 1.1 2006/11/21 16:44:39 rlt Exp $
      subroutine pprint(string)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Print prompt for interactive inputs

c      Inputs:
c         string - Character string to output

c      Outputs:
c         Writes prompt string to screen
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "include/finclude/petsc.h"
#     include   "pfeapa.h"
#     include   "pfeapb.h"

      character  string*(*)
      integer    l,le

      le = len(string)
      do l = le,1,-1
        if(string(l:l).ne.' ') exit
      end do ! l
      if(pfeap_on) then
        if(rank.eq.0) then
          write(*,2000) string(1:l),' '
        endif
      else
        write(*,2000) string(1:l)
      endif

c     Format

2000  format(a,' ')
2001  format(a,' ',$)

      end
