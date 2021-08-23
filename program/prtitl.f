c$Id: prtitl.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine prtitl(prtf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output title to file/screen

c      Inputs:
c         prtf      - Flag, output title if true

c      Outputs:
c         None      - Outputs to file/screen
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'iofile.h'

      logical   prtf

      save

      if(prtf) then
        write(iow,2000) o,head
        if(ior.lt.0) then
          write(*,2000) o,head
        endif
      endif

2000  format(/a1,19a4,a2/1x)

      end
