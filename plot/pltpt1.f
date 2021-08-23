c$Id: pltpt1.f,v 1.1 2006/11/20 20:33:12 rlt Exp $
      subroutine pltpt1(iel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set 3-d Plot Sequence for 1-node point elements

c      Inputs:
c         iel       - Element type number

c      Outputs:
c         none      - Output through common block data
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata5.h'
      include  'pdata6.h'

      integer   iel

      save

c     Set number of points

      if(iel.gt.0) then

        inord(iel)  = -1

      elseif(iel.lt.0) then

        exord(-iel) = -1

      endif

      end
