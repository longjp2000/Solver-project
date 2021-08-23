c$Id: pltq16.f,v 1.1 2006/11/20 20:33:12 rlt Exp $
      subroutine pltq16(iel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set 3-d Plot Sequence for 16-node cubic quadrilaterals

c      Inputs:
c         iel       - Element type number

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata5.h'
      include  'pdata6.h'

      integer   iel

      save

c     Set number of points

      if(iel.gt.0) then

        inord(iel)    = 13

c       Set plot sequence

        ipord( 1,iel) = 1
        ipord( 2,iel) = 5
        ipord( 3,iel) = 6
        ipord( 4,iel) = 2
        ipord( 5,iel) = 7
        ipord( 6,iel) = 8
        ipord( 7,iel) = 3
        ipord( 8,iel) = 9
        ipord( 9,iel) = 10
        ipord(10,iel) = 4
        ipord(11,iel) = 11
        ipord(12,iel) = 12
        ipord(13,iel) = 1

      elseif(iel.lt.0) then

        exord(-iel)    = 13

c       Set plot sequence

        epord( 1,-iel) = 1
        epord( 2,-iel) = 5
        epord( 3,-iel) = 6
        epord( 4,-iel) = 2
        epord( 5,-iel) = 7
        epord( 6,-iel) = 8
        epord( 7,-iel) = 3
        epord( 8,-iel) = 9
        epord( 9,-iel) = 10
        epord(10,-iel) = 4
        epord(11,-iel) = 11
        epord(12,-iel) = 12
        epord(13,-iel) = 1

      endif

      end
