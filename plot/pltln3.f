c$Id: pltln3.f,v 1.1 2006/11/20 20:33:12 rlt Exp $
      subroutine pltln3(iel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set 3-d Plot Sequence for 3-node line elements

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

        inord(iel)    = 5

c       Set plot sequence

        ipord( 1,iel) = 1
        ipord( 2,iel) = 3
        ipord( 3,iel) = 2
        ipord( 4,iel) = 3
        ipord( 5,iel) = 1

      elseif(iel.lt.0) then

        exord(-iel)    = 5

c       Set plot sequence

        epord( 1,-iel) = 1
        epord( 2,-iel) = 3
        epord( 3,-iel) = 2
        epord( 4,-iel) = 3
        epord( 5,-iel) = 1

      endif

      end
