c$Id: pvalues.f,v 1.4 2006/11/21 23:32:12 rlt Exp $
      subroutine pvalues()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute values of constants used in program

c      Inputs:
c         none

c      Outputs:
c         Values are output through common 'pcommon.h'
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      one3  = 1.0d0
      one3  = one3/3.0d0              ! One  third
      two3  = one3*2.0d0              ! Two  thirds
      four3 = two3*2.0d0              ! Four thirds

      one6  = one3*0.5d0              ! One  sixth

      one9  = one3/3.0d0              ! One  ninth
      four9 = one9*4.0d0              ! Four ninths

      pi    = acos(-1.d0)             ! pi
      pi23  = (pi*2.0d0)/3.0d0        ! pi*2/3

      sqrt2 = 2.0d0
      sqrt2 = sqrt(sqrt2)             ! Square root of 2

      sqt13 = 1.0d0
      sqt13 = sqt13/3.0d0
      sqt13 = sqrt(sqt13)             ! Square root of 1/3

      sqt23 = 2.0d0
      sqt23 = sqt23/3.0d0
      sqt23 = sqrt(sqt23)             ! Square root of 2/3

      sqtp6 = 0.6d0
      sqtp6 = sqrt(sqtp6)             ! Sqare root of 0.6

      sqt48 = 4.8d0
      sqt48 = sqrt(sqt48)             ! Sqare root of 4.8

      five9 = 5.0d0
      five9 = five9/9.0d0             ! 5/9

      eight9 = 8.0d0
      eight9 = eight9/9.0d0           ! 8/9

      thty29 = 32.0d0
      thty29 = thty29/9.0d0           ! 32/9

      end
