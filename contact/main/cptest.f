c$Id: cptest.f,v 1.1 2006/11/20 20:32:46 rlt Exp $
      subroutine cptest(cp0,npair,surfl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact Pair TEST

c      Purpose: Test for contact types which require surface generations

c      Inputs:
c         cp0(*)  - Contact pair parameters
c         npair   - Pair set

c      Outputs:
c         surfl   - Surface generation flag (true for surfaces)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'c_0.h'
      include 'c_comnd.h'

      logical  surfl
      integer  npair
      real*8   cp0(nr0,n0c3:nc03,*)

      if(.not.surfl) then
        surfl = nint(cp0(12,-1,npair)).eq.1
     &    .and. nint(cp0(13,-1,npair)).ne.2
      end if

      end
