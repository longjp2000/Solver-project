c$Id: umshlib.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine umshlib(i,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Interface for user mesh commands

c      Inputs:
c         i      - Command number
c         prt    - Flag, output if true

c      Outputs:
c         None   - Users are responsible for providing outputs in
c                  umeshi routines
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   prt
      integer   i

      save

      if(i.eq.1) then
        call umesh1(prt)
      elseif(i.eq.2) then
        call umesh2(prt)
      elseif(i.eq.3) then
        call umesh3(prt)
      elseif(i.eq.4) then
        call umesh4(prt)
      elseif(i.eq.5) then
        call umesh5(prt)
      elseif(i.eq.6) then
        call umesh6(prt)
      elseif(i.eq.7) then
        call umesh7(prt)
      elseif(i.eq.8) then
        call umesh8(prt)
      elseif(i.eq.9) then
        call umesh9(prt)
      elseif(i.eq.10) then
        call umesh0(prt)
      endif

      end
