c$Id: upltlib.f,v 1.1 2006/11/20 20:33:12 rlt Exp $
      subroutine upltlib(i,ct,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Interface for user plot commands

c      Inputs:
c         i      - Command number
c         prt    - Flag, output if true

c      Outputs:
c         None   - Users are responsible for providing outputs in
c                  uploti routines
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   prt
      integer   i
      real*8    ct(3)

      save

      if(i.eq.1) then
        call uplot1(ct,prt)
      elseif(i.eq.2) then
        call uplot2(ct,prt)
      elseif(i.eq.3) then
        call uplot3(ct,prt)
      elseif(i.eq.4) then
        call uplot4(ct,prt)
      elseif(i.eq.5) then
        call uplot5(ct,prt)
      elseif(i.eq.6) then
        call uplot6(ct,prt)
      elseif(i.eq.7) then
        call uplot7(ct,prt)
      elseif(i.eq.8) then
        call uplot8(ct,prt)
      elseif(i.eq.9) then
        call uplot9(ct,prt)
      elseif(i.eq.10) then
        call uplot0(ct,prt)
      endif

      end
