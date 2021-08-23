c$Id: umaclib.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine umaclib(i,lct,ct,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Interface for user command language instructions

c      Inputs:
c         i      - Command number
c         lct    - Character array describing option
c         ct(3)  - Command parameters
c         prt    - Output if true

c      Outputs:
c         N.B.  Users are responsible for generating command options
c               See programmer manual for example.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   prt
      character lct*(*)
      integer   i
      real*8    ct(3)

      save

      if(    i.eq.1) then
        call umacr1(lct,ct,prt)
      elseif(i.eq.2) then
        call umacr2(lct,ct,prt)
      elseif(i.eq.3) then
        call umacr3(lct,ct,prt)
      elseif(i.eq.4) then
        call umacr4(lct,ct,prt)
      elseif(i.eq.5) then
        call umacr5(lct,ct,prt)
      elseif(i.eq.6) then
        call umacr6(lct,ct,prt)
      elseif(i.eq.7) then
        call umacr7(lct,ct,prt)
      elseif(i.eq.8) then
        call umacr8(lct,ct,prt)
      elseif(i.eq.9) then
        call umacr9(lct,ct,prt)
      elseif(i.eq.10) then
        call umacr0(lct,ct,prt)
      endif

      end
