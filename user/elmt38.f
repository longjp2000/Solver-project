c$Id: elmt38.f,v 1.1 2006/11/20 20:33:27 rlt Exp $
      subroutine elmt38(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      real*8  d(*), ul(*), xl(*), tl(*), s(*), p(*)
      integer ix(*), ndf , ndm , nst , isw

      if(isw.gt.0) write(*,2000)
2000  format('    Elmt 38: *WARNING* Dummy subroutine called')
      end
