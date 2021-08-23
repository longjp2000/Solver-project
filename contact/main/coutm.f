c$Id: coutm.f,v 1.1 2006/11/20 20:32:46 rlt Exp $
      subroutine coutm (cs0,cm0,cp0,ics,hic)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Control for output of contact data to file

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_comnd.h'
      include   'c_contac.h'

      integer    n, ofs, neps,dnope,nope, ics(*),hic((c_lp1+c_lp3),*)
      real*8     cs0(nr0,n0c1:nc01,*),cm0(nr0,n0c2:nc02,*)
      real*8     cp0(nr0,n0c3:nc03,*)

c     Surface outputs

      do n = 1,numcs
        ofs   = abs(cs0(2,-1,n))
        neps  = abs(cs0(3,-1,n))
        dnope = abs(cs0(4,-1,n))
        nope  = abs(cs0(2, 0,n))

        call coutms(n,ics(ofs),dnope,nope,neps)
      end do ! n

c     Pair outputs

      do n = 1,numcp
        call setcomp(n,cs0,cm0,cp0,hic)
        call coutmp(n,cp0)
      end do

      end
