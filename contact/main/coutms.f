c$Id: coutms.f,v 1.1 2006/11/20 20:32:46 rlt Exp $
      subroutine coutms(n,ics,dnope,nope,neps)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Output surface data to file

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'iodata.h'
      include   'comblk.h'

      character  type(4)*5
      integer    n,dnope,nope,neps
      integer    ics(dnope,neps), i

      data       type / 'POINT','LINE','TRIA','QUAD' /

      write(ios,2000) n,type(nope)
      do n = 1,neps
        write(ios,2001) n,' 0 ',(mr(ics(i,n)+1),i=1,nope)
      end do ! n

c     Formats

2000  format(/'SURFACE',i5/'  ',a/'    FACET')
2001  format(i8,a,8i8)

      end
