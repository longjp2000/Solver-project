c$Id: pmovec.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pmovec(id,a,b,nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Move compressed array a into uncompressed array b

c      Inputs:
c         a(*)      - Compressed array to move
c         nn        - Length of uncompressed array

c      Outputs:
c         b(*)    - Uncompressed move of a (zero undefined values)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n,nn,j,id(*)
      real*8    a(nn),b(nn)

      save

      do n = 1,nn
        j = id(n)
        if (j.gt.0) then
          b(n) = a(j)
        else
          b(n) = 0.0d0
        endif
      end do ! n

      end
