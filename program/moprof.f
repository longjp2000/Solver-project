c$Id: moprof.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine moprof(ic,hr,neq, name)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Print profile stored array

c      Inputs:
c         ic(*)  - Column pointers
c         hr(*)  - Array of real values to print
c         neq    - Number of row/colums
c         name   - Character identification array

c      Outputs:
c         To file and screen
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iofile.h'

      character  name*(*)
      integer    neq,ic(*)
      integer    n,i, istart
      real*8     hr(*)

      save

c     Print profile stored array by equation number

      write(iow,2000) name
      if(ior.lt.0) then
        write(*,2000) name
      endif
      istart = 1
      do n = 1,neq
        write(iow,2001) n, (hr(i),i=istart,ic(n))
        if(ior.lt.0) then
          write(*,2001) n, (hr(i),i=istart,ic(n))
        endif
        istart = ic(n) + 1
      end do ! n

c     Formats

2000  format(/'  S p a r s e    M a t r i x  : ',a/)
2001  format(i8,':',1p,4e15.5:/(9x,1p,4e15.5))

      end
