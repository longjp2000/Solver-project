c$Id: datric.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine datric(iops,mops,jp,neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Symbolic decomposition of matrix stored in profile form

c      Inputs:
c         jp(*)  - Pointer to ends of rows/columns in A
c         neq    - Number of equations in A

c      Outputs:
c         iops   - Number of operations in excess of mega-op counter
c         mops   - Mega-operation counter
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer   i, id, ih, is, iops, mops, j, jd, jh, neq, jp(*)
      real*4    etime, tt(2), tar

      save

c     Loop through columns to count operations

      jd   = 1
      iops = 0
      mops = 0
      do j = 2,neq
        jh = jp(j) - jp(j-1) - 1
        if(jh.gt.0) then
          is = j - jh

c         Reduce each column

          do i = is,j-1
            id = jp(i)
            ih = min(id-jp(i-1),i-is+1)
            if(ih.gt.0) then
              iops = iops + ih + ih + 1
            endif
          end do ! i
        endif

c       Reduce diagonal

        if(jh.ge.0) then
          iops = iops +3*jh + 2
        endif

c       Accumulate megaflops

        do while(iops.gt.1000000)
          iops = iops - 1000000
          mops = mops + 1
        end do ! while
      end do ! j

c     Total number of operations

      tar = etime(tt)
      write(iow,2000) iops,mops,tt
      if(ior.lt.0) then
        write(  *,2000) iops,mops,tt
      end if

2000  format('   Number of operations =',i8,' plus ',i8,' Mega-ops'/
     &       '  Time: CPU = ',f12.2,' , System = ',f12.2)

      end
