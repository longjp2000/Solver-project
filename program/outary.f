c$Id: outary.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine outary(array,ct)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output contents of named "array"

c      Inputs:
c         array     - name of array to print
c         ct(3)     - range of array to print

c      Outputs:
c         none      - To screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'allotn.h'
      include  'allotd.h'
      include  'cdata.h'
      include  'comblk.h'
      include  'pointer.h'

      character array*(*)
      logical   pcomp
      integer   lengt, i, ilor,iupr,ict(3), ipa, irp(2,2)
      real*8    ct(3)

      save

      data      irp  / 3*1, 2 /

c     Output contact dictionary array values

      ict(1) = ct(1)
      if (((pcomp(array,'C0   ', 5)) .or.
     &     (pcomp(array,'CM   ', 5)) .or.
     &     (pcomp(array,'ICS  ', 5)) .or.
     &     (pcomp(array,'HIC  ', 5)) .or.
     &     (pcomp(array,'CH   ', 5))).and. ict(1).eq.0 ) then

c       Output of contact arrays

        call contact (303)

      else

c       Set range of print

        ilor   = max(0,min(nint(ct(1)),nint(ct(2))))
        iupr   = max(nint(ct(1)),nint(ct(2)))

        do i = 1,ndict
          if(pcomp(array,dict(i),5)) then

c           Assign pointer, length, and precision

            ipa   = irp(iprec(i),ipr)
            lengt = ipoint(i)

c           Set range

            if(iupr.eq.0) then
              iupr = lengt - ilor
            else
              ilor = min(iupr,ilor)  - 1
              iupr = min(iupr,lengt) - ilor
            endif

c           Output array values

            if(iprec(i).eq.1) then
              call iprint(mr(np(dlist(i))+ilor),1,iupr,1,dict(i))
            else
              call mprint(hr(np(dlist(i))+ilor),1,iupr,1,dict(i))
            endif
          endif
        end do ! i

      endif

      end
