c$Id: setlagm.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine setlagm(ilagm,ix, ie)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set element lagrange multiplier equations

c      Inputs:
c        ix(nen1,*) - Element connection array
c        ie(nie,*)  - Element control data

c      Outputs:
c        ilagm(*)   - Lagrange multiplier equation numbers
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'sdata.h'

      integer    ilagm(*),ix(nen1,*), ie(nie,*), ma, n

      save

      do n = 1,numel
        ilagm(n) = 0
      end do ! n

      ndl = 0
      do n = 1,numel
        ma = ix(nen1,n)
        if(ma.gt.0) then
          if(ie(nie-8,ma).gt.0) then
            ix(nen+4,n) = n
            ix(nen+5,n) = ilagm(n)
            ilagm(n)    = ilagm(n) + ie(nie-8,ma)
            ndl         = max(ndl,ilagm(n))
          endif
        endif
      end do ! n

      end

      logical function setlagf(ix,ie,nen1,numel)

      implicit   none

      include   'cdat1.h'

      integer    n,nen1,numel, ix(nen1,*),ie(nie,*)

      save

c     Test for lagrange multipliers in elements

      setlagf = .false.
      do n = 1,numel
        if(ix(nen1,n).gt.0) then
          if(ie(nie-8,ix(nen1,n)).gt.0) then
            setlagf = .true.
            return
          endif
        endif
      end do ! n

      end
