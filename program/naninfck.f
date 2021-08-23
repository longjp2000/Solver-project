c$Id: naninfck.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      logical function naninfck(y,sizey,isw)

c     * * F E A P * * A Finite Element Analysis Program      

c.... Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]      
c     Purpose: Check if a double precision number is NAN or -INF or +INF
c              This is a Fortran  replacement for C: isinf
c                                                    isnan
c       Inputs:
c          y(*)     : The array of values to be checked
c          sizey    : Length of list
c          isw      : Switch - 0: Initialize values; >0: Check values

c       Outputs:
c         naninfck  : A logical variable which is true or false
c-----[--.----+----.----+----.-----------------------------------------]
      use, intrinsic :: iso_fortran_env
      use, intrinsic :: iso_fortran_env
      use, intrinsic :: ieee_arithmetic
      implicit none

      real*8   y(*)
      integer  sizey,isw,ic

c     Set values for checks

      if(isw.eq.0) then

        naninfck = .true.

c     Check values

      else

        naninfck = .false.
        do ic = 1,sizey
           if ( ieee_is_nan ( y(ic) ) ) then
            naninfck = .true.
            return
           else if ( .not. ieee_is_finite ( y(ic) ) )then
            naninfck = .true.
            return
          end if
        end do ! ic
      
      endif ! isw 

      end 
