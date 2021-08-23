c$Id: prwext.F,v 1.1 2006/11/21 16:44:39 rlt Exp $
      subroutine prwext(fname)

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
      implicit   none

#     include  "include/finclude/petsc.h"
#     include  "pfeapa.h"

      character  fname*20
      integer    i

      if(ntasks.gt.1) then
        do i = 1,15
          if(fname(i:i).eq.' ') then
            go to 100
          endif
        end do ! i
        i = 16
100     if    (rank.ge.999) then
          write(fname(i:i+4),'(a,i4)') '_'   ,rank+1
        elseif(rank.ge. 99) then
          write(fname(i:i+4),'(a,i3)') '_0'  ,rank+1
        elseif(rank.gt.  9) then
          write(fname(i:i+4),'(a,i2)') '_00' ,rank+1
        else
          write(fname(i:i+4),'(a,i1)') '_000',rank+1
        endif
      endif

      end
