c$Id: pfeapsr.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pfeapsr(array,tdatabuf,pmax)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Send and receive real data using MPI
c               Not used in serial version

c     Inputs:
c        array(pmax) - Real array from this processor 
c        pmax        - Number items in array

c     Outputs:
c        array(pmax) -  Real array accumulated from all processors
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    pmax
      real*8     array(pmax),tdatabuf(pmax)

      save

c     Dummy function

      end
