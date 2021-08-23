c$Id: pfeapmi.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pfeapmi(intsr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Send and receive data using MPI:
c               Not used in serial version

c     Inputs:
c        intsr - Integer value from current processor

c     Outputs:
c        intsr - Maximum value from all processors
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    intsr

      save

c     Dummy function

      end
