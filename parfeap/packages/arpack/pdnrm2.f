c$Id: pdnrm2.f,v 1.1 2006/11/20 20:43:46 rlt Exp $
      function pdnrm2(n, v1, n1)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Compute Euclidian-norm of v1

c      Inputs:
c         nn    - Length of vectors in current domain
c         v1(*) - Vector 1
c         n1    - Step size

c      Output:
c         pdnrm2 - real double precision
c-----[--+---------+---------+---------+---------+---------+---------+-]
       implicit   none

       include   'pfeapb.h'

       integer    n,  n1
       real*8     pdnrm2,v1(*), pddot

       pdnrm2 = pddot(numpeq, v1,n1, v1,n1 )
       pdnrm2 = sqrt(pdnrm2)

       end
