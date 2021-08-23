c$Id: pddot.f,v 1.1 2006/11/20 20:43:46 rlt Exp $
      function pddot(n,v1,n1,v2,n2)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Compute parallel vector dot product

c      Inputs:
c         n     - Length of vectors in current domain
c         v1(*) - Vector 1
c         n1    - Step size
c         v2(*) - Vector 1
c         n2    - Step size

c      Output:
c         pddot - real double precision
c-----[--+---------+---------+---------+---------+---------+---------+-]
       implicit   none

       include   'pfeapb.h'

       integer    n, n1,n2
       real*8     pddot,v1(*),v2(*), ddot, tdatabuf(2)

       pddot  = ddot(numpeq,v1,n1,v2,n2)

       call pfeapsr(pddot,tdatabuf,1)

       end
