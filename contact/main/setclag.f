c$Id: setclag.f,v 1.1 2006/11/20 20:32:46 rlt Exp $
      subroutine setclag(ixl,nx,ilm,nl,nlag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set Contact LAGrange multiplier equation numbers

c      Inputs:
c        ixl(*)   - List of nodes of contact 'element'
c        nx       - Number of nodes in 'ixl'
c        ilm(*)   - List of nodes with Lagrange multipliers
c        nl       - Numver of nodes in 'ilm'
c        nlag     - Number of multipliers at each node (same at each)

c      Outputs:
c        Passed back in global array: IAD(numnp,3)
c          IAD(i,2) - stores maximum node number for each multiplier node
c          IAD(i,3) - stores Maximum no. multipliers attached to node.
c-----[--+---------+---------+---------+---------+---------+---------+-]
       implicit  none

       include  'cdata.h'     ! Has numnp
       include  'pointer.h'   ! Has np(*) array
       include  'comblk.h'    ! Has mr(*) array

       integer   nx,nl,nlag
       integer   ixl(*),ilm(*)

       call csetleq(ixl,nx,ilm,nl,nlag,mr(np(89)+numnp),mr(np(224)),
     &              numnp)

       end
