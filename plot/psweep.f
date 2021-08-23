c$Id: psweep.f,v 1.1 2006/11/20 20:33:12 rlt Exp $
      subroutine psweep(swang, nsinc,nix,nxd,nxn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Plot sweep of 2-d mesh through angle 'swang'

c      Inputs:
c        swang   - Sweep angle in degrees
c        nsinc   - Sweep increments
c        nix     - Address of connections
c        nxd     - Dimension of nix array
c        nxn     - Number of nodes on array

c      Outputs:
c        to screen
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'pointer.h'
      include   'comblk.h'

      integer    nsinc,nix,nxd,nxn
      real*8     swang

      save

      call pppcol(1,1) ! White

      call pswsub(mr(np(nix)), mr(np(78)), mr(np(62)), hr(np(53)),
     &            nxd,nxn, swang, nsinc )

      end
