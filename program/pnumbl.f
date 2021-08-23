c$Id: pnumbl.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pnumbl(ndm,nr,ns,nt,ntyp, nf,ng, flag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generage blend blocks
c      Inputs:
c        ndm       - Spatial dimension of mesh
c        nr        - Number of 1-direction increments to generate
c        ns        - Number of 2-direction increments to generate
c        nt        - Number of 3-direction increments to generate
c        ntyp      - Element type to generate
c        flag      - 1-d generation

c      Outputs:
c        nf        - Number last element
c        ng        - Number last node
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   flag
      integer   ndm,nr,ns,nt,ntyp,nf,ng

      save

c     Check the 2-D types

      if(ntyp.lt.10) then
        if(flag) then
          ng = nr + 1
          if(ns.eq.1) then
            nf = nr
          else
            nf = nr/2
          endif
          nr = nr + 1
        else
          if (ntyp.eq.0) then
            nf = nr*ns
          elseif (abs(ntyp).eq.7) then
            nf = (nr*ns)/2
          elseif (ntyp.ge.8) then
            nf = (nr*ns)/4
          elseif (ntyp.lt.0) then
            nf = 4*nr*ns
          else
            nf = 2*nr*ns
          endif

c         Determine last node number to be generated

          nr = nr + 1
          ns = ns + 1
          if(ndm.eq.1) ns = 1
          ng = nr*ns
          if(ntyp.eq. -7) then
            ng = ng + (nr-1)*(ns-1)/2
          elseif(ntyp .eq. -1) then
            ng = ng + (nr-1)*(ns-1)
          elseif(ntyp .eq.  8) then
            ng = ng - ((nr-1)*(ns-1))/4
          endif
        endif

c     3-d generations

      elseif(ntyp.lt.20) then
        if(    ntyp.eq.11) then                 ! Linear    tetrahedron
          nf = nr*ns*nt*6
          ng = (nr+1)*(ns+1)*(nt+1)
        elseif(ntyp.eq.12 .or. ntyp.eq.14) then ! Quadratic hexahedron
          nf = nr*ns*nt/8
          ng = (nr+1)*(ns+1)*(nt+1)
        elseif(ntyp.eq.13) then                 ! Quadratic tetrahedron
          nf = 6*nr*ns*nt/8
          ng = (nr+1)*(ns+1)*(nt+1)
        else                                    ! Linear    hexahedron
          nf = nr*ns*nt
          ng = (nr+1)*(ns+1)*(nt+1)
        endif

c     Shell:

      elseif(ntyp.lt.30) then
        if (ntyp.eq.20) then
          nf = nr*ns
        elseif (ntyp.eq.27) then
          nf = (nr*ns)/2
        elseif (ntyp.ge.28) then
          nf = (nr*ns)/4
        else
          nf = 2*nr*ns
        endif

c       Determine last node number to be generated

        ng = (nr+1)*(ns+1)

c     Line:

      elseif(ntyp.lt.40) then
        if(ntyp.eq.33) then
          nf = nr/ns
        else
          nf = nr
        endif

        ng = nr + 1

      endif

      end
