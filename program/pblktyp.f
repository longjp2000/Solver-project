c$Id: pblktyp.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      logical function pblktyp(layer,td, ntyp,ns)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Describe element type and number of nodes/element. 

c      Inputs:
c        layer       - Type of element
c        td(*)       - Number nodes on element

c      Outputs:
c        ntyp        - Element type
c        ns          - Generation type
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      logical    pcomp
      character  layer*15
      integer    ntyp,ns, n,j
      real*8     td(4)

c     Line elements

      if    (pcomp(layer,'line',4)) then
        pblktyp = .false.
        n       = nint(td(1))
        if(n.eq.0 .or. n.eq.2) then
          ns = 1
        elseif(n.eq.3) then
          ns = 2
        endif

c     Triangular elements

      elseif(pcomp(layer,'tria',4)) then
        pblktyp = .false.
        n       = nint(td(1))
        j       = nint(td(2))
        if(n.eq.0 .or. n.eq.3) then  ! 3-node linear triangle
          if(j.ge.0) then
            ntyp = min(6,max(1,j))
          elseif(j.eq.-1) then
            ntyp = -1
          endif
        elseif(n.eq.6) then          ! 6-node linear triangle
          ntyp = 7
        elseif(n.eq.7) then          ! 7-node linear triangle
          ntyp = -7
        endif

c     Quadrilateral elements

      elseif(pcomp(layer,'quad',4)) then
        pblktyp = .false.
        n       = nint(td(1))
        if(n.eq.0 .or. n.eq.4) then ! 4-node Linear quadrilateral
          ntyp = 1
        elseif(n.eq. 8) then  !  8-node Serendipity quadratic quad
          ntyp = 8
        elseif(n.eq. 9) then  !  9-node Lagrangian  quadratic quad
          ntyp = 9
        elseif(n.eq.16) then  ! 16-node Lagrangian  cubic     quad
          ntyp = 16
        endif

c     Tetrahedral elements

      elseif(pcomp(layer,'tetr',4)) then
        pblktyp = .false.
        n       = nint(td(1))
        if(n.eq.0 .or. n.eq.4) then  !  4-node linear tetrahedron
          ntyp = 11
        elseif(n.eq.10) then         ! 10-node quadratic tetrahedron
          ntyp = 13
        elseif(n.eq.11) then         ! 11-node quadratic tetrahedron
          ntyp = 15
        endif

c     Brick elements

      elseif(pcomp(layer,'bric',4)) then
        pblktyp = .false.
        n       = nint(td(1))
        if(n.eq.0 .or. n.eq.8) then
          ntyp = 10
        elseif(n.eq.20) then  ! 20-node Serendipity quadratic brick
          ntyp = 14
        elseif(n.eq.27) then  ! 27-node Lagrangian  quadratic brick
          ntyp = 12
        elseif(n.eq.64) then  ! 64-node Lagrangian  cubic     brick
          write(*,*) ' *ERROR* 64-node cubic element not available'
        endif

      else
        pblktyp = .true.
        n       =  0
      endif

      end
