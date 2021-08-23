c$Id: pelnum.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pelnum(tx,iel,errck)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Get material set numbers for FEAP elements

c              Current FEAP Element Types
c                 Name     |     iel
c              ------------+-------------
c               Solid      |     -1
c               Truss      |     -2
c               Frame      |     -3
c               Plate      |     -4
c               Shell      |     -5
c               Membrane   |     -6
c               Gap        |     -7
c               Thermal    |     -8
c               Convection |     -9
c               Point      |     -10
c               Pressure   |     -11

c      Inputs:
c         tx     - Name of element type requested

c      Outputs:
c         iel    - Element type for request
c         errck  - Flag, true if request found
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   pcomp, errck
      character tx*15
      integer   iel

      save

      errck = .false.

c     Check based on coordinate dimension

      if    (pcomp(tx,'soli',4)) then
        iel = -1
      elseif(pcomp(tx,'trus',4)) then
        iel = -2
      elseif(pcomp(tx,'fram',4)) then
        iel = -3
      elseif(pcomp(tx,'plat',4)) then
        iel = -4
      elseif(pcomp(tx,'shel',4)) then
        iel = -5
      elseif(pcomp(tx,'memb',4)) then
        iel = -6
      elseif(pcomp(tx,'gap',3)) then
        iel = -7
      elseif(pcomp(tx,'ther',4)) then
        iel = -8
      elseif(pcomp(tx,'conv',4)) then
        iel = -9
      elseif(pcomp(tx,'poin',4)) then
        iel = -10
      elseif(pcomp(tx,'pres',4)) then
        iel = -11
      else
        errck = .true.
      endif

      end
