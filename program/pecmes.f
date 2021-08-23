c$Id: pecmes.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pecmes()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Edge and coordinate values

c      Inputs:
c        none

c      Outputs:
c        none
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'edgdat.h'
      include   'corset.h'

      save

c     Set edge boundary codes, forces, displacements, and angles

      if(eanfl.or.ebcfl.or.edifl.or.efcfl.or.eprfl.or.ebsfl) then
        call pedgin()
        eanfl = .false.
        ebcfl = .false.
        edifl = .false.
        efcfl = .false.
        eprfl = .false.
        ebsfl = .false.
      endif

c     Set cordinate angles, boundary codes, forces, displacements,
c         proportional load types and surface loads

      if(boufl .or. surfl .or. angfl .or. disfl .or. cprfl .or.
     &   forfl .or. damfl .or. masfl .or. stifl .or. basfl .or.
     &   lfrfl                                            ) then
        call ploadc()
        boufl = .false.
        surfl = .false.
        angfl = .false.
        disfl = .false.
        cprfl = .false.
        forfl = .false.
        lfrfl = .false.
        damfl = .false.
        masfl = .false.
        stifl = .false.
        basfl = .false.
      endif

c     Set body forces and reactions

      if(reafl .or. intfl ) then
        call pintec()
        reafl = .false.
        intfl = .false.
      endif

      end
