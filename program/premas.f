c$Id: premas.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine premas(fl,cmtyp,unsfl)

c     * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Initialize data for mass computations

c      Inputs:
c         fl(2)  -  Flag for mass type: 1 = lump, 2 = consistent
c         cmtyp  -  Mass type: symmetric or unsymmetric

c      Outputs:
c         unsfl  -  Symmetry flag
c         np(*)  -  Pointers to allocated array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none
      logical    fl(2), unsfl, setvar,palloc, pcomp
      character  cmtyp*4, tname*5
      integer    k1

      include   'cdata.h'
      include   'compas.h'
      include   'eqsym.h'
      include   'ndata.h'
      include   'part0.h'
      include   'pfeapb.h'
      include   'pointer.h'
      include   'comblk.h'

c     Consistent mass initialization

      if(pfeap_on) then
        call upremas(fl)
        unsfl = .false.
      elseif(fl(1)) then
        if(compms) then
          k1 = 0
          call iters(k1,2)
          compms = .false.
        else
          compfl = .true.
        endif
        nx = npart+8
        if(pcomp(cmtyp,'unsy',4)) then
          unsfl  = .true.
          neqs   = 1
          write(tname,'(4hCMAS,i1)') npart
          setvar = palloc(npart+8,tname, nnm+nnm-neq, 2)
          nm     = np(npart+8)
          nmu    = nm  + neq
          nml    = nmu + nnm - neq
          call pzero (hr(nm),nnm+nnm-neq)
        else
          unsfl  = .false.
          neqs   = neq
          nm     = np(npart+8)
          nmu    = nm  + neq
          nml    = nmu
          call pzero (hr(nm),nnm)
        endif
      elseif(fl(2)) then
        write(tname,'(4hLMAS,i1)') npart
        setvar = palloc(npart+12,tname,neq,2)
        nx = npart + 12
        nl = np(nx)
        call pzero(hr(nl),neq)
        unsfl = .false.
      endif

      end
