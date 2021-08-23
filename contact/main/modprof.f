c$Id: modprof.f,v 1.1 2006/11/20 20:32:46 rlt Exp $
      subroutine modprof (ixl,ida,nnod,ndof)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0
c               Robert L. Taylor           May   30, 1996            1.1

c      Acronym: MODify PRoFile

c      Purpose: Modify profile for active contacts

c      Inputs :
c        ixl(*)     - List nodes on contact element
c        ida(*)     - List degree of freedoms for contact
c        nnod       - Number of nodes on contact element ixl(*)
c        ndof       - Number of degree of freedoms in ida(*)

c      Outputs:
c        jp         - Profile returned via pointer np(20+npart)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'allotd.h'
      include  'cdata.h'
      include  'compac.h'
      include  'compas.h'
      include  'part0.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   setvar,palloc
      integer   ixl(*),ida(*),nnod,ndof

      save

      call cdebug0 ('      modprof ->',-1)

c     Compute contact element size parameters

      if(iccom.eq.1) then
        ncen    = max(ncen,nnod)
        numcels = numcels + 1

c     Transfer active contact element connections to single array

      elseif(iccom.eq.2) then
        numcels = numcels + 1
        call storec(mr(np(168)),ncen1,ixl,nnod,numcels)

c     Compute column heights for contact elements

      elseif(iccom.eq.4) then
        setvar = palloc(120,'TEMP0',nnod*ndof,1)
        call conprof(mr(np(20+npart)),mr(np(120)),mr(np(31)),ixl,ida,
     &               nnod,ndof)
        setvar = palloc(120,'TEMP0',        0,1)

      endif

      end
