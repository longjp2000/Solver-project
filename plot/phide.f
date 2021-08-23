c$Id: phide.f,v 1.1 2006/11/20 20:33:12 rlt Exp $
      subroutine phide(ct,nix,nxd,nxn,nne,nface,iln)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: do hidden line removal

c      Inputs:
c         ct        - Plot negative faces if positive
c         iln(2)    - Line type data
c         nface     - Number of faces on surfaces

c      Outputs:
c         nix       - Face connection location
c         nxd       - Face connection dimension
c         nxn       - Number of nodes/face
c         nne       - Number of faces
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comblk.h'
      include  'cdata.h'
      include  'pdata3.h'
      include  'pdatay.h'
      include  'plflag.h'
      include  'pointer.h'
      include  'ppers.h'
      include  'sdata.h'
      include  'iofile.h'

      integer   nix,nxd,nxn,nne,nface, iln(2)
      real*8    ct

      save

c     Plot visible mesh

      call pzeroi(mr(np(66)),numnp)
      nix     = 54
      nxd     = 7
      nxn     = 4
      nne     = nface
      nfac(1) = nne
      call plface(mr(np(nix)),mr(np(62)),hr(np(53)),
     &            3,nxd,numnp,nne,iln,ct)
      if(ndm.eq.3 .and. nen.gt.3) then
        call p3edge(mr(np(nix)),hr(np(43)),ndm,nface,nxd)
        edgfl = .true.
      endif

c     Set plot sequence for z-sort

      if(kpers.ne.0) then
        call perspz(hr(np(53)),mr(np(nix)), mr(np(61)),mr(np(55)),
     &              mr(np(62)),nxd,nxn,  3,numnp,nne)
      endif

      end
