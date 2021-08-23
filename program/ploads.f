c$Id: ploads.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine ploads(u,f1,prop,flg,afl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Surface loading routine: Calls ploadl

c      Inputs:
c         u(*)     - Current solution state
c         prop     - Total proportional load paramenter
c         flg      - Flag, form reaction (uncompressed) if true.
c         afl      - Flag, form tangent if true

c      Outputs:
c         f1(*)    - Residual/reaction array
c                    Tangent returned through blank common address
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'ndata.h'
      include  'part0.h'
      include  'pointer.h'
      include  'sdata.h'
      include  'comblk.h'

      logical   flg,afl
      real*8    prop, u(*),f1(*)

      save

c     Form surface type load and tangent load arrays

      call ploadl(mr(np(31)),mr(np(20+npart)),mr(np(34)),f1,
     &            hr(na),hr(nal),hr(nau),hr(np(35)),hr(np(36)),
     &            hr(np(43)),hr(np(44)),u,hr(np(41)),
     &            prop,ndf,ndm,flg,afl)

      end
