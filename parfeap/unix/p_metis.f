c$Id: p_metis.f,v 1.1 2006/11/20 20:44:10 rlt Exp $
      logical function p_metis(domains,numedg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Partition mesh using METIS

c     Inputs:
c        domains - Number of domains to create

c     Outputs:
c        numedg  - Unused; Output for Windows version only
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'pointer.h'
      include   'comblk.h'

      integer    domains, numedg, edgecut
      integer    options(5), moptions(3)

      save

      data       moptions / 0, 0, 0 /
      data       options  / 0, 0, 0, 0, 0 /

      if(domains.le.8) then
        call Metis_PartGraphRecursive(numnp,mr(np(252)),mr(np(253)),
     &                       moptions(1),moptions(2),moptions(3),1,
     &                       domains,options,edgecut,mr(np(254)))
      else
        call Metis_PartGraphKway(numnp,mr(np(252)),mr(np(253)),
     &                       moptions(1),moptions(2),moptions(3),1,
     &                       domains,options,edgecut,mr(np(254)))
      endif

      p_metis = .true.

      end
