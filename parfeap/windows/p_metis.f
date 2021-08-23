c$Id: p_metis.f,v 1.1 2006/11/20 20:44:16 rlt Exp $
      logical function p_metis(domains,numedg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Partition mesh using METIS

c     Inputs:
c        domains  - Number of domains to create

c     Outputs:
c        numedg   - Number of edges in graph
c-----[--.----+----.----+----.-----------------------------------------]
      use        MSFLIB                ! Windows version

      implicit   none

      include   'cdata.h'
      include   'iodata.h'
      include   'iofile.h'
      include   'pointer.h'
      include   'comblk.h'

      logical    exec
      character  runmetis*40
      integer    domains,numedg, i

c     Write graph to file for processing by metis

      call node_graph_out(domains,numnp,numedg,
     &                    mr(np(252)),mr(np(253)))

      runmetis       = ' '
      runmetis(1:14) = 'kmetis metisgr'
      write(runmetis(16:25),'(i10)') domains

      if(domains.le.8) then
        runmetis(1:1) = 'p'
        exec = systemqq(runmetis)  ! P_metis
      else
        exec = systemqq(runmetis)  ! K_metis
      endif
      if(exec) then
        write( *,*) ' Successful METIS execution'
        runmetis       = ' '
        runmetis(1:13) = 'metisgr.part.'
        if(domains.lt.10) then
          write(runmetis(14:14),'(i1)') domains
        elseif(domains.lt.100) then
          write(runmetis(14:15),'(i2)') domains
        elseif(domains.lt.1000) then
          write(runmetis(14:15),'(i3)') domains
        endif

        open(unit = ios,file = runmetis)
        rewind ios
        do i = 1,numnp
          read(ios,*) mr(np(254)+i-1)
          mr(np(254)+i-1) = mr(np(254)+i-1) + 1  ! Make fortran numbering
        end do ! i

c       Delete scratch files

        open (unit = ios,file = 'metisgr',status = 'old')
        close(unit = ios,status = 'delete')
        open (unit = ios,file = runmetis)
        close(unit = ios,status = 'delete')
      else
        write( *,*) ' Unsuccessful METIS execution'
        call plstop()
      endif

      p_metis = exec

      end
