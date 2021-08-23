c$Id: psetb.F,v 1.2 2006/11/25 02:13:57 rlt Exp $
      subroutine psetb(b,getp,getv,senp,senv,id, ndf,
     &                 rdatabuf, sdatabuf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Revise send/receive data add barrier             24/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Transfer PETSC vector to local arrays including ghost
c               node data via MPI messaging

c     Inputs:
c        getp(ntasks)    - Pointer array for getting ghost data
c        getv(sum(getp)) - Local node numbers for getting ghost data
c        senp(ntasks)    - Pointer array for sending ghost data
c        senv(sum(senp)) - Local node numbers to send out as ghost data
c        id(ndf,numnp)   - Local equation numbers
c        ndf             - dofs per node
c        rdatabuf(*)     - receive communication array
c        sdatabuf(*)     - send    communication array

c     Outputs:

c        b(neq) - Local solution vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "include/finclude/petsc.h"
#     include   "pfeapa.h"

#     include   "cdata.h"
#     include   "iofile.h"
#     include   "pfeapb.h"

      integer    ndf
      integer    i, j, k, lnodei, eqi, roff,soff, rbuf,sbuf,tbuf, idesp
      integer    getp(*),getv(*),senp(*),senv(*)
      integer    id(ndf,*)
      real*8     b(*), rdatabuf(*), sdatabuf(*)

      integer    usolve_msg, req
      parameter (usolve_msg=10)

c     Petsc values

      PetscErrorCode ierr

      integer        msg_stat(MPI_STATUS_SIZE)

c     Sending Data Asynchronously

      soff  = 0
      idesp = 0
      do i = 1, ntasks

        if (senp(i) .gt. 0) then
          sbuf   = soff
          do j = 1, senp(i)
            lnodei = senv(j + idesp)
            do k = 1, ndf
              eqi = id(k,lnodei)
              if (eqi .gt. 0) then
                sbuf = sbuf + 1
                sdatabuf( sbuf ) = b(eqi)
              endif
            end do ! k
          end do ! j
          idesp = idesp + senp(i)
          call MPI_Isend( sdatabuf(soff+1), sbuf, MPI_DOUBLE_PRECISION,
     &                    i-1, usolve_msg, MPI_COMM_WORLD, req, ierr )
          soff = soff + sbuf
        endif

c       Report send error

        if(ierr.ne.0) then
          write(iow,*) ' -->> Send Error[',rank+1,'->',i,']',ierr
        endif

      end do ! i

c     Receiving Data in blocking mode

      idesp = 0
      roff  = 0
      do i = 1, ntasks
        if (getp(i) .gt. 0) then
          tbuf = getp(i)*ndf
          call MPI_Recv( rdatabuf(roff+1), tbuf, MPI_DOUBLE_PRECISION,
     &                   i-1, usolve_msg, MPI_COMM_WORLD,msg_stat, ierr)
          rbuf = roff
          do j = 1, getp(i)
            lnodei = getv(j + idesp)
            do k = 1, ndf
              eqi = id(k,lnodei)
              if (eqi .gt. 0 ) then
                rbuf  = rbuf + 1
                b( eqi ) = rdatabuf( rbuf )
              endif
            end do ! k
          end do ! j
          idesp = idesp + getp(i)
          roff  = roff + rbuf
        endif
      end do ! i

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      end
