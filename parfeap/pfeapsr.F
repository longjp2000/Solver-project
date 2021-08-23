c$Id: pfeapsr.F,v 1.1 2006/11/21 16:44:39 rlt Exp $
      subroutine pfeapsr(array,tdatabuf,pmax)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Send and receive real data using MPI: PETSc version

c     Inputs:
c        array(pmax) - Real array from this processor 
c        pmax        - Number items in array

c     Outputs:
c        array(pmax) -  Real array accumulated from all processors
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "include/finclude/petsc.h"
#     include   "pfeapa.h"

      integer    i, j, pmax
      real*8     array(pmax),tdatabuf(pmax,2)

      integer    usr_msg, req
      parameter (usr_msg=12)

c     Petsc values

      PetscErrorCode ierr

      integer        msg_stat(MPI_STATUS_SIZE)

      save

c     Sending Data Asynchronously

      do j = 1, pmax
        tdatabuf(j,1) = array(j)
      end do ! j

      do i = 0,ntasks-1
        if(rank.ne.i) then
          call MPI_Isend( tdatabuf(1,1), pmax, MPI_DOUBLE_PRECISION, i,
     &                    usr_msg, MPI_COMM_WORLD, req, ierr )
        endif
      end do ! i

c     Receiving Data in blocking mode

      do i = 0,ntasks-1
        if(rank.ne.i) then
          call MPI_Recv( tdatabuf(1,2), pmax, MPI_DOUBLE_PRECISION, i,
     &                     usr_msg, MPI_COMM_WORLD, msg_stat, ierr )
          do j = 1, pmax
            array(j) = array(j) + tdatabuf(j,2)
          end do ! j
        endif
      end do ! i

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      end
