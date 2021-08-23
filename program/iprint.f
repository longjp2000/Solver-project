c$Id: iprint.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine iprint(ia,ii,jj,mm,name)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output array of integer values

c      Inputs:
c         ia(mm,*) - Array to output
c         ii       - Number of rows to output
c         jj       - Number of columns to output
c         mm       - Dimension of array
c         name     - Name of array to appear with outputs

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      character name*(*), aname*30
      integer   ii,jj,mm,nn, ja,jb,j,n,i,icol,irow
      integer   ia(mm,*)

      save

c     Print an ii x jj array whose dimension is mm for first subscript

      aname = name
      icol  = 6
      nn    = (jj+icol-1)/icol
      jb    = 0

      do n = 1,nn

        ja = jb + 1
        jb = min(jj,ja+icol-1)

        if(n.eq.1 .or. ii.gt.1) then
          write(iow,2000) aname,(j,j=ja,jb)
        endif

        do i=1,ii
          if(ii.eq.1) then
            irow = n
          else
            irow = i
          endif
          write(iow,2001)irow,(ia(i,j),j=ja,jb)
        end do ! i

        if(ior.lt.0) then
          if(n.eq.1 .or. ii.gt.1) then
            write(*,2000) aname,(j,j=ja,jb)
          endif
          do i=1,ii
            if(ii.eq.1) then
              irow = n
            else
              irow = i
            endif
            write(*,2001)irow,(ia(i,j),j=ja,jb)
          end do ! i
        endif

      end do ! n

c     Formats

 2000 format(/4x,'Matrix: ',a30/4x,'row/col',6i10)

 2001 format(1i8,1i13,5i10)

      end
