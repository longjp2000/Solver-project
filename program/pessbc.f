c$Id: pessbc.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pessbc(idl,id,ndf,numnp,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Input restraint conditions for each node (alternative)

c      Inputs:
c         ndf        - Number dof/node
c         numnp      - Number nodes in mesh
c         prt        - Output results if true
c         prth       - Output title/header data if true

c      Scratch:
c         idl(*)     - Integer storage for inputs

c      Outputs:
c         id(ndf,*)  - Boundary condition restraint conditsions
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'dstars.h'
      include  'iofile.h'

      logical   prt,prth,errck,pinput,ckno0i
      integer   ndf,numnp,n0,n1,n2,ng
      integer   i,j, ii,il,is, idl(ndf),id(ndf,numnp)
      real*8    td(16)

      save

      n1 = 0
      n2 = 0
      ng = 0

100   if(ior.lt.0) then
        write(*,3000)
        call pprint('   >')
      endif

c     Input restraint records - limit is 16 nos. / record

      il = min(ndf+3,16)
      errck = pinput(td,il)
      if(errck) go to 100
      if(nint(td(1)).gt.0) then
        n0 = min(td(1),td(2))
        n1 = n0 + starnd
        n2 = max(td(1),td(2)) + starnd
        if(n0.eq.0) n1 = n2
        if(n2.gt.numnp) then
          write(ilg,4000) n1,n2
          write(  *,4000) n1,n2
          return
        endif
        ng = max(abs(nint(td(3))),1)
        do i = 1,min(ndf,13)
          idl(i) = td(i+3)
        end do ! i
        if(ndf.gt.13) then
          do ii = 1,(ndf+3)/16
            is = il+1
            il = min(is+15,ndf+3)
101         errck = pinput(td,il-is+1)
            if(errck) go to 101
            do i = 1,il-is+1
              idl(i+is-4) = td(i)
            end do ! i
          end do ! ii
        endif
        do j = n1,n2,ng
          do i = 1,ndf
            id(i,j) = idl(i)
          end do ! i
        end do ! j
        go to 100
      else
        if(prt) then
          call prtitl(prth)
          write(iow,2000) (i,i=1,ndf)
          if(ior.lt.0) then
            write(*,2000) (i,i=1,ndf)
          endif
          do j = 1,numnp
            if(ckno0i(id(1,j),ndf)) then
              write(iow,2001) j,(id(i,j),i=1,ndf)
              if(ior.lt.0) then
                write(*,2001) j,(id(i,j),i=1,ndf)
              endif
            endif
          end do ! j
        endif
      endif

c     Formats

2000  format('  N o d a l   B. C.'//
     1       6x,'Node',9(i2,'-b.c.')/(10x,9(i2,'-b.c.')))
2001  format(i10,9i7/(10x,9i7))

3000  format(' Input: node1, node2, inc., (b. codes, i=1,ndf)')

4000  format(' *ERROR* PESSBC: Attempt to set b.c. for nodes',i8,' to',
     &       i8)

      end
