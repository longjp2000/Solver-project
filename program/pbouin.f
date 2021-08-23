c$Id: pbouin.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pbouin(idl,id,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Data input routine for boundary condition codes

c      Inputs:
c         prt        - Flag, print input data if true
c         prth       - Flag, print title/header if true

c      Scratch:
c         idl(*)     - Local degree of freedom integer data

c      Outputs:
c         id(*)      - B.C./Equation numbers identifiers
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'dstars.h'
      include  'iofile.h'
      include  'sdata.h'

      logical   ckno0i,errck,pinput,prt,prth
      integer   i,ii,is,il,k,l,lg,n,ng,nn
      integer   id(ndf,*),idl(*)
      real*8    td(16)

      save

c     [boun]dary codes - read in restraint conditions for each node

      n = 0
      ng = 0
401   l = n
      lg = ng

c     Input restraint records - limit is 16 nos. / record

402   if(ior.lt.0) then
        write(*,3000)
        call pprint('   >')
      endif
      il = min(ndf+2,16)
      errck = pinput(td,il)
      if(errck) go to 402
      nn = nint(td(1))
      n  = nn + starnd
      ng = nint(td(2))
      do k = 1,min(ndf,14)
        idl(k) = td(k+2)
      end do ! k
      if(ndf.gt.14 .and. nn.gt.0) then
        do ii = 1,(ndf+2)/16
          is = il+1
          il = min(is+15,ndf+2)
403       errck = pinput(td,il-is+1)
          if(errck) go to 403
          do k = 1,il-is+1
            idl(k+is-3) = td(k)
          end do ! k
        end do ! ii
      endif
      if(nn.gt.0.and.n.le.numnp) then
        do i = 1,ndf
          id(i,n) = idl(i)
          if(l.ne.0.and.idl(i).eq.0.and.id(i,l).lt.0) id(i,n) = -1
        end do ! i
        lg = sign(lg,n-l)
404     l = l + lg
        if((n-l)*lg.le.0) go to 401
        do i = 1,ndf
          if(id(i,l-lg).lt.0) id(i,l) = -1
        end do ! i
        go to 404
      end if
      if(prt) then
        call prtitl(prth)
        write(iow,2000) (i,i=1,ndf)
        if(ior.lt.0) then
          write(*,2000) (i,i=1,ndf)
        endif
        do n = 1,numnp
          if(ckno0i(id(1,n),ndf)) then
            write(iow,2001) n,(id(i,n),i=1,ndf)
            if(ior.lt.0) then
              write(*,2001) n,(id(i,n),i=1,ndf)
            endif
          endif
        end do ! n
      endif

c     Formats

2000  format('  N o d a l   B. C.'//
     &       6x,'Node',9(i2,'-b.c.')/(10x,9(i2,'-b.c.')))

2001  format(i10,9i7/(10x,9i7))

3000  format(' Input: node#, inc., (b. codes, i=1,ndf)')

      end
