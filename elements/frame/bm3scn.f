c$Id: bm3scn.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine bm3scn(yy,zz,sig,ii,nqudr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pointer.h'
      include  'comblk.h'

      integer   ii,nqudr
      real*8    yy,zz,sig

      save

      if(np(199).ne.0) then
        call pbm3sc(mr(np(199)),hr(np(200)),hr(np(201)),hr(np(202)),
     &              yy,zz,sig,ii,nqudr)
      else
        write(  *,*) '  *WARNING*: No arrays to project stresses!'
        write(iow,*) '  *WARNING*: No arrays to project stresses!'
      endif

      end
