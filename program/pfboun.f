c$Id: pfboun.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pfboun(td,x,f,ntyp,ndm,ndf,numnp,numprt,type,prt,prth,
     c                  aname)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Input nodal force or displacement values using
c               specified coordinates

c      Inputs:
c         td(*)    - Coordinates and values for point
c         x(ndm,*) - Nodal coordinates of mesh
c         ntyp(*)  - Node type (negative for inactive)
c         ndm      - Spatial dimension of mesh
c         ndf      - Number dof/node
c         numnp    - Number of nodes in mesh
c         numprt   - Print counter
c         type     - 'add' to accumlate, otherwise set value
c         prt      - Output generated values if true
c         prth     - Output title/header data if true
c         aname    - Header type name

c      Outputs:
c         f(ndf,*) - Generated force or displacements
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      character type*4,aname*(*),pname*6
      logical   prt,prth, clflg, pcomp
      integer   ndm,ndf,numnp,numprt, n,nbc, ntyp(*)
      real*8    dotx, xmn, tmn, x(ndm,numnp),td(*),f(ndf,numnp)

      save

c     Find closest node to input coordinates

      pname = aname
      if(prt .and. numprt.le.0) then
        call prtitl(prth)
        write(iow,2000) (n,pname,n=1,ndf)
        if(ior.lt.0) write(*,2000) (n,pname,n=1,ndf)
        numprt = 50
      endif

      clflg = .false.
      do n = 1,numnp
        if(ntyp(n).ge.0) then
          tmn = dotx(td(1),x(1,n),ndm)
          if(clflg) then
            if(tmn.lt.xmn) then
              xmn = tmn
              nbc = n
            endif
          else
            xmn   =  tmn
            nbc   =  n
            clflg = .true.
          endif
        endif
      end do ! n

c     Set b.c. forced values

      if(clflg) then
        if(pcomp(type,'add',3)) then
          do n = 1,ndf
            f(n,nbc) = f(n,nbc) + td(ndm+n)
          end do ! n
        else
          do n = 1,ndf
            f(n,nbc) = td(ndm+n)
          end do ! n
        endif

c       Output current restraint codes set

        if(prt) then
          write(iow,2001) nbc,(td(ndm+n),n=1,ndf)
          if(ior.lt.0) then
            write(*,2001) nbc,(td(ndm+n),n=1,ndf)
          endif
          numprt = numprt - 1
        endif
      endif

c     Format

2000  format('  C o o r d i n a t e    N o d a l    V a l u e s'/
     &       /3x,'Node',6(i6,a6:)/(7x,6(i6,a6:)))

2001  format(i7,1p,6e12.4/(7x,1p,6e12.4))

      end
