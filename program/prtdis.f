c$Id: prtdis.f,v 1.2 2006/12/18 17:00:16 rlt Exp $
      subroutine prtdis(x,b,ttim,prop,ndm,ndf,n1,n2,n3,ii,
     &                  ndf1,ndf2,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Remove fmt4 and print correctly using fmt3         13/11/2006
c     2. Add 'force' output option                          17/12/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output nodal values for real solutions

c      Inputs:
c         x(ndm,*)  - Nodal coordinates of mesh
c         b(*)      - Current value of solution
c         ttim      - Value of solution time
c         prop      - Value of total proportional load
c         ndm       - Spatial dimension of mesh
c         ndf       - Number dof/node
c         n1        - First node to output
c         n2        - Last noed to output
c         n3        - Increment to n1
c         ii        - Type of output: 1 = displacement; 2 = velocity;
c                                     3 = acceleration; 4 = forces;
c                                     5 = eigenvector;
c         ndf1      - First component to output
c         ndf2      - Last  component to output
c         prth      - Output title/header data if true

c      Outputs:
c         None      - Outputs to file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'fdata.h'
      include  'iofile.h'
      include  'xtout.h'
      include  'pfeapb.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   prth,tot, cknonzero
      character nd*4,cd*6,di(5)*6,fmt1*30,fmt2*30,fmt3*42
      integer   ndm,ndf,n1,n2,n3,ii, ndf1,ndf2
      integer   i,n, count, nxt1, lnod,gnod, bserchi
      real*8    ttim,prop, x(ndm,*),b(ndf,*)

      save

      data      nd   /'Node'/,cd /' Coord'/
      data      di   /' Displ',' Veloc',' Accel',' Force',' EigV.'/
      data      fmt1 /'(4x,a4,3(i6,a6)/(8x,6(i6,a6)))'/
      data      fmt2 /'(i8,1p,3e12.4/(8x,1p,6e12.4))'/
      data      fmt3 /'(i8,1p,3e12.4/i8,1p,6e12.4/(i8,1p,6e12.4))'/

      tot   = (ndf2-ndf1+ndm).le.5
      if(.not.tot) then
        write(fmt1(8:8),'(i1)') ndm
        write(fmt2(8:8),'(i1)') ndm
      endif
      write(fmt3(8:8),'(i1)') ndm

      count = 0
      nxt1  = max(1,nxt)
      fp(1) = np(190) - 1
      do n = n1,n2,n3
        if( (mr(fp(1)+n).ge.0) .and. ( nxt.eq.0  .or.
     &     ( abs(x(nxt1,n)-xt).le.xtol ) ) ) then

c         Check non-zero

          if(ii.eq.4) then
            cknonzero = .false.
            do i = ndf1,ndf2
              if(b(i,n).ne.0.0d0) then
                cknonzero = .true.
                exit
              endif
            end do ! i
          else
            cknonzero = .true.
          endif
          if(cknonzero) then
            count = count - 1
            if(count.le.0) then
              call prtitl(prth)
              if(ii.le.4) then
                write(iow,2000) ttim,prop
                if(ior.lt.0.and.pfr) then
                  write(*,2000) ttim,prop
                endif
              else
                write(iow,2001) ttim,prop
                if(ior.lt.0.and.pfr) then
                  write(*,2001) ttim,prop
                endif
              endif
              if(tot) then
                write(iow,2002) (i,cd,i=1,ndm),(i,di(ii),i=ndf1,ndf2)
                if(ior.lt.0.and.pfr) then
                  write(*,2002) (i,cd,i=1,ndm),(i,di(ii),i=ndf1,ndf2)
                endif
              else
                write(iow,fmt1) nd,(i,cd,i=1,ndm),(i,di(ii),i=ndf1,ndf2)
                if(ior.lt.0.and.pfr) then
                  write(*,fmt1) nd,(i,cd,i=1,ndm),(i,di(ii),i=ndf1,ndf2)
                endif
              endif
              count = 48000000
            endif
            if(.not.pfeap_on) then
              if(tot) then
                write(iow,2003) n,(x(i,n),i=1,ndm),(b(i,n),i=ndf1,ndf2)
                if(ior.lt.0.and.pfr) then
                  write(*,2003) n,(x(i,n),i=1,ndm),(b(i,n),i=ndf1,ndf2)
                endif
              else
                write(iow,fmt2) n,(x(i,n),i=1,ndm),(b(i,n),i=ndf1,ndf2)
                if(ior.lt.0.and.pfr) then
                  write(*,fmt2) n,(x(i,n),i=1,ndm),(b(i,n),i=ndf1,ndf2)
                endif
              endif
            else
              if(pfeap_gnod) then
                lnod = bserchi(mr(np(244)),numpn, n)
                gnod = mr(np(244)+lnod-1)
              else
                lnod = n
                gnod = mr(np(244)+n-1)
              endif
              if(gnod.gt.0) then
                write(iow,fmt3) lnod,(x(i,lnod),i=1,ndm),
     &                          gnod,(b(i,lnod),i=ndf1,ndf2)
                if(ior.lt.0.and.pfr) then
                  write(*,fmt3) lnod,(x(i,lnod),i=1,ndm),
     &                          gnod,(b(i,lnod),i=ndf1,ndf2)
                endif
              endif
            endif
          endif
        endif
      end do ! n

c     Format

2000  format('  N o d a l   D i s p l a c e m e n t s',5x,
     &  'Time',e18.5/44x,'Prop. Ld.',1pe13.5/1x)

2001  format('  N o d a l   D i s p l a c e m e n t s',5x,
     &  'Time',e18.5/43x,'Eigenvalue',1pe13.5/1x)

2002  format('    Node',6(i6,a6):/(8x,6(i6,a6)))

2003  format(i8,1p,6e12.4/(8x,1p,6e12.4))

      end
