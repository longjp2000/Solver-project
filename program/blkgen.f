c$Id: blkgen.f,v 1.2 2006/12/14 00:48:06 rlt Exp $
      subroutine blkgen(ndm,nen1,x,ix,rben,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Delete label 100                                 13/12/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate a block of elements and nodes for mesh
c               descriptions.

c      Inputs:
c         ndm        - Dimension of 'x' array
c         nen1       - Dimension for 'ix' array
c         prt        - Print generated data if true
c         prth       - Print headers if true

c      Outputs:
c         x(ndm,*)   - Block of nodal coordinates
c         ix(nen1,*) - Block of elements
c         rben(numel) - Rigid body number associated with element
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cblktr.h'
      include  'cdata.h'
      include  'crotas.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'region.h'
      include  'rigid2.h'
      include  'comblk.h'

      logical   prt,prth,pcomp,errck,pinput,tinput,palloc,bernfl,err
      logical   eltype,pblktyp
      character xh*6,ctype*15,layer*15
      integer   i,j,k,l,n,nm,nn,nr,ns,nt,nf,ng,ni,ntyp,nodinc
      integer   ndm,nen1,ne,ma, dlayer,nlay
      integer   ix(nen1,*),ixl(27),rben(*)
      real*8    dr,ds,dt
      real*8    x(ndm,*),xl(3,27),shp(3,16),tb(7),tc(4),td(16)

      save

      data      xh/' coord'/

c     Block mesh generation routine

      if(ior.lt.0) then
        if(ndm.le.2) then
          write(*,5000)
          call pprint('   >')
        elseif(ndm.ge.3) then
          write(*,5001)
          call pprint('   >')
         endif
      endif

c     Start input of data

      bernfl = .false.
101   errck = tinput(ctype,1,tb,7)
      if(errck) go to 101

c     Check for form of coordinate system

      if(pcomp(ctype,'bern',4)) then
        bernfl = .true.
      elseif(.not.pcomp(ctype,'pola',4) .and.
     &       .not.pcomp(ctype,'cyli',4) .and.
     &       .not.pcomp(ctype,'sphe',4)) then
        ctype = 'cartesian'
      endif

c     Set parameters for block

      nr     = tb(1)
      ns     = tb(2)
      nt     = tb(3)
      ntyp   = tb(7)
      do n = 1,27
        ixl(n) = 0
        do j = 1,ndm
          xl(j,n) = 0.0
        end do ! j
      end do ! n
      nn     = 0      ! Number of block nodes
      nm     = 0      ! Number of block nodes 
      nlay   = 0      ! Number of material layers in block
      dlayer = 0      ! Direction of layering

c     Check for layered material properties

110   errck = tinput(layer,1,td,3)
      if(errck) go to 110
c     Layer data

      if(pcomp(layer,'laye',4)) then
        dlayer = nint(td(1))
        if(dlayer.eq.1) then
          nlay = nr
        elseif(dlayer.eq.2) then
          nlay = ns
        elseif(dlayer.eq.3) then
          nlay = nt
        endif
        errck = palloc(111,'TEMP1',nlay,1)

        j = 1
        do while(j.le.nlay)
120       errck = tinput(layer,0,td,16)
          if(errck) go to 120
          do i = j,min(j+15,nlay)
            mr(np(111)+i-1) = nint(td(i-j+1))
          end do ! i
          j = j + 16
        end do ! while
        go to 110

c     Data is element type or node of block coordinate

      else
        eltype = pblktyp(layer,td, ntyp,ns) 
        if(eltype) then
          call setval(layer,15,tc(1))
          l  = nint(tc(1))
          if(l.eq.0 ) go to 22
          if(l.gt.27) go to 402
          nm = max(nm,l)
          nn = nn + 1
          ixl(l)  = l
          xl(1,l) = td(1)
          xl(2,l) = td(2)
          xl(3,l) = td(3)
        else
          go to 110  ! Input another record
        endif
      end if ! no layers

c     Input block coordinates

21    if(ior.lt.0) then
        write(*,5002)
        call pprint('   >')
      endif
        errck = pinput(tc,4)
        if(errck) go to 21
        l = tc(1)
        if(l.eq.0 ) go to 22
        if(l.gt.27) go to 402
        nm = max(nm,l)
        nn = nn + 1
        ixl(l)  = l
        xl(1,l) = tc(2)
        xl(2,l) = tc(3)
        xl(3,l) = tc(4)
      go to 21

c     Set a default value for unspecified block type

22    if(ntyp.eq.0 .and. nt.gt.0 .and. ndm.eq.3 .and. ixl(8).ne.0) then
        ntyp = 10
      endif

c     2-d generations

      if(ntyp.lt.10) then

        nt     = 1
        ni     = tb(3)
        ne     = tb(4)
        ma     = tb(5)
        nodinc = tb(6)

c     3-d generations

      elseif(ntyp.lt.20) then
        nt     = tb(3)
        ni     = tb(4)
        ne     = tb(5)
        ma     = tb(6)
        nodinc = 0

c     Shell generations

      elseif(ntyp.lt.30) then

        nt     = 1
        ni     = tb(3)
        ne     = tb(4)
        ma     = tb(5)
        nodinc = tb(6)

c     User generations

      else
        if(ndm.le.2) then
          nt     = 1
          ni     = tb(3)
          ne     = tb(4)
          ma     = tb(5)
          nodinc = tb(6)
        else
          nt     = tb(3)
          ni     = tb(4)
          ne     = tb(5)
          ma     = tb(6)
          nodinc = 0
        endif
      endif

c     Reset to default values if necessary

      if(ni.eq.0) ni = nio + 1
      if(ne.eq.0) ne = neo + 1
      if(ma.eq.0) ma = mao

      ma     = max(ma,1)
      nodinc = max(nodinc,0)
      nr     = max(nr,1)
      ns     = max(ns,1)
      nt     = max(nt,1)
      ni     = max(ni,1)

c     Output block data

      if(prt) then
        call prtitl(prth)
        write(iow,2000) nr,ns,nt,ni,ne,ma,nodinc,ntyp
        if(ne.le.0) write(iow,3000)
        if(ne.le.0) write(ilg,3000)
        if(ior.lt.0) then
          write(*,2000) nr,ns,nt,ni,ne,ma,nodinc,ntyp
          if(ne.le.0) write(*,3000)
        endif
        if(nlay.gt.0) then
          write(iow,2005) (j,j=1,min(5,nlay))
          write(iow,2006) (j,mr(np(111)+j-1),j=1,nlay)
          if(ior.lt.0) then
            write(*,2005) (j,j=1,min(5,nlay))
            write(*,2006) (j,mr(np(111)+j-1),j=1,nlay)
          endif
        end if ! nlay > 0
        write(iow,2002) ctype,(i,xh,i=1,ndm)
        if(ior.lt.0) then
          write(*,2002) ctype,(i,xh,i=1,ndm)
        endif
        do l = 1,27
          if(ixl(l).gt.0) then
            write(iow,2001) l,(xl(i,l),i=1,ndm)
            if(ior.lt.0) then
              write(*,2001) l,(xl(i,l),i=1,ndm)
            endif
          endif
        end do ! l
      endif

c     Set generation increments of natural coordinates

      dr = 2.d0/nr
      ds = 2.d0/ns

c     Determine last element number to be generated

      if(ntyp.lt.10) then
        if(bernfl .or. nn.le.3) then
          ng = ni + nr
          if(ns.eq.1) then
            nf = ne + nr - 1
          elseif(ns.le.9) then
            nf = ne + nr/ns - 1
          else
            nf = ne + nr/2 - 1
          endif
          nr = nr + 1
        else

c         Check block for input errors

          if(pcomp(ctype,'cart',4)) then
            call ck2dblk(ixl,xl,shp,4,3, err)
            if(err) call plstop()
          endif

          if (ntyp.eq.0) then             ! 4-node quadrilateral
            nf = ne + nr*ns - 1
          elseif (abs(ntyp).eq.16) then   ! 16-node quadrilateral
            nf = ne + (nr*ns)/9 - 1
            ntyp = 16
          elseif (abs(ntyp).eq.7) then    ! 6 & 7-node triangle
            nf = ne + (nr*ns)/2 - 1
          elseif (ntyp.ge.8) then         ! 8 & 9 -node quadrilateral
            nf = ne + (nr*ns)/4 - 1
          elseif (ntyp.lt.0) then         ! Crossed triangles
            nf = ne + 4*nr*ns - 1
          else                            ! 3-node triangle
            nf = ne + 2*nr*ns - 1
          endif

c         Determine last node number to be generated

          nr = nr + 1
          ns = ns + 1
          if(ndm.eq.1) ns = 1
          ng = nr*ns + nodinc*(ns-1) + ni -1
          if(ntyp.eq. -7) then
            ng = ng + ((nr-1)*(ns-1))/2
          elseif(ntyp .eq. -1) then
            ng = ng + (nr-1)*(ns-1)
          elseif(ntyp .eq.  8) then
            ng = ng - ((nr-1)*(ns-1))/4
          endif
        endif
        if(nf.gt.numel.and.ne.gt.0) go to 401
        if(ng.gt.numnp) go to 400

c       Generate nodes

        call sblkn(nr,ns,xl,ixl,shp,x,dr,ds,ni,n,ndm,nodinc,ntyp,
     &             nm,ctype,prt,prth)

c       Generate elements

        call sblke(nr,ns,x,ix,ni,ne,n,ndm,nen1,nodinc,ntyp,nm,ma,
     &             dlayer,mr(np(111)),ctype)

c     3-d generations

      elseif(ntyp.lt.20) then
        if(pcomp(ctype,'cart',4)) then
          call ck3dblk(ixl,xl,shp,8,3,err)
          if(err) call plstop()
        endif
        dt = 2.d0/nt
        if(ntyp.eq.10) then                     ! 8-node hexahedron
          nf = ne + nr*ns*nt - 1
        elseif(ntyp.eq.11) then                 ! 4-node tetrahedron
          nf = ne + 6*nr*ns*nt - 1
        elseif(ntyp.eq.12 .or. ntyp.eq.14) then ! 20 & 27-node hexahedron
          nf = ne + (nr*ns*nt)/8 - 1
        elseif(ntyp.eq.13) then                 ! 10-node tetrahedron
          nf = ne + 6*(nr*ns*nt)/8 - 1
        elseif(ntyp.eq.15) then                 ! 11-node tetrahedron
          nf = ne + 6*(nr*ns*nt)/8 - 1
        else
          write(ilg,4003) ntyp
          write(iow,4003) ntyp
          call plstop()
        endif
        if(nf.gt.numel.and.ne.gt.0) go to 401
        nr = nr + 1
        ns = ns + 1
        nt = nt + 1
        ng = nr*ns*nt + ni -1
        if(ntyp.eq.15) then                 ! 11-node tetrahedron
          ng = ng + 6*nr*ns*nt
        endif
        if(ng.gt.numnp) go to 400

c       Compute node locations

        call vblkn(nr,ns,nt,xl,x,ixl,dr,ds,dt,
     &             ni,ndm,ctype,prt,prth)

c       Compute element connections

        if(ne.gt.0) then
          call vblke(nr,ns,nt,x,ix,ni,ne,nf,ndm,nen1,ma,ntyp,
     &               dlayer,mr(np(111)))
        endif

c     Shell:

      elseif(ntyp.lt.30) then
        if (ntyp.eq.20) then
          nf = ne + nr*ns - 1
        elseif (ntyp.eq.27) then
          nf = ne + (nr*ns)/2 - 1
        elseif (ntyp.ge.28) then
          nf = ne + (nr*ns)/4 - 1
        else
          nf = ne + 2*nr*ns - 1
        endif
        if(nf.gt.numel.and.ne.gt.0) go to 401

c       Determine last node number to be generated

        nr = nr + 1
        ns = ns + 1
        ng = nr*ns + ni -1
        if(ng.gt.numnp) go to 400

c       If rotational updates have not been set yet:

        if (.not.frotas) then
          frotas = .true.
          errck = palloc(81,'MO   ',numnp*2 ,1)
          errck = palloc(83,'MT   ',numnp   ,2)
          errck = palloc(82,'MR   ',numnp*54,2)
        endif

c       Form block of elements

        call hblk(nr,ns,xl,ixl,shp,x,ix,dr,ds,ni,ne,ndm,
     &            nen1,nodinc,ntyp,nm,ma,hr(np(82)),ctype,prt,prth)

c     User generations

      else
        dt = 2.d0/nt
        nf = ne
        ng = ni
        call ublk(td,nn,nr,ns,nt,xl,x,ixl,ix,dr,ds,dt,
     &            ng,nf,ndm,nen1,ma,ctype,prt, 2)
        if(nf.gt.numel.and.ne.gt.0) go to 401
        if(ng.gt.numnp) go to 400
      endif

c     Set old numbers

      if(ne.gt.0) neo = nf
      nio = ng
      mao = ma

c     Set node type to active

      do n = ni,ng
        mr(np(190)+n-1) = 0
      end do ! n

c     Set region and rigid body number

      do n = ne,nf
        ix(nen1-1,n) = nreg
        if(nrigid.ge.1) then
          rben(n) = nrigid
        else
          rben(n) = -1
        end if
      end do ! n

      if(nlay.gt.0) then
        errck = palloc(111,'TEMP1',0,1)
      endif

c     Print lists if wanted

      if(prt.and.ne.gt.0) then
        do n = ne,nf,50
          call prtitl(prth)
          write(iow,2003) (i,i=1,nen)
          if(ior.lt.0) then
              write(  *,2003) (i,i=1,nen)
          endif
          j = min(nf,n+49)
          do i = n,j
            ma = ix(nen1,i)
            write(iow,2004) i,ma,nreg,(ix(k,i),k=1,nen)
            if(ior.lt.0) then
              write(  *,2004) i,ma,nreg,(ix(k,i),k=1,nen)
            endif
          end do ! i
        end do ! n
      endif

      return

c     Error messages

400   write(ilg,4000) ng,numnp
      write(iow,4000) ng,numnp
      if(ior.lt.0) write(  *,4000) ng,numnp
      return

401   write(ilg,4001) nf,numel
      write(iow,4001) nf,numel
      if(ior.lt.0) write(  *,4001) nf,numel
      return

402   write(ilg,4002) l
      write(iow,4002) l
      if(ior.lt.0) write(*,4002) l

c     Formats

2000  format('   N o d e   G e n e r a t i o n s'//
     &   10x,'Number of xi_1-increments ',i5/
     &   10x,'Number of xi_2-increments ',i5/
     &   10x,'Number of xi_3-increments ',i5/
     &   10x,'First node number         ',i5/
     &   10x,'First element number      ',i5/
     &   10x,'Element material number   ',i5/
     &   10x,'Node line increment       ',i5/
     &   10x,'Block type (0-40)         ',i5/1x)

2001  format(i9,1p,3e12.3)

2002  format(/5x,'Block coordinates specified as: ',a,//,
     &       5x,'node',3(i6,a6))

2003  format('   E l e m e n t   C o n n e c t i o n s'//
     &   '   Elmt Mat Reg',8(i3,' node'):/(15x,8(i3,' node')))

2004  format(i7,2i4,8i8:/(15x,8i8))

2005  format(/5x,'Layered Material Properties'/
     &      /(7x,4(i2,'-Layer Matl')))

2006  format(7x,4(i8,i5))

3000  format(' *WARNING* BLKGEN: No elements are generated ')

4000  format(' *ERROR* BLKGEN: Insufficient storage for nodes'/
     &   10x,'final node =',i5,5x,'numnp =',i5)

4001  format(' *ERROR* BLKGEN: Insufficient storage for elements'/
     &   10x,'final element =',i5,5x,'numel =',i5)

4002  format(' *ERROR* BLKGEN: Block node has number > 27.  No. =',i8)

4003  format(' *ERROR* BLKGEN: Block type incorrect: Ntype =',i8)

5000  format(' Input: type,nr,ns,ni,ne,ma,nodinc,ntyp')

5001  format(' Input: type,nr,ns,nt,ni,ne,ma,ntyp')

5002  format(' Input: node, x-1, x-2, x-3')

      end
