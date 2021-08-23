c$Id: pmesh.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pmesh(iii,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 244 to 256 in palloc of 'NFORC' and
c          change np(244) to np(256) in call genvec         02/11/2006
c       2. Change call pprint('  Enter ...') to write(*,*) 'Enter ...
c          Set xxx(1:13) correctly and remove format 2000.  13/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Data input routine for mesh description

c      Inputs:
c         iii        - Initialization indicator
c         prt        - Flag, print input data if true
c         prth       - Flag, print title/header if true

c      Outputs:
c         Depends on commands specified
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cblend.h'
      include  'cblktr.h'
      include  'crotas.h'
      include  'cdata.h'
      include  'cdat1.h'
      include  'cdat2.h'
      include  'chdata.h'
      include  'codat.h'
      include  'corset.h'
      include  'corfil.h'
      include  'cornum.h'
      include  'comfil.h'
      include  'debugs.h'
      include  'dstars.h'
      include  'edgdat.h'
      include  'eldata.h'
      include  'eqslv.h'
      include  'hlpdat.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'iosave.h'
      include  'linka.h'
      include  'mdata.h'
      include  'modreg.h'
      include  'pdata3.h'
      include  'pointer.h'
      include  'prflag.h'
      include  'prld1.h'
      include  'region.h'
      include  'rigid1.h'
      include  'rigid2.h'
      include  'sdata.h'
      include  'sldata.h'
      include  'trdata.h'
      include  'umac1.h'
      include  'comblk.h'

      include  'p_int.h'

      integer   list   , numesh
      parameter(list=90, numesh=10)

      logical   setvar,palloc, lp_in, lopen
      logical   prt,error,pcomp,lmesh,errck,pinput,tinput,vinput
      logical   prth, umesh, readfl, savefl, snofl,flgco
      character lp_file*10
      character wd(list)*4, cc*4,c2*8,fext*8,tx(2)*15
      integer   ed(list)  , i,j,ii,jj, iii, isd, ibn, lmr, side, face
      integer   ll,llo,lp_lun
      integer   nsblk,nblend, starsvn,starsve
      real*8    td(6)

      save

c     List of command names

      data wd/'coor','elem','mate','boun','forc','temp','end ','prin',
     1        'nopr','titl','bloc','pola','ebou','angl','sloa','cons',
     2        'sphe','btem','icon','pars','nopa','trib','para','efor',
     3        'eang','cbou','cfor','cang','foll','slav','rese','sblo',
     4        'esse','rota','setn','setr','btra','fpro','cpro','regi',
     5        'tran','damp','mass','stif','csur','ereg','reac','manu',
     6        'body','glob','shif','disp','edis','cdis','debu','side',
     7        'face','snod','blen','move','rigi','moda','flex','base',
     8        'epro','mpro','loop','next','file','cdam','cmas','csti',
     9        'ebas','cbas','eule','ceul','rfor','lfor','*nod','*ele',
     u        'mes1','mes2','mes3','mes4','mes5','mes6','mes7','mes8',
     u        'mes9','mes0'/

      data ed/    0,     0,     0,     0,     0,     1,     0,     0,
     1            0,     2,     0,     0,     0,     0,     6,     6,
     2            0,     1,     6,     3,     3,     2,     0,     0,
     3            0,     0,     0,     0,     1,     1,     1,     2,
     4            1,     2,     2,     2,     3,     1,     1,     1,
     5            1,     1,     1,     1,     0,     1,     2,     4,
     6            1,     0,     0,     0,     0,     0,     0,     0,
     7            6,     0,     0,     2,     1,     1,     1,     0,
     8            0,     1,     3,     3,     3,     0,     0,     0,
     9            2,     2,     1,     1,     1,     1,     2,     2,
     u            5,     5,     5,     5,     5,     5,     5,     5,
     u            5,     5/

c     Write message to log file

      if(iii.ge.0) write(ilg,4000) finp

c     Initialize arrays and set error detection values

      lp_in = .true.
      error = .false.
      lmesh = .false.
      ll    = 1
      nneq  = ndf*numnp
      if(iii.ge.0) then
        frotas      = .false.
        nffl        = .false.
        nmfl        = .false.
        snofl       = .false.
        anglefl     = .false.
        eulerfl     = .false.
        nreg        = 0
        mxreg       = 0
        nblend      = 0
        npstr       = 0
        nrigid      = 0
        numg        = 0
        numsl       = 0
        nei         = 0
        ner         = 0
        nio         = 0
        neo         = 0
        mao         = 0
        ibn         = 20
        isd         = 16
        side        = 0
        face        = 0

c       Set angles, boundary code/forced values to zero

        do j = 1,3
          x0(j)    = 0.0d0
          xr(j)    = 0.0d0
          do i = 1,3
            tr(i,j) = 0.0d0
          end do ! i
          tr(j,j) = 1.d0
        end do ! j
        trdet = 1.d0
        if(iii.eq.0) then
          prt = .true.

c         Set node type to undefined

          do j = 0,numnp-1
            mr(np(190)+j) = -1
          end do ! j

          do j = 1,50
            prldv(j) = 1.d0
          end do ! j

c         Insert user functions into list

          do j = 1,10
            i = list-numesh+j
            if(.not.pcomp(umshc(j),wd(i),4)) then
              wd(i) = umshc(j)
              ed(i) = 0
            endif
          end do ! j

        endif

      elseif(iii.eq.-2) then
        nio   = 0
        neo   = 0
        mao   = 0
        nreg  = 0
      endif
100   if(ior.lt.0) then
        write(*,*) ' Enter "help" for list of commands, "end" to exit'
        xxx       = ' '
        xxx(1:13) = '    Mesh    >'
        write(xxx(10:12),'(i3)') ll
        call pprint(xxx)
      endif
      errck = tinput(tx,2,td,0)
      if(errck) go to 100
      utx(1) = tx(1)
      utx(2) = tx(2)
      cc     = tx(1)
      c2     = tx(2)
      if( pcomp(cc,'read',4) ) then
        lmesh = readfl(tx(2))
        if(lmesh) then
          llo = ll
        else
          ll  = llo
        endif
        go to 100
      endif
      if(pcomp(cc,'save',4)) then
        lsave = savefl(tx(2))
        go to 100
      endif
      if(ior.lt.0.and.pcomp(cc,'help',4)) then
        call phelp(c2,wd,ed,list,'MESH')
        go to 100
      endif
      go to 120
110   call  errclr ('PMESH ')
      go to 100
120   do i = 1,list
        if(pcomp(cc,wd(i),4)) go to 130
      end do ! i

c     User mesh commands

      if(.not. pcomp( cc, ' ', 1 ) ) errck = umesh(cc,prt)
      if(.not. errck .and. ior.lt.0) call errclr('PMESH ')
      go to 100
130   ll = ll + 1
      if(iii.ge.0) write(ilg,4001) tx(1),tx(2)

c     [coor]dinates  - Nodal coordinate data input
c     [coor]<all>    - Nodal coordinate data input (no parsing)
c     [coor]<add>    - Accumulate the numbers starting at nio
c     [coor]<bina>ry - Input coordinates in binary mode

      if(i.eq.1) then
        if(pcomp(c2,'all',3)) then
          call pcrdrd(hr(np(43)),mr(np(190)),prt)
        elseif(pcomp(c2,'bina',4)) then
          call pcrbrd(hr(np(43)),mr(np(190)),prt)
        else
          starsvn = starnd
          if(pcomp(c2,'add',3)) then
            starnd  = nio
          endif
          call genvec(ndm,ndm,hr(np(43)),' Coordinates',
     &                prt,prth,error,.true.)
          if(nio.gt.starnd) then
            numnp = max(numnp,nio)
          endif
          starnd = starsvn
        endif

c     [elem]ent      - Data input
c     [elem,all]     - Input all connection data
c     [elem,add]     - Input all connection data
c     [elem<bina>ry  - Input element connections in binary mode
c     [elem,gene,xxx]- Set generation array to xxx

      elseif(i.eq.2) then
        if(pcomp(c2,'all',3)) then
          call pelmrd(mr(np(33)),mr(np(181)),prt)
        elseif(pcomp(c2,'bina',4)) then
          call pelbrd(mr(np(33)),mr(np(181)),prt)
        else
          starsve = starel
          if(pcomp(c2,'add',3)) then
            starel  = neo
          endif
          call pelmin(tx(2),mr(np(34)),mr(np(33)),mr(np(181)),nen1,
     &                prt,prth,error)
          if(neo.gt.starel) then
            numel = max(numel,neo)
          endif
          starel = starsve
        endif

c     [mate]rial,ma: Data input for material set ma

      elseif(i.eq.3) then
        setvar = palloc(151,'USER1', ndf*(nen+1),1)
        call pmatin(tx,hr(np(25)),hr(np(41)),hr(np(44)),hr(np(39)),
     &                 hr(np(36)),hr(np(35)),mr(np(34)),mr(np(32)),
     &                 mr(np(240)),mr(np(151)),prt,prth)
        setvar = palloc(151,'USER1', 0          ,1)

c     [boun]dary codes,<set,add> - input nodal restraint conditions

      elseif(i.eq.4) then
        if(pcomp(c2,'add',3)) then
          setvar = palloc(151,'USER1', ndf*numnp,1)
          call pzeroi(mr(np(151)),ndf*numnp)
          fp(1) = np(151)
          fp(2) = np(31) + nneq
        else
          fp(1) = np(31) + nneq
        endif
        call pbouin(mr(np(34)),mr(fp(1)),prt,prth)
        if(pcomp(c2,'add',3)) then
          do ii = 0,ndf*numnp-1
            mr(ii+fp(2)) = mr(ii+fp(2)) + mr(ii+fp(1))
          end do ! ii
          setvar = palloc(151,'USER1', 0,1)
        endif

c     [forc]e,<set,add>: Data input for nodal generalized forces

      elseif(i.eq.5) then
        if(pcomp(c2,'add',3)) then
          setvar = palloc(151,'USER1', ndf*numnp,2)
          call pzero(hr(np(151)),ndf*numnp)
          fp(1) = np(151)
          fp(2) = np(27)
        else
          fp(1) = np(27)
        endif
        call genvec(ndf,ndf,hr(fp(1)),' Forces',prt,prth,error,.false.)
        if(pcomp(c2,'add',3)) then
          do ii = 0,ndf*numnp-1
            hr(ii+fp(2)) = hr(ii+fp(2)) + hr(ii+fp(1))
          end do ! ii
          setvar = palloc(151,'USER1', 0,2)
        endif

c     [temp]erature data input

      elseif(i.eq.6) then
        call genvec(1,1,hr(np(38)),' Temperatures',prt,prth,error,
     &             .false.)

c     [end] of mesh data inputs

      elseif(i.eq.7) then
        if(lsave) then
          if(iii.ge.0) write(ilg,3006)
          write(iow,3006)
          if(ior.lt.0) then
            write(iow,3006)
          endif
          error = .true.
        endif

        close(unit=iwd)

c       Complete surface load sets

        call ploadi(prt,prth,error,2)

c       Exit if any errors occurred in mesh data

        if(error) then
          call plstop()

        elseif(iii.ge.0) then

c         Perform delayed mesh generation steps

          if(numbd.gt.0 .and. iii.ge.0) then
            if(snofl) then
              call pblendm(isd,ibn,ndm,nen1,prt,prth,.true.,.true.)
            else
              if(ior.le.0) then
                write(*,3007)
              endif
              write(ilg,3007)
              write(iow,3007)
              call plstop()
            endif
          endif

c         Initialize transformation matrices for all configurations

          call chkrot(mr(np(32)),mr(np(33)),nie,nen1)
          if (frotas) then
            setvar = palloc(111,'TEMP1', 3*numnp, 2)
            call pshsurf(hr(np(43)),hr(np(82)),hr(np(111)),
     &                   mr(np(33)),mr(np(32)),mr(np(81)))
            setvar = palloc(111,'TEMP1', 0, 2)

            call setrot(hr(np(43)),mr(np(81)),hr(np(82)),hr(np(83)),
     &                  numnp)

c           Compress rotational data storage

            call cmprot(mr(np(81)),hr(np(82)), numnp, lmr)
            setvar = palloc(82,'MR   ', lmr*54,2)

c           Initialize rotational parameters

            call updrot(hr,ndf,hr(np(82)),mr(np(81)),numnp,0)

          endif
        endif
        return

c     [prin]t/[nopr]int of input data

      elseif(i.eq.8 .or. i.eq.9) then
        prt = i.eq.8

c     [titl] - set title prints on/off

      elseif(i.eq.10) then
        if(pcomp(c2,'off',3)) then
          prth = .false.
        else
          prth = .true.
        endif

c     [bloc]k - generate block of nodes and elements

      elseif(i.eq.11) then
        if(iii.lt.0) write(iow,3004)
        call blkgen(ndm,nen1,hr(np(43)),mr(np(33)),mr(np(181)),prt,prth)

c     [pola]r - convert polar to cartesian coordinates

      elseif(i.eq.12) then
        call polar(mr(np(190)),hr(np(43)),ndm,prt,prth)

c     [ebou] - set edge boundary constraints

      elseif(i.eq.13) then
        ebcfl = .true.
        fext  = 'co0'
        if(nebcs.le.9) then
          write(fext(3:3),'(i1)') nebcs
        elseif(nebcs.le.99) then
          write(fext(2:3),'(i2)') nebcs
        endif
        nebcs = nebcs + 1
        if(.not.pcomp(c2,'set',3)) then
          c2 = 'add'                 ! default mode 'add' for edge b.c.
        endif
        call plinka(fext,c2,'   ')

c     [angl]e - set boundary angles

      elseif(i.eq.14) then
        call genvec(1,1,hr(np(45)),' Angles',prt,prth,error,.false.)
        anglefl = .true.

c     [sloa]ds - set surface loadings

      elseif(i.eq.15) then
        call ploadi(prt,prth,error,1)

c     [para]meter - set parameter variables

      elseif(i.eq.16 .or. i.eq.23) then
        coflg = .true.
        call pconst(prt)

c     [sphe]re - convert spherical to cartesian coordinates

      elseif(i.eq.17) then
        call sphere(mr(np(190)),hr(np(43)),ndm,prt,prth)

c     [btem] - input block of interpolated temperatures

      elseif(i.eq.18) then
        call blktem(ndm,hr(np(38)),prt,prth)

c     [icon] - input and check contact data

      elseif(i.eq.19) then
c       call incon(hr(np(43)),ndm,prt,error)

c     [pars]ing/[nopa]rsing of statements

      elseif(i.eq.20 .or. i.eq.21) then
        coflg = i.eq.20

c     [trib] - triangular block generator

      elseif(i.eq.22) then
        if(iii.lt.0) write(iow,3004)
        call blktri(ndm,nen1,hr(np(43)),mr(np(33)),prt,prth)

c     [efor] set edge force constraints

      elseif(i.eq.24) then
        efcfl = .true.
        fext  = 'ld0'
        if(nefrc.le.9) then
          write(fext(3:3),'(i1)') nefrc
        elseif(nefrc.le.99) then
          write(fext(2:3),'(i2)') nefrc
        endif
        nefrc = nefrc + 1
        if(.not.pcomp(c2,'set',3)) then
          c2 = 'add'                 ! default mode 'add' for edge force
        endif
        call plinka(fext,c2,'   ')

c     [eang] set edge angle constraints

      elseif(i.eq.25) then
        anglefl = .true.
        eanfl   = .true.
        fext    = 'ww0'
        if(neang.le.9) then
          write(fext(3:3),'(i1)') neang
        elseif(neang.le.99) then
          write(fext(2:3),'(i2)') neang
        endif
        neang = neang + 1
        if(.not.pcomp(c2,'set',3)) then
          c2 = 'add'                 ! default mode 'add' for edge angle
        endif
        call plinka(fext,c2,'   ')

c     [cbou] set coordinate boundary constraints - based on coordinates

      elseif(i.eq.26) then
        boufl = .true.
        fext  = 'bn0'
        if(nbouf.le.9) then
          write(fext(3:3),'(i1)') nbouf
        elseif(nbouf.le.99) then
          write(fext(2:3),'(i2)') nbouf
        endif
        nbouf = nbouf + 1
        call plinka(fext,c2,'   ')

c     [cfor] set coordinate nodal forces - based on coordinates

      elseif(i.eq.27) then
        forfl = .true.
        fext  = 'fr0'
        if(nforf.le.9) then
          write(fext(3:3),'(i1)') nforf
        elseif(nforf.le.99) then
          write(fext(2:3),'(i2)') nforf
        endif
        nforf = nforf + 1
        call plinka(fext,c2,'   ')

c     [cang] set coordinate angle - based on coordinates

      elseif(i.eq.28) then
        anglefl = .true.
        angfl   = .true.
        fext    = 'an0'
        if(nangf.le.9) then
          write(fext(3:3),'(i1)') nangf
        elseif(nangf.le.99) then
          write(fext(2:3),'(i2)') nangf
        endif
        nangf = nangf + 1
        call plinka(fext,c2,'   ')

c     [foll]ower loads

      elseif(i.eq.29) then
        fext = 'tem'
        call plinka(fext,c2,'   ')
        nfol = iclink

        setvar = palloc(129,'FOLLI ',nfol*2,1)
        setvar = palloc(130,'FOLLR ',nfol,  2)
c       call pfollo(nfol,hr(np(130)),mr(np(129)),prt,prth)

c     [slav]ed nodes to subspace matrix

      elseif(i.eq.30) then
30      errck = pinput(td,1)
        j     = td(1)
        if(j.gt.0 .and. j.le.numnp) then
          errck = palloc(205,'NSLAV',numg+1,1)
          mr(np(205)+numg) = j
          numg             = numg + 1
        else
          call iprint(mr(np(205)),1,numg,1,'Slaved Nodes')
          go to 100
        endif
        go to 30

c     [rese]t boundary condition codes to zero - permits releases

      elseif(i.eq.31) then
        call pzeroi(mr(np(31)+nneq),ndf*numnp)

c     [sblo,nsblk] - input surface blocks to generate 3-d meshes

      elseif(i.eq.32) then
        cc    = yyy(1:4)
        errck = vinput(yyy(16:30),15,td,1)
        nsblk = max(1,nint(td(1)))

c       Allocate arrays for 3-d generation from surface mesh

        setvar = palloc(111,'TEMP1', 9*nsblk,  1)
        setvar = palloc(112,'TEMP2', 9*nsblk,  1)
        setvar = palloc(113,'TEMP3', 5*nsblk,  1)
        setvar = palloc(114,'TEMP4',27*nsblk,  2)
        setvar = palloc(115,'TEMP5', 9*nsblk,  2)
        setvar = palloc(116,'TEMP6',27*nsblk,  2)

        call sblkgn(nsblk,ndm,nen,nen1,hr(np(43)),mr(np(33)),
     &              mr(np(111)),mr(np(112)),mr(np(113)),hr(np(114)),
     &              hr(np(115)),hr(np(116)),prt,prth)

c       Delete arrays used to generate mesh

        setvar = palloc(116,'TEMP6', 0,  2)
        setvar = palloc(115,'TEMP5', 0,  2)
        setvar = palloc(114,'TEMP4', 0,  2)
        setvar = palloc(113,'TEMP3', 0,  1)
        setvar = palloc(112,'TEMP2', 0,  1)
        setvar = palloc(111,'TEMP1', 0,  1)

c     [esse]ntial - b.c. set option
c     [node1,node2,inc,idl(i),i=1,ndf]

      elseif(i.eq.33) then
        call pessbc(mr(np(34)),mr(np(31)+nneq),ndf,numnp,prt,prth)

c     [rota] - set flags for rotational transformations

      elseif(i.eq.34) then
        if (.not.frotas) then
          frotas = .true.
          setvar = palloc(81,'MO   ',numnp*2 ,1)
          setvar = palloc(83,'MT   ',numnp   ,2)
          setvar = palloc(82,'MR   ',numnp*54,2)
        endif

        call protin(prt,prth)

c     [setn] - define second mesh for director input
c            & set initial normal quantities

      elseif(i.eq.35) then
        call gendir(hr(np(82)),' Director Nodes',prt,prth,error,.true.)

c     [setr] -set initial rotation transformations from nodal directors

      elseif(i.eq.36) then
        call setrot(hr(np(43)),mr(np(81)),hr(np(82)),hr(np(83)),numnp)

c       Read input

        call genvec(9,54,hr(np(82)),' Transforms ',prt,prth,error,
     &             .true.)

c       Set default rotation matrices

        call defrot(hr(np(82)),hr(np(27)),mr(np(31)+nneq),numnp,ndf)

c       Initialize all configurations

        call pmove(hr(np(82)),hr(np(82)+numnp*9 ),numnp*9)
        call pmove(hr(np(82)),hr(np(82)+numnp*18),numnp*9)
        call pmove(hr(np(82)),hr(np(82)+numnp*45),numnp*9)

c     [btra] - transition mesh patch

      elseif(i.eq.37) then
        if(iii.lt.0) write(iow,3004)
        call blktra(ndm,nen1,hr(np(43)),mr(np(33)),prt,prth)

c     [fpro],<set,add> - Proportional load number specification

      elseif(i.eq.38) then
        if(pcomp(c2,'add',3)) then
          setvar = palloc(151,'USER1', ndf*numnp,1)
          call pzeroi(mr(np(151)),ndf*numnp)
          fp(1) = np(151)
          fp(2) = np(29)
        else
          fp(1) = np(29)
        endif
        call genint(ndf,mr(fp(1)),ndf,numnp,'P r o p.  L o a d  N o s.',
     &              '-dof',prt,prth,error,1)
        if(pcomp(c2,'add',3)) then
          do ii = 0,ndf*numnp-1
            mr(ii+fp(2)) = mr(ii+fp(2)) + mr(ii+fp(1))
          end do ! ii
          setvar = palloc(151,'USER1', 0,1)
        endif

c     [cpro]  coordinate proportional load number specification

      elseif(i.eq.39) then
        cprfl = .true.
        fext  = 'yp0'
        if(nprof.le.9) then
          write(fext(3:3),'(i1)') nprof
        elseif(nprof.le.99) then
          write(fext(2:3),'(i2)') nprof
        endif
        nprof = nprof + 1
        call plinka(fext,c2,'   ')

c     [regi,nreg]  set region number: all

      elseif(i.eq.40) then
        cc    = yyy(1:4)
        errck = vinput(yyy(16:30),15,td,1)
        if(errck) go to 110
        nreg  = nint(td(1))
        mxreg = max(mxreg,nreg)
        write(iow,2001) nreg,mxreg
        if(ior.lt.0) then
          write(*,2001) nreg,mxreg
        endif

c     [tran],<inc> - Specify coordinate transformation array

c     Incremental rotation update: tr_new = tinc * tr_old
c                                  xr_new = tinc * xr_old + xr_inc
      elseif(i.eq.41) then
        call ptranf(tx(2),prt)

c     [damp],<set,add> - specify nodal damping quantities
c     [mass],<set,add> - specify nodal mass quantities
c     [stif],<set,add> - specify nodal stiffness quantities

      elseif(i.eq.42 .or. i.eq.43 .or. i.eq.44) then
        if(.not.nmfl) then
          setvar = palloc(88,'NSTI  ',ndf*numnp,2)
          setvar = palloc(87,'NMAS  ',ndf*numnp,2)
          setvar = palloc(86,'NDAM  ',ndf*numnp,2)
          nmfl = .true.
        endif
        if(pcomp(c2,'add',3)) then
          setvar = palloc(151,'USER1', ndf*numnp,2)
          call pzero(hr(np(151)),ndf*numnp)
          fp(1) = np(151)
          fp(2) = np(i+44)
        else
          fp(1) = np(i+44)
        endif
        if(i.eq.42) then
          call genvec(ndf,ndf,hr(fp(1)),' Nodal Dampers',prt,prth,
     &                error,.false.)
        elseif(i.eq.43) then
          call genvec(ndf,ndf,hr(fp(1)),' Nodal Masses' ,prt,prth,
     &                error,.false.)
        elseif(i.eq.44) then
          call genvec(ndf,ndf,hr(fp(1)),' Nodal Springs',prt,prth,
     &                error,.false.)
        endif
        if(pcomp(c2,'add',3)) then
          do ii = 0,ndf*numnp-1
            hr(ii+fp(2)) = hr(ii+fp(2)) + hr(ii+fp(1))
          end do ! ii
          setvar = palloc(151,'USER1', 0,2)
        endif

c     [csur] - surface loading by coordinates

      elseif(i.eq.45) then
        surfl = .true.
        fext  = 'sl0'
        if(nsurf.le.9) then
          write(fext(3:3),'(i1)') nsurf
        elseif(nsurf.le.99) then
          write(fext(2:3),'(i2)') nsurf
        endif
        nsurf = nsurf + 1
        call plinka(fext,c2,'   ')

c     [ereg] - set element regions

      elseif(i.eq.46) then
        call genint(1,mr(np(33)+nen1-2),nen1,numel,'R e g i o n  N o s',
     &              '-regn',prt,prth,error,2)
        do j = nen1-2,numel*nen1,nen1
          mxreg = max(mxreg,mr(np(33)+j))
        end do ! j

c     [reac,file] - input reactions from 'file'

      elseif(i.eq.47) then
        inquire(file=tx(2),exist=reafl)
        if(.not.reafl) then
          if(ior.lt.0) then
            write(*,3005) tx(2)
          else
            if(iii.ge.0) write(ilg,3005) tx(2)
            write(iow,3005) tx(2)
            call plstop()
          endif
        else
          reafi = tx(2)
        endif

c     [manu],hlplev - set Manual help options level

      elseif(i.eq.48) then
        cc    = yyy(1:4)
        errck = vinput(yyy(16:30),15,td,1)
        if(errck) go to 110
        hlplev = max(-1,min(3,nint(td(1))))

c     [body] - compute internal forces, will be stored in nodal loads

      elseif(i.eq.49) then
        intfl = .true.
        fext  = 'in0'
        if(nintf.le.9) then
          write(fext(3:3),'(i1)') nintf
        elseif(nintf.le.99) then
          write(fext(2:3),'(i2)') nintf
        endif
        nintf = nintf + 1
        call plinka(fext,c2,'   ')

c     [glob]al - set global parameters

      elseif(i.eq.50) then
        call global()

c     [shif]t:<x0,y0,z0> - origin for polar/spherical conversions

      elseif(i.eq.51) then
51      errck = pinput(x0,3)
        if(errck) go to 51
        if(prt) then
          write(iow,2002) (x0(i),i=1,ndm)
          if(ior.lt.0) then
            write(*,2002) (x0(i),i=1,ndm)
          endif
        endif

c     [disp]l,<set,add> - data input

      elseif(i.eq.52) then
        if(pcomp(c2,'add',3)) then
          setvar = palloc(151,'USER1', ndf*numnp,2)
          call pzero(hr(np(151)),ndf*numnp)
          fp(1) = np(151)
          fp(2) = np(27) + ndf*numnp
        else
          fp(1) = np(27) + ndf*numnp
        endif
        call genvec(ndf,ndf,hr(fp(1)),' Displacements',
     &              prt,prth,error,.false.)
        if(pcomp(c2,'add',3)) then
          do ii = 0,ndf*numnp-1
            hr(ii+fp(2)) = hr(ii+fp(2)) + hr(ii+fp(1))
          end do ! ii
          setvar = palloc(151,'USER1', 0,2)
        endif

c     [edis] set edge displacement constraints

      elseif(i.eq.53) then
        edifl = .true.
        fext  = 'ud0'
        if(nedis.le.9) then
          write(fext(3:3),'(i1)') nedis
        elseif(nedis.le.99) then
          write(fext(2:3),'(i2)') nedis
        endif
        nedis = nedis + 1
        if(.not.pcomp(c2,'set',3)) then
          c2 = 'add'                 ! default mode 'add' for edge displ
        endif
        call plinka(fext,c2,'   ')

c     [cdis] set coordinate nodal forces - based on coordinates

      elseif(i.eq.54) then
        disfl = .true.
        fext  = 'ds0'
        if(ndisf.le.9) then
          write(fext(3:3),'(i1)') ndisf
        elseif(ndisf.le.99) then
          write(fext(2:3),'(i2)') ndisf
        endif
        ndisf = ndisf + 1
        call plinka(fext,c2,'   ')

c     [debu]g,<on,off>  Activate,deactivate debug option

      elseif(i.eq.55) then
        if(pcomp(c2,'off',3)) then
          debug = .false.
        else
          debug = .true.
        endif

c     [side] - for blending interpolations

      elseif(i.eq.56) then
        if(iii.ge.0) then
          numsd  = max(numsd,1)
          setvar = palloc ( 162, 'BSIDE', numsd*isd, 1)
          call psides(mr(np(162)),side,isd,prt,prth,1)
          setvar = palloc ( 162, 'BSIDE', numsd*isd, 1)
        else
          write(*,3001)
        endif

c     [face] - for blending interpolations

      elseif(i.eq.57) then
        if(iii.ge.0) then
          setvar = palloc ( 165, 'BFACE', numfc*isd, 1)
          call psides(mr(np(165)),face,isd,prt,prth,2)
        else
          write(*,3001)
        endif

c     [snod]e - for blending interpolations

      elseif(i.eq.58) then
        if(iii.ge.0) then
          numsn  = max(numsn,1)
          setvar = palloc ( 161, 'BNODE', numsn*3, 2)
          call pnodes(hr(np(161)),ndm,prt,prth)
          setvar = palloc ( 161, 'BNODE', numsn*3, 2)
          snofl  = .true.
        else
          write(*,3002)
        endif

c     [blen]ding interpolations (Delayed mesh generation feature)

      elseif(i.eq.59) then
        if(iii.ge.0) then
          nblend = nblend + 1
          numbd  = max(numbd,nblend)
          setvar = palloc ( 163, 'BTRAN', numbd*12          , 2)
          setvar = palloc ( 164, 'BLEND', numbd*ibn         , 1)
          setvar = palloc ( 166, 'BNILR', numbd*max(1,mxilr), 1)
          call pblend(hr(np(163)),mr(np(164)),mr(np(166)),nblend,ibn,
     &                ndm,prt,prth)
        else
          write(*,3003)
        endif

c     [move]  relocate coordinates based on some specified locations

      elseif(i.eq.60) then
        call pcormv(hr(np(43)),ndm,numnp,prt,prth)

c     [rigi,nrigid,rbcen(nrigid),rbx0(i,nrigid)]
c     [moda,nrigid,rbcen(nrigid),rbx0(i,nrigid)]

c     Set rigid body number/type
c     rbtype(nrigid)  = 0   :  rigid-flex (default)
c                     = 1   :  rigid-modal
c     rbcen(nrigid)   = 0   :  At center of mass
c                     = >0  :  At rbx0
c     rbx0(ii,nrigid) = x(i):  reference location for dofs

      elseif(i.eq.61 .or. i.eq.62) then
        cc     = yyy(1:4)
        errck  = vinput(yyy(16:105),90,td,6)
        if(errck) go to 110
        nrigid = nint(td(1))

c       Return to non-rigid mode
        if(nrigid.eq.-1) then
          nrigid = 0
          write(iow,2006)
          if(ior.lt.0) write(iow,2006)
          go to 100

c       Set rigid body number if none is given by user

        elseif(nrigid.eq.0) then
          nrigid = nrbody + 1
        end if
        nrbody = max(nrbody,nrigid)

c       Set rigid/modal indicator

        if(i.eq.61) then
          rbtype(nrigid) = 0
        else
          rbtype(nrigid) = 1
        endif

c       Set location for dofs

        rbcen(nrigid)  = nint(td(2))
        if(rbcen(nrigid).gt.0) then
          do ii = 1,ndm
            rbx0(ii,nrigid) = td(ii+2)
          end do ! ii
        endif

c       Set modal body number

        if(rbtype(nrigid).eq.1) then
          nmbody = nmbody + 1
          modbod(nmbody) = nrigid
          write(iow,2003) nrigid
          if(ior.lt.0) then
            write(*,2003) nrigid
          endif
        else
          write(iow,2004) nrigid
          if(ior.lt.0) then
            write(*,2004) nrigid
          endif
        end if
        if(rbcen(nrigid).gt.0) then
          write(iow,2005) (rbx0(ii,nrigid),ii=1,ndm)
          if(ior.lt.0) then
            write(*,2005) (rbx0(ii,nrigid),ii=1,ndm)
          endif
        end if

c     [flex]ible - set element types to deformable (non-rigid)

      elseif(i.eq.63) then
        nrigid = 0
        write(iow,2006)
        if(ior.lt.0) write(iow,2006)

c     [base]  proportional load number specification

      elseif(i.eq.64) then
        if(np(125).eq.0) then
          setvar = palloc(125,'NUBAS',ndf*numnp,1)
        endif
        call genint(ndf,mr(np(125)),ndf,numnp,
     &             'B a s e   P a t t e r n s',
     &            '-dof',prt,prth,setvar,1)

c     [epro] set edge proportional load numbers

      elseif(i.eq.65) then
        eprfl = .true.
        fext  = 'ep0'
        if(nepro.le.9) then
          write(fext(3:3),'(i1)') nepro
        elseif(nepro.le.99) then
          write(fext(2:3),'(i2)') nepro
        endif
        nepro = nepro + 1
        if(.not.pcomp(c2,'set',3)) then
          c2 = 'add'                 ! default mode 'add' for edge prop
        endif
        call plinka(fext,c2,'   ')

c     [mpro],<set,add> -  mass proportional load number specification

      elseif(i.eq.66) then
        if(pcomp(c2,'add',3)) then
          setvar = palloc(151,'USER1', ndf*numnp,1)
          call pzeroi(mr(np(151)),ndf*numnp)
          fp(1) = np(151)
          fp(2) = np(29) + nneq
        else
          fp(1) = np(29) + nneq
        endif
        call genint(ndf,mr(fp(1)),ndf,numnp,
     &             'M a s s  P r o p.  L o a d','-dof',
     &              prt,prth,error,1)
        if(pcomp(c2,'add',3)) then
          do ii = 0,ndf*numnp-1
            fp(3) = ii + ii
            mr(ii+fp(2)) = mr(ii+fp(2)) + mr(fp(3))
          end do ! ii
          setvar = palloc(151,'USER1', 0,1)
        endif

c     [loop,#] - Loop start

      elseif(i.eq.67) then
        call ploops(lp_in,tx(2), 1)

c     [next] - Loop end

      elseif(i.eq.68) then
        call ploops(lp_in,tx(2), 2)

c     [file] - File inputs: Open new input file: lp_file

      elseif(i.eq.69) then
        lp_file = tx(2)
        lopen   = .true.
        lp_lun  = icl
        do while(lopen)
          lp_lun       = lp_lun + 1
          inquire(unit = lp_lun, opened = lopen)
        end do ! while
        ior       = lp_lun
        open(unit = ior, file = lp_file, status = 'old')

c     [cdam] set coordinate nodal dampers   - based on coordinates
c     [cmas] set coordinate nodal masses    - based on coordinates
c     [csti] set coordinate nodal stiffness - based on coordinates

      elseif(i.eq.70 .or. i.eq.71 .or. i.eq.72) then
        if(.not.nmfl) then
          setvar = palloc(88,'NSTI  ',ndf*numnp,2)
          setvar = palloc(87,'NMAS  ',ndf*numnp,2)
          setvar = palloc(86,'NDAM  ',ndf*numnp,2)
          nmfl = .true.
        endif

c       [cdam] set coordinate nodal dampers - based on coordinates

        if(i.eq.70) then
          damfl = .true.
          fext  = 'vd0'
          if(ndamf.le.9) then
            write(fext(3:3),'(i1)') ndamf
          elseif(ndamf.le.99) then
            write(fext(2:3),'(i2)') ndamf
          endif
          ndamf = ndamf + 1
          call plinka(fext,c2,'   ')

c       [cmas] set coordinate nodal masses - based on coordinates

        elseif(i.eq.71) then
          masfl = .true.
          fext  = 'ms0'
          if(nmasf.le.9) then
            write(fext(3:3),'(i1)') nmasf
          elseif(nmasf.le.99) then
            write(fext(2:3),'(i2)') nmasf
          endif
          nmasf = nmasf + 1
          call plinka(fext,c2,'   ')

c       [csti] set coordinate nodal stiffness - based on coordinates

        elseif(i.eq.72) then
          stifl = .true.
          fext  = 'ks0'
          if(nstif.le.9) then
            write(fext(3:3),'(i1)') nstif
          elseif(nstif.le.99) then
            write(fext(2:3),'(i2)') nstif
          endif
          nstif = nstif + 1
          call plinka(fext,c2,'   ')
        endif

c     [ebas] - set edge base proportional loads

      elseif(i.eq.73) then
        if(np(125).eq.0) then
          setvar = palloc(125,'NUBAS',ndf*numnp,1)
        endif
        ebsfl = .true.
        fext  = 'gb0'
        if(nebas.le.9) then
          write(fext(3:3),'(i1)') nebas
        elseif(nebas.le.99) then
          write(fext(2:3),'(i2)') nebas
        endif
        nebas = nebas + 1
        if(.not.pcomp(c2,'set',3)) then
          c2 = 'add'                 ! default mode 'add' for edge b.c.
        endif
        call plinka(fext,c2,'   ')

c     [cbas] set coordinate base proportional load - based on coordinates

      elseif(i.eq.74) then
        if(np(125).eq.0) then
          setvar = palloc(125,'NUBAS',ndf*numnp,1)
        endif
        basfl = .true.
        fext  = 'hb0'
        if(nbasf.le.9) then
          write(fext(3:3),'(i1)') nbasf
        elseif(nbasf.le.99) then
          write(fext(2:3),'(i2)') nbasf
        endif
        nbasf = nbasf + 1
        call plinka(fext,c2,'   ')

c     [eule]r  - angles for nodal transformations

      elseif(i.eq.75) then

        if(np(242).eq.0) then
          setvar = palloc(242,'EULER',numnp*3 , 2)
          setvar = palloc(243,'LEULR',nen*3   , 2)
          eulerfl = .true.
        endif
        call genvec(3,3,hr(np(242)),' Euler',prt,prth,error,.false.)

c     [ceul]er  - coordinate angles for nodal transformations

      elseif(i.eq.76) then

        if(np(242).eq.0) then
          setvar = palloc(242,'EULER',numnp*3 , 2)
          setvar = palloc(243,'LEULR',nen*3   , 2)
          eulerfl = .true.
          fext  = 'qr0'
          if(neule.le.9) then
            write(fext(3:3),'(i1)') neule
          elseif(neule.le.99) then
            write(fext(2:3),'(i2)') neule
          endif
          neule = neule + 1
        call plinka(fext,c2,'   ')
        endif

c     [rfor]ce -- Specify Radial follower nodal FORces

      elseif(i.eq.77) then
        if(.not.nffl) then
          setvar = palloc(256,'NFORC',ndf*numnp,2)
          nffl = .true.
        endif
        call genvec(ndf,ndf,hr(np(256)),' Nodal Follower Forces',
     &              prt,prth,error,.false.)

c     [lfor]ce -- Specify Line nodal FORces

      elseif(i.eq.78) then

        lfrfl = .true.
        fext  = 'pf0'
        if(nforl.le.9) then
          write(fext(3:3),'(i1)') nforl
        elseif(nforl.le.99) then
          write(fext(2:3),'(i2)') nforl
        endif
        nforl = nforl + 1
        call plinka(fext,c2,'   ')

c     [*nod] - Number to add to all input and generated nodes
c     [*ele] - Number to add to all input and generated elements

      elseif(i.eq.79 .or. i.eq.80) then
        j = index(record,'=')
        if(j.eq.0) then
          j = index(record,',')
        endif
        xxx(1:75) = ' '
        ii = 0
        do jj = j+1,75
          if(record(jj:jj).ne.' ') then
            ii = ii + 1
            xxx(ii:ii) = record(jj:jj)
          endif
        end do ! jj

c       Set flag to parse expression

        flgco = coflg
        coflg = .true.
        call setval(xxx,75,td(1))
        coflg = flgco

c       Set *node

        if(i.eq.79) then
          starnd = nint(td(1))
          write(iow,2007) starnd
          if(ior.lt.0) then
            write(*,2007) starnd
          endif

c       Set *element

        else
          starel = nint(td(1))
          write(iow,2008) starel
          if(ior.lt.0) then
            write(*,2008) starel
          endif
        endif

c     [mesn] -> user defined mesh inputs

      elseif(i.gt.list-numesh) then
        n   = i + numesh - list
        uct = wd(i)
        call umshlib(n,prt)

      endif

      go to 100

c     Formats

2001  format(/' -> Region Number:',i4,' Maximum:',i4)
2002  format(' Coordinate transformation origin set to:'/
     &       '   x0 =',1p,e12.4:,'  y0 =',1p,e12.4:,'  z0 =',1p,e12.4)
2003  format(' -> Rigid body:',i4,' with modal deformations')
2004  format(' -> Rigid body:',i4,' initiated')
2005  format('    Location for Rigid Body Unknowns'/
     &   10x,'x0 =',1p,1e12.4:,'  y0 =',1p,1e12.4:,'  y0 =',1p,1e12.4)
2006  format(' -> Flexible elements')
2007  format(' -> Number to be added to all nodal   inputs =',i8)
2008  format(' -> Number to be added to all element inputs =',i8)

3001  format(' *ERROR* PMESH: Cannot regenerate SIDEs')
3002  format(' *ERROR* PMESH: Cannot regenerate SNODes')
3003  format(' *ERROR* PMESH: Cannot regenerate BLENds')
3004  format(' *WARNING* Initial node/element numbers necessary to'
     &      ,' use BLOCk in solution mode.')
3005  format(' *ERROR* PMESH: File:',a,' does not exist.')
3006  format(' *ERROR* PMESH: No SAVE,END statement for input data.')
3007  format(' *ERROR* PMESH: No SNODes were input for blending.')

4000  format('  INPUT FILE NAME: ',a//
     &       '  MESH INPUT RECORDS'/'  ------------------')
4001  format(5x,a,' ',a)

      end
