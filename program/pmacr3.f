c$Id: pmacr3.f,v 1.2 2006/12/18 17:00:16 rlt Exp $
      subroutine pmacr3(lct,ct,j)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'force' prints                               17/12/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Command language instruction subprogram mostly for
c               FEM arrays

c      Inputs:
c         lct(*)     - Command option
c         ct(3,*)    - Command parameters
c         j          - Command number in this routine

c      Outputs:
c         Depends on command number j
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'allotd.h'
      include  'arclel.h'
      include  'arclei.h'
      include  'arcler.h'
      include  'augdat.h'
      include  'auto2.h'
      include  'cdata.h'
      include  'cdat1.h'
      include  'chdata.h'
      include  'comfil.h'
      include  'compac.h'
      include  'compas.h'
      include  'complx.h'
      include  'corset.h'
      include  'corfil.h'
      include  'cornum.h'
      include  'counts.h'
      include  'ddata.h'
      include  'debugs.h'
      include  'edgdat.h'
      include  'elacts.h'
      include  'endata.h'
      include  'eqsym.h'
      include  'evdata.h'
      include  'fdata.h'
      include  'hdatam.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'ldata.h'
      include  'modcon.h'
      include  'modreg.h'
      include  'mxsiz.h'
      include  'ndata.h'
      include  'pdata0.h'
      include  'part0.h'
      include  'part1.h'
      include  'part3.h'
      include  'pfeapb.h'
      include  'plist.h'
      include  'pointer.h'
      include  'prflag.h'
      include  'print.h'
      include  'prlod.h'
      include  'pscal.h'
      include  'ptest.h'
      include  'pview.h'
      include  'rdata.h'
      include  'rdat0.h'
      include  'rdat1.h'
      include  'region.h'
      include  'rjoint.h'
      include  'sdata.h'
      include  'setups.h'
      include  'tdata.h'
      include  'xtout.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   pcomp,pfro,sfl,accrcy,tfl,setfl,setvar,palloc,pinput
      logical   prntsv
      character fint*132,y*1, lct(*)*15
      integer   i,imas, j, kk, k1,k2,k3,k4, lflag
      integer   ml1,leng,prec
      integer   mad, nbfgs,nnn, nparto, ndfeig(20)
      integer   ndf1,ndf2
      real*8    stol,etol, ee,phi2,vphi, chec, tau
      real*8    dotid, dot, ct(3,*),td(20)

      save

c     Transfer to correct process

      go to (1,2,3,4,5,6,7,8,9,1,1,12,13,14,15,16,17,18,19,19,
     &       21,22,23,3,25,26,27,28,1,30,31,32,33), j

c     Print displacements

c     [disp,all]             - print all displ.
c     [disp,,k1,k2,k3]       - print displ. k1 to k2 step k3
c     [disp,gnod,k1,k2,k3]   - print displ. for global nodes k1 to k2 step k3
c     [disp,coor,k1,xt,xtol] - print displ. for all nodes where x-k1=xt
c                                                                +-xtol
c     [disp,node,x1,x2,x3]   - print displ. for node closest to x1,x2,x3
c     [disp,imag,k1,k2,k3]   - print all imaginary k1 to k2 step k3
c     [disp,cmpl,k1,k2,k3]   - print all real/imag k1 to k2 step k3
c     [disp,list,k1]         - print all displacements in list k1

c     [velo ... same options as above
c     [acce ... same options as above
c     [comp ... same options as above

1     k1   = ct(1,l)
      pfro = pfr
      ndf1 = 1
      ndf2 = ndf

      if(j.eq.1) then
        ml1 = 1
        fp(1) = np(40)
        setfl = .true.
      elseif(j.eq.10) then
        ml1 = 2
        fp(1) = np(42)
        setfl = fl(9)
      elseif(j.eq.11) then
        ml1 = 3
        fp(1) = np(42) + nneq
        setfl = fl(9)
      elseif(j.eq.29) then
        ml1  = 1
        fp(1)  = np(40)
        ndf1 = nocomp
        ndf2 = nocomp
        setfl = .true.
      endif
      fp(5) = np(43)

c     Set type of coordinate outputs

      if(iview.eq.0) then ! Cartesian
        fp(3) = fp(1)
        fp(4) = fp(5)
      elseif(iview.gt.0) then ! Cylindrical or Spherical
        setvar = palloc(111,'TEMP1',nneq,  2)
        fp(3) = np(111)
        if(cview) then
          setvar = palloc(112,'TEMP2',ndm*numnp,  2)
          fp(4) = np(112)
          do i = 0,ndm*numnp-1
            hr(fp(4)+i) = hr(fp(5)+i) ! Move nodal coordinates to temp location
          end do ! i
        endif
        if(iview.eq.1) then ! Cylindrical
          call prcylc(hr(fp(4)),hr(fp(1)),hr(fp(3)))
        else                ! Spherical
c         call prcylc(hr(fp(4)),hr(fp(1)),hr(fp(3)))
          write(*,*) ' *ERROR* Not implemented'
        endif
      endif

c     List outputs

      if(pcomp(lct(l),'list',4)) then
        k1 = max(1,min(3,k1))
        if(setfl) then
          call prtlis(hr(fp(4)),hr(fp(3)),ttim,prop,ndm,ndf,
     &                niols(k1),iolist(1,k1),ml1,ndf1,ndf2,prth)
        else
          write(ilg,3001)
          write(iow,3001)
          if(ior.lt.0) then
            write(*,3001)
          endif
        endif

c     Set output limits

      else

        k2 = ct(2,l)
        if(k2.eq.0) k2 = k1
        k3 = ct(3,l)
        if(k3.eq.0) k3 = 1
        nxt = 0

c       Set for imaginary or real part of complex quantities

        if(pcomp(lct(l),'imag',4) .or. pcomp(lct(l),'cmpl',4)) then
          if(cplxfl) then
            k4 = nneq*3
          else
            write(iow,3005)
            if(ior.lt.0) write(*,3005)
            if(iview.gt.0) then
              setvar = palloc(111,'TEMP1',0,  2)
              if(np(112).ne.0) then
                setvar = palloc(112,'TEMP2',0,  2)
              endif
            endif
            return
          endif
        else
          k4 = 0
        endif

c       Set flag for local node prints

        pfeap_gnod = .false.

c       Set for all

        if(pcomp(lct(l),'all ',4)) then
          k1 = 1
          k2 = numnp
          k3 = 1

c       Set for specified coordinate

        elseif(pcomp(lct(l),'coor',4)) then
          nxt = max(1,min(k1,ndm))
          xt  = ct(2,l)
          if(ct(3,l).eq.0.0d0) then
            xtol = 1.0d-6*dxmsh(nxt)
          else
            xtol = abs(ct(3,l))
          endif
          k1 = 1
          k2 = numnp
          k3 = 1

c       Set for specified node

        elseif(pcomp(lct(l),'node',4)) then

          call pclnod(ct(1,l),hr(fp(4)),mr(np(190)),ndm,numnp, k1)
          k2 = k1
          k3 = 1

c       Set for global node prints

        elseif(pcomp(lct(l),'gnod',4) .and. processor.gt.1) then

          pfeap_gnod = .true.

c       Do range

        else
          k1 = max(1,min(k1,numnp))
          k2 = max(1,min(numnp,k2))
          if(k2-k1.ne.0) k3 = sign(k3,k2-k1)
          pfr = .true.
        endif

c       Perform displacement, velocity, acceleration outputs:

        if(setfl) then
          if( .not.pcomp(lct(l),'cmpl',4)) then
            call prtdis(hr(fp(4)),hr(fp(3)+k4),ttim,prop,ndm,ndf,
     &                  k1,k2,k3,ml1,ndf1,ndf2,prth)
          else
            call prtcmp(hr(fp(4)),hr(np(40)),hr(fp(3)+k4),ttim,
     &                  prop,ndm,ndf,k1,k2,k3,ndf1,ndf2,prth)
          endif
        else
          write(ilg,3001)
          write(iow,3001)
          if(ior.lt.0) then
            write(*,3001)
          endif
        endif

      endif
      pfr = pfro

      if(iview.gt.0) then
        setvar = palloc(111,'TEMP1',0,  2)
        if(np(112).ne.0) then
          setvar = palloc(112,'TEMP2',0,  2)
        endif
      endif
      return

c     [test],<on,off>,value - Test step on the solution

2     if(pcomp(lct(l),'off',3)) then
        testfl = .false.
        if(ior.gt.0) then
          write(iow,2032)
        else
          write(  *,2032)
        endif
      else
        testfl = .true.
        if(ct(1,l).eq.0.0d0) then
          testva = 0.001d0
        else
          testva = ct(1,l)
        endif
        if(ior.gt.0) then
          write(iow,2033) testva
        else
          write(  *,2033) testva
        endif
      endif
      return

c     Modify mesh data
c     [mesh]           -  Reenter mesh generation phase
c     [mesh,filename]  -  Read the data from 'filename'

c     Reset boundary codes to represent correct b.c. - can change

3     nparto = npart

c     Set surface flags

      if(j.eq.3) then

c       Set edge input flags all false

        eanfl = .false.
        ebcfl = .false.
        edifl = .false.
        efcfl = .false.
        eprfl = .false.

        neang = 0
        nebcs = 0
        nedis = 0
        nefrc = 0
        nepro = 0

c       Set coordinate input flags all false

        surfl = .false.
        boufl = .false.
        disfl = .false.
        forfl = .false.
        angfl = .false.
        reafl = .false.
        intfl = .false.
        basfl = .false.

        nsurf = 0
        nbouf = 0
        ndisf = 0
        nforf = 0
        nangf = 0
        nintf = 0
        nbasf = 0
        neule = 0
        i     = -1
        if(pcomp(lct(l),'new',1)) then
          i = -2

c       Regenerate mesh data

        elseif(.not.pcomp(lct(l),' ',1)) then

c         Set inputs from file specified in second field

          fint = lct(l)
          if(ior.lt.0) write(*,2014) fint
          call pincld(fint)
        endif

        chflg = .true.
        call pmesh(i,prt,prth)
        chflg = .false.

        if(.not.pcomp(lct(l),' ',1)) then
          call pincld('end')
        endif

c       Do the edge/coord boundary conditions

        call pecmes()

c     [opti]mize node numbering then reset profile
c     [opti,off]  - Turn off profile optimization
c     [opti,cont] - Do optimization during contact solution
c     [opti,hoit] - Set to use old optimizer by Wilson/Hoit
c     [opti,sloa]n- Set to use Scott Sloan algorithm
c                   N.B. Default is Sloan algorithm

      elseif(j.eq.24) then

        i      = -1
        nparto =  npart

        if(pcomp(lct(l),'off',3)) then

c         Set renumber array to input values

          optflg  = .false.
          optmsh  = .false.
          opthoit = .false.
          do k1 = 0,numnp-1
            mr(np(89)+k1) = k1 + 1
          end do ! k1

        elseif(pcomp(lct(l),'hoit',4)) then

          opthoit = .true.
          if(ior.lt.0) write(*,2034)
          return

        elseif(pcomp(lct(l),'sloa',4)) then

          opthoit = .false.
          if(ior.lt.0) write(*,2035)
          return

        elseif(pcomp(lct(l),'cont',4)) then

c       Set renumber array to "optimized" values

          if(ior.lt.0) write(*,2036)
          optflg = .true.
          optmsh = .true.
          if(ct(1,l).eq.0.0d0) return

        else

          optflg = .false.
          optmsh = .true.
          call contact (203)
          call optid()
          if(debug) call iprint(mr(np(89)),1,numnp,1,'Node Renumbers')

        endif

      endif

c     Set ID to current boundary condition array


c     Reset ID equation numbers

      setvar = palloc(111,'TEMP1',nneq,  1)
      call psetid(mr(np(31)),mr(np(33)),mr(np(111)),mr(np(32)),
     &            mr(np(240)),nie,ndf,nen,nen1,numel,numnp,nummat)
      setvar = palloc(111,'TEMP1',   0,  1)
c     do k2 = 0,nneq-1
c       mr(np(31)+k2) = mr(np(31)+k2+nneq)
c     end do ! k2

c     Set region indicator so all active region elements are assembled

301   nreg = -1

c     Compute new equation numbers and profile

      nparto =  npart
      do k1 = 1,4
        if(.not.tflp(k1)) then
          call partpt(k1,tflp(k1),.false.)

c         Set new current profile

          call profil(mr(np(20+k1)),mr(np(34)),mr(np(31)),
     &                mr(np(33)),1,pfr)
          call profil(mr(np(20+k1)),mr(np(34)),mr(np(31)),
     &                mr(np(33)),2,pfr)
          nqp(k1)    =  neq
          nqr(k1)    =  neqr
          mxprop(k1) = (mr(np(20+k1)+neq-1))*ipc
          mxneqp(k1) =  neq*ipc

c         Compute new fill for sparse storage

          if(ittyp.gt.-3) then
            if(ittyp.eq.-1 .and. np(92).ne.0) then
              setvar = palloc(92,'OINB',0,1)
            endif
            k2 = 0
            call iters(k2,1)
          endif
        endif
      end do ! k1
      npart = nparto
      call partpt(npart,tflp(npart),.true.)

      return

c     Plot outputs
c     [plot] - enter plot mode for interactive
c     [plot,optn] - see plot manual for optn's

4     rfl = .false.
      call pplotf(lct(l),ct(1,l),prop)
      return

c     Subspace eigencomputations (for: mass,iden,geom)
c     [subs,,k1,k2]     - subspace for k1 eigenpairs
c     [subs,prin,k1,k2] - subspace for k1 eigenpairs - print matrices
c                       - k2 used to overwrite default no. guard vects.

5     if(fl(4)) then
        write(ilg,3002) 'SUBS'
        write(iow,3002) 'SUBS'
        if(ior.lt.0) write(*,3002) 'SUBS'
        return
      endif
      if(processor.gt.1) then
        write(ilg,3011) 'SUBS'
        write(iow,3011) 'SUBS'
        if(ior.lt.0) write(*,3011) 'SUBS'
        return
      end if
      if(fl(5).and.fl(6)) then
        write(ilg,3003)
        write(iow,3003)
        if(ior.lt.0) then
          write(*,3003)
          return
        endif
        call plstop()
      endif
      if(fl(1)) then
        fp(1) = np(npart+8)
        imas  = 1
      else
        fp(1) = np(npart+12)
        imas  = 2
      endif
      mf   = ct(1,l)
      mad  = ct(2,l)
      stol = ct(3,l)
      vneq = neq

c     Mask tolerance if it is too small

      if(stol.eq.0.d0) then
        stol = max(tol, 1.d-12)
      endif

      if(mad.le.0) mad = max(mf,8)
      mf = min(neq,max(1,mf))
      mq = min(mf+mad,neq)
      call numass(hr(fp(1)),neq,mq)

c     Number of possible terms positive

      if(mq.gt.0) then
        prntsv = prnt
        if(ittyp.eq.-2) then
          prnt = .false.
        endif
        if(mq.lt.mf) then
          write(iow,2001) mq
          if(ior.lt.0) then
            write(*,2001) mq
          endif
        endif
        mf = min(mf,mq)
        setvar = palloc( 76,'EVAL',mq,    2)
        setvar = palloc( 77,'EVEC',mq*neq,2)

        setvar = palloc(111,'TEMP1',neq        ,2)
        setvar = palloc(112,'TEMP2',mq*(mq+1)/2,2)
        setvar = palloc(113,'TEMP3',mq*(mq+1)/2,2)
        setvar = palloc(114,'TEMP4',mq         ,2)
        setvar = palloc(115,'TEMP5',mq         ,2)
        setvar = palloc(116,'TEMP6',mq*mq      ,2)

        sfl    = pcomp(lct(l),'prin',4)
        setvar = pfr
        pfr    = sfl
        call subsp(hr(np(npart)),hr(fp(1))  ,hr(np(77)),hr(np(111)),
     &             hr(np(112))  ,hr(np(113)),hr(np(76)),hr(np(114)),
     &             hr(np(115))  ,hr(np(116)),mf,mq,neq,imas,shift,
     &             stol,sfl,25)
        pfr    = setvar

        setvar = palloc(116,'TEMP6', 0,2)
        setvar = palloc(115,'TEMP5', 0,2)
        setvar = palloc(114,'TEMP4', 0,2)
        setvar = palloc(113,'TEMP3', 0,2)
        setvar = palloc(112,'TEMP2', 0,2)
        setvar = palloc(111,'TEMP1', 0,2)

c       Set IMF array for flexible dof's and eliminate rigid body modes

        if(nmbody.gt.0) then

          setvar = palloc(176,'IMODF', nneq,1)
          call pmovei(mr(np(31)),mr(np(176)),nneq)

          write(*,*) 'No. FLEX EQ, MODES',neq,mf
          call elimrbm(hr(np(76)),hr(np(77)),vneq,mf)

          neqmf = neq

          write(*,*) 'No. FLEX EQ, MODES',neqmf,mf
          setvar = palloc(177,'AFD  ', mf  ,2)
          setvar = palloc(178,'AFL  ', mf*3,2)
          setvar = palloc(179,'AFU  ', mf*3,2)
          setvar = palloc(180,'BFORC', mf*8,2)

        endif

c       Restore print flag

        prnt = prntsv

c     There are no possible modes

      else

        if(imtyp.eq.1) then
          write(ilg,3009)
          write(iow,3009)
          if(ior.lt.0) then
            write(*,3009)
          endif
        elseif(imtyp.eq.2) then
          write(ilg,3010)
          write(iow,3010)
          if(ior.lt.0) then
            write(*,3010)
          endif
        endif
        if(ior.gt.0) call plstop()

      endif

      return

c     Write a file
c     [writ,file]  - open write file named 'file'
c     [writ,disp]  - write displacements to 'file'
c     [writ,stre]  - write nodal streses to 'file'
c     [writ,eige]  - write eigenpairs to 'file'
c     [writ,wind]  - rewind 'file'
c     [writ,clos]  - close 'file'

6     call writer(lct(l),hr(np(40)),nneq)
      return

c     Read a file
c     [read,file]  - open read file named 'file'
c     [read,disp]  - read displacements from 'file'
c     [read,stre]  - read nodal streses from 'file'
c     [read,eige]  - read eigenpairs from 'file'
c     [read,wind]  - rewind 'file'
c     [read,clos]  - close 'file'

7     call reader(lct(l),hr(np(40)),nneq)
      if(pcomp(lct(l),'disp',4)) fl(11) = .false.
      return

c     Set contact flag
c     [cont,chec]k     : Perform Geometry check now, but not at iters.
c     [cont,noch]eck   : Perform Geometry check at each iteration
c     [cont,on]        : Reenable contact
c     [cont,off]       : Disable  contact
c     [cont,pena,n,pen]: n = pair number ; pen = penalty value
c     [cont,fric]tion  : Contact friction ON
c     [cont,nofr]iction: Contact friction OFF

8     if(pcomp(lct(l),'chec',4) .or. pcomp(lct(l),'noch',4)) then
        call contact (304)
        if(ior.lt.0 .and. pcomp(lct(l),'chec',4)) write(*,2021)
        if(ior.lt.0 .and. pcomp(lct(l),'noch',4)) write(*,2022)
      else if(pcomp(lct(l),'off ',4) .or.
     &        pcomp(lct(l),'on  ',4)) then
        call contact (309)
        if(ior.lt.0) write(*,2023) lct(l)
      else if(pcomp(lct(l),'pena',4)
     &   .or. pcomp(lct(l),'fric',4) .or. pcomp(lct(l),'nofr',4)) then
        call contact (310)
        if(ior.lt.0 .and. pcomp(lct(l),'fric',4)) write(*,2024) ct(1,l)
        if(ior.lt.0 .and. pcomp(lct(l),'nofr',4)) write(*,2025)
      endif
      return

c     Restart previously run problem
c     [rest,ext_name]

9     fint = fres
      if(.not.pcomp(lct(l),'    ',4)) then
        call addext(fint,lct(l),128,4)
      endif
      if(ior.lt.0) write(*,2012) fint
      write(iow,2012) fint
      call restrt(fint,ndm,ndf,nneq,1)
      return

c     BFGS algorithm
c     [bfgs,,k1,k2,k3] - BFGS soln k1 = no. steps; k2 = line search tol;
c                                  k3 = bfgs energy tol

c     Allocate memory for arrays
c     mbo : first memory address for oldrsd in bfgs algo. (iterat)
c     mbd : first memory address for d      in bfgs algo. (iterat)
c     mbv : first memory address for v      in bfgs algo. (iterat)
c     mbw : first memory address for w      in bfgs algo. (iterat)
c     mbt : first memory address for t      in bfgs algo. (iterat)
c     mst : first memory address for v      in bfgs algo. (store )

12    if (fl(12)) then
        setvar = palloc( 69,'BFGD',neq,   2)
        setvar = palloc( 70,'BFGO',neq,   2)
        setvar = palloc( 73,'BFGV',neq,   2)
        setvar = palloc( 74,'BFGW',neq,   2)

c     Reserve storage for temporary vectors
c            u du dyn
c            | |   |
c            V V   V
        i  = 3+1 +nrt

        setvar = palloc( 72,'BFGT',nneq*i, 2)

c     15 max. vectors for v, w in BfGS (store)

c       setvar = palloc( 71,'BFGS',neq*30,2)
        fl(12) = .false.
      endif

c     Call BFGS routine
c        Max. of of iterations (nbfgs) = 15
c        Default stol for line search  = .8
c        Default etol for bfgs energy  = tol

      nbfgs = ct(1,l)
      stol  = ct(2,l)
      etol  = ct(3,l)
      if (nbfgs.eq.0)    nbfgs = 15
      if (stol.eq.0.0d0) stol  = 0.8d0
      if (etol.le.0.0d0) etol  = tol
c     nbfgs = min(nbfgs,15)
      stol  = min(stol,0.9d0)

c     Assign storage for BFGS V and W vectors

      setvar = palloc( 71,'BFGS',neq*nbfgs*2,2)

      call iterat(np(40),np(26),hr(np(70)),hr(np(69)),hr(np(72)),
     &            accrcy,hr(np(73)),hr(np(74)),prt,mr(np(31)),
     &            nbfgs,stol,etol)
      return

c     Arc-length method
c     [arcl,,kflag,lflag]   - set arc length parameters
c     [arcl,add,k1,tau]     - add eigvenvector k1, amount = tau
c     [arcl,chec,k1]        - check for bifurcation using eigv. k1
c     [arcl,off]            - set arclength to off

13    if(pcomp(lct(l),'off' ,3)) then
        if(ior.lt.0) write(*,2010)
        write(iow,2010)
        arcf  = .false.
        rlnew = 1.d0
        return
      else
        if(ior.lt.0) write(*,2011)
        write(iow,2011)
      endif
      if(pcomp(lct(l),'add' ,3)) go to 133
      if(pcomp(lct(l),'chec',4)) go to 134
      if(kflag.eq.0) kflag = ct(1,l)
      if(kflag.eq.0) kflag = 2
      lflag  = ct(2,l)
      if(lflag.eq.0) then
        k1 = max(neq,nneq)
        setvar = palloc( 84,'MU1  ',k1,2)
        if(kflag .eq. 2 .or. kflag .eq. 3 .or. kflag .eq. 5) then
          setvar = palloc( 85,'MU2   ',k1,2)
        endif
        arcf = .true.
      endif
      call dicont(mr(np(31)),numnp,ndf,lflag)
      return

c     Add scaled eigenvector to displacement vector

133   tau = ct(2,l)
      k1  = ct(1,l)
      k1 = max(min(mf,k1),1)
      if(np(77).eq.0.and.ior.lt.0) then
        write(*,3004)
        return
      endif
      kk = (k1 - 1) * vneq
      if(tau.eq.0.0d0) then
        vphi = dotid(hr(np(40)),hr(np(77)+kk),mr(np(31)),nneq)
        phi2 = dot(hr(np(77)+kk),hr(np(77)+kk),neq)
        ee   = dot(hr(np(40)),hr(np(40)),nneq)
        tau  = 100.d0 * vphi / sqrt(ee*phi2) + 1.d0
        if(ior.lt.0) write(*,2003) tau,k1
        write(iow,2003) tau,k1
      endif
      call paddv(hr(np(40)),hr(np(77)+kk),nneq,tau,mr(np(31)))
      return

c     Check for bifurcation or limit point

134   k1 = ct(1,l)
      if(k1.eq.0) k1 = 1
      if(np(77).eq.0.and.ior.lt.0) then
        write(*,3004)
        return
      endif
      kk = (k1 - 1)*vneq
      chec = dotid(hr(np(27)),hr(np(77)+kk),mr(np(31)),nneq)
      if(ior.lt.0) write(*,2004) k1,chec
      write(iow,2004) k1,chec
      return

c     Save restart information for intermediate points
c     [save,ext_name]

14    fint = fres
      if(.not.pcomp(lct(l),'    ',4)) then
        call addext(fint,lct(l),128,4)
      endif
      if(ior.lt.0) write(*,2013) fint
      write(iow,2013) fint
      call restrt(fint,ndm,ndf,nneq,2)
      return

c     [paus] Pause on no convergence

15    if(pcomp(lct(l),'inte',4)) then
        write(*,2017)
        read(*,1000) y
      elseif(abs(aengy).ge.100.d0*rnmax) then
        if(ior.lt.0) then
          call pprint('   Solution is diverging. Continue (y or n)? >')
          read(*,1000) y
          if(y.eq.'y' .or. y.eq.'Y') return
          l = lve(1) - 1
        else
          l = lve(1) - 1
        endif
      endif
      return

c     [eige,<vect,k1,k2>] Compute element eigenpairs for element 'k1'.
c      k1 <= 0 for last element; k2 <= 0 for stiffness; k2 > 0 for mass.

16    k1 = nint(ct(1,l))
      if(k1.le.0) then
        k1 = numel
      else
        k1 = min(k1,numel)
      endif
      k2 = nint(ct(2,l))
      if(k2.gt.0) then
        k3 = 5
        if(ior.lt.0) then
          write(*,2026) 'MASS',k1
        endif
        write(iow,2026) 'MASS',k1
      else
        k3 = 3
        if(ior.lt.0) then
          write(*,2026) 'STIFFNESS',k1
        endif
        write(iow,2026) 'STIFFNESS',k1
      endif
      call formfe(np(40),np(26),np(26),np(26),.false.,.false.,.false.,
     &            .false.,k3,k1,k1,1)
      tfl = pcomp(lct(l),'vect',4)
      setvar = palloc( 75,'EIGE ',  nst*nst+  nst, 2)
      setvar = palloc(111,'TEMP1',2*nst*nst+5*nst, 2)
      call peigc(hr(np(36)),mr(np(34)),nst)
      call peige(hr(np(36)),nst,hr(np(111)),tfl)
      setvar = palloc(111,'TEMP1', 0, 2)
      return

c     [expl]icit solutions

c     Initialize for lumped mass

17    if(fl(2).and.fl(9)) then
        if(rnmax.eq.0.0d0) then
          rnmax = rnorm
        endif
        call piacel(hr(nl),hr(np(26)),hr(np(26)),neq)
        call update(mr(np(31)),hr(np(30)),hr(np(40)),hr(np(42)),
     &              hr(np(26)),fl(9),2)
        fl(8) = .false.
        if(abs(rnorm).le.sqrt(tol)*rnmax .and. lv.gt.1) then
          ct(1,lve(lv)) = ct(1,lvs(lv))
          l             = lve(lv) - 1
          floop(1)      = .false.
          autcnv        = .true.
        endif

c     Else write warning

      else
        write(iow,3000)
        if(ior.lt.0) write(*,3000)
      endif
      return

c     [memo]ry usage

18    if(ior.lt.0) write(*,2006)
      return

c     Activation/deactivation of regions

c     [acti,all]       activate all regions
c     [acti,,k1,k2,k3] activate regions: do n = k1,k2,k3
c     [acti,init,k1,k2,k3] (also initialize strains)

c     [deac,all]       deactivate all regions
c     [deac] show current parameters on screen
c     [deac,,k1,k2,k3] deactivate regions: do n = k1,k2,k3

19    if(j.eq.19) k4 =  1 ! activate
      if(j.eq.20) k4 = -1 ! deactivate
      pltact = .true.

c     Set all regions to be active/deactive

      if(pcomp(lct(l),'all',3)) then
        k1 = 1
        k2 = mxreg
        k3 = 1

c     Set specified region only

      else
        k1 = nint(ct(1,l))

c       Output status of active/deactive retions

        if(k1.le.0) then
          do i = 1,mxreg
            do nnn = nen1-2,numel*nen1,nen1
              if(mr(np(33)+nnn).eq.i) then
                write(iow,2029) i
                if(ior.lt.0) write(*,2029) i
                go to 191
              elseif(-mr(np(33)+nnn).eq.i) then
                write(iow,2030) i
                if(ior.lt.0) write(*,2030) i
                go to 191
              end if
            end do ! nnn
            write(iow,2031) i
            if(ior.lt.0) write(*,2031) i
191         continue
          end do ! i
          return

c       Error: Specified region too large

        elseif(k1.gt.mxreg) then
          write(ilg,3008) k1
          write(iow,3008) k1
          if(ior.lt.0) write(*,3008) k1
          return

c       Set region range to active/deactivate

        else
          k2 = max(k1,min(nint(ct(2,l)),mxreg))
          k3 = max(1,nint(ct(3,l)))
        endif
      end if

c     Set 'first' activation indicators in 'ix' array to zero

      do nnn = nen1-3,numel*nen1,nen1
        mr(np(33)+nnn) = 0
      end do ! nnn

      do i = k1,k2,k3
        if(k4.eq.1) then
          if(ior.lt.0) write(*,2007) i
          write(iow,2007) i
          numact = i
          itract = 0
        else
          if(ior.lt.0) write(*,2008) i
          write(iow,2008) i
          numdea = i
          itrdea = 0
        endif

c       Set region indicators on IX array

        do nnn = nen1-2,numel*nen1,nen1
          if(abs(mr(np(33)+nnn)).eq.i) then
            if(k4.eq.1 .and. mr(np(33)+nnn).lt.0) then
              mr(np(33)+nnn-1) =  i                        ! new activation
            endif
            mr(np(33)+nnn) = k4*abs(mr(np(33)+nnn))        ! active/deactive
          end if
        end do ! nnn

c       Initialize strains in activated elements

        if(pcomp(lct(l),'init',4)) then
          if(ior.lt.0) write(*,2016) i
          write(iow,2016) i
          nreg   = i
          hflgu  = .false.
          h3flgu = .true.
          if(k4.eq.1) then
            call formfe(np(40),np(26),np(26),np(26),
     &                 .false.,.false.,.false.,.false.,17,1,numel,1)
          else
            call formfe(np(40),np(26),np(26),np(26),
     &                 .false.,.false.,.false.,.false.,18,1,numel,1)
          endif
        endif

      end do ! i

      nreg = -1

c     Reset ID equation numbers

      setvar = palloc(111,'TEMP1',nneq,  1)
      call psetid(mr(np(31)),mr(np(33)),mr(np(111)),mr(np(32)),
     &            mr(np(240)),nie,ndf,nen,nen1,numel,numnp,nummat)
      setvar = palloc(111,'TEMP1',   0,  1)

      nparto = npart

c     Set projection flag false (forces recompute of projections)

      fl(11) = .false.

      i     = -1
      go to 301

c     [zero]         : Zero solution for entire problem to start over
c     [zero],rate    : Zero transient vectors (VEL)
c     [zero],regi,n1 : Zero displacements in region 'n1' only

21    if(pcomp(lct(l),'regi',4)) then

c       Zero displacements for a region, but not for adjacent nodes

        k1 = abs(ct(1,l))
        if(ior.lt.0) write(*,*) 'Zero region =',k1

        setvar = palloc(111,'TEMP1',numnp,  1)

c       Set region values on the IX array

        fp(2) = np(111) - 1
        do nnn = 0,numel*nen1-1,nen1
          fp(1) = np(33) + nnn - 1
          if(mr(fp(1)+nen1-1).ge.0 .and. mr(fp(1)+nen1-1).ne.k1) then
            do kk = 1,nen
              k2 = mr(fp(1)+kk)
              if(k2.gt.0) mr(fp(2)+k2) = mr(fp(2)+k2) + 1
            end do ! kk
          end if
        end do ! nnn

        do nnn = 0,numel*nen1-1,nen1
          fp(1) = np(33) + nnn - 1

          if(mr(fp(1)+nen1-1).eq.k1) then

c           Zero unconnected node displacements

            do kk = 1,nen
              k2 = mr(fp(1)+kk)
              if(k2.gt.0 .and. mr(fp(2)+k2).eq.0) then
                k2 = ndf*(k2-1) + np(40) - 1
                k3 = k2 + nneq
                k4 = k3 + nneq
                do i = 1,ndf
                  hr(i+k2) = 0.0d0
                  hr(i+k3) = 0.0d0
                  hr(i+k4) = 0.0d0
                end do ! i
              end if
            end do ! kk
          end if
        end do ! nnn

c       Zero incremental displacements

        call pzero(hr(np(40)+nneq),nneq)

c       Copy history t_n+1 variables to t_n location

        call reshis(mr(np(33)+nen),nen1, numel, 2, 1)

        setvar = palloc(111,'TEMP1',0,  1)

c     Zero rate vector only

      elseif(pcomp(lct(l),'rate',4)) then

        if(np(42).ne.0) then
          call pgetd('VEL  ',fp(1), leng, prec, tfl)
          if(tfl) call pzero(hr(fp(1)),nrt*nneq)
        endif

      else
        if(pcomp(lct(l),' ',1)) then
          if(ior.lt.0) then
            call pprint(' *WARNING* Reset solution to zero? (y or n)>')
            read (*,1000) y
            if(y.ne.'y' .and. y.ne.'Y') return
          endif
          ttim   =  0.0d0
          timold = -1.0d0
          rlnew  =  0.0d0
        endif

c       Set displacements: U_fem, URATE_fem, U_rb, U_lam, U_jts

        if(pcomp(lct(l),' ',1) .or. pcomp(lct(l),'node',4)) then

          call pzero( hr(np(40)), nneq*3 )

          if(np(42).ne.0) then
            call pgetd('VEL  ',fp(1), leng, prec, tfl)
            if(tfl) call pzero(hr(fp(1)),nrt*nneq)
          endif

          if(np(95).ne.0) then
            call pgetd('RCG  ',fp(1), leng, prec, tfl)
            call pgetd('RLAMB',fp(2),   k2,   k3, tfl)
            if(tfl) call rcglam(hr(fp(1)),hr(fp(2)),k2/54)
          endif

          if(np(103).ne.0) then
            call pgetd('RJTU ',fp(1), leng, prec, tfl)
            if(tfl) call pzero(hr(fp(1)),leng)
          endif
        elseif(pcomp(lct(l),'dofs',4)) then
          k1 = abs(nint(ct(1,l)))
          k2 = abs(nint(ct(2,l)))
          k3 = nint(ct(3,l))
          if(k1+k2.eq.0) then
            k1 = 1
            k2 = ndf
            k3 = 1
          endif
          k1 = max(1,k1)
          k2 = min(ndf,max(k1,k2))
          k3 = max(1,k3)

c         Zero displacement values

          fp(2) = np(40) - 1
          do i = k1,k2,k3
            do kk = i,3*nneq,ndf
              hr(fp(2)+kk) = 0.0d0
            end do ! kk
          end do ! i

c         Zero time dependent values

          if(np(42).ne.0) then
            call pgetd('VEL  ',fp(1), leng, prec, tfl)
            if(tfl) then
              do i = k1-1,k2-1,k3
                do kk = i,leng-1,ndf
                  hr(fp(1)+kk) = 0.0d0
                end do ! kk
              end do ! i
            endif
          endif

        endif

c       Zero history

        if(pcomp(lct(l),' ',1) .or. pcomp(lct(l),'hist',4)) then
          if(np(49).ne.0) then
            call pgetd('H    ',fp(1), leng, prec, tfl)
            if(tfl) then
              call pzero(hr(fp(1)),leng)
              hflgu  = .true.
              h3flgu = .true.
              call formfe(np(40),np(26),np(26),np(26),
     &                   .false.,.false.,.false.,.false.,14,1,numel,1)
            endif
          endif
        endif
      endif

      return

c     [epri]nt - Output last element stiffness and residual

22    call mprint(hr(np(36)),nst,nst,nst,'Last Element S-Matrix')
      call mprint(hr(np(35)),  1,nst,  1,'Last Element P-Vector')
      if(cplxfl) then
        call mprint(hr(np(36)+nst*nst),nst,nst,nst,
     &             'Last Element S-Matrix (Imag)')
        call mprint(hr(np(35)+nst    ),  1,nst,  1,
     &             'Last Element P-Vector (Imag)')
      endif
      return

c     [moda]l - Integration of linear problems by modal methods

23    if(np(77).eq.0) then
        if(ior.lt.0) then
          write(*,3007) 'MODA'
        else
          write(ilg,3007) 'MODA'
          write(iow,3007) 'MODA'
          call plstop()
        endif
      else
        if(modfl) then
          if(pfeap_on) then
            setvar = palloc(187,'FSMOD', mf*3   ,2)
          else
            setvar = palloc(187,'FSMOD', mf*2   ,2)
          endif
          setvar = palloc(188,'YYMOD', mf*3   ,2)
          setvar = palloc( 42,'VEL  ',nneq*nrt,2)
          setvar = palloc(111,'TEMP1',nneq*2  ,2)
          fl(9) = .true.

          call pmodin(hr(np(77)),hr(np(188)),hr(np(nx)),hr(np(nx)+neq),
     &                mr(np(31)),hr(np(40)),hr(np(42)),hr(np(111)),
     &                hr(np(42)+2*nneq),mf,neq,nneq,fl(1))

          setvar = palloc(111,'TEMP1', 0, 2)
          modfl = .false.
        endif
        call ploa1(ttim,dt)
        call pload(mr(np(31)),hr(np(40)),hr(np(30)),hr(np(26)),
     &             prop*rlnew,.true.,.false.)
        call pmodal(dt,hr(np(76)),hr(np(77)),hr(np(187)),hr(np(188)),
     &              mr(np(31)),hr(np(30)),hr(np(40)),hr(np(42)),mf)
      endif
      return

c     [eigv],dofs/dof-list  : Set DOFS for eigen computation (1= active; 0=not)
c     [eigv],all,nnn        : Output eigenvector nnn (all)
c     [eigv],coor,k1,xt,nnn : Output eigenvector nnn for x_k1 = xt
c     [eigv],list,k1,nnn    : Output eigenvector nnn using list k1
c     [eigv],nnn,k1,k2,k3   : Output eigenvector nnn nodes k1 - k2 @ inc k3


25    if(pcomp(lct(l),'dofs',4)) then

        if(ior.lt.0) then
          write(*,4000)
          call pprint('   >')
        endif
        setvar = pinput(td,ndf)
        do i = 1,ndf
          ndfeig(i) = nint(td(i))
        end do ! i
        write(iow,2018) (i,ndfeig(i),i=1,ndf)
        if(ior.lt.0) then
          write(*,2018) (i,ndfeig(i),i=1,ndf)
        endif

      elseif(np(77).ne.0) then

        call pzero(hr(np(26)),nneq)
        tfl = pfr

c       List outputs

        if(pcomp(lct(l),'list',4)) then
          k1  = max(1,min( 3,nint(ct(1,l))))
          nnn = max(1,min(mq,nint(ct(2,l)))) - 1
          call pmovec(mr(np(31)),hr(np(77)+nnn*vneq),hr(np(26)),nneq)
          if(np(167).ne.0) then
            call ruplnk(hr(np(26)),hr(np(26)),hr(np(43)),mr(np(100)),
     &                  mr(np(167)),1,ndm,ndf,numnp,.false.)
          endif
          call prtlis(hr(np(43)),hr(np(26)),ttim,hr(np(76)+nnn),
     &                ndm,ndf,niols(k1),iolist(1,k1),4,1,ndf,prth)
          return


c       Set for all

        elseif(pcomp(lct(l),'all ',4)) then

          nxt = 0
          nnn = max(1,min(mq,nint(ct(1,l)))) - 1
          k1  = 1
          k2  = numnp
          k3  = 1

c       Set for specified coordinate

        elseif(pcomp(lct(l),'coor',4)) then
          nxt = max(1,min(nint(ct(1,l)),ndm))
          xt  = ct(2,l)
          xtol= 1.0d-6*dxmsh(nxt)
          k1  = 1
          k2  = numnp
          k3  = 1
          nnn = max(1,min(mq,nint(ct(3,l)))) - 1

c       Output vector nnn nodal values

        else
          call setval(lct(l),15,tau)
          nnn = max(1,min(mq,nint(tau))) - 1
          nxt = 0
          k1  = nint(ct(1,l))
          k2  = nint(ct(2,l))
          k3  = nint(ct(3,l))
          if(k2.eq.0) k2 = k1
          if(k3.eq.0) k3 = 1
          k1  = max(1,min(k1,numnp))
          k2  = max(1,min(numnp,k2))
          if(k2-k1.ne.0) k3 = sign(k3,k2-k1)
          pfr = .true.
        endif

c       Print eigenvectors

        call pmovec(mr(np(31)),hr(np(77)+nnn*vneq),hr(np(26)),nneq)
        if(np(167).ne.0) then
          call ruplnk(hr(np(26)),hr(np(26)),hr(np(43)),mr(np(100)),
     &                mr(np(167)),1,ndm,ndf,numnp,.false.)
        endif
        call prtdis(hr(np(43)),hr(np(26)),ttim,hr(np(76)+nnn),ndm,ndf,
     &              k1,k2,k3,5,1,ndf,prth)
        pfr = tfl

c     Warn must compute eigen problem first

      else
        write(ilg,3007) 'EIGV'
        write(iow,3007) 'EIGV'
        if(ior.lt.0) then
          write(*,3007) 'EIGV'
        endif
      endif
      return

c     [rayl]eigh,,a0,a1      - Rayleigh damping specification
c     [rayl,freq,zeta,w1,w2] - Input as two frequencies/damping ratio

26    if(pcomp(lct(l),'freq',4)) then
        rayla1 = 2.d0*ct(1,l)/(ct(2,l) + ct(3,l))
        rayla0 = rayla1*ct(2,l)*ct(3,l)
        if(ior.lt.0) then
          write(*,2020) ct(1,l),ct(2,l),ct(3,l),rayla0,rayla1
        endif
        write(iow,2020) ct(1,l),ct(2,l),ct(3,l),rayla0,rayla1
      else
        rayla0 = ct(1,l)
        rayla1 = ct(2,l)
        if(ior.lt.0) then
          write(*,2019) rayla0, rayla1
        endif
        write(iow,2019) rayla0, rayla1
      endif
      return

c     [cxso]lve,,omega - complex solution for frequence = omega rad/t

27    call pcxsol(ct(1,l),prt)
      return

c     [broy],iter,number  - broyden update for unsymmetric matrices

28    if(niter.gt.0) then
        k1 = max(1,nint(ct(1,l)))
        write(iow,2027) k1
        if(ior.lt.0) then
          write(*,2027) k1
        endif
        setvar = palloc(208,'VTILD',neq*(k1+1),2)
        setvar = palloc(209,'DELTX',neq*k1,2)
        call pzero(hr(np(208)),neq*k1+neq)
        call pzero(hr(np(209)),neq*k1)

c       Set first Delta-x vector from solution

        do i = 0,neq-1
          hr(np(209)+i) = hr(np(26)+i)
        end do ! i
        call broyden(hr(np(208)+neq),hr(np(209)),hr(np(208)),k1)
      else
        write(iow,2028)
      endif
      return

c     [rect]angular form for displacement, velocity, and acceleration outputs

30    iview = 0

      zview0(1) = ct(1,l)
      zview0(2) = ct(2,l)
      zview0(3) = ct(3,l)

      if(ior.lt.0) then
        write(*,2037) (i,zview0(i),i=1,ndm)
      endif

      cview = .false.
      return

c     [cyli],<pola> form for displacement, velocity, and acceleration outputs

31    iview = 1

      zview0(1) = ct(1,l)
      zview0(2) = ct(2,l)
      zview0(3) = ct(3,l)
      if(pcomp(lct(l),'pola',4)) then
        cview = .true.
      else
        cview = .false.
      endif

      if(ior.lt.0) then
        write(*,2038) (i,zview0(i),i=1,ndm)
      endif
      return

c     [sphe]rical form for displacement, velocity, and acceleration outputs

32    iview = 2

      zview0(1) = ct(1,l)
      zview0(2) = ct(2,l)
      zview0(3) = ct(3,l)

      if(ior.lt.0) then
        write(*,2039) (i,zview0(i),i=1,ndm)
      endif
      return

c     [forc]e prints

33    k1  = nint(ct(1,l))
      k2  = nint(ct(2,l))
      k3  = nint(ct(3,l))
      if(k2.eq.0) k2 = k1
      if(k3.eq.0) k3 = 1
      pfro = pfr
      pfr  = .true.

c     Set flag for local node prints

      pfeap_gnod = .false.

c     Set for all

      if(pcomp(lct(l),'all ',4)) then
        k1 = 1
        k2 = numnp
        k3 = 1

c     Set for specified coordinate

      elseif(pcomp(lct(l),'coor',4)) then
        nxt = max(1,min(k1,ndm))
        xt  = ct(2,l)
        if(ct(3,l).eq.0.0d0) then
          xtol = 1.0d-6*dxmsh(nxt)
        else
          xtol = abs(ct(3,l))
        endif
        k1 = 1
        k2 = numnp
        k3 = 1

c     Set for specified node

      elseif(pcomp(lct(l),'node',4)) then

        call pclnod(ct(1,l),hr(fp(4)),mr(np(190)),ndm,numnp, k1)
        k2 = k1
        k3 = 1

c     Set for global node prints

      elseif(pcomp(lct(l),'gnod',4) .and. processor.gt.1) then

        pfeap_gnod = .true.

c     Do range

      else
        k1 = max(1,min(k1,numnp))
        k2 = max(1,min(numnp,k2))
        if(k2-k1.ne.0) k3 = sign(k3,k2-k1)
      endif
      call prtdis(hr(np(43)),hr(np(27)),ttim,prop,ndm,ndf,
     &            k1,k2,k3,4,1,ndf,prth)
      pfr = pfro
      return

c     Formats

1000  format(a)

2001  format('   Number of eigenpairs reduced to',i4,' by number of',
     &       ' nonzero matrix diagonal terms')

2003  format('   Scaling factor tau = ',1p,1e15.5,' using phi',i2)

2004  format('   Bifurcation check: f*phi',i2,' = ',1p,1e15.5)


2006  format('   Use SHOW DICTionary to check memory useage'/)

2007  format('   Activate   region',i4)
2008  format('   Deactivate region',i4)

2010  format('   Arc length set to OFF.')
2011  format('   Arc length set to ON.')

2012  format('   Restart from : ',a)
2013  format('   Save to file : ',a)

2014  format('   Read new mesh data from file : ',a)

2016  format('   Setting initial strain for region =',i5)

2017  format(/'   Press <ENTER> to continue.')

2018  format( '   Eigenpair active DOF (1 = active; 0 = inactive)'/
     &      (7x,6(i3,'-dof =',i3)))

2019  format( '   Rayleigh Damping Parameters'/
     &        '               Mass      Factor (a0) =',1p,1e12.5/
     &        '               Stiffness Factor (a1) =',1p,1e12.5/)

2020  format( '  Rayleigh Damping Frequency Parameters'/
     &        '               Damping value  (zeta) =',1p,1e12.5/
     &        '               Frequency-1    (omg1) =',1p,1e12.5/
     &        '               Frequency-2    (omg2) =',1p,1e12.5/
     &        '               Mass      Factor (a0) =',1p,1e12.5/
     &        '               Stiffness Factor (a1) =',1p,1e12.5/)

2021  format( '  Contact surface check performed. No check during',
     &        ' iterations.'/)
2022  format( '  Contact surface check not performed. Check during',
     &        ' iterations.'/)

2023  format( '  Contact surface check ',a/)

2024  format( '  Contact friction coefficient =',1p,1e12.5/)
2025  format( '  Contact friction coefficient = 0.00000e+00 (OFF)'/)

2026  format( '  Eigenvalues for ',a,' of element',i8/)

2027  format( '  Broyden Solution: Number iterations =',i4/)
2028  format( '  Must use TANG or UTAN before BROYden'/)

2029  format(10x,'Region',i4,' activated')
2030  format(10x,'Region',i4,' deactivated')
2031  format(10x,'Region',i4,' does not exist')

2032  format( 10x,'No test on state imposed'/)
2033  format( 10x,'Test on state imposed: Size =',1p,1e12.4/)

2034  format( '  Optimize profile using Wilson/Hoit algorithm'/)
2035  format( '  Optimize profile using Scott Sloan algorithm'/)

2036  format( '  Optimize profile during contact solution.'/)

2037  format( '  Rectangular Cartesian coordinate outputs with center',
     &        ' at:',/,3('    x(',i1,') = ',1p,1e12.5:))

2038  format( '  Cylindrical coordinate outputs with center at:',/,
     &      3('    x(',i1,') = ',1p,1e12.5:))

2039  format( '  Spherical coordinate outputs with center at:',/,
     &      3('    x(',i1,') = ',1p,1e12.5:))

3000  format(/' *WARNING* Unable to compute incremental acceleration,'
     &       ,'           No mass matrix defined'/)

3001  format(/' *ERROR* Problem not dynamic - no output'/)

3002  format(/' *ERROR* ',a,': No stiffness matrix, use TANG or UTAN'/)

3003  format(/' *ERROR* SUBS: No mass matrix, use MASS or IDEN'/)

3004  format(/' *ERROR* ARCL: Compute eigenvectors first'/)

3005  format(/' *WARNING* Complex not active, no imaginary output'/)


3007  format(/' *ERROR* ',a,': Must use SUBSpace command first'/)

3008  format(/' *ERROR* Cannot ACTIVATE/DEACTIVTE region number',i4/)

3009  format(/' *ERROR* SUBS: Mass has no terms. May need to specify',
     &        ' a density.'/)

3010  format(/' *ERROR* SUBS: Geometric tangent has no terms.  Compute'
     &       /'               solution before GEOM to get values.'/)

3011  format(/' *ERROR* ',a,': Use PSUBspace for parallel eigensolve')

4000  format(' Input DOF for eigen computations')

      end
