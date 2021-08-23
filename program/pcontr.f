c$Id: pcontr.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pcontr()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Control program for FEAP problem input and solution.

c      Inputs:
c        none

c      Outputs:
c        none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'allotd.h'
      include   'allotn.h'
      include   'augdat.h'
      include   'bdata.h'
      include   'cblend.h'
      include   'cdata.h'
      include   'cdat1.h'
      include   'cdat2.h'
      include   'chdata.h'
      include   'codat.h'
      include   'corset.h'
      include   'cornum.h'
      include   'comfil.h'
      include   'compac.h'
      include   'complx.h'
      include   'conval.h'
      include   'counts.h'
      include   'crotas.h'
      include   'c_tanfl.h'
      include   'debugs.h'
      include   'dstars.h'
      include   'edgdat.h'
      include   'eltran.h'
      include   'endata.h'
      include   'errchk.h'
      include   'gltran.h'
      include   'hdatam.h'
      include   'hlpdat.h'
      include   'idata1.h'
      include   'iodata.h'
      include   'iofile.h'
      include   'ioincl.h'
      include   'iosave.h'
      include   'linka.h'
      include   'mdata.h'
      include   'mxsiz.h'
      include   'part0.h'
      include   'part1.h'
      include   'part3.h'
      include   'part7.h'
      include   'pdata5.h'
      include   'pdata6.h'
      include   'pdatps.h'
      include   'pfeapb.h'
      include   'pglob1.h'
      include   'pointer.h'
      include   'plflag.h'
      include   'prflag.h'
      include   'print.h'
      include   'pscal.h'
      include   'rdata.h'
      include   'rdat1.h'
      include   'refng.h'
      include   'region.h'
      include   'rigid1.h'
      include   'rjoint.h'
      include   'sdata.h'
      include   'setups.h'
      include   'tdata.h'
      include   'trdata.h'
      include   'umac1.h'
      include   'vdata.h'
      include   'comblk.h'

      include   'p_point.h'

      logical    bflg,errs,incf,initf,intr,intx,nopart,tfl,tief,tiefl
      logical    setvar,palloc,errc,pinput,tinput,vinput,nocount,pcomp
      logical    evint,gapfl,cprt,mchflg,oprt,oprth,newprob,usetfl(10)
      logical    cxifl,lp_in, naninfck
      character  titl*80,dnam*15, cdate*24, uset(10)*4, type*4, usub*4
      character  fnamr*132, fext*4, fnamp*128, flog*128,wfac(2)*8
      integer    pln,pre, iorsv
      integer    i, iii,irdef,iif, j, l1,l2,l3,l4,l5, nad,npd, nprob
      integer    usetno(10)
      real*4     etime, tary(2)
      real*8     gaptol,td(12)

      save

c     Default names for manipulation sets

      data      uset / 'man1', 'man2', 'man3', 'man4', 'man5',
     &                 'man6', 'man7', 'man8', 'man9', 'man0' /

      data      wfac / 'Material','Region' /

c     Destroy old output and log files if they exist

      inquire(file=fout,exist=initf)
      if(initf) then
        open (unit=iow,file=fout,status='old')
        close(unit=iow,          status='delete')
      endif

c     Set file for log

      call fileset(fout,flog,'O', 'L')

c     Open files for input, output and log

      open(unit=ior,file=finp,status='old')
      open(unit=iow,file=fout,status='new')
      open(unit=ilg,file=flog,status='new')
      call uscreen(1,iow)
      write(ilg,2018)

c     Initial values for include options

      chflg   = .false.
      cprt    = .true.
      everon  = .false.
      evint   = .false.
      hdcpy   = .false.
      incf    = .false.
      intr    = .false.
      intx    = .false.
      lp_in   = .true.
      mchflg  = .false.
      newprob = .false.
      nocount = .true.
      debug   = .false.
      lread   = .false.
      lsave   = .false.
      lmate   = .false.
      eofile  = .false.
      eralloc = .false.
      ucount  = .false.
      lfile   = ios
      icf     = icl
      isf     = 1
      iif     = 1
      irdef   = ior
      fincld(1) = finp
      irecrd(1) = 0

c     Set program constants

      call pvalues()

c     Set problem counter and initialize clock

      nprob    = 0
      call stime()

c     Initial value for variable storage

      cplxfl = .false.

c     Set default to print headers

      prth   = .true.

c     Flags for user manipulation commands

      do j = 1,10
        usetfl(j) = .false.
        usetno(j) = 0
      end do ! j

c     Install user functions

c     Set umesh input names

      do j = 1,10
        if(j.lt.10) then
          write(usub,'(a3,i1)') 'mes',j
        else
          write(usub,'(a4)') 'mes0'
        endif
        uct = usub
        call umshlib(j,prt)
        umshc(j) = uct
      end do ! j

c     Set umacro input names

      do j = 1,10
        if(j.lt.10) then
          write(usub,'(a3,i1)') 'mac',j
        else
          write(usub,'(a4)') 'mac0'
        endif
        uct   = usub
        fnamp = ' '
        call umaclib(j,fnamp,td,.false.)
        umacc(j) = uct
      end do ! j

c     Set umati model names

      do j = 1,10
        if(j.lt.10) then
          write(usub,'(a3,i1)') 'mat',j
        else
          write(usub,'(a4)') 'mat0'
        endif
        uct = 'mate'
        call uconst(usub,td,td,td,l1,l2,l3)
      end do ! j

c     Set uplot input names

      do j = 1,10
        if(j.lt.10) then
          write(usub,'(a3,i1)') 'plt',j
        else
          write(usub,'(a4)') 'plt0'
        endif
        uct = usub
        call upltlib(j,td,.false.)
        upltc(j) = uct
      end do ! j

c     Set umanipulation names

      do j = 1,10
        uct     = uset(j)
        call usetlib(j,prt)
        uset(j) = uct
      end do ! j

c     Input with interactive interactive statements

1     if(intx) then
        if(cprt) then
          call pprint(
     &        ' Continue with interactive input options for control?')
          call pprint('<y or n> :')
        endif
        errck = tinput(dnam,1,td,0)

c       Read command interactively

        if(pcomp(dnam,'y',1)) then
          call pprint(
     &        ' Specify command (INTEractive, INCLude, STOP, etc.) >')
          read (*,1001,err=900,end=910) yyy
          cprt = .true.

c       Read command from current file and turn off intx flag

        else
          cprt = .false.
          intx = .false.
          ior  =  abs(ior)
          read(ior,1001,err=900,end=910) yyy
        endif

c     Input from current file

      else
        ior = abs(ior)
        read(ior,1001,err=900,end=910) yyy
      endif

c     Compare with command list

      call pstrip(xxx,yyy,1)
      l1   = min(80,len(xxx))
      titl = xxx(1:l1)

c     Start solution of new problem

      if(pcomp(titl(1:4),'feap',4) .or.
     &   pcomp(titl(1:4),'geof',4)) then
        go to 100

c     Set count/nocount or input mode to binary

      elseif(pcomp(titl(1:4),'noco',4)) then
        nocount = .false.

      elseif(pcomp(titl(1:4),'coun',4)) then
        nocount = .true.

c     Set flag for constructing interface mesh

      elseif(pcomp(titl(1:4),'matc',4) .or.
     &       pcomp(titl(1:4),'clus',4)) then

c       Increase storage for local solution arrays

        iif    = 2
        setvar = palloc( 34,'LD   ',max(nen+1,nst*4,21),1)
        setvar = palloc( 35,'P    ',nst*6*ipc          ,2)
        setvar = palloc( 36,'S    ',nst*nst*4*ipc      ,2)
        setvar = palloc( 39,'TL   ',nen*2              ,2)
        setvar = palloc( 41,'UL   ',nst*14*ipc         ,2)
        setvar = palloc( 46,'ANGL ',nen*6              ,2)
        setvar = palloc( 44,'XL   ',max(4,nen)*6       ,2)

        if(pcomp(titl(1:4),'matc',4)) then
          mchflg = .true.
        else
          setvar = palloc(210,'MATCH', 8*numel, 1)
          call faclset(hr(np(25)),mr(np(32)),mr(np(33)),
     &                 mr(np(210)),l2, l4)
          setvar = palloc(210,'MATCH', 8*l2   , 1)

c         Allocate storage for interface history

          setvar = palloc(214,'HINTE',max(1,l4),2)
          setvar = palloc(215,'HINT1',max(1,hnimax),2)
          setvar = palloc(216,'HINT2',max(1,hnimax),2)
          setvar = palloc(217,'HINT3',max(1,hni3max),2)
        endif

c     Set input for a binary mode

      elseif(pcomp(titl(1:4),'bina',4)) then
        newprob = .true.
        bflg    = .true.
        call acheck(titl,yyy,18,80,80)
        fnamr   = ' '
        fnamr   = yyy(19:36)
        fext    = titl(1:4)
        call opnfil(fext,fnamr, 3,ios,lread)
        rewind ios
        read(ios) head,numnp,numel,nummat,ndm,ndf,nen,ndd,nud
        npd = ndd - nud - 1
        nad = 0
        go to 200

c     User command sets

      elseif(pcomp(titl(1:4),uset(1),4)) then
        usetno(1) = usetno(1) + 1
        usetfl(1) = .true.
        fext      = 'u1a'
        go to 300
      elseif(pcomp(titl(1:4),uset(2),4)) then
        usetno(2) = usetno(2) + 1
        usetfl(2) = .true.
        fext      = 'u2a'
        go to 300
      elseif(pcomp(titl(1:4),uset(3),4)) then
        usetno(3) = usetno(3) + 1
        usetfl(3) = .true.
        fext      = 'u3a'
        go to 300
      elseif(pcomp(titl(1:4),uset(4),4)) then
        usetno(4) = usetno(4) + 1
        usetfl(4) = .true.
        fext      = 'u4a'
        go to 300
      elseif(pcomp(titl(1:4),uset(5),4)) then
        usetno(5) = usetno(5) + 1
        usetfl(5) = .true.
        fext      = 'u5a'
        go to 300
      elseif(pcomp(titl(1:4),uset(6),4)) then
        usetno(6) = usetno(6) + 1
        usetfl(6) = .true.
        fext      = 'u6a'
        go to 300
      elseif(pcomp(titl(1:4),uset(7),4)) then
        usetno(7) = usetno(7) + 1
        usetfl(7) = .true.
        fext      = 'u7a'
        go to 300
      elseif(pcomp(titl(1:4),uset(8),4)) then
        usetno(8) = usetno(8) + 1
        usetfl(8) = .true.
        fext      = 'u8a'
        go to 300
      elseif(pcomp(titl(1:4),uset(9),4)) then
        usetno(9) = usetno(9) + 1
        usetfl(9) = .true.
        fext      = 'u9a'
        go to 300
      elseif(pcomp(titl(1:4),uset(10),4)) then
        usetno(10) = usetno(10) + 1
        usetfl(10) = .true.
        fext       = 'u0a'
        go to 300

c     Perform inputs from an include file

      elseif(pcomp(titl(1:4),'incl',4)) then
        call acheck(titl,yyy,15,80,80)
        read(yyy,1002,err=900,end=900) titl(1:4),dnam
        if(pcomp(dnam,'end',3)) then
          call pincld(dnam)
          if(evint) then
            write(*,2005) fnamr
          endif
          write(iow,2005) fnamr
        else
          fnamr =  dnam
          call pincld(dnam)
        endif
        incf = .true.
        cprt = .false.

c     Perform inputs for initial conditions

      elseif(pcomp(titl(1:4),'init',4)) then
        call acheck(titl,yyy,15,80,80)
        titl( 1: 4) = yyy( 1: 4)
        titl(16:19) = yyy(16:19)
        setvar = vinput(yyy(31:75),45,td(1),3)
        call pinitl(dnam,td,errs)
        if(errs) call plstop()

c     Set title prints on/off

      elseif(pcomp(titl(1:4),'titl',4)) then
        call acheck(titl,yyy,15,80,80)
        read(yyy,1002,err=900,end=911) titl(1:4),titl(16:19)
        if(pcomp(titl(16:19),'off',3)) then
          prth = .false.
        else
          prth = .true.
        endif

c     Solution mode

      elseif(pcomp(titl(1:4),'inte',4)) then
        ior   = -abs(ior)
        evint = .true.
        intr  = .true.
        intx  = .true.
        cprt  = .true.
        call pltcur()
        go to 400

      elseif(pcomp(titl(1:4),'batc',4) .or.
     &       pcomp(titl(1:4),'macr',4)) then
        intr = .false.
        go to 400

c     Manual level set: 0 = basic; 1 = advanced; 2 = expert

      elseif(pcomp(titl(1:4),'manu',4)) then
        call acheck(titl,yyy,15,80,80)
        read(yyy,1003,err=900,end=911) titl(1:4),hlplev
        hlplev = max(-1,min(3,hlplev))

c     Mesh manipulations: Link and tie

c     Reset id list to link dof's on different nodes - set by node #

      elseif(pcomp(titl(1:4),'link',4)) then
        call plinka('lnk ','set','   ')
        lkflg = .true.

c     Reset id list to link dof's on different nodes - set by coord.

      elseif(pcomp(titl(1:4),'elin',4)) then
        call plinka('eln ','set','   ')
        leflg = .true.

      elseif(pcomp(titl(1:4),'clin',4)) then
        call acheck(titl,yyy,15,80,255)
        read(yyy,1003,err=900,end=911) titl(1:4),(clnk(i),i=1,ndf)
        if(prt) then
          write(iow,2016) (i,clnk(i),i=1,ndf)
        endif
        lcflg = .true.
        if(np(79).eq.0) then
          setvar = palloc( 79,'IPOS ',numnp,  1)
          call pseqn(mr(np(79)),numnp)
        endif
        setvar = palloc(111,'TEMP1',numnp, 1)
        setvar = palloc(112,'TEMP2',numnp, 1)

        l1 = 1
        l2 = numnp
        l3 = 0
        l4 = 0
        l5 = 0
        td(2) = 0.0d0
        if(.not.gapfl) then
          gaptol = 1.d-3/sqrt(dble(max(1,numnp)))
        endif
        call tienod(mr(np(33)),hr(np(43)),mr(np(79)),mr(np(111)),
     &              mr(np(112)),mr(np(78)),ndm,nen,nen1,
     &              numnp,numel,l1,l2,l3,l4,l5,gaptol,td(2))

        setvar = palloc(112,'TEMP2',0, 1)
        setvar = palloc(111,'TEMP1',0, 1)

      elseif(pcomp(titl(1:4),'tie' ,3)) then
        go to 500

c     Parameter sets

      elseif(pcomp(titl(1:4),'para',4) .or.
     &       pcomp(titl(1:4),'cons',4)) then
        coflg = .true.
        call pconst(prt)

c     Rigid body and master/slave sets

c     Create rigid bodies

      elseif(pcomp(titl(1:4),'rigi',4)) then
        call rigidb(1,1,errs)
        if(errs) call plstop()
        tfl = .true.

c     Create list of joints

      elseif(pcomp(titl(1:4),'join',4)) then
        call rigidb(1,2,errs)
        if(errs) call plstop()

c     Rigid body loads

      elseif(pcomp(titl(1:4),'rloa',4)) then
        call rigidb(1,3,errs)
        if(errs) call plstop()

c     Rigid body boundary conditions

      elseif(pcomp(titl(1:4),'rbou',4)) then
        call rigidb(1,4,errs)
        if(errs) call plstop()

c     Rigid body displacements

      elseif(pcomp(titl(1:4),'rdis',4)) then
        call rigidb(1,5,errs)
        if(errs) call plstop()

c     Master/slave

      elseif(pcomp(titl(1:4),'mast',4)) then
        setvar = palloc( 167, 'RLINK', ndf*numnp, 1)
        call pmastr(mr(np(100)),mr(np(167)),mr(np(190)),hr(np(43)),prt)

c     Optimize profile

      elseif(pcomp(titl(1:4),'opti',4)) then
        call acheck(titl,yyy,15,80,80)
        read(yyy,1002,err=900,end=911) titl(1:4),dnam
        if(pcomp(dnam,'hoit',4)) then
          opthoit = .true.
        else
          opthoit = .false.
        endif
        call optid()
        tfl = .true.

c     Dictionary search for name of array

      elseif(pcomp(titl(1:4),'dict',4)) then
        call acheck(titl,yyy,15,80,80)
        read(yyy,1002,err=900,end=911) titl(1:4),dnam
        write(*,*) 'Looking for ',dnam
        call pgetd(dnam,point,pln,pre,errs)
        write(*,*) point,pln,pre,errs

c     Loop start

      elseif(pcomp(titl(1:4),'loop',4)) then
        call acheck(titl,yyy,15,80,80)
        read(yyy,1002,err=900,end=911) titl(1:4),dnam
        call ploops(lp_in,dnam,1)

c     Loop end

      elseif(pcomp(titl(1:4),'next',4)) then
        call acheck(titl,yyy,15,80,80)
        read(yyy,1002,err=900,end=911) titl(1:4),dnam
        call ploops(lp_in,dnam,2)

c     Mesh generator outputs

      elseif(pcomp(titl(1:4),'tri2',4)) then
        call tridat(hr(np(43)),mr(np(33)),numnp,numel,nen,ndm,nen1,
     &              evint)

c     Unit multipliers for length, force and time

      elseif(pcomp(titl(1:4),'unit',4)) then
        call punits()

c     Data storage modes: Real/complex

      elseif(pcomp(titl(1:4),'*rea',4)) then
        cplxfl = .false.
        write(iow,2006)

      elseif(pcomp(titl(1:4),'*com',4)) then
        cplxfl = .true.
        write(iow,2007)

c     Contact data inputs

      elseif(pcomp(titl(1:4),'cont',4)) then
        backspace (ior)
        call plinka('cxi ','set','end')
        cxifl = .true.

c     Partition/ODE order sets

      elseif(pcomp(titl(1:4),'part',4)) then
        errc = pinput(td,ndf)
        call pzeroi(nqp , 5)
        do i = 1,ndf
          ndfp(i)      = td(i)
          ndfg(i)      = ndfp(i)
          nqp(ndfp(i)) = ndfp(i)
        end do ! i
        write(iow,2002) (i,ndfp(i),i=1,ndf)
        do i = 1,4
          npart = nqp(i)
          if(npart.ge.1 .and. npart.le.4) then
            call partpt(npart,tflp(npart),.false.)
          elseif(npart.ne.0) then
            write(  *,3000) npart
            write(ilg,3000) npart
          endif
        end do ! i

c       Set initial allocation for monolithic case

        npart = 5
        call partpt(npart,tflp(npart),.false.)

      elseif(pcomp(titl(1:4),'orde',4)) then
        errc = pinput(td,ndf)
        do i = 1,ndf
          ndfo(i) = td(i)
          ndog(i) = ndfo(i)
        end do ! i
        write(iow,2003) (i,ndfo(i),i=1,ndf)

c     Remarks to output file

      elseif(pcomp(titl(1:4),'rema',4)) then
        write(*,2008) titl(1:78)

c     Debug set

      elseif(pcomp(titl(1:4),'debu',4)) then
        debug = .true.

c     Domain inputs for parallel solutions

      elseif(pcomp(titl(1:4),'doma',4)) then

       errck = pinput(td,3)

       numpn  = nint(td(1))
       numtn  = nint(td(2))
       numteq = nint(td(3))

       write(iow,2023) numpn, numtn, numteq

       call pdomain(prt)

c     Stop execution

      elseif(pcomp(titl(1:4),'stop',4)) then
        if(evint) write(*,2004) fout
        if(abs(ior).eq.irdef) then
          call plstop()
        endif

c     Return without a call to plstop

      elseif(pcomp(titl(1:4),'cntn',4)) then
        return

      endif

c     Read again

      go to 1

c     Start Problem: Read and print control information

100   newprob = .true.
      do i = 1,20
        l2      = 4*i
        l1      = l2 - 3
        head(i) = titl(l1:l2)
      end do ! i
      bflg   = .false.
      gapfl  = .false.
      incf   = .false.
      intr   = .false.
      intx   = .false.
      call fdate( cdate )
      errc   = pinput(td,9)
      numnp  = td(1)
      numel  = td(2)
      nummat = td(3)
      ndm    = td(4)
      ndf    = td(5)
      nen    = td(6)
      nad    = td(7)
      npd    = td(8)
      nud    = td(9)
      nnn    = 0

c     Adjust storage for material parameters

      npd    = max(npd,250)
      nud    = max(nud,150)
      ndd    = npd + nud + 1

c     Star node/element initialization

      starnd = 0
      starel = 0

c     Blending function initialization

      numsn  = 0
      numsd  = 0
      numbd  = 0

c     Contact array initialization

      numcels = 0
      ncen    = 0
      optflg  = .false.
      optmsh  = .false.
      opthoit = .false.

c     Serial & parallel solution by unblocked equations

      pfeap_blk = .false.

c     Set filenames for multiple problem case

      if(irdef.ne.ior) then

        inquire(unit=ior,name=fnamp,exist=errs)

        if(errs) then

c         Clear plot files, delete scratch files and close output file

          if(hdcpy) call fpplcl()
          if(niter.gt.1) then
            write(ilg,3004) nstep,niter,nform,ttim,dt,rel0,rnorm,rnmax,
     &                      aengy,etime(tary)
          else
            write(ilg,3005) nstep,niter,nform,ttim,dt,rel0,rnorm,rnmax,
     &                      etime(tary)
          endif
          titer = titer + niter
          tform = tform + nform
          write(ilg,3006) titer,tform
          rfl = .false.
          call ptimpl()
          call pdelfl()

c         Delete memory use

          do i = ndict,1,-1
            setvar = palloc(dlist(i),dict(i),0,iprec(i))
          end do ! i

          i = index(flog,' ')
          if(nprob.eq.0) write(iow,2017) flog(1:i-1)
          close(unit=iow)

c         Extract file name

          i = index(fnamp,' ')
          if(i.eq.0) i = 128
          do j = i,1,-1
            if(pcomp(fnamp(j:j),char(92),1)) go to 110
          end do ! j
          j = 0
110       fnamr = fnamp(j+1:j+21)

c         Set new plot file name

          fnamr(1:1) = 'P'
          fplt(1:128) = ' '
          fplt(1: 17) = fnamr
          i = index(fplt,'.')
          if(i.gt.0) then
            fplt(i: 21) = ' '
          endif
          i = min(index(fplt,' '), 16)
          if(i.eq.0) then
            i = 16
          endif

c         Add problem counter to name

          nprob = nprob + 1
          write(fplt(i:i+2),'(a)') '000'
          if(nprob.lt.10) then
            write(fplt(i+2:i+2),'(i1)') nprob
          elseif(nprob.lt.100) then
            write(fplt(i+1:i+2),'(i2)') nprob
          elseif(nprob.lt.1001) then
            write(fplt(  i:i+2),'(i3)') nprob
          else
            write(*,*) 'Exceeded limit of multiple files (PCONTR)'
          endif

c         Set file names for new problem

          fout = fplt
          fout(1:1) = 'O'
          fres = fplt
          fres(1:1) = 'R'
          fsav = fplt
          fsav(1:1) = 'S'

c         Create clean output file

          inquire(file=fout,exist=initf)
          if(initf) then
            open (unit=iow,file=fout,status='old')
            close(unit=iow,          status='delete')
          endif
          open(unit=iow,file=fout,status='new')
          if(nprob.gt.1) write(ilg,2019)
          write(ilg,2020) nprob,fout

c       Error in file structure

        else
          write(  *,3003)
          write(ilg,3003)
          call plstop()
        endif
      endif

c     If number of nodes, or elements is zero compute number from data

      if(nocount) then
        ucount      = .true.
        call pnums()
        irecrd(isf) = 2
        ucount      = .false.

c       Star node/element re-initialization

        starnd = 0
        starel = 0
      endif

c     Output problem size data

200   write(iow,2000) head,cdate,versn,fincld(isf),
     &               numnp,numel, ndm,ndf,nad,nen, nummat,npd,nud
      call stime()

c     Set parameters for page eject and rotation dof

      o   = '    '
      errck = .false.
      lsave = .false.
      lkflg = .false.
      leflg = .false.
      lcflg = .false.
      initf = .true.
      cxifl = .false.
      eanfl = .false.
      ebcfl = .false.
      ebsfl = .false.
      edifl = .false.
      efcfl = .false.
      eprfl = .false.
      surfl = .false.
      basfl = .false.
      boufl = .false.
      cprfl = .false.
      disfl = .false.
      forfl = .false.
      lfrfl = .false.
      angfl = .false.
      reafl = .false.
      intfl = .false.
      tiefl = .true.
      tief  = .false.
      damfl = .false.
      masfl = .false.
      stifl = .false.

c     Dynamic contact flags

      rattlfl = .false.
      shakefl = .false.

      do i = 1,50
        ia(1,i)  = 1
        ia(2,i)  = 2
        ir(1,i)  = 0
        ir(2,i)  = 0
        ea(1,i)  = 1
        ea(2,i)  = 2
        er(1,i)  = 0
        er(2,i)  = 0
        ia3(1,i) = 1
        ia3(2,i) = 2
        ia3(3,i) = 3
        ir3(1,i) = 0
        ir3(2,i) = 0
        ir3(3,i) = 0
        ea3(1,i) = 1
        ea3(2,i) = 2
        ea3(3,i) = 3
        er3(1,i) = 0
        er3(2,i) = 0
        er3(3,i) = 0
        inord(i) = 0
        exord(i) = 0
        do j = 1,30
          ipord(j,i) = 0
          epord(j,i) = 0
        end do ! j
      end do ! i
      nadd   = 0
      nprof  = 0
      nsurf  = 0
      nbouf  = 0
      ndisf  = 0
      nforf  = 0
      nforl  = 0
      nangf  = 0
      nintf  = 0
      neang  = 0
      nebcs  = 0
      nedis  = 0
      nefrc  = 0
      nepro  = 0
      ndamf  = 0
      nmasf  = 0
      nstif  = 0
      nbasf  = 0
      neule  = 0

c     Zero global parameters

      gtypfl = .false.
      gdeffl = .false.
      gtdofl = .false.
      gomgfl = .false.
      grayfl = .false.
      groufl = .false.
      if(    ndm.eq.2) then
        g2type = 2           ! default plane strain
      elseif(ndm.eq.3) then
        g2type = 7           ! default 3-d
      else
        g2type = 8           ! other cases
      endif
      gdtype = 1
      gtdof  = 0
      gref   = 0
      do i = 1,3
        units(i)  = 1.0d0
        grefx(i)  = 0.0d0
        gtref(i)  = 0.0d0
        gomega(i) = 0.0d0
        gomex(i)  = 0.0d0
        gomev(i)  = 0.0d0
      end do ! i
      do i = 1,2
        gray(i) = 0.0d0
      end do ! i
      do i = 1,14
        gfac(i) = 0.0d0
      end do ! i
      gquadn =  0.0d0        ! Default quadrature: Gauss
      augf   =  1.0d0
      augg   =  0.0d0
      auggfl = .false.
      augmfl = .false.
      lvaug  =  0
      lvsol  =  0

c     Zero pointer array

      setvar = palloc( 0 ,'START', 0 , 0 )

c     Set initial values for nan and inf checks

      setvar = naninfck(td(1),1,0)

c     Set contact flags for a new problem

      call contact (300)

c     Set rigid body flags for a new problem

      call rigidb (0,0,errs)

c     Set pointers for allocation of mesh arrays

      nen1      = nen + 8
      nie       = 11
      nst       = max(nen*ndf + nad,1)
      ndict     = 0
      nneq      = ndf*numnp
      do i = 1,ndf
        ndfo(i) = 0
        ndog(i) = 10
      end do ! i

c     Set pointers for allocation of mesh arrays

      if(cplxfl) then
        ipc = 2
      else
        ipc = 1
      end if

c     Allocate size for arrays for mesh and solution vecors

      l1   = ndm*numnp
      l2   = max(ndf*numnp,1)
      l3   = max(nen+1,7*nst,21)
      l4   = numnp*max(ndf,ndm)*ipc
      l5   = ndf*nen

c     Allocate and zero arrays

      setvar = palloc( 26,'DR   ',l4          ,  2)
      setvar = palloc( 34,'LD   ',l3          ,  1)
      setvar = palloc( 35,'P    ',nst*3       ,  2)
      setvar = palloc( 36,'S    ',nst*nst*2   ,  2)
      setvar = palloc( 39,'TL   ',nen         ,  2)
      setvar = palloc( 41,'UL   ',nst*14      ,  2)
      setvar = palloc( 44,'XL   ',max(4,nen)*3,  2)
      setvar = palloc( 25,'D    ',nummat*ndd  ,  2)
      setvar = palloc( 32,'IE   ',nummat*nie  ,  1)
      setvar = palloc(240,'IEDOF',nummat*l5   ,  1)
      setvar = palloc( 31,'ID   ',l2*2        ,  1)
      setvar = palloc( 33,'IX   ',nen1*numel  ,  1)
      setvar = palloc(190,'NDTYP',numnp       ,  1)
      setvar = palloc(100,'RIXT ',numnp       ,  1)
      setvar = palloc(181,'RBEN ',numel       ,  1)
      setvar = palloc( 43,'X    ',l1          ,  2)
      setvar = palloc( 45,'ANG  ',numnp       ,  2)
      setvar = palloc( 46,'ANGL ',nen         ,  2)
      setvar = palloc( 27,'F    ',2*l2        ,  2)
      setvar = palloc( 28,'F0   ',4*l2        ,  2)
      setvar = palloc( 29,'FPRO ',2*l2        ,  1)
      setvar = palloc( 30,'FTN  ',4*l2        ,  2)
      setvar = palloc( 38,'T    ',numnp       ,  2)
      setvar = palloc( 40,'U    ',3*l2*ipc    ,  2)
      setvar = palloc( 89,'NREN ',numnp*2     ,  1)

c     Set initial numbering in renumber vector and mark nodes as unused.

      do i = 0,numnp-1
        mr(np( 89)+i      ) = i+1  ! Remap list
        mr(np( 89)+i+numnp) = i+1  ! Reverse list
        mr(np(190)+i      ) = 0
      end do ! i

c     Set pointers for allocation of loads/partition data

      do j = 1,5
        tflp(j)   = .true.
        flp(9,j)  = .false.
        scale(j)  = .false.
        nittyp(j) = -3
        nsolver(j)= solver
        nitolp(j) =  1.d-08
        natolp(j) =  1.d-16
      end do ! j
      do j = 1,ndf
        ndfp(j) = 1
        ndfg(j) = ndfp(j)
      end do ! j
      nqp(1)    =  1
      npart     =  1
      nopart    = .true.

c     Open file to store material data

      fmtl      = finp
      fmtl(1:1) = 'M'
      open(unit=iwd,file=fmtl,status='unknown')
      rewind iwd

c     Input a mesh from binary file (if it exists)

      if(bflg) then
        call bmesh(mr(np(32)),hr(np(25)),mr(np(31)+nneq),hr(np(43)),
     &             mr(np(33)),hr(np(27)),hr(np(38)),hr(np(45)),
     &             ndd,nie,ndm,ndf,nen,nen1,numnp,numel,nummat)
        close(ios)
        prt   = .true.
        iii   = -2
      else
        iii   =  0
      endif

c     Input mesh data from file

      call pmesh(iii,prt,prth)

c     Perform simple check on mesh to ensure basic data was input

      setvar = palloc(111,'TEMP1',max(numnp*ndf,1), 1)
      call meshck(mr(np(111)),mr(np(32)),mr(np(240)),mr(np(31)+nneq),
     &            mr(np(190)),mr(np(33)),nie,nen,nen1,ndf,numnp,numel,
     &            nummat,errs)
      setvar = palloc(111,'TEMP1',0, 1)
      if(errs) then
        call plstop()
      endif

c     Compute boundary nodes (before ties)

      if(tiefl) then
        setvar = palloc( 78,'EXTND',numnp   ,1)
        call pextnd()
        tiefl  = .false.
      endif

      tfl = .true.
      go to 1

c     [mani] - Perform user manipulation commands

300   errs  = .true.
      j     = 0
      do while(errs .and. j.lt.26)
        j     = j + 1
        fnamr = fsav
        write(fext(3:3),'(a1)') char(96+j)
        call addext(fnamr,fext,128,4)
        inquire(file = fnamr, exist = errs)
      end do !
      call plinka(fext,'set','   ')
      go to 1

c     Establish profile of resulting equations for stiffness, mass, etc
c     [batc]h and [inte]ractive execution

400   if(.not.newprob) then
        write(  *,3001)
        write(ilg,3001)
        call plstop()
      elseif(intx .and. .not.intr .and. .not.incf) then
        write(  *,3002)
        write(ilg,3002)
        go to 1
      endif

      if(tfl) then

c       Check if nadd greater than nad

        if(nadd.gt.nad) then

          nst    = max(ndf*nen + nadd,1)
          setvar = palloc( 34,'LD   ',max(nen+1,7*nst*iif,21), 1)
          setvar = palloc( 35,'P    ',nst*3*iif*ipc          , 2)
          setvar = palloc( 36,'S    ',nst*nst*2*iif*ipc      , 2)
          setvar = palloc( 41,'UL   ',nst*14*ipc             , 2)

        endif

c       If ties have occurred merge boundary conditions, forces & contact

        if(tief) then
          call pmerbc(nopart,prt)
c         tief = .false.
        endif

c       Compute active boundary loading after ties

        call pecmes()

c       Compute boundary nodes (after ties)

        call pextnd()

c       Rotational 5/6 degree-of-freedom checks

        if(frotas) then
          call rotred(mr(np(81)),mr(np(31)+nneq),mr(np(190)),ndf,numnp)
        endif

c       Allocate memory to store all possible equations
c       (include: rigid body dof x rigid body)
c       (include: rigid joints   x 6         )

        do j = 1,4
          neq = 0
          do i = 1,ndf
            if(ndfp(i).eq.j) then
              neq = neq + numnp
            endif
          end do ! i
          if(nrbprt.eq.j) then
            neq = neq + nrbdof*nrbody + numjts*6
          elseif(numjts.gt.0) then
            neq = neq + numjts*6
          endif
          if(neq.gt.0) then
            write(dnam,'(a2,i1)') 'JP',j
            setvar = palloc( 20+j, dnam, neq, 1)
          endif
        end do ! j

c       Set user commands

        do j = 1,10
          fext = 'u1a'
          write(fext(2:2),'(i1)') j
            if(usetfl(j)) then
            do l3 = 1,26
              write(fext(3:3),'(a1)') char(96+l3)
              fnamr =  fsav
              call addext(fnamr,fext,128,4)
              inquire(file = fnamr, exist = errs)
              if(errs) then
                call opnfil(fext,fnamr,-1,ios,prt)

c               Read data from file

                iorsv = ior
                ior   = ios

                do l1 = 0,36
                  do l2 = 1,26
                    vvsave(l2,l1) = vvv(l2,l1)
                  end do ! l2
                end do ! l1
                do l1 = 1,3
                  x0sav(l1) = x0(l1)
                end do ! l1
                oprt  = prt
                oprth = prth

                read(ior,1004) type,fincld(isf),irecrd(isf),prt,prth
                read(ior,1005) vvv
                read(ior,1005) tr,xr,trdet,x0

                call usetlib(j,prt)

                close(ior,status='delete')
                ior   = iorsv

                do l1 = 0,36
                  do l2 = 1,26
                    vvv(l2,l1) = vvsave(l2,l1)
                  end do ! l2
                end do ! l1
                do l1 = 1,3
                  x0(l1) = x0sav(l1)
                end do ! l1
                prt  = oprt
                prth = oprth

              endif
            end do ! l3
          endif
        end do ! j

c       Set default partition data

        if(nopart) then
          nopart = .false.
          npart  = 1
          call partpt(npart,tflp(npart),.false.)
        endif

c       Construct interface mesh

        if(mchflg) then

c         Compute sparse structure for matrix

          setvar = palloc(111,'TEMP1', numnp+1, 1)         ! IC
          call optic(numnp,numel,numcels,nen,nen1,ncen,ncen1,
     &               mr(np(33)),mr(np(168)),mr(np(111)), l1)

          setvar = palloc(112,'TEMP2', l1, 1)              ! IP
          call pzeroi(mr(np(112)),l1)

          setvar = palloc(113,'TEMP3', numel+numcels, 1)   ! NNEL
          call opcon(numel,numcels,nen,nen1,ncen,ncen1,mr(np(33)),
     &                mr(np(168)),mr(np(111)),mr(np(112)),mr(np(113)))
          setvar = palloc(113,'TEMP3', 0, 1)

c         Set an upper bound for number of faces

          setvar = palloc(210,'MATCH', 50*numel, 1)

          call faceset(hr(np(25)),mr(np(32)),mr(np(33)),mr(np(111)),
     &                 mr(np(112)),mr(np(210)),l2, l4)

          setvar = palloc(210,'MATCH', 8*l2, 1)

          setvar = palloc(112,'TEMP2', 0, 1)
          setvar = palloc(111,'TEMP1', 0, 1)

c         Allocate storage for interface history

          setvar = palloc(214,'HINTE',max(1,l4),2)
          setvar = palloc(215,'HINT1',max(1,hnimax),2)
          setvar = palloc(216,'HINT2',max(1,hnimax),2)
          setvar = palloc(217,'HINT3',max(1,hni3max),2)

        endif
      endif

c     Check if contact has been specified

      if(cxifl) then

c       Generate the contact surfaces, pairs and material sets

        call pextnd()
        fext = 'cxi'
        call pinpfl('PCONTR',fext, type, 1)
        call contact (1)
        call pinpfl('PCONTR',fext, type, 2)

c       Check contact surfaces for eliminated tied nodes

        if(tief) then
          call contact (312)
          tief = .false.
        endif

        call contact (313)

        cxifl = .false.
        tfl   = .true.

      endif

c     Set initial history in elements and contact

      if(tfl) then

c       Set up stress history addresses

        call sethis(mr(np(32)),mr(np(33)),mr(np(181)),nie,nen,nen1,
     &              numel,nummat,prt)

c       Initialize history database items: Include 2 x for DG interface

        if(nhmax.gt.0) then
          setvar = palloc( 50,'NH1  ', nhmax*2,  2)
          setvar = palloc( 51,'NH2  ', nhmax*2,  2)
        endif
        if(nh3max.gt.0) then
          setvar = palloc( 52,'NH3  ', nh3max*2, 2)
        endif

        hflgu  = .true.
        h3flgu = .true.
        ctan(1) = 1.d0
        ctan(2) = 0.d0
        ctan(3) = 0.d0

c       call formfe(   u  ,   dr ,   dr ,   dr
        call formfe(np(40),np(26),np(26),np(26),
     &             .false.,.false.,.false.,.false.,14,1,numel,1)
      endif

c     Determine current profile

      if(tfl) then
        do j = 0,nneq-1
          mr(np(31)+j) = mr(np(31)+j+nneq)
        end do ! j

        mxpro = 0
        mxneq = 0
        do j = 1,4
          if(.not.tflp(j)) then
            if(prt.and.ior.lt.0) write(  *,2001) j
            if(prt)              write(iow,2001) j
            call partpt(j,tflp(j),.false.)
            mxprop(j) = 0
            mxneqp(j) = 0

c           Set current profile

            call profil(mr(np(20+j)),mr(np(34)),mr(np(31)),
     &                  mr(np(33)),1,prt)
            call profil(mr(np(20+j)),mr(np(34)),mr(np(31)),
     &                  mr(np(33)),2,prt)
            nqp(j)    = neq
            nqr(j)    = neqr
            mxprop(j) = max(mxprop(j),(mr(np(20+j)+neq-1))*ipc)
            mxneqp(j) = max(mxneqp(j),neq*ipc)
            mxpro     = max(mxpro,mxprop(j))
            mxneq     = max(mxneq,mxneqp(j))
          endif
        end do ! j
        if(pfeap_blk) then  ! Form when all equations are included
          call pidreset(mr(np(31)))
        endif

        tfl = .false.

      endif

c     Macro module for establishing solution algorithm (partition 1)

      call partpt(1,tflp(1),.false.)
      call pmacr(initf)
      go to 1

c     Tie nodes within tolerance of one another
c     [tie ] - merge regions with common coordinates

500   call acheck(titl,yyy,15,90,90)
      titl( 1: 4) = yyy( 1: 4)
      titl(16:19) = yyy(16:19)
      setvar = vinput(yyy(31:90),60,td(1),4)

c     Retrieve current boundary connection status

      if(.not.tief) then
        setvar = palloc( 79,'IPOS ',numnp,  1)
        call pseqn(mr(np(79)),numnp)
        tief = .true.
      endif

c     Tie line elements to regions

      if(pcomp(titl(16:19),'line',4)) then
        l2 = max(    1,min(nummat,int(td(1))))
        call ptiend(mr(np(32)),mr(np(33)),mr(np(78)),mr(np(79)),
     &              hr(np(43)),l2,nie,nen,nen1,ndm,numel)
      else

        if(pcomp(titl(16:19),'node',4)) then
          l1 = max(    1,int(td(1)))
          l2 = min(numnp,int(td(2)))
          l5    = 0
          td(2) = 0.0d0
          if(.not.gapfl) then
            gaptol = 1.d-3/sqrt(dble(max(1,numnp)))
          endif
          write(iow,2011) l1,l2
        elseif(pcomp(titl(16:19),'regi',4)) then
          l1 = 1
          l2 = numnp
          l3 = max(    0,int(td(1)))
          l4 = min(mxreg,int(td(2)))
          l5 =-1
          if(.not.gapfl) then
            gaptol = 1.d-3/sqrt(dble(max(1,numnp)))
          endif
          write(iow,2012) l3,l4
        elseif(pcomp(titl(16:19),'mate',4)) then
          l1 = 1
          l2 = numnp
          l3 = max(     1,int(td(1)))
          l4 = min(nummat,max(1,int(td(2))))
          l5 =-2
          if(.not.gapfl) then
            gaptol = 1.d-3/sqrt(dble(max(1,numnp)))
          endif
          write(iow,2013) l3,l4
        elseif(pcomp(titl(16:19),'coor',4)) then
          l1 = 1
          l2 = numnp
          l3 = 1
          l5 =-3
          do j = 4,1,-1
            td(j+1) = td(j)
          end do ! j
          if(.not.gapfl) then
            gaptol = 1.d-3/sqrt(dble(max(1,numnp)))
          endif
          write(iow,2021) gaptol,(td(j),j=2,ndm+1)
        elseif(pcomp(titl(16:19),'tol',3)  .or.
     &         pcomp(titl(16:19),'gap',3)) then
          if(td(1).eq.0.0d0) then
            gaptol = 1.d-3/sqrt(dble(max(1,numnp)))
            gapfl  = .false.
          else
            gaptol = td(1)
            gapfl  = .true.
          endif
          go to 1
        elseif(pcomp(titl(16:19),'face',4)) then
          l1 = nint(td(1))
          l2 = nint(td(2))
          l3 = max(0,min(1,nint(td(3))))
          write(iow,2022) wfac(l3+1),l1, wfac(l3+1),l2
          if(.not.gapfl) then
            gaptol = 1.d-3/sqrt(dble(max(1,numnp)))
          endif
          call ptiefac(hr(np(43)),mr(np(33)),mr(np(79)),
     &                 l1,l2,l3,gaptol)
          tfl = .true.
          go to 1
        else
          l1 = 1
          l2 = numnp
          l3 = 0
          l4 = 0
          l5 = td(1)
          if(.not.gapfl) then
            gaptol = 1.d-3/sqrt(dble(max(1,numnp)))
          endif
          if(l5.gt.0) then
            write(iow,2014) l5,td(2)
          else
            write(iow,2015)
          endif
        endif
        setvar = palloc(111,'TEMP1',numnp, 1)
        setvar = palloc(112,'TEMP2',numnp, 1)

        call tienod(mr(np(33)),hr(np(43)),mr(np(79)),mr(np(111)),
     &              mr(np(112)),mr(np(78)),ndm,nen,nen1,
     &              numnp,numel,l1,l2,l3,l4,l5,gaptol,td(2))

        setvar = palloc(112,'TEMP2',0, 1)
        setvar = palloc(111,'TEMP1',0, 1)
      endif
      tfl = .true.
      go to 1

c     Error treatments

900   call  errclr ('PCONTR')
      call plstop()

910   if(ior.eq.icf) then
        call pincld('end')
        incf = .false.
        intx = evint
        cprt = evint
        go to 1
      endif

911   call  endclr ('PCONTR',titl)
      call plstop()

c     Input/output formats

1001  format(a)

1002  format(a4,11x,a)

1003  format(a4,11x,15i15)

1004  format(a4,2x,a12,i8,2l5)

1005  format(4f20.0)

2000  format(1x,19a4,a3//5x,'Solution date: ',a//14x,a/14x,a/
     &                /5x,'Input Data Filename: ',a/
     &                /5x,'Number of Nodal Points  - - - - - - :',i9
     &                /5x,'Number of Elements  - - - - - - - - :',i9/
     &                /5x,'Spatial Dimension of Mesh - - - - - :',i9
     &                /5x,'Degrees-of-Freedom/Node (Maximum) - :',i9
     &                /5x,'Equations/Element       (Maximum) - :',i9
     &                /5x,'Number Element Nodes    (Maximum) - :',i9/
     &                /5x,'Number of Material Sets - - - - - - :',i9
     &                /5x,'Number Parameters/Set   (Program) - :',i9
     &                /5x,'Number Parameters/Set   (Users  ) - :',i9)

2001  format(/5x,'P a r t i t i o n',i4)

2002  format(/'   N o d a l   P a r t i o n   D a t a'/
     &        '      ndf       Partition '/ (i8,i12))

2003  format(/'   N o d a l   M a x i m u m   PDE   O r d e r'/
     &        '      ndf       Order'/ (i8,i12))

2004  format(/' *End of <FEAP> solution,  File: ',a/1x)

2005  format(/' *End of INCLUDE solution, File: ',a/1x)

2006  format(//'  **Variable Storage REAL**'///)

2007  format(//'  **Variable Storage COMPLEX**'///)

2008  format(/' ',a/)

2011  format(/5x,'Tie nodes from',i8,' to ',i8/1x)

2012  format(/5x,'Tie from region',i4,' to region',i4/1x)

2013  format(/5x,'Tie from material',i4,' to material',i4/1x)

2014  format(/5x,'Tie: direction =',i3,' X =',1p,1e12.5/1x)

2015  format(/5x,'Tie all nodes with common coordinates'/1x)

2016  format(/'   C o o r d i n a t e     L i n k s'/
     &        '      ndf        Link'/ (i8,i12))

2017  format(/'  Problem definitions are specified by include files.'
     &      //'  Output for each problem is written to separate files.'
     &      //'  Check file ',a,' for problem list and errors.')

2018  format(/'  ',33('-'),' START OF FEAP LOG ',34('-')/)

2019  format(/'  ',70('-'))

2020  format(/'  --> Problem',i4,': Output in file: ',a)

2021  format(/5x,'Tie all nodes with common coordinates'/
     &       10x,'Tol  = ',1p,1e12.4 /
     &       10x,'x_1  = ',1p,1e12.4:/
     &       10x,'x_2  = ',1p,1e12.4:/
     &       10x,'x_3  = ',1p,1e12.4:)

2022  format(/5x,'Tie 4-node elements with common faces'/
     &       10x,a,' 1 =',i4/ 10x,a,' 2 =',i4)

2023  format(//5x,'DOMAIN Information',//
     &        10x,'Number of partition nodes    : ',i7/
     &        10x,'Number of total problem nodes: ',i7,/
     &        10x,'Number of total problem eqns.: ',i7,/)

3000  format(/' *ERROR* PCONTR: Partition wrong: Npart =',i3)

3001  format(/' *ERROR* PCONTR: Attempt to solve problem before mesh'/
     &        '         input.  Check for error on FEAP or BINAry',
     &        ' record.'/1x)

3002  format(/' *ERROR* PCONTR: Can not do BATCH execution from this',
     &        ' mode.'/
     &        '         Do INTERACTIVE or put in INCLUDE file.'/1x)

3003  format(/' *ERROR* PCONTR: File name error')

3004  format(2i6,i5,1p,1e11.4,1p,5e10.2,    0p,1f10.2)
3005  format(2i6,i5,1p,1e11.4,1p,4e10.2,10x,0p,1f10.2)

3006  format(/'Total',i7,i5)

      end
