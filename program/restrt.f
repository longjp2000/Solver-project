c$Id: restrt.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine restrt(fres,ndm,ndf,nneq,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Read/save restart files for resolutions

c      Inputs:
c         fres    - Name of restart file to read/save
c         ndm     - Spatial dimension of mesh
c         ndf     - Number dof/node
c         nneq    - Total dumber of parameters in solutions
c         isw     - Switch: = 1 for read; =2 for save.

c      Outputs:
c         u(*)    - Solution state read
c         none    - from/to blank common
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'arclel.h'
      include  'arcler.h'
      include  'cdata.h'
      include  'counts.h'
      include  'ddata.h'
      include  'dyndat.h'
      include  'evdata.h'
      include  'fdata.h'
      include  'gltran.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'ndata.h'
      include  'part7.h'
      include  'pointer.h'
      include  'print.h'
      include  'prlod.h'
      include  'rdata.h'
      include  'tdata.h'
      include  'tdato.h'
      include  'umac1.h'
      include  'comblk.h'

      include  'p_point.h'

      logical   exst,sfl,fl9,setvar,ralloc,walloc
      character fres*17,yorn*1,lct*15
      integer   i,ndm,ndf,nneq,isw, nnpo,nnlo,nnmo,ndmo,ndfo,nh2
      real*8    ctl(3)

      save

      data      ctl /3*0.0d0/

c     Check file status

1     inquire(file=fres,exist=exst)
      if(.not.exst.and.isw.eq.1) then
        write(iow,3002) fres
        if(ior.lt.0) then
          write(*,3002) fres
          call pprint(
     &       '           Specify new name for restart file? (y or n) >')
          read (*,1000) yorn
          if(yorn.eq.'y' .or. yorn.eq.'Y') then
            call pprint('           New Restart File Name >')
            read (*,1000) fres
            goto  1
          endif
        endif
        return
      endif

c     Open file

      if(exst) then
        open (ios,file=fres,form='unformatted',status='old')
      else
        open (ios,file=fres,form='unformatted',status='new')
      endif
      rewind ios

c     Read restart files

      if(isw.eq.1) then

c       Control information

        read(ios) nnpo,nnlo,nnmo,ndmo,ndfo,fl(9)
        if((nnpo.eq.numnp).and.(nnlo.eq.numel).and.(nnmo.eq.nummat)
     &          .and.(ndmo.eq.ndm).and.(ndfo.eq.ndf)) then

c         Solution parameters

          read(ios) theta,nrk,nrc,nrm,nrt,noi,numint,alpha,gtan,
     &              nstep,niter,naugm,titer,taugm,iaugm,iform,
     &              ttim,dt,dtold,rnmax,prop,rlnew,c0,cs01,cs02,
     &              ds0,r,det0,xn,fl9,mf,mq

c         Input displacement solution

          setvar = ralloc( 40,'U    ',ios)
          write(iow,2000) 'I n p u t',nstep,ttim,dt,'input'
          if(ior.lt.0) then
            write(*,2000) 'I n p u t',nstep,ttim,dt,'input'
          endif

c         Eigenpairs

          if(mq.gt.0) then
            setvar = ralloc( 76,'EVAL ',ios)
            setvar = ralloc( 77,'EVEC ',ios)
            write(iow,2001) 'input',mf,mq
            if(ior.lt.0) then
              write(*,2001) 'input',mf,mq
            endif
          endif

c         Transient data

          write(iow,2002) prop,rlnew
          if(ior.lt.0) then
            write(*,2002) prop,rlnew
          endif
          if(fl9) then
            setvar = ralloc( 42,'VEL  ',ios)
            write(iow,2003) 'input',noi
            if(ior.lt.0) then
              write(*,2003) 'input',noi
            endif
          endif

c         Current load state

          setvar = ralloc( 30,'FTN  ',ios)
          write(iow,2004) 'input'
          if(ior.lt.0) then
            write(*,2004) 'input'
          endif

c         History data

          setvar = ralloc( 49,'H    ',ios)
          write(iow,2005) 'input'
          if(ior.lt.0) then
            write(*,2005) 'input'
          endif
          refl   = .true.
          fl(11) = .false.

c         Contact data

          call contact (306)

c         User data

          urest = 1
          lct   = 'restart'
          do i = 1,10
            call umshlib(i,prt)
          end do ! i
          do i = 1,10
            call umaclib(i,lct,ctl,prt)
          end do ! i
          urest = 0

        else
          write(iow,3001)
          if(ior.lt.0) then
            write(*,3001)
          endif
        endif

c     Save information for restart

      elseif(isw.eq.2) then

c       Control information

        write(ios) numnp,numel,nummat,ndm,ndf,fl(9)

c       Solution parameters

        fl9 = flp(9,1) .or. flp(9,2) .or. flp(9,3) .or. flp(9,4)
     &                 .or. fl(9)
        write(ios) theta,nrk,nrc,nrm,nrt,noi,numint,alpha,gtan,
     &             nstep,niter,naugm,titer,taugm,iaugm,iform,
     &             ttim,dt,dtold,rnmax,prop,rlnew,c0,cs01,cs02,
     &             ds0,r,det0,xn,fl9,mf,mq

c       Solution state

        setvar = walloc( 40,'U    ',nneq*3,2, ios)
        write(iow,2000) 'O u t p u t',nstep,ttim,dt,'output'
        if(ior.lt.0) then
          write(*,2000) 'O u t p u t',nstep,ttim,dt,'output'
        endif

c       Eigenpairs

        if(mq.gt.0) then
          setvar = walloc( 76,'EVAL ', mq    , 2, ios)
          setvar = walloc( 77,'EVEC ', mq*neq, 2, ios)
          write(iow,2001) 'output',mf,mq
          if(ior.lt.0) then
            write(*,2001) 'output',mf,mq
          endif
        endif

c       Transient data

        if(fl9) then
          setvar = walloc( 42,'VEL  ',nrt*nneq,2, ios)
          write(iow,2003) 'output',noi
          if(ior.lt.0) then
            write(*,2003) 'output',noi
          endif
        endif

c       Current load state

        setvar = walloc( 30,'FTN  ', 4*nneq, 2, ios)
        write(iow,2004) 'output'
        if(ior.lt.0) then
          write(*,2004) 'output'
        endif

c       History data

        call pgetd('H   ',point,nh2, i,sfl)
        if(.not.sfl) then
          nh2 = 1
        endif
        setvar = walloc( 49,'H    ', nh2, 2, ios)
        write(iow,2005) 'output'
        if(ior.lt.0) then
          write(*,2005) 'output'
        endif

c       Contact data

        call contact (307)

c       User data

        urest = 2
        lct   = 'restart'
        do i = 1,10
          call umshlib(i,prt)
        end do ! i
        do i = 1,10
          call umaclib(i,lct,ctl,prt)
        end do ! i
        urest = 0

      endif

c     Close file

      close(ios)

c     Formats

1000  format(a)

2000  Format('   R e s t a r t   ',a,'   D a t a'/
     &       10x,'Time step number  =',i8/
     &       10x,'Time at restart   =',1p,1e12.5/
     &       10x,'Time increment    =',1p,1e12.5/
     &       10x,'Displacements ',a)

2001  Format(10x,'Eigenpairs ',a,' for',i4,' modes',i4,' total values')

2002  Format(10x,'Proportional load =',1p,1e12.5/
     &       10x,'Arc-length   load =',1p,1e12.5)

2003  Format(10x,'Transient states ',a,' (noi =',i2,')')

2004  Format(10x,'Force vector ',a)

2005  Format(10x,'History data ',a)

3001  format(' *ERROR* Incorrect information in a restart')

3002  format(' *ERROR* Restart file ',a17,' does not exist')

      end
