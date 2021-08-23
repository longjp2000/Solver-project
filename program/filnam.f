c$Id: filnam.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine filnam

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set up files for input/output, etc.

c      Inputs:
c         None

c      Outputs:
c         File names returned in common /comfil/
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'chdata.h'
      include  'comfil.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'machnc.h'
      include  'prmptd.h'
      include  'vdata.h'

      logical   pcomp,lfil,linp,lout,lres,lsav,lplt
      character dinp*128,dout*128,dres*128,dplt*128, wd(2)*6,y*1
      character pinp*128,pout*128,pres*128,psav*128,pplt*128,filein*128
      integer   iorsav, nargs, ierr, iinp, ioup,ires,isav,iplt

      save

      data      dinp,dout,dres,dplt/'NONE','NONE','NONE','NONE'/
      data      wd/'New','Exists'/

c     Save input unit number set command line arguments to 0

      iorsav = ior
      nargs  = 0

c     Look to see if any problem has been run

101   inquire(file='feapname',exist=lfil)
      if(lfil) then
        open(ios,file='feapname',status='unknown')
        read(ios,2010,err=900,end=900) pinp,pout,pres,psav,pplt,epmac
        close(ios)
        finp = pinp
        fout = pout
        fres = pres
        fsav = psav
        fplt = pplt

c       Add extenders or get filenames from command line

        call filargs(nargs)

c       Compute solution timings

        if(nargs.gt.0) then
          call getcon(epmac)
        else
          ior = -abs(ior)
          write(*,2000) versn
          call uscreen(1,ior)
        endif
        go to 200
      else

c       Set default file names

        ior = -abs(ior)
        write(*,2000) versn
        call uscreen(1,ior)
        pinp = dinp
        pout = dout
        pres = dres
        psav = dres
        pplt = dplt

        call getcon(epmac)
      endif

c     Name file for input data

100   ierr = 1
      finp = pinp
1     write(*,2001)
      write(*,2002) pinp
      call pprint('                   Enter Name -->')
      read (*,1000,err=901,end=902) filein
      if (pcomp(filein,'stop',4)) then
       stop
c      call plstop()
      endif
      if (pcomp(filein,'exit',4)) call plstop()
      if (pcomp(filein,'quit',4)) call plstop()
      if (.not.pcomp(filein,' ',1)) finp = filein

c     Check if input files exists

10    inquire(file=finp,exist=linp)
      if(.not.linp) then
        write(*,3000)
        go to 1
      endif
      pinp = finp

c     Set default files for a filname beginning with 'I' or 'i'

      if(pcomp(pinp,'I',1)) then
        pinp(1:1) = 'I'
        pout      = pinp
        pout(1:1) = 'O'
        pres      = pinp
        pres(1:1) = 'R'
        psav      = pinp
        psav(1:1) = 'R'
      endif

c     Name file for output data

      ierr = 2
      fout = pout
2     write(*,2003) pout
      call pprint('                   Enter Name -->')
      read (*,1000,err=901,end=902) filein
      if (.not.pcomp(filein,' ',1)) fout = filein
20    pout = fout

c     Name file for restart read data

      ierr = 3
      fres = pres
3     write(*,2004) pres
      call pprint('                   Enter Name -->')
      read (*,1000,err=901,end=902) filein
      if (.not.pcomp(filein,' ',1)) fres = filein
30    pres = fres

c     Name file for restart save data

      ierr = 4
      fsav = psav
4     write(*,2005) psav
      call pprint('                   Enter Name -->')
      read (*,1000,err=901,end=902) filein
      if (.not.pcomp(filein,' ',1)) fsav = filein
40    psav = fsav

c     Set default plot window

      if(pcomp(pinp,'I',1)) then
        pplt      = pinp
        pplt(1:1) = 'P'
      endif

c     Name file for plot data

      ierr = 5
      fplt = pplt
5     write(*,2006) pplt
      call pprint('                   Enter Name -->')
      read (*,1000,err=901,end=902) filein
      if (.not.pcomp(filein,' ',1)) fplt = filein
50    pplt = fplt

c     Check file status and input if necessary

200   inquire(file=finp,exist=linp)
      if(.not.linp.and.nargs.gt.0) call plstop()
      if(.not.linp) go to 100
      iinp = 2
      inquire(file=fout,exist=lout)
      ioup = 1
      if(lout) ioup = 2
      inquire(file=fres,exist=lres)
      ires = 1
      if(lres) ires = 2
      inquire(file=fsav,exist=lsav)
      isav = 1
      if(lsav) isav = 2
      inquire(file=fplt,exist=lplt)
      iplt = 1
      if(lplt) iplt = 2
      if(nargs.gt.0) go to 300
      write(*,2007) wd(iinp),finp,wd(ioup),fout,wd(ires),fres,
     &              wd(isav),fsav,wd(iplt),fplt

c     Start execution

      if(.not.fileck) go to 300
      xxx = ' '
      write(xxx,2008)
      call pprint(xxx)
      ierr = 6
6     read(*,1001,err=901,end=901) y
      if(y.eq.'S' .or. y.eq.'s') call plstop()
      if(y.ne.'Y' .and. y.ne.'y') go to 100

c     Save a copy of current filenames

300   open(ios,file='feapname',status='unknown')
      rewind ios
      write(ios,2010) finp,fout,fres,fsav,fplt,epmac
      close(ios)
      write(*,2009)
      call uscreen(2,ior)

c     Reset ior unit number

      ior = iorsav
      return

c     Error in form of 'feapname' file

900   lfil = .false.
      write(*,'(/a)') ' ERROR IN feapname FILE -- Read filenames'
      close(ios,status = 'delete')
      go to 101

c     Error trap

901   write(*,3001)
      if(    ierr.eq.1) then
        go to 1
      elseif(ierr.eq.2) then
        go to 2
      elseif(ierr.eq.3) then
        go to 3
      elseif(ierr.eq.4) then
        go to 4
      elseif(ierr.eq.5) then
        go to 5
      elseif(ierr.eq.6) then
        go to 6
      else
        call plstop()
      endif

c     EOF encountered

902   call  endclr ('FILNAM',filein)
      if(    ierr.eq.1) then
        go to 10
      elseif(ierr.eq.2) then
        go to 20
      elseif(ierr.eq.3) then
        go to 30
      elseif(ierr.eq.4) then
        go to 40
      elseif(ierr.eq.5) then
        go to 50
      elseif(ierr.eq.6) then
        go to 6
      else
        call plstop()
      endif

c     Formats

1000  format(a)

1001  format(75a1)

2000  format(//
     & '    F I N I T E   E L E M E N T   A N A L Y S I S',
     & '   P R O G R A M'
     &  //11x,'FEAP (C) Regents of the University of California'
     &   /25x,'All Rights Reserved.'/23x,'VERSION: ',a/26x,'DATE: ',a)

2001  format(/8x,' I n p u t    F i l e n a m e s',
     & //8x,' Specify filenames:'/)
2002  format(11x,'Input   Data (default: ',a32,') :')
2003  format(11x,'Output  Data (default: ',a32,') :')
2004  format(11x,'Restart Read (default: ',a32,') :')
2005  format(11x,'Restart Save (default: ',a32,') :')
2006  format(11x,'Plot File    (default: ',a32,') :')

2007  format(/8x,' Files are set as:   Status    Filename'//
     &       11x,'Input   (read ) : ',a6,2x,a32/
     &       11x,'Output  (write) : ',a6,2x,a32/
     &       11x,'Restart (read ) : ',a6,2x,a32/
     &       11x,'Restart (write) : ',a6,2x,a32/
     &       11x,'Plots   (write) : ',a6,2x,a32//
     &  8x,' Caution, existing write files will be overwritten.'/)

2008  format( 8x,' Are filenames correct? ( y or n; s = stop) :')
2009  format(/8x,' R U N N I N G    F E A P    P R O B L E M    N O W')

2010  format(a/a/a/a/a/1p,1e25.15)

3000  format(/' *ERROR* FILNAM: Specified input file does not exist, ',
     &        'reinput name.'/)
3001  format(/' *ERROR* FILNAM: Reinput data')

      end
