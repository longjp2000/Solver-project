c$Id: feap81.f,v 1.2 2006/12/05 22:42:53 rlt Exp $
      program feap

c-----[--.----+----.----+----.-----------------------------------------]

c      * * F E A P * * A Finite Element Analysis Program
c                        -      -       -        -
c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c     The University of California  does not guaranteed this program to
c     be error free. Every effort has been made to ensure proper coding
c     and declarations, but testing of all program options has not been
c     performed.

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1.  Version 8_1.a2                                  17/11/2006
c       2.  Version 8_1.a3                                  05/12/2006
c-----[--.----+----.----+----.-----------------------------------------]

c     Finite Element Analysis Program - (FEAP)  for solution of general
c     problem classes using the finite element method.     Problem size
c     is controlled by the dimension of blank common and value of  MAXM
c     as set in the parameter statement below.

c     Programmed by:

c                R. L. Taylor
c                Department of Civil and Environmental Engineering
c                University of California at Berkeley
c                Berkeley, California 94720-1710
c     E-mail:
c                rlt@ce.berkeley.edu
c     Web address:
c                www.ce.berkeley.edu/~rlt/feap
c-----[--.----+----.----+----.-----------------------------------------]

c     Notes:

c     1. Precision is controlled by ipr:

c        Set ipr = 1 for 8-byte integers; = 2 for 4-byte integers.

c     2. User written subprograms should include type specification
c        for all variables in each subprogram.

c        e.g.    real*8    a(12)
c                integer   name
c                character word*6
c                logical   flag
c                etc.

c     3. FEAP may create temporary input files during use.
c        Users should periodically check and delete files
c        which are no longer needed.  File names are normally
c        either the name of the data input file with an extender
c        or the name of the plot save file with an extender.

c     4. System dependent routines are included in:

c           a. File for subroutine doargs.

c           b. Plot routines include special characters which
c              must have case preserved (i.e., upper and lower
c              case letters in formats, etc.).


c     5. Input/Output is performed to files during execution of FEAP.
c        In general, the following files are used during executions:

c           a.  User logical unit numbers should be from 1 to 8.
c           b.  iot =  9 : Used for dynamic memory to expand arrays
c           c.  ilg = 10 : Used for write of log file.
c           d.  iop = 11 : Used for read/write delayed inputs.
c           e.  ios = 12 : Used for read/write scratch files.
c           f.  ird = 13 : Used to read results data from disk.
c           g.  iwd = 14 : Used to write results data to disk.
c           h.  ior = 15 : Use to read from the input data file.
c                          (specified when a problem is initiated).
c           i.  iow = 16 : Use to write output result data to file.
c                          (specified when a problem is initiated).
c           j.  lun = 17 : For PostScript file outputs.
c           k.  icl = 18+: Used for include file inputs.
c                          (additional include files may be opened).
c           l.  24+      : Used to save time history data.

c     End of Notes

c-----[--.----+----.----+----.-----------------------------------------]

c     Set variable types

      implicit none

      integer iii, lll

      include 'cdata.h'
      include 'codat.h'
      include 'hlpdat.h'
      include 'iodata.h'
      include 'iofile.h'
      include 'prmptd.h'
      include 'psize.h'
      include 'pathn.h'
      include 'setups.h'
      include 'vdata.h'
      
      iii=1

      do while( iii .ge. 0)
         lll = 0
      end do
      
c-----[--.----+----.----+----.-----------------------------------------]

c     Set version header for output to file and screen

      versn(1) = 'Release 8.1.a3'
      versn(2) = '05 December 2006'

c-----[--.----+----.----+----.-----------------------------------------]

c     Set ratio for real to integer variables: Set ipr = 1 or 2
                          ! ipr = 1 for equal  length real to integers
      ipr = 2             ! ipr = 2 for double length real to integers

c-----[--.----+----.----+----.-----------------------------------------]

c     Set path names to where manual files are stored (path name should
c     be 37 characters or less.  Feap will append a '/xxxx.t' where xxxx
c     is the manual page name to be displayed)

      file(1) = '/home/jp/proj/feap/manual/report/'
      file(2) = '/home/jp/proj/feap/manual/report/'
      file(3) = '/home/jp/proj/feap/manual/report/'
      file(4) = '/home/jp/proj/feap/manual/report/'
      file(5) = '/home/jp/proj/feap/manual/report/'

c-----[--.----+----.----+----.-----------------------------------------]

c     Set default logical unit numbers for files

      ilg = 10
      iop = 11
      ios = 12
      ird = 13
      iwd = 14
      ior = 15
      iow = 16
      lun = 17
      icl = 18

c-----[--.----+----.----+----.-----------------------------------------]

c     Set data input parsing flag

      coflg = .true.  ! Parse all input as expressions   (slower mode)
                      ! N.B. Use of 'parse' and 'noparse' mesh commands
                      !      permit change of this default setting.
      ciflg = .false. ! Get input file from menu window if .true.
                      ! or from keyboard entry if .false.
                      ! N.B. Used for windows version only

c-----[--.----+----.----+----.-----------------------------------------]

c     Set graphics default options

      defalt = .true. ! Graphics runs with default contour intervals, etc.
      prompt = .true. ! Prompt for graphics inputs when defalt = .false.

c-----[--.----+----.----+----.-----------------------------------------]

c     Set PostScript default mode

      pscolr = .true. ! PostScript outputs are in color
      psrevs = .false.! Color order is normal

c-----[--.----+----.----+----.-----------------------------------------]

c     Set help display level:

      hlplev = 0      ! Basic

c-----[--.----+----.----+----.-----------------------------------------]

c     Set increment for reducing array size

      incred = 2      ! No reduction unless array is less by 'incred'

c-----[--.----+----.----+----.-----------------------------------------]

c     Set solver flag: Program - solver = .true.; User - solver = .false.
c                      processor = Number of processors used by compiler
c                      soltyp    = Profile solver type: 1 = by 1 column
c                                                       2 = by 2 column
      processor = 1
      soltyp    = 2
      solver    = .true.

c-----START SOLUTION---------------------------------------------------]

c     Initialize solution system

      call pstart()

c     Display messages

      call pmessage()

c     Solve problem

      call pcontr()

c-----END SOLUTION-----------------------------------------------------]

c     N.B. No additional calls permitted after PCONTR

      end ! FEAP8_1

      subroutine pmessage()

      implicit   none

c     Place any messages for display here

      write(*,2000)

2000  format(/10x,'--> Please report errors by e-mail to:'/
     &        10x,'    feap-help@ce.berkeley.edu '/)

      end
