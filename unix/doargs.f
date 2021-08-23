c$Id: doargs.f,v 1.1 2006/11/20 20:33:21 rlt Exp $
      subroutine doargs(inp,outp,res,sav,plt,narg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Get names of files from command line input

c      Inputs:
c         None

c      Outputs:
c         inp   - Name of input file
c         outp  - Name of output file
c         res   - Name of restart read file
c         sav   - Name of restart save file
c         plt   - Name of plot file
c         narg  - Number of arguments read
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      character*12 inp, outp, res, sav, plt, argv*14
      integer      i, narg, nchars
      integer      iargc, lnblnk

      save

      narg = iargc()

      if(narg.gt.0) then

c       Set files to blank

        inp = ' '
        outp= ' '
        plt = ' '
        res = ' '
        sav = ' '

c       Check arguments set on command line

        do i = 1, narg

          call getarg(i,argv)
          nchars = lnblnk(argv)

c         Device specification

          if (argv(1:1) .eq. '-') then

c           Input file specification

            if      (argv(2:2).eq.'i') then
              inp = argv(3:nchars)

c           Output file specification

            else if (argv(2:2).eq.'o') then
              outp = argv(3:nchars)

c           Plot file specification

            else if (argv(2:2).eq.'p') then
              plt = argv(3:nchars)

c           Restart read file specification

            else if (argv(2:2).eq.'r') then
              res = argv(3:nchars)

c           Restart save file specification

            else if (argv(2:2).eq.'s') then
              sav = argv(3:nchars)

c           Error on command line

            else
              write( *, 2000) argv(2:nchars)
              call plstop()
            endif
          else
            write( *, 2001)  argv(1:nchars)
            call plstop()
          endif
        end do ! i

c     Check that files are correct if narg > 0

        if(inp.ne.' ') then
        call dosets(inp,outp,res,sav,plt)
        else
          call plstop()
        endif
      endif

2000  format( ' Unknown command line option   -> ',a12)
2001  format( ' Unknown command line argument -> ',a12)

      end
