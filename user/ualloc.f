c$Id: ualloc.f,v 1.1 2006/11/20 20:33:27 rlt Exp $
      logical function ualloc(num,name,length,precis)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Define, delete, or resize a user dictionary entry.
c               Pointer defined for integer (single) and real
c               (double precision) arrays.

c               N.B. Currently limited to 200 names by dimension of
c               common blocks 'allotd','allotn','pointer'

c      Inputs:
c         num        - Entry number for array (see below)
c         name       - Name of array          (see below)
c         length     - Length of array defined: =0 for delete
c         precis     - Precision of array: 1 = integers; 2 = reals

c      Output:
c         up(num)    - Pointer to first word of array in blank common
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'allotd.h'

      logical   usetmem
      character name*(*)
      integer   i, num,length,precis

c     Storage definitions for UALLOC variables

      integer   list
      parameter (list = 1)

      character names(list)*5

      save

c     Define and store names

      data   (names(i),i=1,list)/

     &         'DUMMY'/

c     Short description of variables

c              'DUMMY',     !     1: Start here with user defined names

c     Do memory management operations

      ualloc = usetmem(list,names,num,name,length,precis)

      end
