c$Id: dosets.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine dosets(inp,otp,res,sav,plt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Add initial character to default files using input
c               filename

c      Inputs:
c         inp  - Input   filename without inital character
c         otp  - Output  filename without inital character
c         res  - Restart filename without inital character
c         sav  - Save    filename without inital character
c         plt  - Plot    filename without inital character

c      Outputs:
c         inp  - Input   filename with inital character
c         otp  - Output  filename with inital character
c         res  - Restart filename with inital character
c         sav  - Save    filename with inital character
c         plt  - Plot    filename with inital character
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      character inp*12,otp*12,res*12,sav*12,plt*12

      save

      if(otp.eq.' ') then
        otp      = inp
        otp(1:1) = 'O'
      endif
      if(res.eq.' ') then
        res      = inp
        res(1:1) = 'R'
      endif
      if(sav.eq.' ') then
        sav      = inp
        sav(1:1) = 'R'
      endif
      if(plt.eq.' ') then
        plt = inp
        plt(1:1) = 'P'
      endif

      end
