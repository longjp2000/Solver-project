
      integer          num_nps       , num_ups
      parameter       (num_nps = 400 , num_ups = 200)

      integer*8        np         , up                             ! int8
      common /pointer/ np(num_nps), up(num_ups)
