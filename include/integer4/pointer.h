
      integer          num_nps       , num_ups
      parameter       (num_nps = 400 , num_ups = 200)

      integer          np         , up                             ! int4
      common /pointer/ np(num_nps), up(num_ups)
