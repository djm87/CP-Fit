500   2   0.001   298.         nsteps  ictrl  eqincr  temp
* boundary conditions
    0       1       1           iudot    |    flag for vel.grad.
    1       1       1                    |    (0:unknown-1:known)
    1       1       0                    |
                                         |
   0.0005    0.      0.           udot   |    vel.grad
    0.    -0.001   0.                    |
    0.      0.    0.0005                 |
                                         |
    1       0       0           iscau    |    flag for Cauchy
            0       0                    |
                    1                    |
                                         |
    0.      0.      0.          scauchy  |    Cauchy stress
            0.      0.                   |
                    0.                   |
