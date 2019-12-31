       subroutine getvphi (rho, vphi)

C///////////////////////////////////////////////////////////////////////
C//
C// Function returns the rotational velocity of the cloud (vs. r).
C//
C///////////////////////////////////////////////////////////////////////

#include "prob.h"

        implicit none

        REAL_T rho (0:NR+1, 0:NZ+1)
        REAL_T vphi (0 : NR + 1)
        REAL_T r, dr, dz, mass, jangd
        REAL_T Jang, beta0, M, rhoc, rhoext, ciso, specmom
        integer i

        common / block1 / dr, dz
        common / block2 / Jang, beta0, M, rhoc, rhoext, ciso, specmom

        do i = 1, NR

          r = i * dr

          call getstripmass (rho, r, mass)

          vphi (i) = (Jang / M) * jangd (mass / M) / r

       end do

       return
       end
