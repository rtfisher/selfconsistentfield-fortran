        function getrotener (rho)

C///////////////////////////////////////////////////////////////////////
C//
C// Function returns the rotational energy of the cloud.
C//
C///////////////////////////////////////////////////////////////////////  

#include "const.h"
#include "prob.h"

        implicit none

        REAL_T getrotener, energy, rho (0:NR+1, 0:NZ+1)
        REAL_T r, dr, dz, omega, mass, jangd, dm
        REAL_T Jang, beta0, M, rhoc, rhoext, ciso, specmom
        integer i, j

        common / block1 / dr, dz
        common / block2 / Jang, beta0, M, rhoc, rhoext, ciso, specmom

        energy = 0.0

        do i = 1, NR 

          r = i * dr

          do j = 1, NZ

            if (rho (i, j) .gt. rhoext) then

              call getstripmass (rho, r, mass)

              dm = 4.0 * pi * r * dr * dz * rho (i, j)

              omega = (Jang / M) * jangd (mass / M) / r**2

              energy = energy + 4.0 * pi *
     &          r * dr * dz * rho (i, j) *
     &          (0.5 * r**2 * omega**2)

            end if

         end do

       end do

       getrotener = energy

       return
       end
