        function getthermener (rho)

C//////////////////////////////////////////////////////////////////////
C//
C// Function returns the total thermal energy content of the cloud,
C// defined as all densities above the external ambient medium density.
C//
C//////////////////////////////////////////////////////////////////////

#include "const.h"
#include "prob.h"

        implicit none

        REAL_T getthermener, rho (0:NR+1, 0:NZ+1), mass
        REAL_T r, dr, dz
        REAL_T Jang, beta0, M, rhoc, rhoext, ciso, specmom

        integer i, j

        common / block1 / dr, dz
        common / block2 / Jang, beta0, M, rhoc, rhoext, ciso, specmom

        mass = 0.0

        do i = 1, NR

          r = i * dr

          do j = 1, NZ
            if (rho (i, j) .gt. rhoext) then

              mass = mass + 4.0 * pi *
     &          r * dr * dz * rho (i, j)

            end if

          end do
        end do

        getthermener = 3.0 / 2.0 * ciso**2 * mass

        return
        end
