        function getbindener (rho, phig)

C//////////////////////////////////////////////////////////////////////
C//
C// Function returns the gravitational binding energy content of the cloud,
C// defined as all densities above the external ambient medium density.
C//
C//////////////////////////////////////////////////////////////////////

#include "const.h"
#include "prob.h"

        implicit none

        REAL_T getbindener, energy, rho (0:NR+1, 0:NZ+1)
        REAL_T phig (0:NR+1, 0:NZ+1), gradphir, gradphiz
        REAL_T r, z, dr, dz
        REAL_T Jang, beta0, M, rhoc, rhoext, ciso, specmom

        integer i, j

        common / block1 / dr, dz
        common / block2 / Jang, beta0, M, rhoc, rhoext, ciso, specmom

        energy = 0.0

        do i = 1, NR

          r = i * dr

          do j = 1, NZ

            z = j * dz

            if (rho (i, j) .gt. rhoext) then

              gradphir = (phig (i, j) - phig (i - 1, j) ) / dr
              gradphiz = (phig (i, j) - phig (i, j - 1) ) / dz

              energy = energy + 4.0 * pi *
     &          r * dr * dz * rho (i, j) * (r * gradphir +
     &          z * gradphiz)

            end if

          end do
        end do

        getbindener = energy

        return
        end
