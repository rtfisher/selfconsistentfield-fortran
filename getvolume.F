       function getvolume (rho, rhoext, dr, dz)

C///////////////////////////////////////////////////////////////////////
C//
C// Function returns the volume of the cloud.
C//
C///////////////////////////////////////////////////////////////////////

#include "const.h"
#include "prob.h"

        implicit none

        REAL_T getvolume, volume, rho (0:NR+1, 0:NZ+1)
        REAL_T r, dr, dz
        REAL_T Jang, beta0, M, rhoc, rhoext, ciso, specmom
        integer i, j

!        common / block1 / dr, dz
!        common / block2 / Jang, beta0, M, rhoc, rhoext, ciso, specmom

        volume = 0.0

        do i = 1, NR

          r = i * dr

          do j = 1, NZ

            if (rho (i, j) .gt. rhoext) then

              volume = volume + 4.0 * pi *
     &          r * dr * dz

            end if

         end do
       end do

       getvolume =  volume

       return
       end


