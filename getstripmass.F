        subroutine getstripmass (rho, omega, mass)

C////////////////////////////////////////////////////////////////////////
C//
C// Subroutine integrates the mass along a strip in z, up to 
C// a value of radius omega. Used in angular momentum
C// distribution j (x), to evaluate dimensionless Lagrangian
C// coordinate x.
C//
C///////////////////////////////////////////////////////////////////////

#include "const.h"
#include "prob.h"

        implicit none

        REAL_T rho (0: NR+1, 0: NZ + 1), omega, mass, radius
        REAL_T dr, dz
        REAL_T Jang, beta0, M, rhoc, rhoext, ciso, specmom
        integer i, j

        common / block1 / dr, dz
        common / block2 / Jang, beta0, M, rhoc, rhoext, ciso, specmom

        mass = 0.0

        i = 1

 40     radius = i * dr

          do j = 1, NZ

C...Additional factor of 2 for mass below the plane.

            if (rho (i, j)  .gt. rhoext) then
              mass = mass + dz * 4.0 * pi
     &          * dr * radius * rho (i, j)
            end if

          end do

        if (radius .lt. omega) then
          i = i + 1

          if (i .gt. NR) then
            print *, 'getstripmass ERROR : i > NR!'
            goto 50
           end if

          goto 40
        end if

 50     end
