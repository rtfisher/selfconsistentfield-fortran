        subroutine getphic (rho, phic)

C//////////////////////////////////////////////////////////////////////
C//
C// Subroutine determines the centripetal potential phic by using the
C// distribution of specific angular momentum (which is preserved in
C// Lagrangian coordinates) and the current mass
C// density distribution rho.
C//
C//////////////////////////////////////////////////////////////////////

#include "const.h"
#include "prob.h"

        implicit none

        REAL_T rho (0:NR+1, 0:NZ +1), phic (NR)
        REAL_T radius, sum, mass, M, Jang, jangd, rhoc, rhoext, ciso
        REAL_T dr, dz, specmom, beta0
        integer i

        common /block1 / dr, dz
        common /block2 / Jang, beta0, M, rhoc, rhoext, ciso, specmom

C...Update total angular momentum, holding beta_0 fixed, accounting
C...for changes of the cloud mass from last iteration.

        Jang = sqrt (12. / 25.) *
     &         (3. / (4. * pi ) )**(1. / 6.) *
     &         G**(1. / 2.) *
     &         M**(10. / 6.) / rhoext**(1. / 6.) *
     &         beta0**(1. / 2.)

        sum = 0.0

        i = 1

20      do i = 1, NR
          radius = i * dr

          call getstripmass (rho, radius, mass)

          sum = sum + ( (jangd (mass / M) )**2.0 / radius**3) * dr

c          phic (i) = (Jang / M)**2 * sum
          phic (i) = sum
        end do

        do i = 1, NR
          phic (i) = (Jang / M)**2 * phic (i)
c           phic (i) = specmom**2 * phic (i)
        end do

        end
