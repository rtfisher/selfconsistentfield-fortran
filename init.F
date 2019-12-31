        subroutine init (rho, phig)

C///////////////////////////////////////////////////////////////////////
C//
C// Subroutine to initialize program.
C//
C///////////////////////////////////////////////////////////////////////

#include "const.h"
#include "prob.h"

        implicit none

        REAL_T rho (0 : NR + 1, 0 : NZ + 1),
     &           phig (0 : NR + 1, 0 : NZ + 1)
        REAL_T  dr, dz, radius, totmass, Jang, M, rhoc, rhoext, ciso
        REAL_T specmom, mstar, rstar, rsphere, beta0, r
        integer i, j

        common / block1 / dr, dz
        common / block2 / Jang, beta0, M, rhoc, rhoext, ciso, specmom

        ciso = 1.8e4

        rhoc = RHOSPHERE
        rhoext = (RHOSPHERE) / CONTRAST

        mstar = ciso**3.0 / (sqrt (rhoext * G**3.0) )
        rstar = ciso / (sqrt (G * rhoext) )

C..As a starting guess, assume a uniform cloud of 1 mstar and density
c..rhoc.

        rsphere = (3.0 * MASS * mstar / (4.0 * pi *
     &    rhoc) )**(1.0 / 3.0)

C..Choose a problem domain size, scaled to rstar in each
C...coordinate direction.

c        dr = max ((RDOMAIN * (rsphere) ), 1.0 * rstar)
c     &        / float (NR - 1)
c        dz = max ((ZDOMAIN * (rsphere) ), 1.0 * rstar)
c     &        / float (NZ - 1)

        dr = RDOMAIN * rstar / float (NR - 1)
        dz = ZDOMAIN * rstar / float (NZ - 1)

        totmass = 0.0

        do i = 0, NR + 1

          radius = i * dr

          do j = 0, NZ + 1

            r = sqrt (i**2 * dr**2 + j**2 * dz**2)

C...Additional factor of 2 for mass below the plane.

            if (r**2 .lt. rsphere**2) then
              rho (i, j) = rhoc
              totmass = totmass + rho (i, j)
     &         * 4.0 * pi * radius * dr * dz
            else
              rho (i, j) = rhoext
c               rho (i, j) = smallnum
            end if

          end do
        end do

        M = totmass

        do i= 0, NR + 1
          do j = 0, NZ + 1
            phig (i, j) = 0.0
          end do
        end do

c        Jang = 0.0
        beta0 = BETA

        Jang = sqrt (12. / 25.) *
     &         (3. / (4. * pi ) )**(1. / 6.) *
     &         G**(1. / 2.) *
     &         M**(10. / 6.) / rhoext**(1. / 6.) *
     &         beta0**(1. / 2.)

c       Jang = totmass * ciso * CHI * (rsphere) 

C...Factor of half here because we do not take into account the 
C...material beneath the plane z = 0 as well. Since the problem
C...domain includes only that portion above z = 0, we have only
C...half of that mass, i.e. M = totmass and not M = 2.0 totmass,
C...which is the actual mass of the cloud.

c        specmom = Jang / M

#ifdef VERBOSE
c        print *, 'INIT  dr : ', dr
c        print *, 'INIT  dz : ', dz
c        print *, 'INIT  Mass : ', M
#endif

        end
