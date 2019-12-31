        subroutine evalrho (phig, phic, rho)

C//////////////////////////////////////////////////////////////////////
C//
C// Subroutine determines the mass density distribution implied
C// by the potentials phig and phic.
C//
C//////////////////////////////////////////////////////////////////////

#include "prob.h"

        implicit none

        REAL_T rho (0: NR + 1, 0: NZ + 1),
     &      phig (0 : NR + 1, 0 : NZ + 1),
     &      phic (NR), rhoc, rhoext,
     &      phig0, Jang, beta0, M, specmom, update
        REAL_T ciso

        common / block2 / Jang, beta0, M, rhoc, rhoext, ciso, specmom

        integer i, j

        phig0 = 0.0

        do i = 1, NR

          do j = 1, NZ

c            rho (i, j) = rhoc * exp ( (-phig (i, j) + phic (i)
c     &       + phig (1, 1) ) / ciso**2)

            update = rhoc * exp ( (-phig (i, j) + phic (i)
     &       + phig (1, 1) ) / ciso**2)

C...Choose a relaxation scheme to improve convergence. Without relaxation,
C...solution sometimes oscillates between two extreme states, and never
C...settles down.

            rho (i, j) = (1. - RELAX) * rho (i, j) + RELAX * update

C...Apply outer boundary condition.

            if (rho (i, j) .lt. rhoext ) then
              rho (i, j) = rhoext
            end if

          end do
        end do

c        do i = 1, NR
c        print *, i, ' ', (-phig (i, 1)+ phic(i)+ phig (1, 1) ) /
c     &         ciso**2
c        end do

        end
