       subroutine dump (rho, phig, phic, vphi)

C///////////////////////////////////////////////////////////////////////
C//
C// Subroutine dumps solution to various files.
C//
C///////////////////////////////////////////////////////////////////////

#include "const.h"
#include "prob.h"

        implicit none

        REAL_T rho (0 : NR + 1, 0 : NZ + 1),
     &        phig (0 : NR + 1, 0 : NZ + 1),
     &        phic (0 : NR + 1),
     &        vphi (0 : NR + 1)

        REAL_T Jang, beta0, M, rhoc, rhoext, ciso, specmom
        REAL_T minrho, maxrho, radiuseq, height
        REAL_T beta, volume, rstar, mstar, jstar
        REAL_T getrotener, getbindener, getvolume, getthermener
        REAL_T r, dr, dz
        REAL_T T, U, V, W, virial, pext, remainder, totalmass

        integer i, j, val, d1, d2, d3

        character*1 string1, string2, string3, string4,
     &      string5, string6, string7, string8

        character*26 filen

        common / block1 / dr, dz
        common / block2 / Jang, beta0, M, rhoc, rhoext, ciso, specmom
        common / block3 / radiuseq, height

 100    totalmass = 0.0
        volume = 0.0

        maxrho = 0.0
        minrho = 1.0 / smallnum

        do i = 1, NR

          r = i * dr

          do j = 1, NZ

C...Additional factor of 2 for mass below the xy plane.

             if (rho (i, j) .gt. rhoext) then
               totalmass = totalmass + 4.0 * pi *
     &                     r * dr * dz *
     &                     rho (i, j)
               volume = volume + 4.0 * pi *
     &                     r * dr * dz
             end if

             maxrho = max (maxrho, rho (i, j) )
             minrho = min (minrho, rho (i, j) )

          end do
        end do

         beta = 25.0 / 12.0 *
     &       (4.0 * pi / 3.0)**(1.0/3.0) *
     &       rhoext**(1.0 / 3.0) *
     &       Jang**2.0 / (G *
     &       totalmass**(10.0 / 3.0) )

         rstar = ciso / (sqrt (G * rhoext) )
         mstar = ciso**3.0 / (sqrt (rhoext * G**3.0) )
         jstar = ciso**6.0 / (rhoext * G**2.0)
         pext = rhoext * ciso**2

        T = getrotener (rho)
        U = getthermener (rho)
        V = getvolume (rho, rhoext, dr, dz)
        W = getbindener (rho, phig)

        virial = 2.0 * U + 2.0 * T - 3.0 * pext * V - W
        virial = abs (virial / (2.0 * W) )

C.......find equatorial radius and height along z-axis.


         i = 1

 101     if ( (rho (i, 1) .gt. rhoext) .and.
     &        (i .lt. NR) ) then
           i = i + 1
           goto 101
         else if (i .ge. NR) then
           print *, 'Cloud off grid in radial direction.'
           stop
         end if

         radiuseq = i * dr

c.........find height along z-axis.                                       

         i = 1

 201     if ( (rho (1, i) .gt. rhoext) .and.
     &       (i .lt. NZ) ) then
           i = i + 1
           goto 201
         else if (i .ge. NZ) then
           print *, 'Cloud off grid in z-direction.'
           stop
         end if

         height = i * dz

C...convert the numerical value of the density contast to a
c...string for use in the file name

        val = (rhoc / rhoext)

        remainder = (rhoc / rhoext) - val

        if (remainder .ge. 0.5) then
          val = val + 1
        end if     

c        val = round (val)

        d1 = val / 10
        d2 = (val - 10 * d1)

        string1 = CHAR (d1 + 48)
        string2 = CHAR (d2 + 48)

C...convert the numerical value of beta_0 to a string for
c...use in the file name

        val = beta0 * 1000.0

        d1 = val / 100
        d2 = (val - 100 * d1) / 10
        d3 = (val - 100 * d1 - 10 * d2)

        string3 = CHAR (d1 + 48)
        string4 = CHAR (d2 + 48)
        string5 = CHAR (d3 + 48)

        val = NR

        d1 = val / 100
        d2 = (val - 100 * d1) / 10
        d3 = (val - 100 * d1 - 10 * d2)

        string6 = CHAR (d1 + 48)
        string7 = CHAR (d2 + 48)
        string8 = CHAR (d3 + 48)

c...write a log of parameters used in run

        filen = 'isorot_'//string1//string2//'_'//
     &      string3//string4//string5//'_'//
     &      string6//string7//string8//
     &       '.log'

        open (6, file = filen, status = 'new')

C...Include the most important parameters in the first line of
C...the file, so that they may be easily read back in later.


        write (6, 121) NR, NZ, (NR * dr / rstar),
     &   (NZ * dz / rstar), totalmass / mstar
        write (6, 109) 'Solution converged.'
        write (6, 109) '==================='
        write (6, 102) beta
        write (6, 113) ciso * Jang / (G * M**2.)
        write (6, 103) rhoc / rhoext
        write (6, 104) maxrho / minrho
        write (6, 105) totalmass / mstar
        write (6, 109) '-------------------'
        write (6, 120) log (virial) / log (10.0)
        write (6, 109) '-------------------'
        write (6, 106) U / W
        write (6, 107) T / W
        write (6, 108) 3.0 * pext * V / (2.0 * W)
        write (6, 115) totalmass / volume / (rhoext)
        write (6, 109) '-------------------'
        write (6, 110) rstar
        write (6, 111) mstar
        write (6, 112) jstar
        write (6, 109) '-------------------'
        write (6, 114) totalmass
        write (6, 109) '-------------------'
        write (6, 116) radiuseq / rstar
        write (6, 117) radiuseq / rstar * mstar / totalmass
        write (6, 118) height / rstar
        write (6, 119) radiuseq / height

        close (6)

 118     format ('h / rstar = ', d22.16)
 119     format ('r_e / h = f = ', d22.16)
 116     format ('radius_equator / rstar = ', d22.16)
 117     format ('r_e / rstar / (m / mstar) = ', d22.16)

 120    format ('log virial error = ', d22.16)
 102    format ('beta_0 = ', d22.16)
 103    format ('density contrast rhoc / rhoext = ', d22.16)
 104    format ('actual density contrast = ', d22.16)
 105    format ('m / mstar = ', d22.16)
 106    format ('alpha = ', d22.16)
 107    format ('beta = ', d22.16)
 108    format ('3 P V / (2 W) = ', d22.16)
 109    format (A19)
 110    format ('rstar = ', d22.16)
 111    format ('mstar = ', d22.16)
 112    format ('jstar = ', d22.16)
 113    format ('cs J / (G M^2) = ', d22.16)
 114    format ('physical mass (gm) = ', d22.16)
 115    format ('average density / external density = ', d22.16)             

 121    format (i5, 2x, i5, 2x, d22.16, 2x, d22.16, 2x, d22.16)

c...write density in three-coiumn ASCII form

        filen = 'isorot_'//string1//string2//'_'//
     &      string3//string4//string5//'_'//
     &      string6//string7//string8//
     &      '.rho.dat'

        open (6, file = filen, status = 'new')

        do i = 1, NR
          do j = 1, NZ
            write (6, 11) i, j, rho (i, j) / rhoext
          end do
        end do

 11     format (i5, ' ', i5, ' ', d22.16)

        close (6)

c...write potential in three-column ASCII form  

        filen = 'isorot_'//string1//string2//'_'//
     &      string3//string4//string5//'_'//
     &      string6//string7//string8//
     &      '.phig.dat'

        open (6, file = filen, status = 'new')

        do i = 1, NR
          do j = 1, NZ
            write (6, 11) i, j, phig (i, j) /
     &        (G * mstar / rstar)
          end do
        end do

        close (6)

        call getvphi (rho, vphi)

c...write vphi (r) in two-column ASCII form

        filen = 'isorot_'//string1//string2//'_'//
     &      string3//string4//string5//'_'//
     &      string6//string7//string8//
     &      '.vphi.dat'

        open (6, file = filen, status = 'new')

        do i = 1, NR
          write (6, 12) i, vphi (i) / ciso
        end do

  12    format (i5, ' ', d22.16)

        close (6)

c...write 1-D cut of density along r in two-column ASCII form

        filen = 'isorot_'//string1//string2//'_'//
     &      string3//string4//string5//'_'//
     &      string6//string7//string8//
     &      '.rhor.dat'

        open (6, file = filen, status = 'new')

        do i = 1, NR
          write (6, 12) i, rho (i, 1) / rhoext
        end do

        close (6)

c...write 1-D cut of density along z in two-column ASCII form

        filen = 'isorot_'//string1//string2//'_'//
     &      string3//string4//string5//'_'//
     &      string6//string7//string8//
     &      '.rhoz.dat'

        open (6, file = filen, status = 'new')

        do i = 1, NZ
          write (6, 12) i, rho (1, i) / rhoext
        end do

        close (6)

        end
