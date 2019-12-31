        program selfcons 

#include "const.h"  
#include "prob.h"

#define MAXITERS 1000 
#define VERBOSE

        implicit none

C///////////////////////////////////////////////////////////////////////
C//
C// Self-consistent field program. Determines a self-consistent 
C// equilibrium solution to the axisymmetric rotating isothermal
C// body problem, subject to the external boundary condition
C// P_ext = const.
C//
C// The SCF method implemented here employs a 
C// successive overrelaxation Gauss-Seidel relaxation method to solve
C// the Poisson equation on a uniform cylindrical mesh with NR x NZ
C// cells in the radial and vertical directions, respectively. Vertical 
C// symmetry is assumed, so only the top half plane (z > 0) is modeled. 
C// The successive overrelaxation parameter omega (0 < omega < 2) has been
C// hand-tuned to optimal values on specific mesh sizes.
C//  
C// For papers describing the self-consistent field method employed here,
C// see Stahler, 1983 (doi:10.1086/160941) and
C// Ostriker & Mark, 1968 (doi:10.1086/149506).
C//
C///////////////////////////////////////////////////////////////////////

        REAL_T rho (0:NR+1, 0: NZ + 1),
     &          phig (0 : NR + 1, 0 : NZ + 1),
     &          phic (0 : NR + 1), totalmass,
     &          oldrho (NR, NZ), maxerr, dr, dz, r
     &          , vphi (0 : NR + 1)
        REAL_T Jang, M, beta0, rhoc, rhoext, ciso, beta, specmom
        REAL_T jstar, mstar, rstar, radiuseq, height, volume
        REAL_T getrotener, getbindener, getthermener, minrho, maxrho
        REAL_T getvolume, virial, T, U, V, W, pext

#ifdef COMMANDLINE
        character*8 arg
#endif

        integer i, j, iter, pass

#ifdef COMMANDLINE
        integer nargs, iargc
#endif
        common / block1 / dr, dz
        common / block2 / Jang, beta0, M, rhoc, rhoext, ciso, specmom
        common / block3 / radiuseq, height

        pass = 0

#ifdef COMMANDLINE
C...Read in arguments from command line.

      nargs = iargc ()

      if (nargs .eq. 0) then

        print *, '====================================================='
        print *, 'Rotating Isothermal Equilibrium Generator'
        print *, '     -- written Robert Fisher, 11/1/98'
        print *, '====================================================='
        print *, 'usage : isorotgen nr nz contrast beta0 relaxpar mass'
        print *
        print *, 'Where'
        print *, '   nr = number of cells in the radial direction'
        print *, '   nz = number of cells in the z direction'
        print *, '   contast = ratio of central to external density'
        print *, '   beta0 = beta_0 parameter'
        print *, '   relaxpar = relaxation parameter (0 - 1)'
        print *, '   mass = initial guess for dimensionless mass'

        stop

      else

        print *, 'Argument calls not implemented.'
        stop

        call getarg (1, arg)
        read(arg,*,iostat=stat)  nr
 
        call getarg (2, arg)
        read(arg,*,iostat=stat)  nz

        call getarg (3, arg)
        read(arg,*,iostat=stat) contrast

        call getarg (4, arg)
        read(arg,*,iostat=stat)  beta0

        call getarg (5, arg)
        read(arg,*,iostat=stat)  relaxpar 

        call getarg (6, arg)
        read(arg,*,iostat=stat)  mass 
      end if
#endif

#ifdef SWEEP 
C...Initialize program with initial density distribution and
C...grid properties on first time through (pass = 0). Once we have
C...calculated the first equilibrium structure, use that as
C...an initial guess for succesive sweeps (pass > 0).

 999    if (pass .eq. 0) then
#endif
          call init (rho, phig)

#ifdef SWEEP
        else

C...Set parameter values for successive sweeps.
C...Keep rhoext / rhoc fixed here, but vary beta0.

          if (pass .eq. 1) then
            beta0 = .01
          else if (pass .eq. 2) then
            beta0 = .02
          else if (pass .eq. 3) then
            beta0 = .03
          else if (pass .eq. 4) then
            beta0 = .04
          else if (pass .eq. 5) then
            beta0 = .05
          else if (pass .eq. 6) then
            beta0 = .06
          else if (pass .eq. 7) then
            beta0 = .07
          else if (pass .eq. 8) then
            beta0 = .08
          else if (pass .eq. 9) then
            beta0 = .09
          else if (pass .eq. 10) then
            beta0 = .10
          else if (pass .eq. 11) then
            beta0 = .11
          else if (pass .eq. 12) then
            beta0 = .12  
          end if       
#endif
          Jang = sqrt (12. / 25.) *
     &         (3. / (4. * pi ) )**(1. / 6.) *
     &         G**(1. / 2.) *
     &         M**(10. / 6.) / rhoext**(1. / 6.) *
     &         beta0**(1. / 2.)

C...Rescale grid with previous model's size distribution,
C...to keep it from going off the grid.

c        dr = RDOMAIN * rstar / float (NR - 1)
c        dz = ZDOMAIN * rstar / float (NZ - 1)


#ifdef SWEEP
        end if
#endif

        beta = 25.0 / 12.0 *
     &       (4.0 * pi / 3.0)**(1.0/3.0) *
     &       rhoext**(1.0 / 3.0) *
     &       Jang**2.0 / (G *
     &       M**(10.0 / 3.0) )

        mstar = ciso**3.0 / (sqrt (rhoext * G**3.0) )
        rstar = ciso / (sqrt (G * rhoext) )
        jstar = ciso**6.0 / (rhoext * G**2.0)
        pext = rhoext * ciso**2

        totalmass = 0.0

        do i = 1, NR

          r = i * dr

          do j = 1, NZ
             if (rho (i, j) .gt. rhoext) then

C...Additional factor of 2 for mass below the plane.

               totalmass = totalmass + 4.0 * pi *
     &                     r * dr * dz *
     &                     rho (i, j)
             end if
          end do
        end do


        do i = 0, NR + 1
          phic (i) = 0.0
        end do

c///////////////////////////////////////////////////////////////////
c//
c// Begin main iteration loop.
c//
c///////////////////////////////////////////////////////////////////

        do iter = 1, MAXITERS 

c... Main SCF update steps.

          call getphic (rho, phic)        ! calculate phic from rho
          call sor (rho, phig, dr, dz)    ! calculate phig from rho
          call evalrho (phig, phic, rho) ! update rho 
                                          ! from phig and phic

c... Update cloud mass.
          call getcloudmass (rho, M)

#ifdef STDOUT 
          print *, 'Iteration # ', iter , ' ', maxerr
          print *, 'mass / mstar : ', M / mstar
#endif

c... Calculate the maximum fractional error in density over this iteration.

          maxerr = 0.0

          do i = 1, NR
            do j = 1, NZ
              maxerr = max(maxerr, abs(oldrho (i, j)
     &         -rho (i,j) )/
     &         (oldrho (i, j) + rho (i, j) ) )
              oldrho (i, j) = rho (i, j)
            end do
          end do

C...Check to see if result is converged to desired accuracy.

          if (maxerr .lt. 1.e-04) then
            goto 100
          end if

C...Compute virial error as a check on solution.

        T = getrotener (rho)
        U = getthermener (rho)
        V = getvolume (rho, rhoext, dr, dz)
        W = getbindener (rho, phig)

        virial = 2.0 * U + 2.0 * T - 3.0 * pext * V - W
      
#ifdef STDOUT 
        print *, 'pressure term = ', 3.0 * pext * V / (2.0 * W)
        print *, 'T / W = ', T / W
        print *, 'U / W = ', U / W
        print *, 'log10 virial error = ', iter, ' ', log (abs
     &    (virial / (2.0 * W) ) ) / log (10.0)
#endif    

#ifdef FILEOUT   

!        if (pass .eq. 1) then
!          open (unit = 3, file = 'selfcons.out', status = 'new')
!        else
!          open (unit = 3, file = 'selfcons.out', status = 'old')
!        endif

        open (unit = 3, file = 'selfcons.out')

        write (3, 99) iter, maxerr, M / mstar, 
     &    3.0 * pext * V / (2.0 * W), T / W, V / W,
     &    abs (virial / (2.0 * W) ),
     &    Jang / M / (jstar / mstar)

        close (unit = 99)
        call flush (3)
#endif
       
#ifdef STDOUT 
        print *, iter, ' ', maxerr, ' ', M / mstar, ' ',
     &    3.0 * pext * V / (2.0 * W), ' ', T / W, ' ', U / W,
     &    ' ', abs (virial / (2.0 * W) ), ' ',
     &    Jang / M / (jstar / mstar)
#endif

        end do   ! main iteration loop

#ifdef FILEOUT
        write (3, *) 'Solution not converged.'  ! get here only if not
        write (3, *) '======================='  ! converged.
        stop
#endif

 99     format (i5, ' ', d22.15, ' ', d22.15, ' ', d22.15, ' ',
     &    d22.15, ' ', d22.15, ' ', d22.15, ' ', d22.15)

C...Compute and display properties of the solution.

 100    totalmass = 0.0
        volume = 0.0

        maxrho = 0.0
        minrho = 1.0 / smallnum

        do i = 1, NR
          r = i * dr
          do j = 1, NZ

C...Additional factor of 2 for mass below the plane.

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

#ifdef FILEOUT
        open (unit = 3,file='selfcons.out', status = 'old')
#endif

        T = getrotener (rho)
        U = getthermener (rho)
        V = getvolume (rho, rhoext, dr, dz)
        W = getbindener (rho, phig)
     
#ifdef FILEOUT 
        write (3, *) 'Solution converged.'
        write (3, *) '==================='
        write (3, 102) beta
        write (3, 113) ciso * Jang / (G * M**2.)
        write (3, 103) rhoc / rhoext
        write (3, 104) maxrho / minrho
        write (3, 105) totalmass / mstar
        write (3, *) '--------------------------------------'
        write (3, 106) U / W
        write (3, 107) T / W
        write (3, 108) 3.0 * pext * V / (2.0 * W)
        write (3, 115) totalmass / volume / (rhoext)
        write (3, *) '--------------------------------------'
        write (3, 110) rstar
        write (3, 111) mstar
        write (3, 112) jstar
        write (3, *) '--------------------------------------'
        write (3, 114) totalmass
        write (3, *) '--------------------------------------'
#endif
 
 102    format ('beta_0 = ', d22.16)
 103    format ('density contrast rhoc / rhoext = ', d22.16)
 104    format ('actual density contrast = ', d22.16)
 105    format ('m / mstar = ', d22.16)
 106    format ('alpha = ', d22.16)
 107    format ('beta = ', d22.16)
 108    format ('3 P V / (2 W) = ', d22.16)
c 109    format (1x)
 110    format ('rstar = ', d22.16)
 111    format ('mstar = ', d22.16)
 112    format ('jstar = ', d22.16)  
 113    format ('cs J / (G M^2) = ', d22.16)
 114    format ('physical mass (gm) = ', d22.16)
 115    format ('average density / external density = ', d22.16)

#ifdef STDOUT
        print *, 'Solution converged.'
        print *, '==================='
        print *, 'beta_0 = ', beta
        print *, 'cs J / (G M^2)', ciso * Jang /
     &                              (G * M**2.)
        print *, 'density contast rhoc / rhoext = ', rhoc / rhoext
        print *, 'actual contrast = ', maxrho / minrho
        print *, 'm / mstar = ', totalmass / mstar
        print *
        print *, 'alpha  = ', U / W
        print *, 'beta = ', T / W
        print *, 'surface term = ', 3.0 * pext * V / (2 * W)
        print *, 'length scale rstar = ',
     &     ciso / sqrt (G * rhoext)
        print *, 'mass scale mstar = ', ciso**3.0 /
     &     sqrt (rhoext * G**3.0)
        print *
        print *, 'cloud mass = ', totalmass
        print * 
        print *, 'rhobar / rhoext = ', (totalmass / volume)  /
     &     rhoext
#endif

C.......find equatorial radius and height along z-axis.

         i = 1 

 101     if ( (rho (i, 1) .gt. rhoext) .and.
     &        (i .lt. NR) ) then
           i = i + 1
           goto 101 
         else if (i .ge. NR) then
#ifdef FILEOUT
           write (3, *) 'Cloud off grid in radial direction.'
#endif

#ifdef STDOUT
           print *, 'Cloud off grid in radial direction.'
#endif

           stop

         end if

         radiuseq = i * dr 

#ifdef FILEOUT
         write (3, 116) radiuseq / rstar
         write (3, 117) radiuseq / rstar * mstar / totalmass
#endif

#ifdef STDOUT
          print *, 'radius_equator / rstar = ', radiuseq / rstar
          print *, 'r_e / rstar / (m / mstar) = ',
     &     radiuseq / rstar * mstar / totalmass
#endif

 116     format ('radius_equator / rstar = ', d22.16)
 117     format ('r_e / rstar / (m / mstar) = ', d22.16) 

c.........find height along z-axis.

         i = 1
  
 201     if ( (rho (1, i) .gt. rhoext) .and. 
     &       (i .lt. NZ) ) then
           i = i + 1 
           goto 201
         else if (i .ge. NZ) then
#ifdef FILEOUT
           write (3, *) 'Cloud off grid in z-direction.'
#endif

#ifdef STDOUT
           print *, 'Cloud off grid in z-direction.'
#endif
           stop
         end if

         height = i * dz

#ifdef FILEOUT
         write (3, 118) height / rstar
         write (3, 119) radiuseq / height
#endif

#ifdef STDOUT 
         print *, 'height / rstar = ', height / rstar
         print *, 'radiuseq / height = ', radiuseq / height
#endif

 118     format ('h / rstar = ', d22.16)
 119     format ('r_e / h = f = ', d22.16) 

#ifdef FILEOUT 
         close (unit = 3)
         call flush (3)
#endif
 
#ifdef SWEEP
         if (pass .lt. 13) then
           pass = pass + 1

	   call getvphi (rho, vphi)
           call dump (rho, phig, phic, vphi)

           goto 999
         end if
#endif

#ifndef SWEEP
         call dump (rho, phig, phic, vphi)
#endif

c... Optionally display cuts of vphi and density along r and z.

#ifdef PLOT
         call getvphi (rho, vphi)

         print *, 'vphi'
  
         do i = 1, NR
           print *, i * dr / radiuseq, ' ', vphi (i) / ciso
         end do

          print *, 'rhoz'

         do j = 1, NZ
           print *, log (j * dz / rstar) / log (10.0),
     &        ' ', log (rho (1, j) / rhoext) /
     &         log (10.0)
         end do

         print *, 'rhor'

         do i = 1, NR
           print *, log (i * dr / rstar) / log (10.0),
     &          ' ', log (rho (i, 1) / rhoext) /
     &          log (10.0)
         end do
#endif
        end
