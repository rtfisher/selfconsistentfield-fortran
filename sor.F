        subroutine sor (rho, phi, dr, dz)

*//////////////////////////////////////////////////////////////////////
*/
*/ fixed grid wth spacing h in two dimensions via Gauss-Seidel
*/  with red-black ordering and overrelaxation.
*/
*/ altered to handle axisymmetric grids and overrelaxation,
*/  as well as Chebyshev acceleration. 
*/                                       R.F. 10/24/98
*/
*/ Note that the Cartesian system has rhs = rho, while
*/  axisymmetric system has rhs = 4 pi G rho. 
*/
*//////////////////////////////////////////////////////////////////////

       implicit none

#include "const.h"
#include "prob.h" 

c#define CHEBYSHEV

      REAL_T rho (0 : NR + 1, 0 : NZ + 1),
     &   phi (0 : NR + 1 , 0 : NZ + 1)
      REAL_T sourcesum, residsum, resid0, 
     &        lamb, omega
      integer xmin (2), xmax (2)

      REAL_T rhojac

#ifdef AXI
      REAL_T a, b, r, dr, dz
#endif

#ifdef CARTESIAN
      REAL_T h, hi2
#endif
       integer i, j, niters

* Assume a fixed grid with spacing h in both directions.

#ifdef CARTESIAN
       h = delta (1)
       hi2 = 1. / h**2
       lamb = h**2 / 4.0
#endif

#ifdef AXI
       lamb = 1.0 / 4.0

C...define aspect ratio of grid for use below.

       a = dr / dz
#endif

C...maximum iteration value suitable for straight GSRB.

#define MAXITERS 25 * NR * NZ 
#define ACC 1.e-06

       xmin (1) = 1 
       xmin (2) = 1 
       xmax (1) = NR 
       xmax (2) = NZ

C...Estimate for spectral radius taken from NR. (p. 866, C edition).
C...Different bcs. Doesn't work for self-consistent problem.

#ifdef CHEBYSHEV
c       rhojac = cos (pi / NR) + a**2 * cos (pi / NZ) /
c     &         (1.0 + a**2)
#endif

c       omega = 2.0 / (1.0 + sqrt (1.0 - rhojac**2) ) 

C...SOR omega parameters for 20^2, 40^2, 100^2, and 200^2 grids.

#if (NR == 20 && NZ == 20) 
       omega = 1.75
#elif (NR == 40 && NZ == 40) 
       omega = 1.875 
#elif (NR == 100 && NZ == 100) 
       omega = 1.975
#elif (NR == 200 && NZ == 200)
       omega = 1.9875
#elif (NR == 50 && NZ == 10)
       omega = 1.990
#elif (NR == 100 && NZ == 20)
       omega = 1.990
#elif (NR == 200 && NZ == 40)
       omega = 1.999
#elif (NR == 250 && NZ == 50)
       omega = 1.9999
#elif
       ERROR : omega not defned for nr, nz.
#endif

C...For Chebysheyv acceleration, note that we use the empirical
C... optimal omega to estimate the spectral radius. Omega itself
C... gets reset below.

       rhojac = sqrt (1.0 - (2.0 / omega - 1.0)**2.0 )

#ifndef CHEBYSHEV
       lamb = lamb * omega
#endif

c       print *, 'omega : ', omega
 
       call setupbcphi (phi, dr, dz)

c       print *, 'SOR : Computing rhs norm ...'

c       do i = 0, NZ + 1
c         print *, 'SOR phi : ', i, ' ', phi (1, i)
c       end do

c       do i = 0, NR + 1
c         print *, 'SOR phi : ', i, ' ', phi (i, 1)
c       end do 

       sourcesum = 0.0

       do i = xmin (1), xmax (1)

         r = i * dr

         do j = xmin (2), xmax (2)

c           sourcesum = sourcesum + abs (rho (i, j) )

#ifdef AXI
            sourcesum = sourcesum + abs (4.0 * 
     &        pi * G 
     &        *  rho (i, j) ) * 2.0 * pi * r
     &        *  dr * dz
#endif

#ifdef CARTESIAN
            ERROR : Need to rewrite sourcesum for cartesian coords. 
#endif

         end do  ! j-loop
       end do    ! i-loop

c       print *, 'sourcesum 0 : ', sourcesum

c       print *, 'Entering relaxation loop ...'

#ifdef CHEBYSHEV

C...define initial omega parameter for CHEBYSHEV acceleration

       omega = 1.0
#endif

       do niters = 1, MAXITERS

* Compute residual to determine whether solution is sufficiently
* accurate.

#ifdef CHEBYSHEV
         lamb = 0.25 * omega
#endif

         residsum = 0.0

c         print *, 'iteration : ', niters

         do i = xmin (1), xmax (1)

           b = 1.0 / (2.0 * i)
           r = i * dr

           do j = xmin (2), xmax (2)

#ifdef AXI
             residsum = residsum + abs (1.0 / dr**2 * (
     &        -2.0 * phi (i, j)     * (1.0 + a**2) +
     &               phi (i + 1, j) * (1.0 + b) +
     &               phi (i - 1, j) * (1.0 - b) +
     &               phi (i, j + 1) * a**2 +
     &               phi (i, j - 1) * a**2 ) -
     &               4.0 * pi * G * rho (i, j) )
     &            * 2.0 * pi * r
     &            * dr * dz
#endif

#ifdef CARTESIAN
             residsum = residsum + abs (hi2 * (- 4.0 * phi (i, j) 
     &        + phi (i - 1, j) + phi (i, j - 1) + phi (i + 1, j)
     &        + phi (i, j + 1)) - rho (i, j) )
#endif 
           end do
         end do

         if (niters .eq. 1) then
           resid0 = residsum
         end if

c         print *, niters, ' ', log (residsum) / log (10.0)

c         if (residsum .lt. ACC * sourcesum) then
         if (residsum .lt. ACC * resid0) then  ! converged if a factor ACC less than resid0

#ifdef FILEOUT
!           open (unit = 3, file = 'selfcons.out')
!           write (3, 99) residsum / sourcesum
!           close (3)
#endif

#ifdef STDOUT
           print *, 'SOR : converged ', residsum / sourcesum
#endif

  99       format ('SOR :  converged. ', d22.16)

* Solution sufficiently accurate --> end. 

           goto 40 

         else

* Continue to relax.

           do i = xmin (1), xmax (1)

#ifdef AXI
C..define curvature factor for use below. Note that i != 0 must 
C..hold here to avoid coordinate singularity at r = 0. 

c             b = dr / (2.0 * i * dr)   ! b = dr / (2 r)
              b = 1.0 / (2.0 * i)
#endif

             do j = xmin (2), xmax (2)
               if (mod (niters, 2) .eq. 0) then   !red cycle    

                 if (mod (i + j, 2) .eq. 0) then
                   phi (i, j) = phi (i, j)
                 else

#ifdef AXI
                phi (i, j) = phi (i, j) + lamb * (
     &           2.0 * phi (i, j) * (-1.0 - a**2) +
     &           phi (i + 1, j) * (1.0 + b) +  
     &           phi (i - 1, j) * (1.0 - b) + 
     &           phi (i, j + 1) * a**2 +
     &           phi (i, j - 1) * a**2 -
     &           4.0 * pi * G *
     &                  rho (i, j) * dr**2 ) 
#endif
 
#ifdef CARTESIAN
                phi (i, j) = phi (i, j) + lamb * ( hi2 * (
     &           phi (i + 1,j) + phi (i, j + 1) + phi (i - 1, j) +  
     &           phi (i, j - 1) - 4.0 * phi (i, j) ) - rho (i, j) )
#endif

                 end if    

               else                                !black cycle

                 if (mod (i + j, 2) .eq. 0) then

#ifdef AXI
                phi (i, j) = phi (i, j) + lamb *  (
     &           2.0 * phi (i, j) * (-1.0 - a**2) +
     &           phi (i + 1, j) * (1.0 + b) +
     &           phi (i - 1, j) * (1.0 - b) +
     &           phi (i, j + 1) * a**2 +
     &           phi (i, j - 1) * a**2  -
     &           4.0 * pi * G *
     &                 rho (i, j) * dr**2 ) 
#endif

#ifdef CARTESIAN
                phi (i, j) = phi (i, j) + lamb * (hi2 * (  
     &           phi (i + 1, j) + phi (i, j + 1) + phi (i - 1, j) +
     &           phi (i, j - 1) - 4.0 * phi (i, j) ) - rho (i, j) )
#endif

                 else
                   phi (i, j) = phi (i, j)
                 end if
               end if

             end do ! end j-loop
           end do  ! end i-loop 

         end if

         call setupbcphi (phi, dr, dz) 

#ifdef CHEBYSHEV
         if (niters .eq. 1) then
           omega = 1.0 / (1.0 - rhojac**2 / 2.0) 
         else 
           omega = 1.0 / (1.0 - rhojac**2 * omega /
     &                4.0)              
         end if
#endif

       end do   ! niters loop

* Solution not accurate enough.
#ifdef FILEOUT
           open (unit = 3, file = 'selfcons.out', status = 'unknown')
           write (3, *) 'SOR ERROR : not converged.'
           close (3)
#endif

#ifdef STDOUT
         print *, 'SOR ERROR : not converged.'
#endif
         stop

c         do i = 1, NR
c           print *, i, ' ', phi (i, 1)
c         end do
  
 40    end

