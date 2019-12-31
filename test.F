       program test

C...Mass contained within a cylindrical mass shell of radius r, 
C...determined analytically. Use as check.

       implicit none

#include "const.h"

       real*8 mass, r, rcloud, dr, rho
       integer i

       rcloud = 1.0e17
       dr = rcloud / 100.0
       rho = 5.0e-18

       do i = 1, 100 

         r = i * dr
         mass = 2.0 / 3.0 * pi * rcloud**3.0  * rho *
     &    (1.0 - (1.0 - (r / rcloud)**2.0 )**(1.5) )
 
         print *, r, ' ', mass

       end do

       end 
