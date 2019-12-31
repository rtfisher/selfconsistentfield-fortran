      Subroutine setupbcphi (phi, dr, dz) 

       implicit none

*/////////////////////////////////////////////////////////////////////////
*/
*/  Subroutine establishes the the boundary conditions in the ghost region
*/  which is the order of the scheme thick.
*/
*/  Options include
*/
*/  bcs (_, _)  = 1	==>     periodic
*/  bcs (_, _)  = 2     ==>     inflow
*/  bcs (_, _)  = 3	==>     outflow 
*/  bcs (_, _)  = 4	==>     reflecting
*/
*/  Note that at reflecting bcs, the component of the velocity field 
*/  normal to the interface is reversed, while all other quantities,
*/  including the normal magnetic field, remain the same. This stems
*/  from the pseudo-vector nature of the B-field.
*/
*/  Altered : 10/19/96  -- Errors in inner and outer j outflow corrected
*/                      -- Reflecting bc indexing cleared up  R.F.
*/
*/  Altered : 10/26/98  -- Hardwired for reflecting bc,
*/                          and adjusted to allow only for
*/                          assignment of phi bcs. R.F.
*/
*/////////////////////////////////////////////////////////////////////////

#include "const.h"
#include "prob.h"

      REAL_T phi (0 : NR + 1, 0 : NZ + 1)
      integer xmin (2), xmax (2)

      REAL_T Gcons, mcloud, rhocloud
      REAL_T Jang, beta0, M, rhoc, rhoext, ciso, specmom

      common / block2 / Jang, beta0, M, rhoc, rhoext, ciso, specmom    

#ifdef AXI
      REAL_T dr, dz
#endif

      Integer i, j, order 

      xmin (1) = 1 
      xmax (1) = NR
      xmin (2) = 1 
      xmax (2) = NZ 

      order = 1

      Gcons = G
      rhocloud = RHOSPHERE
      mcloud = M

*..inner i boundary

*....reflecting

        Do j = xmin (2) - order, xmax (2) + order
          Do i = xmin (1) - order, xmin (1) - 1 
            phi (i, j) = phi (2 * xmin (1) - i - 1, j)
          End Do
        End Do              

*..outer i boundary

*.... inhomogeneous Dirichlet

        Do j = xmin (2) - order, xmax (2) + order
          Do i = xmax (1) + 1, xmax (1) + order 
             phi (i, j) = - Gcons
     & * mcloud /
     &         ( j**2 * dz**2 + ( (xmax (1) + 1) * dr)**2 )**(0.5)
          End Do
        End Do      

*..inner j boundary

*....reflecting 

        Do i = xmin (1) - order, xmax (1) + order
          Do j = xmin (2) - order, xmin (2) - 1 
            phi (i, j) = phi (i, 2 * xmin (2) - j - 1)
          End Do
        End Do         

*..outer j boundary

*.... inhomogeneous Dirichlet

        Do i = xmin (1) - order, xmax (1) + order
          Do j = xmax (2) + 1, xmax (2) + order
             phi (i, j) = - Gcons
     & * mcloud /
     &         ( i**2 * dr**2 + ( (xmax (2) + 1) * dz)**2 )**(0.5)
          End Do 
        End Do     
      End    
