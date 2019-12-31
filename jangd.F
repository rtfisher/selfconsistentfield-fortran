        function jangd (x)

C//////////////////////////////////////////////////////////////////////
C//
C// Specific angular momentum distribution as a function of 
C// x = mass coordinate m / total mass M.
C// Chosen as Maclaurin-type, appropriate for a rigidly rotating
C// sphere.
C// 
C//////////////////////////////////////////////////////////////////////

#include "prob.h"

        implicit none

        REAL_T jangd, x

        if (x .gt. 1.0) then


c#ifdef VERBOSE
          print *, 'jangd : x > 1! ', x
c#endif

          jangd = 0.0

        else

          jangd = 5.0 / 2.0 * (1.0 - (1.0 - x)**(2.0 / 3.0) )

C...rigidly rotating disk
c          jangd = 2.0 * x 

        end if


        return
        end

