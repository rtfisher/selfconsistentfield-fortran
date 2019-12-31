c...following declares real type.
c.....note that one might have to redefine function 
c.....calls as well to get full precision.
c.....double precision values are required for
c.....typical stellar quantities in cgs

#define REAL_T real*8

c...following specifies geometry for sor routine

#define AXI

c...following program flow and I/O options

c#define COMMANDLINE
#define SWEEP 
c#define STDOUT
#define FILEOUT
c#define PLOT

c...following defines grid resolution in radial
c.....and z directions (assuming axisymmetry)

#define NR 100 
#define NZ 100 

c...following two parameters specify model

#define CONTRAST 20.0
#define BETA 0.0

c...following two parameters can be varied
c.....to assist convergence; relax specifies
c.....an underrelaxation parameter used in
c.....the update; mass is the initial mass
c.....of the uniform cloud guessed at the 
c.....start of a model sweep in units of
c.....mstar

#define RELAX 0.125 
#define MASS 1.0 

c...following specifies central density for
c.....physical model in cgs units; note, however,
c.....that the output is scaled to dimensionless
c.....units, so the choice is somewhat arbitrary

#define RHOSPHERE 5.0e-18

c...problem domain size in terms of RSTAR; 
c...... 1.5 in each direction seems to work 
c...... reasonably well for low-beta0 models

#define RDOMAIN 1.5 
#define ZDOMAIN 1.5 
