// Retrieved from:
// http://ergodic.ugr.es/cphys/UTILES/aleatorios.h 

// Header file for using dranxor2C.f through C.
// The description of the functions is stored in dranxorInfo.txt

# pragma once

// Initializes the random number generator.
extern void dranini_(int *);

// Returns a random number uniformly distributed in (0,1)
extern double dranu_(void);

// Write infort.entero the state of the random number generator in order to continue a run
extern void dranwrite_(int *);

// Read from fort.entero the state of the random number generator in order to continue a run
extern void dranread_(int *);

// Devuelve entero aleatorio en (1,n)
extern int idran_(int *);

// Devuelve un gaussiano de media 0 y varianza 1, obtenido con el algoritmo de "numerical inversion method"
extern double drang_(void);

// Devuelve en x1 y x2 dos gaussianos de media 0 y varianza 1 obtenido con el m�odo de Box-Mller-Wiener
extern double dranbm_(double *x1,double *x2);

// Idem que  dranwrite y read pero ademas recibe el nombre del fichero. Mejor que las anteriores.
extern void dranwrite2_(int *,char *);
extern void dranread2_(int *,char *);
