///* rkf45d.c
//   Sample driver for the ODE integrator rkf45(). */
//
//#include <math.h>
//#include <stdio.h>
//#include "cmath.h"
//
//main ()
//{
//
//int     f();
//double  h, relerr, abserr, x1, x2;
//int     n, flag, nfe, maxfe, fail, step;
//double  y[2], yp[2];
//
//printf ("\n\n  --- CMATH --- Design Software 1989\n");
//printf ("\nExercise ODE solver rkf45().\n\n");
//
//n      = 2;
//flag   = 1;
//maxfe  = 5000;
//relerr = 1.0e-5;
//abserr = 1.0e-5;
//rkfinit (n, &fail);
//printf ("%s\n\n", cmathmsg(RKFINIT_C, fail));
//
//if (fail == 0)
//   {
//   y[0]   = 1.0;
//   y[1]   = 0.0;
//
//   printf ("       x      y[0]       y[1]\n");
//   printf ("----------------------------------\n");
//
//   for (step = 1; step <= 5; ++step)
//      {
//      x2 = 0.4 * step;
//      x1 = x2 - 0.4;
//      rkf45 (f, n, y, yp, &x1, x2, &relerr, abserr,
//             &h, &nfe, maxfe, &flag);
//      printf ("%10.6f %10.6f %10.6f \n", x1, y[0], y[1]);
//      if (flag != 2)
//         {
//         printf ("%s\n", cmathmsg(RKF45_C, flag));
//         break;
//         }
//      }
//
//   rkfend ();
//   printf ("\n%s\n", cmathmsg(RKF45_C, flag));
//   printf ("nfe         : %d\n", nfe);
//   printf ("step size   : %f\n", h);
//   printf ("\ncorrect answer : y[0] = 0.270671,  y[1] = -0.135335\n");
//   printf ("If rkf45() does not reach the final time, use ...\n");
//   printf ("y[0] = 2 * exp(-x),  y[1] = -exp(-x) \n");
//   }
//
//return (0);
//}
//
//
//
//int f (n, x, y, dydx)
//int    n;
//double x, y[], dydx[];
//{
//dydx[0] = 998.0 * y[0] + 1998.0 * y[1];
//dydx[1] = -999.0 * y[0] - 1999.0 * y[1];
//return (0);
//}
//
