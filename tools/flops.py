from numpy import arange

#!/usr/local/bin/python2.5

#/*****************************/
#/*          flops.c          */
#/* Version 2.0,  18 Dec 1992 */
#/*         Al Aburto         */
#/*      aburto@nosc.mil      */
#/*****************************/

__author__ = 'Brian Olson'

"""
/*
   Flops.c is a 'c' program which attempts to estimate your systems
   floating-point 'MFLOPS' rating for the FADD, FSUB, FMUL, and FDIV
   operations based on specific 'instruction mixes' (discussed below).
   The program provides an estimate of PEAK MFLOPS performance by making
   maximal use of register variables with minimal interaction with main
   memory. The execution loops are all small so that they will fit in
   any cache. Flops.c can be used along with Linpack and the Livermore
   kernels (which exersize memory much more extensively) to gain further
   insight into the limits of system performance. The flops.c execution
   modules also include various percent weightings of FDIV's (from 0% to
   25% FDIV's) so that the range of performance can be obtained when
   using FDIV's. FDIV's, being computationally more intensive than
   FADD's or FMUL's, can impact performance considerably on some systems.
   
   Flops.c consists of 8 independent modules (routines) which, except for
   module 2, conduct numerical integration of various functions. Module
   2, estimates the value of pi based upon the Maclaurin series expansion
   of atan(1). MFLOPS ratings are provided for each module, but the
   programs overall results are summerized by the MFLOPS(1), MFLOPS(2),
   MFLOPS(3), and MFLOPS(4) outputs.

   The MFLOPS(1) result is identical to the result provided by all
   previous versions of flops.c. It is based only upon the results from
   modules 2 and 3. Two problems surfaced in using MFLOPS(1). First, it
   was difficult to completely 'vectorize' the result due to the 
   recurrence of the 's' variable in module 2. This problem is addressed
   in the MFLOPS(2) result which does not use module 2, but maintains
   nearly the same weighting of FDIV's (9.2%) as in MFLOPS(1) (9.6%).
   The second problem with MFLOPS(1) centers around the percentage of
   FDIV's (9.6%) which was viewed as too high for an important class of
   problems. This concern is addressed in the MFLOPS(3) result where NO
   FDIV's are conducted at all. 
   
   The number of floating-point instructions per iteration (loop) is
   given below for each module executed:

   MODULE   FADD   FSUB   FMUL   FDIV   TOTAL  Comment
     1        7      0      6      1      14   7.1%  FDIV's
     2        3      2      1      1       7   difficult to vectorize.
     3        6      2      9      0      17   0.0%  FDIV's
     4        7      0      8      0      15   0.0%  FDIV's
     5       13      0     15      1      29   3.4%  FDIV's
     6       13      0     16      0      29   0.0%  FDIV's
     7        3      3      3      3      12   25.0% FDIV's
     8       13      0     17      0      30   0.0%  FDIV's
   
   A*2+3     21     12     14      5      52   A=5, MFLOPS(1), Same as
	   40.4%  23.1%  26.9%  9.6%          previous versions of the
						flops.c program. Includes
						only Modules 2 and 3, does
						9.6% FDIV's, and is not
						easily vectorizable.
   
   1+3+4     58     14     66     14     152   A=4, MFLOPS(2), New output
   +5+6+    38.2%  9.2%   43.4%  9.2%          does not include Module 2,
   A*7                                         but does 9.2% FDIV's.
   
   1+3+4     62      5     74      5     146   A=0, MFLOPS(3), New output
   +5+6+    42.9%  3.4%   50.7%  3.4%          does not include Module 2,
   7+8                                         but does 3.4% FDIV's.

   3+4+6     39      2     50      0      91   A=0, MFLOPS(4), New output
   +8       42.9%  2.2%   54.9%  0.0%          does not include Module 2,
						and does NO FDIV's.

   NOTE: Various timer routines are included as indicated below. The
	timer routines, with some comments, are attached at the end 
	of the main program.

   NOTE: Please do not remove any of the printouts.

   EXAMPLE COMPILATION:
   UNIX based systems
       cc -DUNIX -O flops.c -o flops
       cc -DUNIX -DROPT flops.c -o flops 
       cc -DUNIX -fast -O4 flops.c -o flops 
       .
       .
       .
     etc.

   Al Aburto
   aburto@nosc.mil
*/
"""

import math
import time

def dtime(p):
   q = p[2]
   p[2] = time.time()
   p[1] = p[2] - q

def flops(quick=True):
   TimeArray = [0.0, 0.0, 0.0]
#double nulltime, TimeArray[3];   /* Variables needed for 'dtime()'.     */
#double TLimit;                   /* Threshold to determine Number of    */
#				 /* Loops to run. Fixed at 15.0 seconds.*/

   T = [0.0 for x in arange(36)]
#double T[36];                    /* Global Array used to hold timing    */
#				 /* results and other information.      */

#double sa,sb,sc,sd,one,two,three;
#double four,five,piref,piprg;
#double scale,pierr;

   A0 = 1.0
   A1 = -0.1666666666671334
   A2 = 0.833333333809067E-2
   A3 = 0.198412715551283E-3
   A4 = 0.27557589750762E-5
   A5 = 0.2507059876207E-7
   A6 = 0.164105986683E-9

   B0 = 1.0
   B1 = -0.4999999999982
   B2 = 0.4166666664651E-1
   B3 = -0.1388888805755E-2
   B4 = 0.24801428034E-4
   B5 = -0.2754213324E-6
   B6 = 0.20189405E-8

   C0 = 1.0
   C1 = 0.99999999668
   C2 = 0.49999995173
   C3 = 0.16666704243
   C4 = 0.4166685027E-1
   C5 = 0.832672635E-2
   C6 = 0.140836136E-2
   C7 = 0.17358267E-3
   C8 = 0.3931683E-4

   D1 = 0.3999999946405E-1
   D2 = 0.96E-3
   D3 = 0.1233153E-5

   E2 = 0.48E-3
   E3 = 0.411051E-6

   print("\n")
   print("   FLOPS Python Program (Double Precision), V2.0 18 Dec 1992\n\n")

# Initial number of loops. Original code claims this is a magic number.
   loops = 15625

#/****************************************************/
#/* Set Variable Values.                             */
#/* T[1] references all timing results relative to   */
#/* one million loops.                               */
#/*                                                  */
#/* The program will execute from 31250 to 512000000 */
#/* loops based on a runtime of Module 1 of at least */
#/* TLimit = 15.0 seconds. That is, a runtime of 15  */
#/* seconds for Module 1 is used to determine the    */
#/* number of loops to execute.                      */
#/*                                                  */
#/* No more than NLimit = 512000000 loops are allowed*/
#/****************************************************/

   T[1] = 1.0E+06/loops

   TLimit = 15.0
   NLimit = 512000000

   piref = 3.14159265358979324
   one   = 1.0
   two   = 2.0
   three = 3.0
   four  = 4.0
   five  = 5.0
   scale = one

   print("   Module     Error        RunTime      MFLOPS\n")
   print("                            (usec)\n")
#/*************************/
#/* Initialize the timer. */
#/*************************/
   
   dtime(TimeArray)
   dtime(TimeArray)
   
#/*******************************************************/
#/* Module 1.  Calculate integral of df(x)/f(x) defined */
#/*            below.  Result is ln(f(1)). There are 14 */
#/*            double precision operations per loop     */
#/*            ( 7 +, 0 -, 6 *, 1 / ) that are included */
#/*            in the timing.                           */
#/*            50.0% +, 00.0% -, 42.9% *, and 07.1% /   */
#/*******************************************************/
   n = loops
   sa = 0.0

   while sa < TLimit:
      n = 2 * n
      x = one / n    #                        /*********************/
      s = 0.0                #                        /*  Loop 1.          */
      v = 0.0                #                        /*********************/
      w = one 

      dtime(TimeArray)
      for i in arange(1,n):
         v = v + w
         u = v * x
         s = s + (D1+u*(D2+u*D3))/(w+u*(D1+u*(E2+u*E3)))
      dtime(TimeArray)
      sa = TimeArray[1]

      if n == NLimit:
         break
      #/* printf(" %10ld  %12.5lf\n",n,sa); */

   scale = 1.0E+06 / n
   T[1]  = scale

#/****************************************/
#/* Estimate nulltime ('for' loop time). */
#/****************************************/
   dtime(TimeArray)
   for i in arange(1, n):
      pass
   dtime(TimeArray)
   nulltime = T[1] * TimeArray[1]
   if nulltime < 0.0:
      nulltime = 0.0

   T[2] = T[1] * sa - nulltime

   sa = (D1+D2+D3)/(one+D1+E2+E3)
   sb = D1

   T[3] = T[2] / 14.0#                             /*********************/
   sa = x * ( sa + sb + two * s ) / two#           /* Module 1 Results. */
   sb = one / sa#                                  /*********************/
   n  = int( ( 40000 * sb ) / scale )
   sc = sb - 25.2
   T[4] = one / T[3]
#						   /********************/
#						   /*  DO NOT REMOVE   */
#						   /*  THIS PRINTOUT!  */
#						   /********************/
#//   printf("     1   %13.4le  %10.4lf  %10.4lf\n",sc,T[2],T[4])
   print("     1   {:13.4e}  {:10.4f}  {:10.4f}\n".format(sc,T[2],T[4]))

   m = n

#/*******************************************************/
#/* Module 2.  Calculate value of PI from Taylor Series */
#/*            expansion of atan(1.0).  There are 7     */
#/*            double precision operations per loop     */
#/*            ( 3 +, 2 -, 1 *, 1 / ) that are included */
#/*            in the timing.                           */
#/*            42.9% +, 28.6% -, 14.3% *, and 14.3% /   */
#/*******************************************************/

   s  = -five#                                      /********************/
   sa = -one#                                       /* Loop 2.          */
#						   /********************/
   dtime(TimeArray)
   for i in arange(1, m+1):
      s  = -s
      sa = sa + s
   dtime(TimeArray)
   T[5] = T[1] * TimeArray[1]
   if T[5] < 0.0:
      T[5] = 0.0

   sc   = m

   u = sa#                                         /*********************/
   v = 0.0#                                        /* Loop 3.           */
   w = 0.0#                                        /*********************/
   x = 0.0

   dtime(TimeArray)
   for i in arange(1, m+1):
      s  = -s
      sa = sa + s
      u  = u + two
      x  = x +(s - u)
      v  = v - s * u
      w  = w + s / u
   dtime(TimeArray)
   T[6] = T[1] * TimeArray[1]

   T[7] = ( T[6] - T[5] ) / 7.0#                   /*********************/
   m  = int( sa * x  / sc )#                    /*  PI Results       */
   sa = four * w / five#                           /*********************/
   sb = sa + five / v
   sc = 31.25
   piprg = sb - sc / (v * v * v)
   pierr = piprg - piref
   T[8]  = one  / T[7]
#						  /*********************/
#						  /*   DO NOT REMOVE   */
#						  /*   THIS PRINTOUT!  */
#						  /*********************/
   print("     2   {:13.4e}  {:10.4f}  {:10.4f}\n".format(pierr,T[6]-T[5],T[8]))

#/*******************************************************/
#/* Module 3.  Calculate integral of sin(x) from 0.0 to */
#/*            PI/3.0 using Trapazoidal Method. Result  */
#/*            is 0.5. There are 17 double precision    */
#/*            operations per loop (6 +, 2 -, 9 *, 0 /) */
#/*            included in the timing.                  */
#/*            35.3% +, 11.8% -, 52.9% *, and 00.0% /   */
#/*******************************************************/

   x = piref / ( three * m )#              /*********************/
   s = 0.0#                                        /*  Loop 4.          */
   v = 0.0#                                        /*********************/

   dtime(TimeArray)
   for i in arange(1, m):
      v = v + one
      u = v * x
      w = u * u
      s = s + u * ((((((A6*w-A5)*w+A4)*w-A3)*w+A2)*w+A1)*w+one)
   dtime(TimeArray)
   T[9]  = T[1] * TimeArray[1] - nulltime

   u  = piref / three
   w  = u * u
   sa = u * ((((((A6*w-A5)*w+A4)*w-A3)*w+A2)*w+A1)*w+one)

   T[10] = T[9] / 17.0#                            /*********************/
   sa = x * ( sa + two * s ) / two#                /* sin(x) Results.   */
   sb = 0.5#                                       /*********************/
   sc = sa - sb
   T[11] = one / T[10]
#						  /*********************/
#						  /*   DO NOT REMOVE   */
#						  /*   THIS PRINTOUT!  */
#						  /*********************/
   print("     3   {:13.4e}  {:10.4f}  {:10.4f}\n".format(sc,T[9],T[11]))


   if quick:
     return  0
   
#/************************************************************/
#/* Module 4.  Calculate Integral of cos(x) from 0.0 to PI/3 */
#/*            using the Trapazoidal Method. Result is       */
#/*            sin(PI/3). There are 15 double precision      */
#/*            operations per loop (7 +, 0 -, 8 *, and 0 / ) */
#/*            included in the timing.                       */
#/*            50.0% +, 00.0% -, 50.0% *, 00.0% /            */
#/************************************************************/
   A3 = -A3
   A5 = -A5
   x = piref / ( three * m )#              /*********************/
   s = 0.0#                                        /*  Loop 5.          */
   v = 0.0#                                        /*********************/

   dtime(TimeArray)
   for i in arange(1, m):
      u = i * x
      w = u * u
      s = s + w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1)+one
   dtime(TimeArray)
   T[12]  = T[1] * TimeArray[1] - nulltime

   u  = piref / three
   w  = u * u
   sa = w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1)+one

   T[13] = T[12] / 15.0#                             /*******************/
   sa = x * ( sa + one + two * s ) / two#            /* Module 4 Result */
   u  = piref / three#                               /*******************/
   w  = u * u
   sb = u * ((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w+A0)
   sc = sa - sb
   T[14] = one / T[13]
#						  /*********************/
#						  /*   DO NOT REMOVE   */
#						  /*   THIS PRINTOUT!  */
#						  /*********************/
   print("     4   {:13.4e}  {:10.4f}  {:10.4f}\n".format(sc,T[12],T[14]))

#/************************************************************/
#/* Module 5.  Calculate Integral of tan(x) from 0.0 to PI/3 */
#/*            using the Trapazoidal Method. Result is       */
#/*            ln(cos(PI/3)). There are 29 double precision  */
#/*            operations per loop (13 +, 0 -, 15 *, and 1 /)*/
#/*            included in the timing.                       */
#/*            46.7% +, 00.0% -, 50.0% *, and 03.3% /        */
#/************************************************************/

   x = piref / ( three * m )#              /*********************/
   s = 0.0#                                        /*  Loop 6.          */
   v = 0.0#                                        /*********************/

   dtime(TimeArray)
   for i in arange(1, m):
      u = i * x
      w = u * u
      v = u * ((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w+one)
      s = s + v / (w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1)+one)
   dtime(TimeArray)
   T[15]  = T[1] * TimeArray[1] - nulltime

   u  = piref / three
   w  = u * u
   sa = u*((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w+one)
   sb = w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1)+one
   sa = sa / sb

   T[16] = T[15] / 29.0#                             /*******************/
   sa = x * ( sa + two * s ) / two#                  /* Module 5 Result */
   sb = 0.6931471805599453#                          /*******************/
   sc = sa - sb
   T[17] = one / T[16]
#						  /*********************/
#						  /*   DO NOT REMOVE   */
#						  /*   THIS PRINTOUT!  */
#						  /*********************/
   print("     5   {:13.4e}  {:10.4f}  {:10.4f}\n".format(sc,T[15],T[17]))

#/************************************************************/
#/* Module 6.  Calculate Integral of sin(x)*cos(x) from 0.0  */
#/*            to PI/4 using the Trapazoidal Method. Result  */
#/*            is sin(PI/4)^2. There are 29 double precision */
#/*            operations per loop (13 +, 0 -, 16 *, and 0 /)*/
#/*            included in the timing.                       */
#/*            46.7% +, 00.0% -, 53.3% *, and 00.0% /        */
#/************************************************************/

   x = piref / ( four * m )#               /*********************/
   s = 0.0#                                        /*  Loop 7.          */
   v = 0.0#                                        /*********************/

   dtime(TimeArray)
   for i in arange(1, m):
      u = i * x
      w = u * u
      v = u * ((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w+one)
      s = s + v*(w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1)+one)
   dtime(TimeArray)
   T[18]  = T[1] * TimeArray[1] - nulltime

   u  = piref / four
   w  = u * u
   sa = u*((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w+one)
   sb = w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1)+one
   sa = sa * sb

   T[19] = T[18] / 29.0#                             /*******************/
   sa = x * ( sa + two * s ) / two#                  /* Module 6 Result */
   sb = 0.25#                                        /*******************/
   sc = sa - sb
   T[20] = one / T[19]
#						  /*********************/
#						  /*   DO NOT REMOVE   */
#						  /*   THIS PRINTOUT!  */
#						  /*********************/
   print("     6   {:13.4e}  {:10.4f} {:10.4f}\n".format(sc,T[18],T[20]))


#/*******************************************************/
#/* Module 7.  Calculate value of the definite integral */
#/*            from 0 to sa of 1/(x+1), x/(x*x+1), and  */
#/*            x*x/(x*x*x+1) using the Trapizoidal Rule.*/
#/*            There are 12 double precision operations */
#/*            per loop ( 3 +, 3 -, 3 *, and 3 / ) that */
#/*            are included in the timing.              */
#/*            25.0% +, 25.0% -, 25.0% *, and 25.0% /   */
#/*******************************************************/

#                                                  /*********************/
   s = 0.0#                                        /* Loop 8.           */
   w = one#                                        /*********************/
   sa = 102.3321513995275
   v = sa / m

   dtime(TimeArray)
   for i in arange(1, m):
      x = i * v
      u = x * x
      s = s - w / ( x + w ) - x / ( u + w ) - u / ( x * u + w )
   dtime(TimeArray)
   T[21] = T[1] * TimeArray[1] - nulltime
#						  /*********************/
#						  /* Module 7 Results  */
#						  /*********************/
   T[22] = T[21] / 12.0                                  
   x  = sa                                      
   u  = x * x
   sa = -w - w / ( x + w ) - x / ( u + w ) - u / ( x * u + w )
   sa = 18.0 * v * (sa + two * s )

   m  = -2000 * int(sa)
   m = int( m / scale )

   sc = sa + 500.2
   T[23] = one / T[22]
#						  /********************/
#						  /*  DO NOT REMOVE   */
#						  /*  THIS PRINTOUT!  */
#						  /********************/
   print("     7   {:13.4e}  {:10.4f]  {:10.4f}\n".format(sc,T[21],T[23]))

#/************************************************************/
#/* Module 8.  Calculate Integral of sin(x)*cos(x)*cos(x)    */
#/*            from 0 to PI/3 using the Trapazoidal Method.  */
#/*            Result is (1-cos(PI/3)^3)/3. There are 30     */
#/*            double precision operations per loop included */
#/*            in the timing:                                */
#/*               13 +,     0 -,    17 *          0 /        */
#/*            46.7% +, 00.0% -, 53.3% *, and 00.0% /        */
#/************************************************************/

   x = piref / ( three * m )#              /*********************/
   s = 0.0#                                        /*  Loop 9.          */
   v = 0.0#                                        /*********************/

   dtime(TimeArray)
   for i in range(1, m):
      u = i * x
      w = u * u
      v = w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1)+one
      s = s + v*v*u*((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w+one)
   dtime(TimeArray)
   T[24]  = T[1] * TimeArray[1] - nulltime

   u  = piref / three
   w  = u * u
   sa = u*((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w+one)
   sb = w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1)+one
   sa = sa * sb * sb

   T[25] = T[24] / 30.0#                             /*******************/
   sa = x * ( sa + two * s ) / two#                  /* Module 8 Result */
   sb = 0.29166666666666667#                         /*******************/
   sc = sa - sb
   T[26] = one / T[25]
#						  /*********************/
#						  /*   DO NOT REMOVE   */
#						  /*   THIS PRINTOUT!  */
#						  /*********************/
   print("     8   {:13.4e}  {:10.4f}  {:10.4f}\n".format(sc,T[24],T[26]))

#/**************************************************/   
#/* MFLOPS(1) output. This is the same weighting   */
#/* used for all previous versions of the flops.c  */
#/* program. Includes Modules 2 and 3 only.        */
#/**************************************************/ 
   T[27] = ( five * (T[6] - T[5]) + T[9] ) / 52.0
   T[28] = one  / T[27]

#/**************************************************/   
#/* MFLOPS(2) output. This output does not include */
#/* Module 2, but it still does 9.2% FDIV's.       */
#/**************************************************/ 
   T[29] = T[2] + T[9] + T[12] + T[15] + T[18]
   T[29] = (T[29] + four * T[21]) / 152.0
   T[30] = one / T[29]

#/**************************************************/   
#/* MFLOPS(3) output. This output does not include */
#/* Module 2, but it still does 3.4% FDIV's.       */
#/**************************************************/ 
   T[31] = T[2] + T[9] + T[12] + T[15] + T[18]
   T[31] = (T[31] + T[21] + T[24]) / 146.0
   T[32] = one / T[31]

#/**************************************************/   
#/* MFLOPS(4) output. This output does not include */
#/* Module 2, and it does NO FDIV's.               */
#/**************************************************/ 
   T[33] = (T[9] + T[12] + T[18] + T[24]) / 91.0
   T[34] = one / T[33]


   print("\n")
   print("   Iterations      = {:10d\n}".format(m))
   print("   NullTime (usec) = {:10.4f}\n".format(nulltime))
   print("   MFLOPS(1)       = {:10.4f}\n".format(T[28]))
   print("   MFLOPS(2)       = {:10.4f}\n".format(T[30]))
   print("   MFLOPS(3)       = {:10.4f}\n".format(T[32]))
   print("   MFLOPS(4)       = {:10.4f}\n\n".format(T[34]))

if __name__ == "__main__":
   flops()
