/* 	$Id: ephem_earth.c,v 1.1 2006/11/22 20:31:50 observe Exp $	 */
/* 	$Id: ephem_earth.c,v 1.1 2006/11/22 20:31:50 observe Exp $	 */
#ifndef lint
static char vcid[] = "$Id: ephem_earth.c,v 1.1 2006/11/22 20:31:50 observe Exp $"; 
#endif /* lint */
/*** ephem_earth.c  - I've changed the below to include a routine explicitly
*** returning the location of earth geocenter relative to SSBARY.  Also
*** have eliminated the nutation & libration routines.
*** 8/9/99 gmb
***
*** Callers pass UTC Julian Dates (MPC times).  DE405 is argumented in
*** Teph, which differs from TT by under 2 ms, so the interpolators
*** evaluate the ephemeris at TT.  TT - UTC = 32.184 s + (TAI - UTC).
*** TAI - UTC is the leap-second offset and is not constant: the old
*** fixed 64 s was the 1999 value (actually 64.184 s) and is 5.184 s
*** short for dates since the 2017 January 1 leap second.
*** Local sidereal time is still computed from the UTC date.
***/
#include "orbfit.h"
#include <string.h>
#include <ctype.h>
#define FNAMESIZE 256
#define BUFFSIZE  512
/* structures used to keep observatory information: */
static int nsites=0;
typedef struct {
  int code;
  int	space;		/*no fixed location: a space observatory*/
  int	warned;		/*already warned that this site has no location*/
  double	lon;		/*observatory west longitude in hours*/
  double	rhocos;		/*rho cos(phi') in Earth equatorial radii*/
  double	rhosin;		/*rho sin(phi') in Earth equatorial radii*/
  char	name[80];
} SITE;
static SITE *sitelist=NULL;

char  observatory_file[FNAMESIZE]="";
void
set_observatory_file(char *fname) {
  strncpy(observatory_file, fname, FNAMESIZE-1);
  observatory_file[FNAMESIZE-1]=0;
}

char  ephem_file[FNAMESIZE]="";
void
set_ephem_file(char *fname) {
  strncpy(ephem_file, fname, FNAMESIZE-1);
  ephem_file[FNAMESIZE-1]=0;
}



/******************************************************************************/
/**                                                                          **/
/**  SOURCE FILE: ephem_read.c                                               **/
/**                                                                          **/
/**     This file contains a set of functions and global variables that      **/
/**     implements an ephemeris server program. Client programs can use      **/
/**     use this server to access ephemeris data by calling one of the       **/
/**     following functions:                                                 **/
/**                                                                          **/
/**        Interpolate_Libration -- returns lunar libration angles           **/
/**        Interpolate_Nutation  -- returns (terrestrial) nutation angles    **/
/**        Interpolate_Position  -- returns position of planet               **/
/**        Interpolate_State     -- returns position and velocity of planet  **/
/**                                                                          **/
/**     Note that client programs must make one, and only one, call to the   **/
/**     function Initialize_Ephemeris before any of these other functions    **/ 
/**     be used. After this call any of the above functions can be called    **/
/**     in any sequence.                                                     **/
/**                                                                          **/
/**  Programmer: David Hoffman/EG5                                           **/
/**              NASA, Johnson Space Center                                  **/
/**              Houston, TX 77058                                           **/
/**              e-mail: david.a.hoffman1@jsc.nasa.gov                       **/
/**                                                                          **/
/******************************************************************************/

#include <stdio.h>
#include <math.h>
#ifndef TYPES_DEFINED
#include "ephem_types.h"
#endif

#ifndef PI
#define PI 3.1415926535
#endif

/**==========================================================================**/
/**  Global Variables                                                        **/
/**==========================================================================**/

   static headOneType  H1;
   static headTwoType  H2;
   static recOneType   R1;
   static FILE        *Ephemeris_File;
   static double       Coeff_Array[ARRAY_SIZE] , T_beg , T_end , T_span;

   static int Debug = FALSE;             /* Generates detailed output if true */

/**==========================================================================**/
/**  Read_Coefficients                                                       **/
/**                                                                          **/
/**     This function is used by the functions below to read an array of     **/
/**     Tchebeychev coefficients from a binary ephemeris data file.          **/
/**                                                                          **/
/**  Input: Desired record time.                                             **/
/**                                                                          **/
/**  Output: None.                                                           **/
/**                                                                          **/
/**==========================================================================**/

void Read_Coefficients( double Time )
{
  double  T_delta = 0.0;
  long     Offset  =  0 ;		/*** ??? change to long 8/9/99 ***/

  /*--------------------------------------------------------------------------*/
  /*  Find ephemeris data that record contains input time. Note that one, and */
  /*  only one, of the following conditional statements will be true (if both */
  /*  were false, this function would not have been called).                  */
  /*--------------------------------------------------------------------------*/

  if ( Time < T_beg )                    /* Compute backwards location offset */
     {
       T_delta = T_beg - Time;
       Offset  = (int) -ceil(T_delta/T_span);	/***Needed negative sign here???*/
     }

  if ( Time > T_end )                    /* Compute forewards location offset */
     {
       T_delta = Time - T_end;
       Offset  = (int) ceil(T_delta/T_span);
     }

  /*--------------------------------------------------------------------------*/
  /*  Retrieve ephemeris data from new record.                                */
  /*--------------------------------------------------------------------------*/

  fseek(Ephemeris_File,(Offset-1)*ARRAY_SIZE*sizeof(double),SEEK_CUR);
  fread(&Coeff_Array,sizeof(double),ARRAY_SIZE,Ephemeris_File);
  
  T_beg  = Coeff_Array[0];
  T_end  = Coeff_Array[1];
  T_span = T_end - T_beg;

  if (Time < T_beg || Time > T_end) {
    fprintf(stderr,"JD %.2f is out of range of ephemeris file\n",Time);
    exit(1);
  }
  /*--------------------------------------------------------------------------*/
  /*  Debug print (optional)                                                  */
  /*--------------------------------------------------------------------------*/

  if ( Debug ) 
     {
       printf("\n  In: Read_Coefficients \n");
       printf("\n      ARRAY_SIZE = %4d",ARRAY_SIZE);
       printf("\n      Offset  = %3ld",Offset);
       printf("\n      T_delta = %7.3f",T_delta);
       printf("\n      T_Beg   = %7.3f",T_beg);
       printf("\n      T_End   = %7.3f",T_end);
       printf("\n      T_Span  = %7.3f\n\n",T_span);
     }

}

/**==========================================================================**/
/**  Initialize_Ephemeris                                                    **/
/**                                                                          **/
/**     This function must be called once by any program that accesses the   **/
/**     ephemeris data. It opens the ephemeris data file, reads the header   **/
/**     data, loads the first coefficient record into a global array, then   **/
/**     returns a status code that indicates whether or not all of this was  **/
/**     done successfully.                                                   **/
/**                                                                          **/
/**  Input: A character string giving the name of an ephemeris data file.    **/
/**                                                                          **/
/**  Returns: An integer status code.                                        **/
/**                                                                          **/
/**==========================================================================**/
/** This routine no longer takes filename as argument, it goes looking
 ** for the file itself, as environment-specified file or as default in
 ** this directory
 **/
int Initialize_Ephemeris()
{
  int headerID;
  char fileName[FNAMESIZE];
  /*** gmb: don't duplicate this call */
  static int init=0;
  if (init) return SUCCESS;
  init = 1;

  /*--------------------------------------------------------------------------*/
  /*  Open ephemeris file.                                                    */
  /*--------------------------------------------------------------------------*/

  /** use previously specified filename,
   ** or environment-specified file, or the default filename, in that order
   */
  if (strlen(ephem_file)>0)
    strncpy(fileName, ephem_file, FNAMESIZE-1);
  else if (getenv(EPHEM_ENVIRON)!=NULL)
    strncpy(fileName, getenv(EPHEM_ENVIRON), FNAMESIZE-1);
  else
    strncpy(fileName, DEFAULT_EPHEM_FILE, FNAMESIZE-1);


  fileName[FNAMESIZE-1]=0;

  
  Ephemeris_File = fopen(fileName,"r");

  /*--------------------------------------------------------------------------*/
  /*  Read header & first coefficient array, then return status code.         */
  /*--------------------------------------------------------------------------*/

  if ( Ephemeris_File == NULL ) /*........................No need to continue */
     {
       printf("\n Unable to open ephemeris file: %s.\n",fileName);
       return FAILURE;
     }
  else 
     { /*.................Read first three header records from ephemeris file */
         
       fread(&H1,sizeof(double),ARRAY_SIZE,Ephemeris_File);
       fread(&H2,sizeof(double),ARRAY_SIZE,Ephemeris_File);
       fread(&Coeff_Array,sizeof(double),ARRAY_SIZE,Ephemeris_File);
       
       /*...............................Store header data in global variables */
       
       R1 = H1.data;
              
       /*..........................................Set current time variables */

       T_beg  = Coeff_Array[0];
       T_end  = Coeff_Array[1];
       T_span = T_end - T_beg;

       /*..............................Convert header ephemeris ID to integer */

       headerID = (int) R1.DENUM;
       
       /*..............................................Debug Print (optional) */

       if ( Debug ) 
          {
            printf("\n  In: Initialize_Ephemeris \n");
            printf("\n      ARRAY_SIZE = %4d",ARRAY_SIZE);
            printf("\n      headerID   = %3d",headerID);
            printf("\n      T_Beg      = %7.3f",T_beg);
            printf("\n      T_End      = %7.3f",T_end);
            printf("\n      T_Span     = %7.3f\n\n",T_span);
          }

       /*..................................................Return status code */
       
       if ( headerID == EPHEMERIS ) 
          {
            return SUCCESS;
          }
       else 
          {
            printf("\n Opened wrong file: %s",fileName);
            printf(" for ephemeris: %d.\n",EPHEMERIS);
            return FAILURE;
          }
     }
}

/**==========================================================================**/
/**  utc_jd_to_tt                                                            **/
/**                                                                          **/
/**  TT = TAI + 32.184 s, and TAI - UTC is the IERS leap-second offset.      **/
/**  Before 1972 the offset also drifts (USNO tai-utc.dat).  After the      **/
/**  last tabulated leap second the latest offset is held.  IERS Bulletin   **/
/**  C 72 (2026 July 6) keeps TAI - UTC = 37 s through 2026 December; a     **/
/**  newly announced leap second has to be added to the table below.        **/
/**  Dates before 1960 (when UTC was not defined) use the 1960 Jan 1        **/
/**  expression.                                                             **/
/**                                                                          **/
/**==========================================================================**/

/* TAI - UTC (seconds) at a UTC Julian Date.  The fraction of jd_utc is the
 * fraction of that UTC day, including a smeared leap second. */
static double tai_minus_utc(double jd_utc)
{
  /* jd0: UTC JD when this expression starts.
     dat0, mjd0, rate: TAI-UTC = dat0 + (MJD - mjd0) * rate seconds.
     Values are the IERS/USNO tai-utc series.  1960-01-01 continues the
     1961 rate back to the start of UTC. */
  static const struct {
    double jd0, dat0, mjd0, rate;
  } leap[] = {
    { 2436934.5,  1.4178180, 37300.0, 0.0012960 },  /* 1960-01-01 */
    { 2437300.5,  1.4228180, 37300.0, 0.0012960 },  /* 1961-01-01 */
    { 2437512.5,  1.3728180, 37300.0, 0.0012960 },  /* 1961-08-01 */
    { 2437665.5,  1.8458580, 37665.0, 0.0011232 },  /* 1962-01-01 */
    { 2438334.5,  1.9458580, 37665.0, 0.0011232 },  /* 1963-11-01 */
    { 2438395.5,  3.2401300, 38761.0, 0.0012960 },  /* 1964-01-01 */
    { 2438486.5,  3.3401300, 38761.0, 0.0012960 },  /* 1964-04-01 */
    { 2438639.5,  3.4401300, 38761.0, 0.0012960 },  /* 1964-09-01 */
    { 2438761.5,  3.5401300, 38761.0, 0.0012960 },  /* 1965-01-01 */
    { 2438820.5,  3.6401300, 38761.0, 0.0012960 },  /* 1965-03-01 */
    { 2438942.5,  3.7401300, 38761.0, 0.0012960 },  /* 1965-07-01 */
    { 2439004.5,  3.8401300, 38761.0, 0.0012960 },  /* 1965-09-01 */
    { 2439126.5,  4.3131700, 39126.0, 0.0025920 },  /* 1966-01-01 */
    { 2439887.5,  4.2131700, 39126.0, 0.0025920 },  /* 1968-02-01 */
    { 2441317.5, 10.0, 0.0, 0.0 },                 /* 1972-01-01 */
    { 2441499.5, 11.0, 0.0, 0.0 },                 /* 1972-07-01 */
    { 2441683.5, 12.0, 0.0, 0.0 },                 /* 1973-01-01 */
    { 2442048.5, 13.0, 0.0, 0.0 },                 /* 1974-01-01 */
    { 2442413.5, 14.0, 0.0, 0.0 },                 /* 1975-01-01 */
    { 2442778.5, 15.0, 0.0, 0.0 },                 /* 1976-01-01 */
    { 2443144.5, 16.0, 0.0, 0.0 },                 /* 1977-01-01 */
    { 2443509.5, 17.0, 0.0, 0.0 },                 /* 1978-01-01 */
    { 2443874.5, 18.0, 0.0, 0.0 },                 /* 1979-01-01 */
    { 2444239.5, 19.0, 0.0, 0.0 },                 /* 1980-01-01 */
    { 2444786.5, 20.0, 0.0, 0.0 },                 /* 1981-07-01 */
    { 2445151.5, 21.0, 0.0, 0.0 },                 /* 1982-07-01 */
    { 2445516.5, 22.0, 0.0, 0.0 },                 /* 1983-07-01 */
    { 2446247.5, 23.0, 0.0, 0.0 },                 /* 1985-07-01 */
    { 2447161.5, 24.0, 0.0, 0.0 },                 /* 1988-01-01 */
    { 2447892.5, 25.0, 0.0, 0.0 },                 /* 1990-01-01 */
    { 2448257.5, 26.0, 0.0, 0.0 },                 /* 1991-01-01 */
    { 2448804.5, 27.0, 0.0, 0.0 },                 /* 1992-07-01 */
    { 2449169.5, 28.0, 0.0, 0.0 },                 /* 1993-07-01 */
    { 2449534.5, 29.0, 0.0, 0.0 },                 /* 1994-07-01 */
    { 2450083.5, 30.0, 0.0, 0.0 },                 /* 1996-01-01 */
    { 2450630.5, 31.0, 0.0, 0.0 },                 /* 1997-07-01 */
    { 2451179.5, 32.0, 0.0, 0.0 },                 /* 1999-01-01 */
    { 2453736.5, 33.0, 0.0, 0.0 },                 /* 2006-01-01 */
    { 2454832.5, 34.0, 0.0, 0.0 },                 /* 2009-01-01 */
    { 2456109.5, 35.0, 0.0, 0.0 },                 /* 2012-07-01 */
    { 2457204.5, 36.0, 0.0, 0.0 },                 /* 2015-07-01 */
    { 2457754.5, 37.0, 0.0, 0.0 }                  /* 2017-01-01 */
  };
  const int nleap = (int)(sizeof leap / sizeof leap[0]);
  int i;
  double mjd;

  if (jd_utc < leap[0].jd0)
    i = 0;
  else {
    for (i = nleap - 1; i > 0; i--)
      if (jd_utc >= leap[i].jd0) break;
  }

  mjd = jd_utc - 2400000.5;
  return leap[i].dat0 + (mjd - leap[i].mjd0) * leap[i].rate;
}

double utc_jd_to_tt(double jd_utc)
{
  /* UTC JD runs one calendar day per JD day, so a leap-second day is
   * squeezed into 86400 JD seconds.  Undo that, then add TAI-UTC at 0h
   * and the 32.184 s TT-TAI offset.  Same construction as SOFA utctai. */
  double jd0, fd, dat0, dat12, dat24, dlod, dleap;
  const double tt_minus_tai = 32.184;
  const double sec_per_day = 86400.0;

  jd0 = floor(jd_utc - 0.5) + 0.5;
  fd = jd_utc - jd0;
  dat0 = tai_minus_utc(jd0);
  dat12 = tai_minus_utc(jd0 + 0.5);
  dat24 = tai_minus_utc(jd0 + 1.0);
  dlod = 2.0 * (dat12 - dat0);
  dleap = dat24 - (dat0 + dlod);
  fd *= (sec_per_day + dleap) / sec_per_day;
  fd *= (sec_per_day + dlod) / sec_per_day;
  return jd0 + fd + (dat0 + tt_minus_tai) / sec_per_day;
}

/**==========================================================================**/
/**  Interpolate_Position                                                    **/
/**                                                                          **/
/**     This function computes a position vector for a selected planetary    **/
/**     body from Chebyshev coefficients read in from an ephemeris data      **/
/**     file. These coefficients are read from the data file by calling      **/
/**     the function Read_Coefficients (when necessary).                     **/
/**                                                                          **/
/**  Inputs:                                                                 **/
/**     Time     -- UTC Julian Date.  Converted to TT before interpolation.  **/
/**     Target   -- Solar system body for which position is desired.         **/
/**     Position -- Pointer to external array to receive the position.       **/
/**                                                                          **/
/**  Returns: Nothing explicitly.                                            **/
/**                                                                          **/
/**==========================================================================**/

void Interpolate_Position( double Time , int Target , double Position[3] )
{
  double    A[50] , Cp[50]  , sum[3] , T_break , T_seg , T_sub , Tc;
  int       i , j;
  long int  C , G , N , offset = 0;

  Time = utc_jd_to_tt(Time);
  /*--------------------------------------------------------------------------*/
  /* This function doesn't "do" nutations or librations.                      */
  /*--------------------------------------------------------------------------*/

  if ( Target >= 11 )             /* Also protects against weird input errors */
     {
       printf("\n This function does not compute nutations or librations.\n");
       return;
     }
 
  /*--------------------------------------------------------------------------*/
  /* Initialize local coefficient array.                                      */
  /*--------------------------------------------------------------------------*/

  for ( i=0 ; i<50 ; i++ )
      {
        A[i] = 0.0;
      }

  /*--------------------------------------------------------------------------*/
  /* Determine if a new record needs to be input (if so, get it).             */
  /*--------------------------------------------------------------------------*/
    
  if (Time < T_beg || Time > T_end)  Read_Coefficients(Time);

  /*--------------------------------------------------------------------------*/
  /* Read the coefficients from the binary record.                            */
  /*--------------------------------------------------------------------------*/
  
  C = R1.coeffPtr[Target][0] - 1;          /*   Coefficient array entry point */
  N = R1.coeffPtr[Target][1];              /* Number of coeff's per component */
  G = R1.coeffPtr[Target][2];              /*      Granules in current record */

  /*...................................................Debug print (optional) */

  if ( Debug )
     {
       printf("\n  In: Interpolate_Position\n");
       printf("\n  Target = %2d",Target);
       printf("\n  C      = %4ld (before)",C);
       printf("\n  N      = %4ld",N);
       printf("\n  G      = %4ld\n",G);
     }

  /*--------------------------------------------------------------------------*/
  /*  Compute the normalized time, then load the Tchebeyshev coefficients     */
  /*  into array A[]. If T_span is covered by a single granule this is easy.  */
  /*  If not, the granule that contains the interpolation time is found, and  */
  /*  an offset from the array entry point for the ephemeris body is used to  */
  /*  load the coefficients.                                                  */
  /*--------------------------------------------------------------------------*/

  if ( G == 1 )
     {
       Tc = 2.0*(Time - T_beg) / T_span - 1.0;
       for (i=C ; i<(C+3*N) ; i++)  A[i-C] = Coeff_Array[i];
     }
  else if ( G > 1 )
     {
       T_sub = T_span / ((double) G);          /* Compute subgranule interval */
       T_seg = T_beg + T_sub ;                 /* Set T_seg to the smallest value that is reasonable */
       for ( j=G ; j>0 ; j-- ) 
           {
             T_break = T_beg + ((double) j-1) * T_sub;
             if ( Time > T_break ) 
                {
                  T_seg  = T_break;
                  offset = j-1;
                  break;
                }
            }
            
       Tc = 2.0*(Time - T_seg) / T_sub - 1.0;
       C  = C + 3 * offset * N;
       
       for (i=C ; i<(C+3*N) ; i++) A[i-C] = Coeff_Array[i];
     }
  else                                   /* Something has gone terribly wrong */
     {
       printf("\n Number of granules must be >= 1: check header data.\n");
     }

  /*...................................................Debug print (optional) */

  if ( Debug )
     {
       printf("\n  C      = %4ld (after)",C);
       printf("\n  offset = %4ld",offset);
       printf("\n  Time   = %12.7f",Time);
       printf("\n  T_sub  = %12.7f",T_sub);
       printf("\n  T_seg  = %12.7f",T_seg);
       printf("\n  Tc     = %12.7f\n",Tc);
       printf("\n  Array Coefficients:\n");
       for ( i=0 ; i<3*N ; i++ )
           {
             printf("\n  A[%2d] = % 22.15e",i,A[i]);
           }
       printf("\n\n");
     }

  /*..........................................................................*/

  /*--------------------------------------------------------------------------*/
  /* Compute interpolated the position.                                       */
  /*--------------------------------------------------------------------------*/
  
  for ( i=0 ; i<3 ; i++ ) 
      {                           
        Cp[0]  = 1.0;                                 /* Begin polynomial sum */
        Cp[1]  = Tc;
        sum[i] = A[i*N] + A[1+i*N]*Tc;

        for ( j=2 ; j<N ; j++ )                                  /* Finish it */
            {
              Cp[j]  = 2.0 * Tc * Cp[j-1] - Cp[j-2];
              sum[i] = sum[i] + A[j+i*N] * Cp[j];
            }
        Position[i] = sum[i];
      }

  return;
}

void Interpolate_State(double Time, 
		       int Target, 
		       double Position[3], 
		       double Velocity[3])
{
  double    A[50]   , B[50] , Cp[50] , P_Sum[3] , V_Sum[3] , Up[50] ,
            T_break , T_seg , T_sub  , Tc;
  int       i , j;
  long int  C , G , N , offset = 0;

  Time = utc_jd_to_tt(Time);

  /*--------------------------------------------------------------------------*/
  /* This function doesn't "do" nutations or librations.                      */
  /*--------------------------------------------------------------------------*/

  if ( Target >= 11 )             /* Also protects against weird input errors */
     {
       printf("\n This function does not compute nutations or librations.\n");
       return;
     }

  /*--------------------------------------------------------------------------*/
  /* Initialize local coefficient array.                                      */
  /*--------------------------------------------------------------------------*/

  for ( i=0 ; i<50 ; i++ )
      {
        A[i] = 0.0;
        B[i] = 0.0;
      }

  /*--------------------------------------------------------------------------*/
  /* Determine if a new record needs to be input.                             */
  /*--------------------------------------------------------------------------*/
  
  if (Time < T_beg || Time > T_end)  Read_Coefficients(Time);

  /*--------------------------------------------------------------------------*/
  /* Read the coefficients from the binary record.                            */
  /*--------------------------------------------------------------------------*/
  
  C = R1.coeffPtr[Target][0] - 1;               /*    Coeff array entry point */
  N = R1.coeffPtr[Target][1];                   /*          Number of coeff's */
  G = R1.coeffPtr[Target][2];                   /* Granules in current record */

  /*...................................................Debug print (optional) */

  if ( Debug )
     {
       printf("\n  In: Interpolate_State\n");
       printf("\n  Target = %2d",Target);
       printf("\n  C      = %4ld (before)",C);
       printf("\n  N      = %4ld",N);
       printf("\n  G      = %4ld\n",G);
     }

  /*--------------------------------------------------------------------------*/
  /*  Compute the normalized time, then load the Tchebeyshev coefficients     */
  /*  into array A[]. If T_span is covered by a single granule this is easy.  */
  /*  If not, the granule that contains the interpolation time is found, and  */
  /*  an offset from the array entry point for the ephemeris body is used to  */
  /*  load the coefficients.                                                  */
  /*--------------------------------------------------------------------------*/

  if ( G == 1 )
     {
       Tc = 2.0*(Time - T_beg) / T_span - 1.0;
       for (i=C ; i<(C+3*N) ; i++)  A[i-C] = Coeff_Array[i];
     }
  else if ( G > 1 )
     {
       T_sub = T_span / ((double) G);          /* Compute subgranule interval */
       T_seg = T_beg + T_sub ;                 /* Set T_seg to the smallest value that is reasonable */
       for ( j=G ; j>0 ; j-- ) 
           {
             T_break = T_beg + ((double) j-1) * T_sub;
             if ( Time > T_break ) 
                {
                  T_seg  = T_break;
                  offset = j-1;
                  break;
                }
            }
            
       Tc = 2.0*(Time - T_seg) / T_sub - 1.0;
       C  = C + 3 * offset * N;
       
       for (i=C ; i<(C+3*N) ; i++) A[i-C] = Coeff_Array[i];
     }
  else                                   /* Something has gone terribly wrong */
     {
       printf("\n Number of granules must be >= 1: check header data.\n");
     }

  /*...................................................Debug print (optional) */

  if ( Debug )
     {
       printf("\n  C      = %4ld (after)",C);
       printf("\n  offset = %4ld",offset);
       printf("\n  Time   = %12.7f",Time);
       printf("\n  T_sub  = %12.7f",T_sub);
       printf("\n  T_seg  = %12.7f",T_seg);
       printf("\n  Tc     = %12.7f\n",Tc);
       printf("\n  Array Coefficients:\n");
       for ( i=0 ; i<3*N ; i++ )
           {
             printf("\n  A[%2d] = % 22.15e",i,A[i]);
           }
       printf("\n\n");
     }

  /*..........................................................................*/

  /*--------------------------------------------------------------------------*/
  /* Compute the interpolated position & velocity                             */
  /*--------------------------------------------------------------------------*/
  
  for ( i=0 ; i<3 ; i++ )                /* Compute interpolating polynomials */
      {
        Cp[0] = 1.0;           
        Cp[1] = Tc;
        Cp[2] = 2.0 * Tc*Tc - 1.0;
        
        Up[0] = 0.0;
        Up[1] = 1.0;
        Up[2] = 4.0 * Tc;

        for ( j=3 ; j<N ; j++ )
            {
              Cp[j] = 2.0 * Tc * Cp[j-1] - Cp[j-2];
              Up[j] = 2.0 * Tc * Up[j-1] + 2.0 * Cp[j-1] - Up[j-2];
            }

        P_Sum[i] = 0.0;           /* Compute interpolated position & velocity */
        V_Sum[i] = 0.0;

        for ( j=N-1 ; j>-1 ; j-- )  P_Sum[i] = P_Sum[i] + A[j+i*N] * Cp[j];
        for ( j=N-1 ; j>0  ; j-- )  V_Sum[i] = V_Sum[i] + A[j+i*N] * Up[j];

        Position[i] = P_Sum[i];
        Velocity[i] = V_Sum[i] * 2.0 * ((double) G) / (T_span * 86400.0);
      }

  /*--------------------------------------------------------------------------*/
  /*  Return computed values.                                                 */
  /*--------------------------------------------------------------------------*/

  return;
}

/********************************************************* END: ephem_read.c **/

/* Give the Earth geocenter wrt SSBary in AU. */
void
geocenter_ssbary(double jd,
		double *xyz)
{
  double embary[3],moon[3];
  static int init=0;
  int i;

  if (!init) {
    if (Initialize_Ephemeris()) exit(1);
    init = 1;
  }

  Interpolate_Position(jd, EARTH, embary);
  Interpolate_Position(jd, MOON, moon);
  for (i=0; i<3; i++) {
    xyz[i] = (embary[i] - moon[i]/(1.+R1.EMRAT)) / R1.AU;
  }

  return;
}

/********** Stealing some routines from skycalc to do the shift
******** from geocenter to observatory.
********/
#define  J2000             2451545.        /* Julian date at standard epoch */
#define  SEC_IN_DAY        86400.
#define  FLATTEN           0.003352813   /* flattening of earth, 1/298.257 */
#define  EQUAT_RAD         6378137.    /* equatorial radius of earth, meters */

double 
lst(double jd,
    double longit)
{
	/* returns the local MEAN sidereal time (dec hrs) at julian date jd
	   at west longitude long (decimal hours).  Follows
	   definitions in 1992 Astronomical Almanac, pp. B7 and L2.
	   Expression for GMST at 0h ut referenced to Aoki et al, A&A 105,
	   p.359, 1982.  On workstations, accuracy (numerical only!)
	   is about a millisecond in the 1990s. */

	double t, ut, jdmid, jdint, jdfrac, sid_g;
	long jdin, sid_int;

	jdin = jd;         /* fossil code from earlier package which
			split jd into integer and fractional parts ... */
	jdint = jdin;
	jdfrac = jd - jdint;
	if(jdfrac < 0.5) {
		jdmid = jdint - 0.5;
		ut = jdfrac + 0.5;
	}
	else {
		jdmid = jdint + 0.5;
		ut = jdfrac - 0.5;
	}
	t = (jdmid - J2000)/36525;
	sid_g = (24110.54841+8640184.812866*t+0.093104*t*t-6.2e-6*t*t*t)/SEC_IN_DAY;
	sid_int = sid_g;
	sid_g = sid_g - (double) sid_int;
	sid_g = sid_g + 1.0027379093 * ut - longit/24.;
	sid_int = sid_g;
	sid_g = (sid_g - (double) sid_int) * 24.;
	if(sid_g < 0.) sid_g = sid_g + 24.;
	return(sid_g);
}

void 
topo(double lmst, double rhocos, double rhosin,
	double *x_geo, double *y_geo, double *z_geo)
/* computes the geocentric equatorial vector (AU) of a site from its
 * local mean sidereal time (radians, positive eastward) and its MPC
 * parallax constants rho cos(phi') and rho sin(phi') (Earth radii).
 */
{
	*x_geo = rhocos * cos(lmst);
	*y_geo = rhocos * sin(lmst);
	*z_geo = rhosin;

	/* EQUAT_RAD is in m and R1.AU is in km */
	*x_geo *= EQUAT_RAD / (1000.*R1.AU);
	*y_geo *= EQUAT_RAD / (1000.*R1.AU);
	*z_geo *= EQUAT_RAD / (1000.*R1.AU);
}

/* OK, now put all of the above together to give the observatory coordinates
 * wrt SSBary in ICRS/AU.
 */
void
earth_ssbary(double jd,
	     int    obscode,
	     double *xtopo, double *ytopo, double *ztopo)
{
  double	xgeo[3], xobs,yobs,zobs;

  /* Return zero if we are "observing" from ssbary */
  if (obscode==OBSCODE_SSBARY) {
    *xtopo = *ytopo = *ztopo = 0.;
    return;
  }

  /* Use Horizons to get the coordinates of geocenter */
  geocenter_ssbary(jd,xgeo);


  /* fprintf(stderr,"geocenter ICRS:  %g %g %g\n",xgeo[0],xgeo[1],xgeo[2]); */
  /* Now get the geocenter->observatory vector */
  observatory_geocenter(jd, obscode, &xobs, &yobs, &zobs);

  /* fprintf(stderr,"observ ICRS:  %g %g %g\n", xobs, yobs, zobs); */

  /* Add in to give observatory coords in ICRS */
  *xtopo = xgeo[0] + xobs;
  *ytopo = xgeo[1] + yobs;
  *ztopo = xgeo[2] + zobs;


  return;
}

/* convert x/y/z geo location to ssbary x/y/z */
void geo_to_ssbary(double jd, double *x, double *y, double *z)
{
  double xgeo[3];

  /* determine where the geocenter is */
  geocenter_ssbary(jd, xgeo);
  /* fprintf(stderr, "XGEO: %lf %lf %lf\n", xgeo[0], xgeo[1], xgeo[2]);
   fprintf(stderr, "OBS: %lf %lf %lf (km)\n", *x, *y, *z); */
  *x /= R1.AU;
  *y /= R1.AU;
  *z /= R1.AU;
  /* fprintf(stderr, "OBS: %lf %lf %lf (au)\n", *x, *y, *z); */
  /* Add in to give observatory coords in ICRS */
  *x = xgeo[0] + *x/R1.AU;
  *y = xgeo[1] + *y/R1.AU;
  *z = xgeo[2] + *z/R1.AU;
  /* fprintf(stderr, "PRO: %lf %lf %lf\n", *x, *y, *z); */
  return ;
}


int
obscode_from_string(const char *code)
{
  int i, n, prefix;

  n = strlen(code);
  if (n<1) return OBSCODE_INVALID;
  for (i=1; i<n; i++)
    if (!isdigit((unsigned char) code[i])) return OBSCODE_INVALID;
  if (isdigit((unsigned char) code[0])) return atoi(code);
  if (n!=3) return OBSCODE_INVALID;
  if (code[0]>='A' && code[0]<='Z') prefix = code[0]-'A'+10;
  else if (code[0]>='a' && code[0]<='z') prefix = code[0]-'a'+36;
  else return OBSCODE_INVALID;
  return prefix*100 + atoi(code+1);
}

static SITE *
find_site(int obscode)
{
  int i;
  if (nsites<=0) read_observatories(NULL);
  for (i=0; i<nsites; i++)
    if (sitelist[i].code==obscode) return &(sitelist[i]);
  return NULL;
}

/* return vector from geocenter->observatory in ICRS coordinates*/
void
observatory_geocenter(double jd,
		      int obscode,
		      double *xobs,
		      double *yobs,
		      double *zobs) {

  static int last_unknown=OBSCODE_INVALID;
  SITE *site;

  *xobs=*yobs=*zobs=0.;
  if (obscode==OBSCODE_GEOCENTER) return;

  site = find_site(obscode);
  if (site==NULL) {
    if (obscode!=last_unknown)
      fprintf(stderr,"Unknown observatory code %d, using geocenter\n",obscode);
    last_unknown = obscode;
    return;
  }
  if (site->space) {
    /* Space observatory positions must come with the observation */
    if (!site->warned)
      fprintf(stderr,"Observatory code %d (%s) has no fixed location, using geocenter\n",
	      obscode, site->name);
    site->warned = 1;
    return;
  }

  /* Get the LMST and calculate to ICRS vector */
  topo(lst(jd,site->lon)*PI/12., site->rhocos, site->rhosin, xobs, yobs, zobs);

  return;
}


/* Copy columns [start, end) of a line, or fewer if the line is short */
static void
column(const char *line, int start, int end, char *out)
{
  int n = strcspn(line, "\r\n");
  if (n<start) {
    out[0]=0;
    return;
  }
  if (end>n) end=n;
  strncpy(out, line+start, end-start);
  out[end-start]=0;
}

/* Read the look-up table for observatories, in the fixed-column format of
 * https://minorplanetcenter.net/iau/lists/ObsCodes.html :
 *   cols 1-3 code, 5-13 east longitude (deg), 14-21 rho cos(phi'),
 *   22-30 rho sin(phi'), 31- name.
 * Space observatories leave the longitude and parallax constants blank.
 * Lines that do not start with a code (headers, <pre>) are skipped.
 * Ground-based longitudes are stored as west longitude in hours. */
void
read_observatories(char *fname)
{

  FILE *sitefile;
  char  inbuff[BUFFSIZE], code[4], lonstring[16], cosstring[16], sinstring[16];
  int obscode, maxsites=0;
  char fileName[FNAMESIZE];

  nsites = 0;
  free(sitelist);
  sitelist = NULL;

  /** use passed filename, or a previously specified filename,
   ** or environment-specified file, or the default filename, in that order
   */
  if (fname != NULL)
    strncpy(fileName, fname, FNAMESIZE-1);
  else if (strlen(observatory_file)>0)
    strncpy(fileName, observatory_file, FNAMESIZE-1);
  else if (getenv(OBS_ENVIRON)!=NULL)
    strncpy(fileName, getenv(OBS_ENVIRON), FNAMESIZE-1);
  else
    strncpy(fileName, DEFAULT_OBSERVATORY_FILE, FNAMESIZE-1);
  
  fileName[FNAMESIZE-1]=0;
  /* fprintf(stderr,"Loading from %s\n", fileName); */

  if ((sitefile = fopen(fileName,"r"))==NULL) {
    fprintf(stderr,"Error opening observatories file %s\n",fileName);
    exit(1);
  }

  while (fgets_nocomment(inbuff, BUFFSIZE-1, sitefile, NULL)!=NULL) {
    SITE *sss;
    double lon;
    int nlon, ncos, nsin;
    char *eol;

    column(inbuff, 0, 3, code);
    if ((obscode = obscode_from_string(code))==OBSCODE_INVALID) continue;
    if (strchr(inbuff, ':')!=NULL && strchr(inbuff, ':') < inbuff+30) {
      fprintf(stderr,"Observatories file %s is not in the MPC ObsCodes format:\n->%s\n",
	      fileName, inbuff);
      exit(1);
    }

    if (nsites>=maxsites) {
      maxsites = maxsites ? 2*maxsites : 1024;
      if ((sitelist = realloc(sitelist, maxsites*sizeof(SITE)))==NULL) {
	fprintf(stderr,"Out of memory reading observatories file %s\n", fileName);
	exit(1);
      }
    }
    sss = &(sitelist[nsites]);
    sss->code = obscode;
    sss->warned = 0;

    column(inbuff, 4, 13, lonstring);
    column(inbuff, 13, 21, cosstring);
    column(inbuff, 21, 30, sinstring);
    nlon = sscanf(lonstring, "%lf", &lon);
    ncos = sscanf(cosstring, "%lf", &(sss->rhocos));
    nsin = sscanf(sinstring, "%lf", &(sss->rhosin));
    if (nlon==1 && ncos==1 && nsin==1) {
      sss->space = 0;
      sss->lon = (360. - lon)/15.;	/*east degrees to west hours*/
    } else if (nlon<=0 && ncos<=0 && nsin<=0) {
      sss->space = 1;
      sss->lon = sss->rhocos = sss->rhosin = 0.;
    } else {
      fprintf(stderr,"Bad line in observatories file %s:\n->%s\n",
	      fileName, inbuff);
      exit(1);
    }

    column(inbuff, 30, 30+79, sss->name);
    for (eol=sss->name+strlen(sss->name); eol>sss->name && isspace((unsigned char) eol[-1]); eol--) ;
    *eol = 0;
    nsites++;
  }
  fclose(sitefile);
  
  if (nsites<1) {
    fprintf(stderr,"Error: no observatory sites found\n");
    exit(1);
  }
}

/* Get arbitrary planetary barycenter position.  Get the
velocity (in AU/YR) too, if desired. */
void
bodycenter_ssbary(double jd,
		  double *xyz,
		  int body,
		  double *vxyz)
{
  double posn[3], vel[3];
  static int init=0;
  int i;

  if (!init) {
    if (Initialize_Ephemeris()) exit(1);
    init = 1;
  }

  if (vxyz==NULL) {
    Interpolate_Position(jd, body, posn);
  } else {
    Interpolate_State(jd, body, posn, vel);
  }
  for (i=0; i<3; i++) {
    xyz[i] = posn[i] / R1.AU;
  }
  if (vxyz!=NULL) /* convert km/s to AU/YR: */
    for (i=0; i<3; i++) {
      vxyz[i] = vel[i] * (86400. * 365.25) / R1.AU ;
    }
  return;
}

/* Return the angle from zenith to the horizon for this observatory.
 * Observatories without a fixed location have no horizon (PI). */
double
zenith_horizon(int obscode) {
  SITE *site;
  if (obscode==OBSCODE_GEOCENTER) return PI;
  site = find_site(obscode);
  if (site==NULL || site->space) return PI;
  return PI/2.;
}
