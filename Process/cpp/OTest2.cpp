
#ifdef WIN32
   #include "stdafx.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef WIN32

const char 
/*
  *trace_txt  = "..\\tmpProcess\\trace.txt",
  *travel_txt = "..\\tmpProcess\\travel.txt",
  *latv_txt   = "..\\tmpProcess\\latv.txt",
  *lonv_txt   = "..\\tmpProcess\\lonv.txt",
  *semb_txt   = "..\\tmpProcess\\semb.txt";
  */
  *trace_txt  = "E:\\Helmut\\develop5\\ARRAYTOOL_2\\tmpProcess\\trace.txt",
  *travel_txt = "E:\\Helmut\\develop5\\ARRAYTOOL_2\\tmpProcess\\travel.txt",
  *latv_txt   = "E:\\Helmut\\develop5\\ARRAYTOOL_2\\tmpProcess\\latv.txt",
  *lonv_txt   = "E:\\Helmut\\develop5\\ARRAYTOOL_2\\tmpProcess\\lonv.txt",
  *semb_txt   = "E:\\Helmut\\develop5\\ARRAYTOOL_2\\tmpProcess\\semb.txt";

#else 

const char 
  *trace_txt  = "../tmpProcess/trace.txt",
  *travel_txt = "../tmpProcess/travel.txt",
  *latv_txt   = "../tmpProcess/latv.txt",
  *lonv_txt   = "../tmpProcess/lonv.txt",
  *semb_txt   = "../tmpProcess/semb.txt";

#endif

// --------------------------------------------------------------------------------------

static void deb (const char *s)
{   //printf ("%s",s);
}


static void deb_int (const char *s, int v)
{   //printf ("%s = %d ", s,v);
}


static double * stringToDouble (char * s, int size)
{
  double     * vect  = new double  [size];
  const char * delim = ",";
  char       * ptr   = strtok (s, delim);

  for (int i = 0; (ptr = strtok (NULL, delim)) != NULL; i++) {
	  sscanf (ptr, "%lf", &vect[i]);
  }

  return vect;
}

static double * readVector (const char *file, int size)
{
	deb ("readVector\n");
	FILE *fp = fopen (file, "r");
	assert (fp != NULL);

	int      N    = 20 * size;
	char   * buf  = new char [N];

	fgets (buf, N-1, fp);
	double *vect = stringToDouble (buf, size);
    delete (buf);

    return vect;
}

static double ** readMatrix (const char *file, int nLines, int nColumns)
{
	deb ("readMatrix\n");
	FILE *fp = fopen (file, "r");
	assert (fp != NULL);

	double ** matrix = new double * [nLines];
	int       N      = 20 * nColumns;
	char    * buf    = new char [N];
	int       i;

	for (i = 0; fgets (buf, N-1, fp) != NULL; i++) {
	  matrix[i] = stringToDouble (buf, nColumns);
	}

   //assert (i == nLines);
    delete (buf);
    return matrix;
}

static void  writeMatrix (const char *file, int nLines, int nColumns, double ** matrix)
{
	FILE *fp = fopen (file, "w");
	assert (fp != NULL);

	for (int i = 0; i < nLines; i++) {
		double * line = matrix[i];

		for (int j = 0; j < nColumns; j++)
			fprintf (fp,"%f,", line[j]);

        fprintf (fp, "\n");
	}

	fclose (fp);
}

// --------------------------------------------------------------------------------------

void otest (int nostat, int nsamp, int ntimes, int nstep, int dimX, int dimY, 
	        double mint, double new_freq, int minSampleCount)
{
	deb ("otest\n");
    double **trace      = readMatrix (trace_txt, nostat, minSampleCount);
    double **traveltime = readMatrix (travel_txt, nostat, dimX * dimY) ;
    double * latv       = readVector (latv_txt, dimX * dimY);
    double * lonv       = readVector (lonv_txt, dimX * dimY);

    double ** backSemb = new double * [ntimes];
	
    deb ("x1\n");

	for (int i = 0; i < ntimes; i++) {
	  backSemb[i] = new double [dimX * dimY];
	}

    deb ("x2\n");

    for (int i = 0; i < ntimes; i++) {
        deb_int ("i ", i);

        //  loop over grid points
        double sembmax  = 0, sembmaxX = 0, sembmaxY = 0;

        for (int j = 0; j < dimX * dimY - 1; j++) {
            deb_int ("j ", j);
	        double semb = 0, nomin = 0, denom = 0;

            for (int l = 0; l < nsamp; l++) {
		 //    deb ("x3 ");
               double sum = 0;

               for (int k = 0; k < nostat; k++) {
		           deb ("x4 ");  // Fehler
                   int relstart = (int) ((traveltime[k][j] - mint) * new_freq + 0.5) + i * nstep;
		           deb ("a ");
		           double val = trace[k][relstart + l];

                   sum   += val;
                   denom += (val * val);   			 
		       } 

               nomin += sum * sum;			
	        }
          
            deb ("x5 ");
            semb = nomin / (float (nostat) * denom);
            deb ("x6 ");
            backSemb[i][j] = semb;
            deb ("x7 ");

            if (semb >= sembmax) {
               deb ("x8\n");

               sembmax  = semb;   // search for maximum and position of maximum on semblance grid for given time step	   
               sembmaxX = latv[j];
               sembmaxY = lonv[j];				
	     }

	   } //endfor dimX *dimY

       printf ("\n%3d : max semblance: %f at lat/lon: %f,%f", i, sembmax, sembmaxX, sembmaxY);
	} // endfor ntimes

    writeMatrix (semb_txt, ntimes, dimX*dimY, backSemb);
}

// ------------------------------------------------------------------------------------------------

int main (int argc, char * argv[])
{
  int    nostat,  nsamp,  ntimes,  nstep,  dimX,  dimY;
  double mint, new_freq;
  int    minSampleCount;
  const char * param = NULL;
  bool  isTest = (argc == 1);  
  
  /*
  #ifdef WIN32
    system ("cd ..\\tmpProcess");
  #else
    chdir ("../tmpProcess");
  #endif
  */

  printf ("\nStart OTest");

  if (! isTest) param = argv[1];
  else          param = "6,300,140,50,40,40,719.473999023,10.0,8156";


  printf ("\nParam = %s", param);

  sscanf (param, "%d,%d,%d,%d,%d,%d,%lf,%lf,%d", &nostat, &nsamp, &ntimes, &nstep, &dimX, &dimY, 
	             &mint, &new_freq, &minSampleCount);

  printf ("\nnostat         : %d",    nostat);
  printf ("\nnsamp          : %d",    nsamp);
  printf ("\nntimes         : %d",    ntimes);
  printf ("\nnstep          : %d",    nstep);
  printf ("\n(x,y)          : %d,%d", dimX, dimY); 
  printf ("\nmint           : %lf",	  mint);
  printf ("\nfreq           : %lf",   new_freq);
  printf ("\nminSampleCount : %d",    minSampleCount);

  otest (nostat, nsamp, ntimes, nstep, dimX, dimY, mint, new_freq, minSampleCount);
  printf ("\nEnd  OTest");

  #ifdef WIN32
    if (isTest) getchar();
  #endif

  return 0;
}
