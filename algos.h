#ifndef algos_h_
#define algos_h_

#include "cv.h"
#include "highgui.h"

#include <stdio.h>

class opencv_algos
{
public:
   static void quickSort(CvMat* numbers, int array_size);
   static void q_sort(CvMat * numbers, int left, int right);

	static int sobel(IplImage *, int , int);

	static void histograma(IplImage *, CvMat *, int, int, int, int = 1);
	static float homo(CvMat *);
   static void homogeneidad(IplImage *, int = 5);
	static double homogeneidad(IplImage *, int, int, int = 5);

	static void printfractal(IplImage *);
	static void printfractalc(IplImage *);
	static void on_change(int);
	static void re_escalar(int);

	static void fractal(IplImage *, IplImage *);
	static double fractal(IplImage *, int, int);
	static double fractal(IplImage *, float*, int, int);
	static double fractalN(IplImage *, int, int);
	static float Hurst(IplImage *, int, int, float *);
	static int MinMax4(int, int, int, int, int *, int *);
	static int MinMax8(int, int, int, int, int, int, int, int, int *, int *);
	static void fit(float *, float *, int, float *, float *, float *);

	static double c_v(IplImage *, int, int);
	static double rough(IplImage *, int, int);
	static void compute_mean_var(IplImage *, int, int, double *, double *);
};

#endif

