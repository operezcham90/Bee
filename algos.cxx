
#include "algos.h"


#include <iostream>
using namespace std;

/*******************************************************
Compute Quick Sort method in order to increasingly sort a vector

functions:
- quickSort       [main]
- qsort

******************************************************/
CvMat * Pdif;
IplImage *buf_temp;

void opencv_algos::quickSort(CvMat *numbers, int array_size)
{
  q_sort(numbers, 0, array_size - 1);
}

void opencv_algos::q_sort(CvMat* numbers, int left, int right)
{
  double pivot, iesimo;
  int l_hold, r_hold;

  l_hold = left;
  r_hold = right;
  pivot = cvmGet(numbers, left, 0);
  iesimo = cvmGet(numbers, left, 1);
  while (left < right)
  {
    while ((cvmGet(numbers, right, 0) >= pivot) && (left < right))
      right--;
    if (left != right)
    {
      cvmSet(numbers, left, 0, cvmGet(numbers, right, 0));
		cvmSet(numbers, left, 1, cvmGet(numbers, right, 1));
      left++;
    }
    while ((cvmGet(numbers, left, 0) <= pivot) && (left < right))
      left++;
    if (left != right)
    {
		cvmSet(numbers, right, 0, cvmGet(numbers, left, 0));
      cvmSet(numbers, right, 1, cvmGet(numbers, left, 1));
      right--;
    }
  }
  cvmSet(numbers, left, 0, pivot);
  cvmSet(numbers, left, 1, iesimo);
  pivot = (double)left;
  left = l_hold;
  right = r_hold;
  if (left < (int)pivot)
    q_sort(numbers, left, (int)pivot-1);
  if (right > pivot)
    q_sort(numbers, (int)pivot+1, right);
}

/*******************************************************
Compute sobel gradient operator

*********************************************************/

int opencv_algos::sobel(IplImage *buf, int u1, int v1)
{
	double Gx, Gy, res, frac, ent;

	Gx = cvGetReal2D(buf, u1-1, v1+1) + 2*cvGetReal2D(buf, u1, v1+1) + cvGetReal2D(buf, u1+1, v1+1) -
		  cvGetReal2D(buf, u1-1, v1-1) - 2*cvGetReal2D(buf, u1, v1-1) - cvGetReal2D(buf, u1+1, v1-1);
   Gy = cvGetReal2D(buf, u1+1, v1-1) + 2*cvGetReal2D(buf, u1+1, v1) + cvGetReal2D(buf, u1+1, v1+1) -			
		  cvGetReal2D(buf, u1-1, v1-1) - 2*cvGetReal2D(buf, u1-1, v1) - cvGetReal2D(buf, u1-1, v1+1);	
	res = sqrt(Gx*Gx + Gy*Gy);
	//res = fabs(Gx) + fabs(Gy);

   frac = modf(res, &ent);
   res = ent;
   if ((res>=0) && (frac>0.5))
   	res ++;

   return (int)res;   
}


/*******************************************************
Compute the difference histogram in order to speeding up homogeneity texture operator.

Extracted from J.R. Parker's book

functions:
- homogeneidad (2)    [main]
- homo
- histograma

********************************************************/

void opencv_algos::histograma(IplImage *buf, CvMat *Pdif, int u1, int v1, int w, int dist)
{
	int index, k = 0;

	cvSetZero(Pdif);

	if ((u1 >= w) && (u1 < buf->height-w) && (v1 >= w) && (v1 < buf->width-w)) 
   {
		for (int i=-w; i <= w; i++)
   		for (int j = -w; j <= w; j++)
      	{
				if (j+dist <= w)		/* Horizontal */
				{
					index = (int)cvGetReal2D(buf, u1+i, v1+j) - (int)cvGetReal2D(buf, u1+i, v1+j + dist) + 256;
					cvmSet(Pdif, index, 0, cvmGet(Pdif, index, 0) + 1);				
					k ++;
				}
			}

		//Normalize
		for(int i = 0; i < 512; i++)
		{
			cvmSet(Pdif, i, 0, cvmGet(Pdif, i, 0) / k);				
		}

	}

}

float opencv_algos::homo(CvMat *Pdif)
{
	float y = 0.0;

	for (int i= -255; i < 255; i++)
	{
	  y += 1.0/(1.0 + i*i) * cvmGet(Pdif, i+255, 0);
	}
	return y;
}

double opencv_algos::homogeneidad(IplImage *buf, int u1, int v1, int ventana)
{
	float h;
	CvMat * Pdif = cvCreateMat(512, 1, CV_32FC1);

	histograma(buf, Pdif, u1, v1, ventana);
	h = homo(Pdif);
	
	cvReleaseMat(&Pdif);
	return h;
}

void opencv_algos::homogeneidad(IplImage *buf, int vent)
{
	float h;
	cvNamedWindow("Homogeneidad", 0);

	CvMat * Pdif = cvCreateMat(512, 1, CV_32FC1);
	IplImage* H = cvCreateImage(cvSize(buf->width,buf->height), IPL_DEPTH_8U, 1);
	IplImage* gx = cvCreateImage(cvSize(buf->width,buf->height), IPL_DEPTH_32F, 1);
	int ventana = vent;
	int desplazamiento = 2;

	cerr<<endl<<"Procesando ... "<<endl;
	for(int v = 10; v < buf->width-10; v++)
		for(int u = 10; u < buf->height-10; u++)
		{
			h = 0.0;
			for(int desp = 1; desp <= desplazamiento; desp++)
			{
				histograma(buf, Pdif, u, v, ventana, desp);
				h += homo(Pdif);
			}
			cvSetReal2D(gx, u, v, h*255);
		}
	cerr<<"sale";
	cvConvertScaleAbs(gx, H);
	cvShowImage("Homogeneidad", H);
	cvWaitKey(400);
	cout<< "Guardar en archivo? (s / n)"<<endl;
	char c = (unsigned char)getchar();
	if( c == 's' || c == 'S')
	{
		cvSaveImage("./homo.tif", H);
		cout<< "archivo guardado"<<endl;
	}
	
	getchar();
	cvReleaseMat(&Pdif);
	cvReleaseImage(&gx);
	cvReleaseImage(&H);
}

/*********************************************************************
Procedures for estimating the fractal dimension in order to compute texture.

Extracted from J.R. Parker's book

functions:
- fractal    [main]
- Hurst
- MinMax4
- MinMax8
- fit
**********************************************************************/
double opencv_algos::fractalN(IplImage *buf, int u1, int v1)
{
	double x = fractal(buf, u1, v1);
	double k = 3.5;
	double c = 200.0;
//sigmoide
	double fN = 1 / (1 + c*exp(-k*x));

//	cout<<c<<","<<k<<","<<x<<", "<<fN; getchar();
	return fN;
}
double opencv_algos::fractal(IplImage *buf, int u1, int v1)
{ 
	float f, err = 10;

	f = Hurst(buf, u1, v1, &err);			
	//if(err <= 0.15)
	{	
		if(f < 0) 
		{	
			cerr<<"NEGATIVE!!\n"; 
			getchar(); 
			return -1; 
		}
		//if(f >= 0.1 && f <= 3.2)
			return f;	
	}
	return -1;
}

double opencv_algos::fractal(IplImage *buf, float *e, int u1, int v1)
{ 
	float f, err = 10;

	f = Hurst(buf, u1, v1, &err);			
	*e = err;
	if(f < 0) 
	{
		cerr<<"NEGATIVE!!\n"; 
		getchar(); 
		return -1; 
	}
	
	return f;	
}


float opencv_algos::Hurst(IplImage *buf, int u1, int v1, float *herr)
{
	int a, b, i, j;
	float cept, slope = -10.0, err = 10;
	float log_dif_pix[7] = {0,0,0,0,0,0,0};
	static float log_dist[7] = {0.0, 0.34657359, 0.69314718, 0.80471896, 1.0397208, 1.0986123, 1.1512925};

	if((u1 >= 5) && (u1 < buf->height-5) && (v1 >= 5) && (v1 < buf->width-5)) 
	{
		a = b = (int)cvGetReal2D(buf, u1, v1+1);
		MinMax4((int)cvGetReal2D(buf, u1, v1+1), (int)cvGetReal2D(buf, u1, v1-1), 
					 (int)cvGetReal2D(buf, u1-1, v1), (int)cvGetReal2D(buf, u1+1, v1), &a, &b);
		if(a-b) 
			log_dif_pix[0] = (float)log((double)(a - b));

		MinMax4((int)cvGetReal2D(buf, u1-1, v1-1), (int)cvGetReal2D(buf, u1-1, v1+1),
					 (int)cvGetReal2D(buf, u1+1, v1-1), (int)cvGetReal2D(buf, u1+1, v1+1), &a, &b);
		if(a-b) 
			log_dif_pix[1] = (float)log((double)(a - b));

		MinMax4((int)cvGetReal2D(buf, u1-2, v1), (int)cvGetReal2D(buf, u1+2, v1),
					 (int)cvGetReal2D(buf, u1, v1-2), (int)cvGetReal2D(buf, u1, v1+2), &a, &b);
		if(a-b) 
			log_dif_pix[2] = (float)log((double)(a - b));

		MinMax8((int)cvGetReal2D(buf, u1-1, v1+2), (int)cvGetReal2D(buf, u1-1, v1-2),
					 (int)cvGetReal2D(buf, u1+1, v1+2), (int)cvGetReal2D(buf, u1+1, v1-2),
		 			 (int)cvGetReal2D(buf, u1-2, v1+1), (int)cvGetReal2D(buf, u1-2, v1-1),
		 			 (int)cvGetReal2D(buf, u1+2, v1+1), (int)cvGetReal2D(buf, u1+2, v1-1), &a, &b);
		if(a-b) 
			log_dif_pix[3] = (float)log((double)(a - b));

		MinMax4((int)cvGetReal2D(buf, u1+1, v1-2), (int)cvGetReal2D(buf, u1+2, v1+2),
		 			(int)cvGetReal2D(buf, u1-2, v1-2), (int)cvGetReal2D(buf, u1-2, v1+2), &a, &b);
		if(a-b) 
			log_dif_pix[4] = (float)log((double)(a - b));

		MinMax4((int)cvGetReal2D(buf, u1+3, v1), (int)cvGetReal2D(buf, u1-3, v1),
		 			(int)cvGetReal2D(buf, u1, v1-3), (int)cvGetReal2D(buf, u1, v1+3), &a, &b);
		if(a-b) 
			log_dif_pix[5] = (float)log((double)(a - b));

		MinMax8((int)cvGetReal2D(buf, u1+3, v1-1), (int)cvGetReal2D(buf, u1+3, v1+1),
					 (int)cvGetReal2D(buf, u1-3, v1+1), (int)cvGetReal2D(buf, u1-3, v1-1),
		 			 (int)cvGetReal2D(buf, u1-1, v1+3), (int)cvGetReal2D(buf, u1-1, v1-3),
		    	 	 (int)cvGetReal2D(buf, u1+1, v1+3), (int)cvGetReal2D(buf, u1+1, v1-3), &a, &b);
		if(a-b) 
			log_dif_pix[6] = (float)log((double)(a - b));	

		/* Now find the LSBF line to log(dG) and log(distance)	*/

		if(log_dif_pix[0] <= log_dif_pix[1] && log_dif_pix[1] <= log_dif_pix[2] && log_dif_pix[2] <= log_dif_pix[3] &&
			log_dif_pix[3] <= log_dif_pix[4] && log_dif_pix[4] <= log_dif_pix[5] && log_dif_pix[5] <= log_dif_pix[6])
			fit(log_dist, log_dif_pix, 7, &cept, &slope, &err);
		else
		{
			cout<<"no se pudo :s"; 
			for(int l = 0; l < 7; l ++)
				cout<<log_dif_pix[l]<<",";
			getchar();
		}

		*herr = err;
		/*for (i=1; i<7; i++)
		  if (log_dif_pix[i] > log_dif_pix[i-1]) 
				err+=1;
			if (err >=6) 
			{
				for (i=0; i<7; i++)
	  				cerr<<"Error :"<<log_dif_pix[i]<<endl;
				*herr = 100;
			}
		*/
	}

	return slope;
}

int opencv_algos::MinMax4(int a, int b, int c, int d, int *dmax, int *dmin)
{
	int x[4], i;

	x[0] = a; x[1] = b; x[2] = c; x[3] = d;
	for (i=0; i<=3; i++)
	{
	  if (x[i] > *dmax) *dmax = x[i];
	  if (x[i] < *dmin) *dmin = x[i];
	}
}

int opencv_algos::MinMax8(int a, int b, int c, int d, int e, int f, int g, int h, int *dmax, int *dmin)
{
	int x[8], i;

	x[0] = a; x[1] = b; x[2] = c; x[3] = d;
	x[4] = e; x[5] = f; x[6] = g; x[7] = h;

	for (i=0; i<=7; i++)
	{
	  if (x[i] > *dmax) *dmax = x[i];
	  if (x[i] < *dmin) *dmin = x[i];
	}
}

void opencv_algos::fit(float *x, float *y, int ndata, float *a, float *m, float *r)
{
	float s1=0.0, s2=0.0, s3=0.0, xbar = 0.0, ybar = 0.0;
	int i;

	for (i=0; i<ndata; i++)
	{
	  xbar += x[i];
	  ybar += y[i];
	}
	xbar /= (float)ndata;
	ybar /= (float)ndata;

	for (i=0; i<ndata; i++)
	{
	  s1 = (x[i] - xbar);
	  s2 += s1*(y[i]-ybar);
	  s3 += s1*s1;
	}
	*m = s2/s3;
	*a = ybar - *m*xbar;

/* Compute error estimate */
	s1 = 0.0;
	for (i=0; i<ndata; i++)
	  s1 += (y[i] - *m * x[i] - *a) * (y[i] - *m * x[i] - *a);
	*r = (float)sqrt((double)s1);

}

void opencv_algos::fractal(IplImage *buf, IplImage *dest)
{
	float f, err = 0;
	cvNamedWindow("gx", 0);

	IplImage* H = cvCreateImage(cvSize(buf->width,buf->height), IPL_DEPTH_8U, 3);
	IplImage* gx = cvCreateImage(cvSize(buf->width,buf->height), IPL_DEPTH_32F, 3);
	
	for(int v = 10; v < buf->width-10; v++)
		for(int u = 10; u < buf->height-10; u++)
		{
			f = Hurst(buf, u, v, &err);			
			/*if(err > 0.15) //|| (f < 2.0))
				continue;//cerr<< err<<endl;*/
			//else
			{	if(f < 0) 
				{	cerr<<"NEGATIVE!!\n";getchar();}
				/*else if(f * 128 > 255) 
				{	cerr<<"Overflow "<< f <<endl; getchar(); }*/
				else
				{	
					if(f >= 0 && f <= 3.2)
					{ 
						//cvSetReal2D(gx, u, v, 255);
						CvPoint pt;
						pt.x=v;
						pt.y=u;
						if(f >= 0.1 && f <= 0.5)
							cvCircle(gx, pt, 1, CV_RGB(0, 0, 255), -1, 8,0); 
						if(f > 0.5 && f <= 1.0)
							cvCircle(gx, pt, 1, CV_RGB(255, 0,0), -1, 8,0); 
						if(f > 1.0 && f <= 1.5)
							cvCircle(gx, pt, 1, CV_RGB(0, 255, 0), -1, 8,0);
						if(f > 1.5)
							cvCircle(gx, pt, 1, CV_RGB(255, 255, 255), -1, 8,0); 
						
					}
					/*else
					{if(f > 2.5)
						{cerr<<"Mayor a 1.5 ->"<<f<<"--"<<u<<","<<v;
						getchar();}
					}*/
				}
			}
		}
	cerr<<"sale";
	cvConvertScaleAbs(gx, H);
	cvShowImage("gx", H);
	//cvWaitKey(0);
	cerr<<"imagen"<< cvSaveImage("./frac.tif", H);
	dest = cvCloneImage(gx);
	cvReleaseImage(&gx);
	cvReleaseImage(&H);
}

void opencv_algos::printfractalc(IplImage *buf)
{
	float f, err = 0, max = 0.0, maxE = 0.0, min = 100000.0, minE = 100000.0;
	cvNamedWindow("Fractales-Todos", 0);
	cvNamedWindow("Fractales 0 = x <= 0.1", 0);
	cvNamedWindow("Fractales 0.1 <= x <= 0.5", 0);
	cvNamedWindow("Fractales 0.5 < x <= 1.5", 0);
	cvNamedWindow("Fractales 1.5 < x <= 2.5", 0);
	cvNamedWindow("Fractales 2.5 < x <= 3.5", 0);
	cvNamedWindow("Fractales 3.5 < x ", 0);

	IplImage* H = cvCreateImage(cvSize(buf->width,buf->height), IPL_DEPTH_8U, 3);
	IplImage* gx = cvCreateImage(cvSize(buf->width,buf->height), IPL_DEPTH_32F, 3);
	IplImage* gx0 = cvCreateImage(cvSize(buf->width,buf->height), IPL_DEPTH_32F, 3);
	IplImage* gx1 = cvCreateImage(cvSize(buf->width,buf->height), IPL_DEPTH_32F, 3);
	IplImage* gx2 = cvCreateImage(cvSize(buf->width,buf->height), IPL_DEPTH_32F, 3);
	IplImage* gx3 = cvCreateImage(cvSize(buf->width,buf->height), IPL_DEPTH_32F, 3);
	IplImage* gx4 = cvCreateImage(cvSize(buf->width,buf->height), IPL_DEPTH_32F, 3);
	IplImage* gx5 = cvCreateImage(cvSize(buf->width,buf->height), IPL_DEPTH_32F, 3);

	for(int v = 10; v < buf->width-10; v++)
		for(int u = 10; u < buf->height-10; u++)
		{
			f = Hurst(buf, u, v, &err);			
			if(err > 0.05) //|| (f < 2.0))
				continue;//cerr<< err<<endl;
			//else
			{	if(f < 0) 
				{	cerr<<"NEGATIVE!!\n";getchar();}
				/*else if(f * 128 > 255) 
				{	cerr<<"Overflow "<< f <<endl; getchar(); }*/
				else
				{	
					if(max < f)
						max = f;
					if(maxE < err)
						maxE = err;
					if(min > f)
						min = f;
					if(minE > err)
						minE = err;
					//if(f >= 0 && f <= 3.2)
					{ 
						CvPoint pt;
						pt.x=v;
						pt.y=u;
						if(f >= 0.0 && f <= 0.1)
						{
							
							cvCircle(gx, pt, 1, CV_RGB(0, 255, 255), -1, CV_AA,0); 
							cvCircle(gx0, pt, 1, CV_RGB(0, 255, 255), -1, CV_AA,0); 
						}
						else if(f >= 0.1 && f <= 0.5)
						{
							cvCircle(gx, pt, 1, CV_RGB(0, 0, 255), -1, 8,0); 
							cvCircle(gx1, pt, 1, CV_RGB(0, 0, 255), -1, 8,0); 
						}
						else if(f > 0.5 && f <= 1.5)
						{
							cvCircle(gx, pt, 1, CV_RGB(255, 0,0), -1, 8,0); 
							cvCircle(gx2, pt, 1, CV_RGB(255, 0, 0), -1, 8,0); 
						}
						else if(f > 1.5 && f <= 2.5)
						{
							cvCircle(gx, pt, 1, CV_RGB(0, 255, 0), -1, CV_AA,0);
							cvCircle(gx3, pt, 1, CV_RGB(0, 255, 0), -1, CV_AA,0); 
						}
						else if(f > 2.5 && f <= 3.5)
						{
							cvCircle(gx, pt, 1, CV_RGB(255, 255, 0), -1, 8,0); 
							cvCircle(gx4, pt, 1, CV_RGB(255, 255, 0), -1, 8,0); 
						}
						else if(f > 3.5)
						{
							cvCircle(gx, pt, 1, CV_RGB(255, 255, 255), -1, 8,0); 
							cvCircle(gx5, pt, 1, CV_RGB(255, 255, 255), -1, 8,0); 
						}
						
					}
					
				}
			}
		}
	cout<<"Máximo valor: "<<max<<" - Mínimo Valor: "<<min<<endl<<"Máximo error: "<<maxE<<" - Mínimo Valor: "<<minE<<endl;
	cvConvertScaleAbs(gx, H);
	cvShowImage("Fractales-Todos", H);
	cvConvertScaleAbs(gx0, H);
	cvShowImage("Fractales 0 = x <= 0.1", H);
	cvConvertScaleAbs(gx1, H);
	cvShowImage("Fractales 0.1 <= x <= 0.5", H);
	cvConvertScaleAbs(gx2, H);
	cvShowImage("Fractales 0.5 < x <= 1.5", H);
   cvConvertScaleAbs(gx3, H);
	cvShowImage("Fractales 1.5 < x <= 2.5", H);
   cvConvertScaleAbs(gx4, H);
	cvShowImage("Fractales 2.5 < x <= 3.5", H);
   cvConvertScaleAbs(gx5, H);
	cvShowImage("Fractales 3.5 < x ", H);
	

	cvReleaseImage(&gx);
	cvReleaseImage(&gx0);
	cvReleaseImage(&gx1);
	cvReleaseImage(&gx2);
	cvReleaseImage(&gx3);
	cvReleaseImage(&gx4);
	cvReleaseImage(&gx5);
	cvReleaseImage(&H);
}

void opencv_algos::printfractal(IplImage *buf)
{
	int slide = 1, d;
	cvNamedWindow("gx", 1);

	cvCreateTrackbar( "bar_gx", "gx", &slide, 10, on_change);

	buf_temp = buf;
	re_escalar(slide);
	
}


void opencv_algos::on_change(int v)
{
	re_escalar(cvGetTrackbarPos("bar_gx", "gx"));
}


void opencv_algos::re_escalar(int pos)
{
	float f, err = 0;
	IplImage* H = cvCreateImage(cvSize(buf_temp->width,buf_temp->height), IPL_DEPTH_8U, 1);
	IplImage* gx = cvCreateImage(cvSize(buf_temp->width,buf_temp->height), IPL_DEPTH_32F, 1);
	
	//cout<<buf_temp->width<<","<<buf_temp->height;getchar();
	
	for(int v = 10; v < buf_temp->width; v++)
		for(int u = 10; u < buf_temp->height; u++)
		{
			//cout<<cvGetReal2D(buf_temp, u, v);getchar();
			f = Hurst(buf_temp, u, v, &err);			
			/*if(err > 0.15) //|| (f < 2.0))
				continue;//cerr<< err<<endl;*/
			//else
			{	if(f < 0) 
				{	//cerr<<"NEGATIVE!!\n";getchar();
				}
				/*else if(f * 128 > 255) 
				{	cerr<<"Overflow "<< f <<endl; getchar(); }*/
				else
				{	
					if(f >= 0 && f <= 3.2)
					{ 
						//CvScalar val;
						CvPoint pt;
						pt.x=v;
						pt.y=u;
						//val = (int) (f * 255 / slide);
						
						cvCircle(gx, pt, 1, cvScalarAll((f * 255.0 / (double)pos)), -1, 8,0); 
						
					}
					/*else
					{if(f > 2.5)
						{cerr<<"Mayor a 1.5 ->"<<f<<"--"<<u<<","<<v;
						getchar();}
					}*/
				}
			}
		}
	cvConvertScaleAbs(gx, H);
	cvShowImage("gx", H);
	//cvWaitKey(0);
	cvReleaseImage(&gx);
	cvReleaseImage(&H);

}


/*********************************************************************
Procedure for estimate a simple analysis of texture by the means of first 
order statistics; using variance estimates

"Coefficient of variation"
- c_v   
"Roughness"
- rough   

auxiliar function

- compute_mean_var
**********************************************************************/

void opencv_algos::compute_mean_var(IplImage *buf_F1, int u, int v, double *m, double *var)
{
	double media = 0.0;
	double varianza = 0.0;
	int paso_vent = 1;

	for (int k=-paso_vent;k<=paso_vent;k++)
   {
		for (int l=-paso_vent;l<=paso_vent;l++)
			media +=(int)cvGetReal2D(buf_F1, u + l, v + k);
   }
   media /= ((paso_vent*2+1) * (paso_vent*2+1));
	for (int k=-paso_vent;k<=paso_vent;k++)
   {
		for (int l=-paso_vent;l<=paso_vent;l++)
			varianza += ((int)cvGetReal2D(buf_F1, u + l, v + k) - media) * ((int)cvGetReal2D(buf_F1, u + l, v + k) - media);
   }
	varianza /= ((paso_vent*2+1) * (paso_vent*2+1));

	//cout<< media << ", "<< varianza <<endl; //getchar();
	*m = media;
	*var = varianza;
}


double opencv_algos::c_v(IplImage *buf_F1, int u, int v)
{
	double var, media, CV;
	
	opencv_algos::compute_mean_var(buf_F1, u, v, &media, & var);
	//////////  Coeficiente de variación ////////////
	CV = sqrt(var) / media;
	
	//cout<< "Coeficiente de variación:" << CV <<endl; //getchar();
	return CV;
}

double opencv_algos::rough(IplImage *buf_F1, int u, int v)
{
	double var, media, R;
	
	opencv_algos::compute_mean_var(buf_F1, u, v, &media, & var);
	//////////  Rugosidad ////////////
	R = 1.0 - 1.0 /(1 + var);
	
	//cout<< "Rugosidad:" << R <<endl; //getchar();
	return R;
}
