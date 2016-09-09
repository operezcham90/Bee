
#include "abejas_foto.h"
#include "algos.h"
//#define _LINUX_TIME_H 1	/* to get things compile on kernel 2.6.x */
//#include <linux/videodev.h>
#include <libv4l1-videodev.h>
#include <sys/ioctl.h>
#include <fcntl.h>
#include <errno.h>

/*************** variables globales *********************/

CvMat *espacios_lower;
CvMat *espacios_upper;
CvMat *sitios;

IplImage *buf_F1 = NULL;
IplImage *buf_F2 = NULL;
IplImage *Im1 = NULL;
IplImage *Im2 = NULL;
CvMat* Mrotacion_izq;
CvMat* Mrotacion_der;
CvMat* Vrotacion_izq;
CvMat* Vrotacion_der;
CvMat* Vtraslacion_izq;
CvMat* Vtraslacion_der;
CvMat* MIntrinsecos_izq;
CvMat* MIntrinsecos_der;
int inicio = 0;
Machine m;
int stage = 0; // modificacion
int return_to_inactive = 0; // modificacion
int bad_bees = 0;
double r_sharing = 0; // modificacion ratio sharing
double rate_dif = 0;
double current_rate_cross = 0;
double current_rate_mut = 0;
int pop_dif = 0; // modificacion
int real_pop_size = 0; // modificacion
ofstream result_file;

CvMat* x1;
CvMat* x2;
CvRNG rng = cvRNG(0xffffffff);

/****************************************************************/

#define PWC_FPS_SHIFT		16
#define PWC_FPS_MASK		0x00FF0000
#define PWC_FPS_FRMASK		0x003F0000


void show_bees()
{
	for (int i = 0; i < 6; i++)
	{
		//cout << m.getStateName(i) << ": " << m.countBeesInState(i) << endl;
	}
}

int main(int argc, char** argv)
{
  string config_param;
  string nombreImg;
  string nombreImgDER;
  string nombreParam;
  string nombreParamDER;
  string nombreParamCEN;
  char* inicio = argc == 2 ? argv[1] : (char*)"./archivos.ini";

  ifstream in(inicio);
  string cadena;

  if (!in)
  {
     cerr<<"Error al leer archivo de parametros de calibracion"<<endl;
     return 0;
  }

  /* * lectura de archivo que contiene nombres de imágenes y archivos de parámetros **/
  in >> nombreImg;
  cout<<nombreImg<<endl;
  in >> nombreImgDER;
  cout<<nombreImgDER<<endl;
  in >> nombreParam;
  cout<<nombreParam<<endl;
  in >> nombreParamDER;
  cout<<nombreParamDER<<endl;
  in >> nombreParamCEN;
  cout<<nombreParamCEN<<endl;
  in >> config_param;
  cout<<config_param<<endl;

  in.close();

//////////////////////// Loading Projection Matrices: ////////////////
  Mrotacion_izq = cvCreateMat(3, 3, CV_32FC1);
  Mrotacion_der = cvCreateMat(3, 3, CV_32FC1);
  Vrotacion_izq = cvCreateMat(1, 3, CV_32FC1);
  Vrotacion_der = cvCreateMat(1, 3, CV_32FC1);
  Vtraslacion_izq = cvCreateMat(1, 3, CV_32FC1);
  Vtraslacion_der = cvCreateMat(1, 3, CV_32FC1);
  MIntrinsecos_izq = cvCreateMat(3, 3, CV_32FC1);
  MIntrinsecos_der = cvCreateMat(3, 3, CV_32FC1);

  ifstream fpr(nombreParam.c_str());
  if(! fpr.is_open())
  {
    cerr<< "El archivo "<< nombreParam << " NO existe, operacion abortada"<<endl; return 0;
  }

/*------------------------------------------------------------------------
   LECTURA DE LOS PARAMETROS EXTRINSECOS Izquierda
-------------------------------------------------------------------------*/
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
    {
      float tmp;
      fpr >> tmp;
      cvmSet(Mrotacion_izq, i, j, tmp);
    }
  cvRodrigues2(Mrotacion_izq, Vrotacion_izq);
  for (int i=0;i<3;i++)
    fpr >> Vtraslacion_izq->data.fl[i];

/*------------------------------------------------------------------------
   LECTURA DE LOS PARAMETROS INTRINSECOS Izquierda
-------------------------------------------------------------------------*/
  float tmp;
  cvSetZero(MIntrinsecos_izq);
   fpr >> tmp;
   cvmSet(MIntrinsecos_izq, 0, 0, tmp);
   fpr >> tmp;
  cvmSet(MIntrinsecos_izq, 1, 1, tmp);
   fpr >> tmp;
   cvmSet(MIntrinsecos_izq, 0, 2, tmp);
   fpr >> tmp;
  cvmSet(MIntrinsecos_izq, 1, 2, tmp);
  cvmSet(MIntrinsecos_izq, 2, 2, 1.0);
   fpr.close();

   fpr.open(nombreParamDER.c_str());
  if(! fpr.is_open())
  {
    cout << "El archivo NO existe, operacion abortada"<<endl; return 0;
  }

/*------------------------------------------------------------------------
   LECTURA DE LOS PARAMETROS EXTRINSECOS Derecha
-------------------------------------------------------------------------*/
  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
    {
      float tmp;
      fpr >> tmp;
      cvmSet(Mrotacion_der, i, j, tmp);
    }
  cvRodrigues2(Mrotacion_der, Vrotacion_der);
  for (int i=0;i<3;i++)
  	 fpr >> Vtraslacion_der->data.fl[i];

/*------------------------------------------------------------------------
   LECTURA DE LOS PARAMETROS INTRINSECOS Derecha
-------------------------------------------------------------------------*/
  cvSetZero(MIntrinsecos_der);
   fpr >> tmp;
  cvmSet(MIntrinsecos_der, 0, 0, tmp);
   fpr >> tmp;
  cvmSet(MIntrinsecos_der, 1, 1, tmp);
   fpr >> tmp;
  cvmSet(MIntrinsecos_der, 0, 2, tmp);
   fpr >> tmp;
  cvmSet(MIntrinsecos_der, 1, 2, tmp);
  cvmSet(MIntrinsecos_der, 2, 2, 1.0);
   fpr.close();

  opencv_abejas::function_main(config_param, nombreImg, nombreImgDER);

  cvReleaseMat(&Mrotacion_izq);
  cvReleaseMat(&Mrotacion_der);
  cvReleaseMat(&Vrotacion_izq);
  cvReleaseMat(&Vrotacion_der);
  cvReleaseMat(&Vtraslacion_izq);
  cvReleaseMat(&Vtraslacion_der);
  cvReleaseMat(&MIntrinsecos_izq);
  cvReleaseMat(&MIntrinsecos_der);

  cvReleaseImage(&buf_F1);
  cvReleaseImage(&buf_F2);

  return 0;
}
/////////////////////////////////////////////////////////////////////////////

void opencv_abejas::function_main(string arch_config, string nombreImg, string nombreImgDER)
{
//se lee el archivo de parametros del algoritmo

  if(!input_param_file(arch_config))
      return;

//se prepara la estructura para la población inicial de abejas
   nextpop = new INDIVIDUAL[pop_size];
   if (nextpop == NULL) nomemory("nextpop in function_main()");

  espacios_lower = cvCreateMat(pop_size, nvar_real, CV_32FC1);
  espacios_upper = cvCreateMat(pop_size, nvar_real, CV_32FC1);
  sitios = cvCreateMat(pop_size, 1, CV_32FC1);


   int R_pop_size = pop_size;
   int R_pop_size_off = pop_size_off;
   double R_rate_cross = rate_cross;
   double R_rate_mut = rate_mut;
   double R_rate_rand = rate_rand;
   double sigma_share_orig = sigma_share;

// add bees
for (int i = 0; i < pop_size + pop_size_for; i++) {
	m.addBee();
}
cout << endl;
show_bees(); // modificacion

   //establecer tamaño del espacio de búsqueda
   for (int k=0; k < nvar_real; k++)
   {
    xreal_lower_orig[k] = xreal_lower[k];
    xreal_upper_orig[k] = xreal_upper[k];
   }

  select_memory(); // revisa que la poblacion total sea mayor a la que se ocupa para el torneo

   time_t tiempo;
   struct tm *tmPtr;
  struct timeb tmb1, tmb2;

  MINM = -1; // use -1 for maximization  ; 1 for minimization

//////// Creación de ventanas ////////////

   cvNamedWindow("izquierda", 1);
   cvResizeWindow("izquierda", 768, 484);
   cvNamedWindow("derecha", 1);
   cvResizeWindow("derecha", 768, 484);
   cvWaitKey(300);

//comienza ciclo global
  do
   {

    ///////  Loading Image / Video ///////////
    buf_F1 = cvLoadImage(nombreImg.c_str(),-1);
    buf_F2 = cvLoadImage(nombreImgDER.c_str(),-1);

    cvShowImage("izquierda",buf_F1);
    cvShowImage("derecha",buf_F2);

    cvWaitKey(2);
    int maxX = (int)buf_F1->width;
    int maxY = (int)buf_F1->height;

    // modificacion
    // calcular proporsion de sigma
    sigma_share = sigma_share_orig;
    double vol = 1;
    for (int k=0; k < nvar_real; k++)
    {
        double lk = abs(xreal_lower_orig[k] - xreal_upper_orig[k]);
        vol *= lk;
    }
    //cout << "-- volume of search space --" << vol << endl;
    double edge = cbrt(vol);
    //cout << "-- edge of search space --" << edge << endl;
    r_sharing = sigma_share / edge;
    //cout << "-- r_sharing --" << sigma_share << "/" << edge << " = " << r_sharing << endl;

/**********************    E T A P A   D E  E X P L O R A C I O N *************************/
    gen_no = 0;
    stage = 0; // modificacion
    return_to_inactive = 0;
    initialize();   //inicializacion de la poblacion y otras var. globales
    //cout << "-- initialize() exploration --" << inicio << endl;
    return_to_inactive = 1;
    show_bees(); // modificacion

    tiempo = time(NULL);
    tmPtr = localtime(&tiempo);
    cerr<< "Inicio: "<< asctime(tmPtr) << endl;       //REGISTRA TIEMPO DE INICIO
    ftime(&tmb1);
    for(int s = 0; s < pop_size; s++)         //evaluación población inicial
      {
         fitness_function_exp(&oldpop[s], maxX, maxY, buf_F1, buf_F2);
      }
    rate_dif = 0; // modificacion
    pop_dif = 0; // modificacion
    real_pop_size = pop_size; // modificacion
    verifica_pop();
    for(gen_no = 1; gen_no <= max_gen; gen_no++)
       {
         if (gen_no > 1) // modificacion
         {
            //cout << "-- exploration new generation --" << endl;
            // en este punto pop_size podria haber cambiado
            for (int i = 0; i < pop_size; i++)
            {
                oldpop[i].state = m.giveFeedback(4, 0); // exp exitosa -> inactividad
                oldpop[i].state = m.giveFeedback(0, 0); // inactividad -> exploracion
            }
            show_bees(); // modificacion
         }

         // modificacion en este punto pop_size podria haber cambiado
         rate_cross = R_rate_cross;
         rate_mut = R_rate_mut;
         rate_rand = R_rate_rand;
         current_rate_cross = rate_cross;
         current_rate_mut = rate_mut;
         verifica_pop();

         numero_cross = numero_mut = numero_otros = 0;
         //GENERATION OF NEW POPULATION through SELECTION, XOVER , MUTATION & RANDOM
	 // and PROJECTION into the images
         generate_new_pop();

         for(int i = 0; i < pop_size_off; i++)     //evaluación nueva población
         {
            fitness_function_exp(&newpop[i], maxX, maxY, buf_F1, buf_F2);
         }
         
/*for (int i=0; i<pop_size_off; i++)
	cout << newpop[i].xreal[0] << " " << newpop[i].xreal[1]  << " " << newpop[i].xreal[2]  << "\t" << newpop[i].obj 
		<< "\t" << newpop[i].u1 << " " << newpop[i].v1 << " " << newpop[i].u2 << " " << newpop[i].v2 
		<< "\t" << newpop[i].type << "\t" << newpop[i].nicho << endl;
getchar();*/
         
         //juntar poblaciones:  mu + lamda (padres e hijos)
         merge_pop();
	//sharing
         if(sigma_share != 0)
         {
            sharing3D();
            //cout << "-- sharing3D explorarion -- " << bad_bees << "/" << (pop_size + pop_size_off) << " = " << ((bad_bees*100)/(pop_size + pop_size_off)) << "%" << endl;
         }
         else
         {
            sharing2D();
         }
         best_mu();             // ordenar poblacion de acuerdo al fitness y obtener los mu mejores
         //cout << "-- best_mu() exploration --" << endl;
         show_bees(); // modificacion
       }
       /* Registrar tiempo final sin considerar tiempo de impresion de abejas en pantalla*/
      ftime(&tmb2);
      int t_diff = (int) (1000.0 * (tmb2.time - tmb1.time) + (tmb2.millitm - tmb1.millitm));
      tiempo = time(NULL);
      tmPtr = localtime(&tiempo);
      cerr<< "Etapa de Exploración. Milisegundos:"<< t_diff<<endl;

      // Escribir los resultados
      //result_file.open("resultados.txt", fstream::app);
      //result_file << "===========RESULTADOS============" << endl;
      //result_file << "SESION: " << asctime(tmPtr);
      //result_file << "=================================" << endl;
      //result_file << "u1, v1, u2, v2, x, y, z, obj" << endl;
      //result_file.close();
      imprime_abejas();
      guarda_exploradoras();

      /**********************  E T A P A   D E    R E C L U T A M I E N T O *************************/
      ftime(&tmb1);

// modificacion
// agregar las abejas que sobran de exploracion
pop_size_for += pop_dif;
//cout << "-- sobraron " << pop_dif << " abejas exploradoras pop_size_for: " << pop_size_for << endl;

      asigna_recolectoras();
      int sit = pop_size;

      stage = 1;
      //cout << "-- recruit --" << endl;
      for (int i = 0; i < sit; i++)
         oldpop[i].state = m.giveFeedback(4, 1); // exp exitosa -> recluta
      show_bees(); // modificacion

      ftime(&tmb2);
      t_diff = (int) (1000.0 * (tmb2.time - tmb1.time) + (tmb2.millitm - tmb1.millitm));
      cerr<< "Etapa de Reclutamiento. Milisegundos:"<< t_diff<<endl;

	  /***********
	  for (int i=0; i<pop_size_off; i++)
		cout << newpop[i].xreal[0] << " " << newpop[i].xreal[1]  << " " << newpop[i].xreal[2]  << "\t" << newpop[i].obj 
			 << "\t" << newpop[i].u1 << " " << newpop[i].v1 << " " << newpop[i].u2 << " " << newpop[i].v2 
			 << "\t" << newpop[i].type << "\t" << newpop[i].nicho << endl;
	  break;
	  /************/
	 
      /**********************  E T A P A   D E     R E C O L E C C I O N *************************/
    inicio = 0;
    stage = 2; // modificacion
    ftime(&tmb1);

    for(int i = sit - 1; i >= 0; i--)   // por cada sitio que indicó la abeja exploradora:
    {
          if(CV_MAT_ELEM(*sitios, float, i, 0) > 0)
          {
	    //prepara las estructuras de datos para la poblacion de abejas
	    free_all();
	    pop_size = (int)CV_MAT_ELEM(*sitios, float, i, 0);
   	    pop_size_off = (int)CV_MAT_ELEM(*sitios, float, i, 0)*2;
	    rate_cross = rate_cross_rec;
	    rate_mut = rate_mut_rec;
	    rate_rand = rate_rand_rec;

// modificacion
// agregar las abejas que sobran de la recoleccion anterior
pop_size += pop_dif;
//cout << "-- sobraron " << pop_dif << " abejas recolectoras pop_size: " << pop_size << endl;

            rate_dif = 0; // modificacion
            pop_dif = 0; // modificacion
            real_pop_size = pop_size; // modificacion
   	    verifica_pop();
	    for(int k = 0; k < nvar_real; k++)
	    {
	      xreal_lower[k] = CV_MAT_ELEM(*espacios_lower, float, i, k);
              xreal_upper[k] = CV_MAT_ELEM(*espacios_upper, float, i, k);
	    }
            // modificacion
            // calcular proporsion de sigma
            sigma_share = sigma_share_orig;
            vol = 1;
            for (int k=0; k < nvar_real; k++)
            {
               double lk = abs(xreal_lower[k] - xreal_upper[k]);
               vol *= lk;
            }
            //cout << "-- volume of search space [foraging] --" << vol << endl;
            edge = cbrt(vol);
            //cout << "-- edge of search space --" << edge << endl;
            sigma_share = edge * r_sharing;
            //cout << "-- sigma_share --" << edge << "*" << r_sharing << " = " << sigma_share << endl;

	    //inicializacion de la poblacion y otras var. globales
	    initialize();
            //cout << "-- initialize() foraging --" << inicio << endl;
            show_bees(); // modificacion


	    for(int s = 0; s < pop_size; s++)    //evaluación población inicial
	    {
	      fitness_function_rec(&oldpop[s], maxX, maxY, buf_F1, buf_F2);
	    }
	    for(gen_no = 1; gen_no <= max_gen/2; gen_no++) //numero de generacion es la mitad de la etapa de exploracion
	    {

	       if (gen_no > 1) // modificacion
	       {
	          //cout << "-- foraging new generation --" << endl;
	          // en este punto pop_size podria haber cambiado
	          for (int i = 0; i < pop_size; i++)
	          {
	             oldpop[i].state = m.giveFeedback(5, 0); // rec exitosa -> inactividad
	             oldpop[i].state = m.giveFeedback(0, 1); // inactividad -> recoleccion
	          }
	          show_bees(); // modificacion
	       }

              // modificacion en este punto pop_size podria haber cambiado
              rate_cross = rate_cross_rec;
              rate_mut = rate_mut_rec;
              rate_rand = rate_rand_rec;
              current_rate_cross = rate_cross;
              current_rate_mut = rate_mut;
              verifica_pop();

	      numero_cross = numero_mut = numero_otros = 0;
	      //GENERATION OF NEW POPULATION through SELECTION, XOVER, MUTATION &RANDOM
	      generate_new_pop();
	      //evaluacion
	      for(int i = 0; i < pop_size_off; i++)
	      {
		fitness_function_rec(&newpop[i], maxX, maxY, buf_F1, buf_F2);
	      }
	      //juntar poblaciones:  mu + lamda
   	      merge_pop();
      	      //sharing
	      if(sigma_share != 0)
	      {
		sharing3D();
                //cout << "-- sharing3D foraging -- " << bad_bees << "/" << (pop_size + pop_size_off) << " = " << ((bad_bees*100)/(pop_size + pop_size_off)) << "%" << endl;
	      }
	      else
   	      {
      	        sharing2D();
	      }
	       // ordenar poblacion de acuerdo al fitness y obtener los mu mejores
   	      best_mu();
              //cout << "-- best_mu() foraging --" << endl;
              show_bees(); // modificacion
	    }
      		//captura_abejas();
	    imprime_abejas();
	  }
       }
       // registro del tiempo en que se realiza la etapa de recolección
      ftime(&tmb2);
      t_diff = (int) (1000.0 * (tmb2.time - tmb1.time) + (tmb2.millitm - tmb1.millitm));
      cerr<< "Etapa de Recolección: "<< "Milisegundos:"<< t_diff<<endl;

      // Se regresan las variables del espacio de búsqueda a sus valores originales,
      // listos para la siguiente iteración
      free_all();
      pop_size = R_pop_size;
      pop_size_off = R_pop_size_off;
      rate_cross = R_rate_cross;
      rate_mut = R_rate_mut;
      rate_rand = R_rate_rand;
      sigma_share = sigma_share_orig;
      for (int k=0; k < nvar_real; k++)
      {
   	xreal_lower[k] = xreal_lower_orig[k];
    	xreal_upper[k] = xreal_upper_orig[k];
      }
      //inicio = 1; // modificacion, es correcto?
  }while(1);
  delete[] nextpop;
  cvReleaseMat(&sitios);
  cvReleaseMat(&espacios_lower);
  cvReleaseMat(&espacios_upper);
  cvDestroyWindow("izquierda");
   cvDestroyWindow("derecha");

}


/*====================================================================
SUBROUTINE FOR INPUTTING GLOBAL PARAMETERS -- FILE:
====================================================================*/
int opencv_abejas::input_param_file(string archivo)
{
   ifstream inI;
   int k, num_arch;
   char ans;
   string una_cad;
  cout<< archivo <<" yiyo"<<endl;
   inI.open(archivo.c_str());
   getchar();
   /***************************************************************
   inI.open(una_cad.c_str());

   **************************************************************/
	if(! inI.is_open())
  	{
		//free(config_param);
   	cout << "El archivo "<< archivo << "NO existe, operacion abortada"<<endl; return 0;
  	}
   parametros = archivo;
   parametros.erase(parametros.find("."));
   parametros.erase(0, parametros.find_last_of("/")+1);
   parametros.insert(0, "./imagenesGeneradas/abejas_");
   cout<<"Leyendo archivo..."<<endl;
   //Cuantas generationes
   inI>>max_gen;
   cout<<"Cuantas generaciones "<<max_gen<<endl;
   // Tamanio de poblacion
   inI>>pop_size;
   cout<<"tam población padres "<<pop_size<<endl;
   if (pop_size > MAXPOPSIZE)
   {
        cerr<<"\n Increase the value of MAXPOPSIZE in program  and re-run the program"<<endl;
        exit(-1);
   }
   //Numero de variables binarias
   inI>>nvar_bin;
   cout<<"bin "<<nvar_bin<<endl;
   //Numero de variables reales
   inI>>nvar_real;
   cout<<"real "<<nvar_real<<endl;
   if (nvar_bin > 0)
     for (k=0, lchrom = 0; k < nvar_bin; k++)
       {
          //String length and Lower and Upper bounds of xbin
          inI>>chr_len[k]>>xbin_lower[k]>>xbin_upper[k];
          cout<<"\nString"<<chr_len[k]<<xbin_lower[k]<<xbin_upper[k]<<endl;
         lchrom += chr_len[k];
       }
   if (nvar_real > 0)
     for (k=0; k < nvar_real; k++)
       {
         inI>>xreal_lower[k]>>xreal_upper[k];
         cout<<"\nxreal "<<xreal_lower[k]<<","<<xreal_upper[k]<<endl;
       }
   if (nvar_real > 0)
   {
      inI>>ans;
      cout<<"rigid ? "<<ans<<endl;
     if (ans == 'y')
         RIGID = TRUE;
     else
         RIGID = FALSE;
   }
   inI>>ans;
   cout<<"3D sharing ? "<<ans<<endl;
   if (ans == 'y')
   {
      SHARING3D = TRUE;
      inI>>sigma_share;
      cout<<"niching val "<<sigma_share<<endl;
   }
   else
      SHARING3D = FALSE;
   inI>>ans;
   cout<<"2D sharing ? "<<ans<<endl;
   if (ans == 'y')
   {
      SHARING2D = TRUE;
      inI>>alpha_share;
      cout<<"alpha val "<<alpha_share<<endl;
   }
   else
      SHARING2D = FALSE;
   inI>>ans;
   cout<<"printed ? "<<ans<<endl;
   if (ans == 'y')
      REPORT = TRUE;
   else
      REPORT = FALSE;
   inI>>maxrun;
   cout<<"runs "<<maxrun<<endl;
   tourneysize=2;

   inI>>p_xover;
   cout<<"CrossOver Prob "<<p_xover<<endl;
   if (nvar_bin > 0)
     {
        inI>>p_mutation_bin;
        cout<<"Mutation Prob binar "<<p_mutation_bin<<endl;
     }
   if (nvar_real > 0)
     {
        inI>>p_mutation_real;
        cout<<"Mutation Prob real "<<p_mutation_real<<endl;
     }
   if (nvar_real > 0)
     {
        inI>>ans;
        cout<<"Poly mut? "<<ans<<endl;
         if (ans == 'y')
         {
            mut_type = 1;
            inI>>n_distribution_m;
            cout<<"eta for Poly"<<n_distribution_m<<endl;
         }
        inI>>ans;

        cout<<"Normal mut? "<<ans<<endl;
         if (ans == 'y')
        {    for (k=0; k < nvar_real; k++)
             {
               inI>>sigma[k];
               cout<<"sigma"<<k<<" "<<sigma[k]<<endl;
             }
             mut_type = 0;
        }
        inI>>ans;
        cout<<"SBX Cross? "<<ans<<endl;
         if (ans == 'y')
        {
           cross_type = 1;
            inI>>n_distribution_c;
           cout<<"eta for SBX"<<n_distribution_c<<endl;
        }
        inI>>ans;
        cout<<"Barycentric Cross? "<<ans<<endl;
         if (ans == 'y')
        {
           cross_type = 0;
           cout<<"Baryc. to be used"<<endl;
        }
     }
   inI>>basic_seed;
   cout<<"rand seed "<<basic_seed<<endl;
   inI>>pop_size_off;
   cout<<"tam poblacion hijos"<<pop_size_off<<endl;
   inI>>rate_cross;
   cout<<"pobl. en cross "<<rate_cross<<endl;
   inI>>rate_mut;
   cout<<"pobl. en mut "<<rate_mut<<endl;
   inI>>rate_rand;
   cout<<"poblacion aleatoria "<<rate_rand<<endl;
   inI>>ans;
   cout<<"gradiente o sobel ? "<<ans<<endl;
   if (ans == 'g')
   {
      inI>>archivos[0];
      cout<<"archivo 1"<<archivos[0].c_str()<<endl;
      inI>>archivos[1];
      cout<<"archivo 2"<<archivos[1].c_str()<<endl;
      grad = 0;
   }
   else if(ans == 's')
   {
      inI>>archivos[0];
      cout<<"archivo 1"<<archivos[0].c_str()<<endl;
      inI>>archivos[1];
      cout<<"archivo 2"<<archivos[1].c_str()<<endl;
      inI>>archivos[2];
      cout<<"archivo 3"<<archivos[2].c_str()<<endl;
      inI>>archivos[3];
      cout<<"archivo 4"<<archivos[3].c_str()<<endl;
      grad = 1;
   }
   else
   {
      cout<<"Error al elegir el tipo de gradiente. Funcion abortada. Revise archivo de parametros"<<endl;
      return 0;
   }
   string correlac;
   inI>>correlac;
   if(correlac == "ZNCC")
      corre = 1;
   else if(correlac == "SSD")
      corre = 0;
   else
   {
      cout<<"Error al elegir el tipo de correlacion. Funcion abortada. Revise archivo de parametros"<<endl;
      return 0;
   }
   inI >> num_abejas;
   inI >> pop_size_for;
   cout<<"Poblacion de abejas recolectoras: "<<pop_size_for<<endl;
   if(num_abejas <= pop_size_for)
      cout<<"Numero de abejas que se imprimiran del total: "<<num_abejas<<endl;
   else
   {
      num_abejas = pop_size;
      cout<<"WARNING: Numero de abejas mayor al de la poblacion inicial"<<endl;
   }
   inI>>rate_cross_rec;
   cout<<"pobl. en cross "<<rate_cross_rec<<endl;
   inI>>rate_mut_rec;
   cout<<"pobl. en mut "<<rate_mut_rec<<endl;
   inI>>rate_rand_rec;
   cout<<"poblacion aleatoria "<<rate_rand_rec<<endl;

	string tipo_M;
	tipo_M.clear();
	inI >> tipo_M;
	cout<<tipo_M;
	if(tipo_M.empty())
		tipo_matriz = 0;
	else
	{
		if(tipo_M == "FT")
			tipo_matriz = 0;
		else
		{
			tipo_matriz = 1;
			inI >> traslacion;
		}
	}

   critical_size = pop_size/4;
   inI.close();
   return 1;

}

void opencv_abejas::select_memory()
{
  unsigned nbytes;

  if(tourneysize > pop_size)
    {
      printf("FATAL: Tournament size (%d) > pop_size (%d)\n",tourneysize,pop_size);
      exit(-1);
    } ;
}

/*====================================================================
Initialses zero'th generation and global parameters
====================================================================*/

void opencv_abejas::initialize()
{
	oldpop = newpop = tempop = NULL;
	oldpop = new INDIVIDUAL[pop_size];
	newpop = new INDIVIDUAL[pop_size_off];
	tempop = new INDIVIDUAL[pop_size + pop_size_off];

	if (oldpop == NULL) nomemory("oldpop in initialize()");
	if (newpop == NULL) nomemory("newpop in initialize()");
	if (tempop == NULL) nomemory("tempop in initialize()");

        if (!return_to_inactive) {
		int sf = m.countBeesInState(5);
		int rec = m.countBeesInState(3);
		for (int k = 0; k < sf; k++) {
			m.giveFeedback(5, 0);
		}
		for (int k = 0; k < rec; k++) {
			m.giveFeedback(3, 0);
		}
        }

	if(!inicio)
	{
		for(int k=0; k < pop_size; k++)
		{
			oldpop[k].obj = 0.0;
			oldpop[k].type = 0;

			for (int j=0; j < nvar_real; j++)
			{
				oldpop[k].xreal[j] = (double)(cvRandReal(&rng)*(xreal_upper[j]-xreal_lower[j]) + xreal_lower[j]);
			}

			if (stage == 0) // modificacion fase exploracion
			{
				oldpop[k].state = m.giveFeedback(0, 0); // inactive -> exploration
			} else if (stage == 2) // fase recoleccion
                        {
				oldpop[k].state = m.giveFeedback(0, 1); // inactive -> foraging
                        }
		}
	}
	else
		for(int k=0; k < pop_size; k++)
			copy_individual(&nextpop[k], &oldpop[k]);
				no_xover = no_mutation = binmut = 0;
				proyecta(pop_size);

}
/**********************************************************
 * Proyecta las coordenadas 3D de las abejas en cada imagen 2D
 * **********************************************************/
void opencv_abejas::proyecta(int popu, int op)
{
   CvMat* X = cvCreateMat(popu, 3, CV_32FC1);
	CvMat* x1 = cvCreateMat(popu, 2, CV_32FC1);
	CvMat* x2 = cvCreateMat(popu, 2, CV_32FC1);

	for(int s = 0; s < popu; s++)
   {
		if(op == 1)
		{
			CV_MAT_ELEM(*X, float, s, 0) = oldpop[s].xreal[0];
			CV_MAT_ELEM(*X, float, s, 1) = oldpop[s].xreal[1];
			CV_MAT_ELEM(*X, float, s, 2) = oldpop[s].xreal[2];
		}
		else if(op == 2)
		{
			CV_MAT_ELEM(*X, float, s, 0) = newpop[s].xreal[0];
			CV_MAT_ELEM(*X, float, s, 1) = newpop[s].xreal[1];
			CV_MAT_ELEM(*X, float, s, 2) = newpop[s].xreal[2];
		}
	}

	cvProjectPoints2(X, Vrotacion_izq, Vtraslacion_izq, MIntrinsecos_izq, NULL, x1);
	if(tipo_matriz == 0)
	{
		cvProjectPoints2(X, Vrotacion_der, Vtraslacion_der, MIntrinsecos_der, NULL, x2);
	}
	else
	{
		CvMat* X1  = cvCreateMat(4, 1, CV_32FC1);
		CvMat* X2  = cvCreateMat(4, 1, CV_32FC1);
		CvMat* As  = cvCreateMat(4, 4, CV_32FC1);
		cvSetZero(As);
		cvmSet(As, 0, 0, 1); cvmSet(As, 1, 1, 1); cvmSet(As, 2, 2, 1); cvmSet(As, 3, 3, 1);
		cvmSet(As, 1, 3, traslacion);
		for(int s = 0; s < popu; s++)
   	{
			CvMat* sub = NULL;
			CvRect rect;
			rect. x = 0; rect.y = 0; rect.width = 1; rect.height = 3;
			if(op == 1)
			{
				CV_MAT_ELEM(*X1, float, s, 0) = oldpop[s].xreal[0];
				CV_MAT_ELEM(*X1, float, s, 1) = oldpop[s].xreal[1];
				CV_MAT_ELEM(*X1, float, s, 2) = oldpop[s].xreal[2];
			}
			else if(op == 2)
			{
				CV_MAT_ELEM(*X1, float, s, 0) = newpop[s].xreal[0];
				CV_MAT_ELEM(*X1, float, s, 1) = newpop[s].xreal[1];
				CV_MAT_ELEM(*X1, float, s, 2) = newpop[s].xreal[2];
			}
			CV_MAT_ELEM(*X1, float, s, 3) = (float) 1;

			cvMatMulAdd(&As, &X1, 0, &X2);
			cvGetSubRect(X2, sub, rect);
			CV_MAT_ELEM(*X, float, s, 0) = CV_MAT_ELEM(*sub, float, 0, 0);
			CV_MAT_ELEM(*X, float, s, 1) = CV_MAT_ELEM(*sub, float, 1, 0);
			CV_MAT_ELEM(*X, float, s, 2) = CV_MAT_ELEM(*sub, float, 2, 0);
		}
		cvProjectPoints2(X, Vrotacion_der, Vtraslacion_der, MIntrinsecos_der, NULL, x2);
		cvReleaseMat(&X1);
		cvReleaseMat(&X2);
		cvReleaseMat(&As);
	}
	for(int s = 0; s < popu; s++)
   {
		if(op == 1)
		{
			oldpop[s].u1 = CV_MAT_ELEM(*x1, float, s, 0);
			oldpop[s].v1 = CV_MAT_ELEM(*x1, float, s, 1);

			oldpop[s].u2 = CV_MAT_ELEM(*x2, float, s, 0);
			oldpop[s].v2 = CV_MAT_ELEM(*x2, float, s, 1);
		}
		else if(op == 2)
		{
			newpop[s].u1 = CV_MAT_ELEM(*x1, float, s, 0);
			newpop[s].v1 = CV_MAT_ELEM(*x1, float, s, 1);

			newpop[s].u2 = CV_MAT_ELEM(*x2, float, s, 0);
			newpop[s].v2 = CV_MAT_ELEM(*x2, float, s, 1);
		}
	}
	cvReleaseMat(&X);
	cvReleaseMat(&x1);
	cvReleaseMat(&x2);
}

/********************************************************************
 * Función de aptitud para las abejas exploradoras
 * *****************************************************************/
void opencv_abejas::fitness_function_exp(INDIVIDUAL *indv, int maxX, int maxY, IplImage *buf_F1, IplImage *buf_F2)
{
   int u1, v1, u2, v2;
   double g_12, f_12, F;

	u1 = (int)indv->u1;
   if(modf(indv->u1, &F) >= 0.5)
      u1++;
   v1 = (int)indv->v1;
  if(modf(indv->v1, &F) >= 0.5)
      v1++;
   u2 = (int)indv->u2;
   if(modf(indv->u2, &F) >= 0.5)
      u2++;
   v2 = (int)indv->v2;
  if(modf(indv->v2, &F) >= 0.5)
      v2++;

  if((u1 > maxY - 10 || u1 < 10 || v1 < 10  || v1 > maxX - 10) ||
      (u2 > maxY - 10 || u2 < 10 || v2 < 10 || v2 > maxX - 10))
  {
     // cout<<"se paso:"<<u1<<","<<v1<<"<-->"<<u2<<","<<v2<<endl;
      indv->obj = -10000.0;
      return;
  }

///////////// Compute Texture (statistics) ///////////////////

	double f1 = opencv_algos::rough(buf_F1, u1, v1);
	double f2 = opencv_algos::rough(buf_F2, u2, v2);

	if(fabs(f2 - f1)  >  0.05)
  {
     indv->obj = 0.0;
     return;
  }
	if((f2 + f1)/2 < 0.70)
  {
     indv->obj = 0.0;
     return;
  }

///////////// Compute sobel ///////////////////
	g_12 = (opencv_algos::sobel(buf_F1, u1, v1)/255.0) * (opencv_algos::sobel(buf_F2, u2, v2)/255.0);

///////////// Compute ZNCC correlation /////////////////
	f_12 = correlacionZNCC(u1, v1, u2, v2, buf_F1, buf_F2, maxX, maxY);

   if(f_12 >= 0.95)
   {
		F = g_12 * f_12;
   }
   else
   {
     F = -10000.0;
   }

  // Write output:
     indv->obj = F;

}

/********************************************************************
 * Función de aptitud para las abejas recolectoras
 * *****************************************************************/
void opencv_abejas::fitness_function_rec(INDIVIDUAL *indv, int maxX, int maxY, IplImage * buf_F1, IplImage * buf_F2)
{
   int u1, v1, u2, v2;
   double g_12, f_12, F;

	u1 = (int)indv->u1;
   if(modf(indv->u1, &F) >= 0.5)
      u1++;
   v1 = (int)indv->v1;
  if(modf(indv->v1, &F) >= 0.5)
      v1++;
   u2 = (int)indv->u2;
   if(modf(indv->u2, &F) >= 0.5)
      u2++;
   v2 = (int)indv->v2;
  if(modf(indv->v2, &F) >= 0.5)
      v2++;

  if((u1 > maxY - 10 || u1 < 10 || v1 < 10  || v1 > maxX - 10) ||
      (u2 > maxY - 10 || u2 < 10 || v2 < 10 || v2 > maxX - 10))
  {
      indv->obj = -10000.0;
      return;
  }

///////////// Compute Texture (statistics) ///////////////////
	double h_izq = opencv_algos::rough(buf_F1, u1, v1);
	double h_der = opencv_algos::rough(buf_F2, u2, v2);

	g_12 = h_der * h_izq;

///////////// Compute ZNCC correlation ///////////////////
   f_12 = correlacionZNCC(u1, v1, u2, v2, buf_F1, buf_F2, maxX, maxY);

	if((f_12 >= 0.95) && (fabs(h_der - h_izq) < 0.05))
   {
	   F = g_12 * f_12;
   }
   else
   {
   	F = -10000.0;
   }
   // Write output:
   indv->obj = F;
}

//////////////////////////////////////////////////////////////////////////////////
//
//  Función que se asegura de que el numero de individuos generados por cada operador
//  esté dentro del tamano de la poblacion
//
//////////////////////////////////////////////////////////////////////////////////

void opencv_abejas::verifica_pop()
{
	int num_cross, num_mut, num_rand, band;
   double temp1, temp2;

//modificacion
if (rate_dif > 0) {
   rate_rand += rate_dif * 2;
   rate_mut -= rate_dif;
   rate_cross -= rate_dif;
}
//cout << "-- rate_rand: " << rate_rand << endl;
//cout << "-- rate_mut: " << rate_mut << endl;
//cout << "-- rate_cross: " << rate_cross << endl;
//cout << "-- pop_size_off begin: " << pop_size_off << endl;

   temp1 = pop_size_off * rate_cross;
   num_cross = (int)temp1;
   if(modf(temp1, &temp2) >= 0.5)
      num_cross++;
   temp1 = pop_size_off * rate_mut;
   num_mut = (int)temp1;
   if(modf(temp1, &temp2) >= 0.5)
      num_mut++;
   temp1 = pop_size_off * rate_rand;
   num_rand = (int)temp1;
   if(modf(temp1, &temp2) >= 0.5)
      num_rand++;
   if((num_mut + num_cross + num_rand) != pop_size_off)
   {
      //cout<<"Advertencia: La Poblacion cambiara para cumplir proporcion de operadores evolutivos, de: "<<pop_size_off<<" a: "<<num_mut + num_cross + num_rand<<endl;
      
      //modificacion lo anterior podria causar error de segmentacion porque podria aumentar el tamaño de pop_size_offset
      if (pop_size_off > num_mut + num_cross + num_rand) {
         //pop_size_off = num_mut + num_cross + num_rand;
         int dif = pop_size_off - (num_mut + num_cross + num_rand);
         num_mut += dif;
      } else {
         int dif = (num_mut + num_cross + num_rand) - pop_size_off;
         num_mut -= dif; // se usa num_mut porque suele ser la mayor subpoblacion
      }
      band = 0;
   }
 	if(num_cross % 2 != 0)
   {
		if(num_rand > 0)
      {
      	num_cross++;
         num_rand --;
      }
      else
      {
      	num_cross --;
         num_rand ++;
       }
   }
   rate_cross = num_cross;
   rate_mut = num_mut;
   rate_rand = num_rand;

//cout << "-- num_rand: " << rate_rand << endl;
//cout << "-- num_mut: " << rate_mut << endl;
//cout << "-- num_cross: " << rate_cross << endl;
//cout << "-- pop_size_off end: " << pop_size_off << endl;
}

/*====================================================================
Copys contents of one individual into another.
====================================================================*/
void opencv_abejas::copy_individual(INDIVIDUAL *indiv1, INDIVIDUAL *indiv2)
{
   int k;

   if (indiv1==NULL) error_ptr_null("indiv1 in copy_individual");
   if (indiv2==NULL) error_ptr_null("indiv2 in copy_individual");

   for (k=0; k < nvar_real; k++)
   {
      indiv2->xreal[k] = indiv1->xreal[k];
   }
      indiv2->obj = indiv1->obj;
      indiv2->u1 = indiv1->u1;
      indiv2->v1 = indiv1->v1;
      indiv2->u2 = indiv1->u2;
      indiv2->v2 = indiv1->v2;
      indiv2->type = indiv1->type;
}


/*====================================================================
Prints an error message and terminates the program
====================================================================*/
void opencv_abejas::nomemory(const char *string)
{
   cout<< "\nmalloc: out of memory making " << string <<"!!"<<endl;
   cout<< "\n Program is halting ....."<<endl;
   exit(-1);
}

/*==============================================================
Gives error message of null pointer  and terminates the program.
==============================================================*/
void opencv_abejas::error_ptr_null(const char *string)
{
   cout << "\n Error !! Pointer "<< string << "found Null !"<<endl;
   cout << "\n Program is halting ....."<< endl;
   exit(-1);
}

/**************************************************************
 * Funcion de correlacion cruzada normalizada en cero
 * ***********************************************************/
double opencv_abejas::correlacionZNCC(int u1, int v1, int u2, int v2, IplImage * buf_F1, IplImage * buf_F2, int maxX, int maxY)
{
   int k, l, int_F1, int_F2, ventana = 11, paso_vent;
   double sumnum, sumden1, sumden2, my_int_F1, my_int_F2;

    paso_vent = (ventana -1)/2;


      if ( (u2 >= paso_vent) && (u2 < maxY-paso_vent) && (v2 >= paso_vent) && (v2 < maxX-paso_vent) &&
           (u1 >= paso_vent) && (u1 < maxY-paso_vent) && (v1 >= paso_vent) && (v1 < maxX-paso_vent) )
      {
         my_int_F1 = 0.0;
         my_int_F2 = 0.0;

         for (k=-paso_vent;k<=paso_vent;k++)
         {
            for (l=-paso_vent;l<=paso_vent;l++)
            {
               int_F1=(int)cvGetReal2D(buf_F1, u1+l, v1+k);
               int_F2=(int)cvGetReal2D(buf_F2, u2+l, v2+k);
               my_int_F1 += int_F1;
               my_int_F2 += int_F2;
            }
         }
         my_int_F1 = my_int_F1/(ventana*ventana);
         my_int_F2 = my_int_F2/(ventana*ventana);
         sumnum=0.0;
         sumden1=0.0;
         sumden2=0.0;

         for (k=-paso_vent; k<= paso_vent;k++)
         {
            for (l=-paso_vent;l<=paso_vent;l++)
            {
               int_F1=(int)cvGetReal2D(buf_F1, u1+l, v1+k);
               int_F2=(int)cvGetReal2D(buf_F2, u2+l, v2+k);
               sumnum += (int_F1 - my_int_F1) * (int_F2 - my_int_F2);
               sumden1 += (int_F1 - my_int_F1) * (int_F1 - my_int_F1);
               sumden2 += (int_F2 - my_int_F2) * (int_F2 - my_int_F2);
            }
         }
         return (sumnum/(sqrt(sumden1*sumden2)));
      }
}

/*====================================================================
GENERATION OF NEW POPULATION through SELECTION, XOVER, RANDOM & MUTATION :
//modified by Cesar Puente.

Note. This function works only with real variables. Not binary
====================================================================*/
void opencv_abejas::generate_new_pop()
{
   int k = 0,mate1,mate2, num_cross, num_mut, num_rand;
   double temp1, temp2;
   preselect_tour();  //reinicia el orden de los individuos para escoger en el torneo

//"Mutacion .";
   for(int i = 0; i < rate_mut; i++)
   {
	    mate1 = tour_select();    //realiza torneo

       copy_individual(&oldpop[mate1], &newpop[k]);
       // mutation polinomial
       if(mut_type)
          mutation(&newpop[k]);
       else
       {
	cerr<< "Error en archivo de parametros: Operador de mutacion no valido.  Abort";
	exit(-1);
       }
       newpop[k].type = 1;
       k++;
   }

//" Crossover: .";
   for(int i = 0; i < rate_cross; i += 2)
   {
		mate1 = tour_select();    //realiza torneo
      mate2 = tour_select();

       if(cross_type)  // crossover SBX
          cross_over(mate1,mate2,k,k+1);
       else
		 {
			cerr<< "Error en archivo de parametros: Operador de cruzamiento no valido.  Abort";
			exit(-1);
		 }
       newpop[k].type = -1;
       newpop[k+1].type = -1;
      k += 2;
   }
//" Random: .";
   for(int i = 0; i < rate_rand; i++)
   {
     newpop[k].obj = 0.0;
     newpop[k].type = 0;
     for (int j=0; j < nvar_real; j++)
     {
         newpop[k].xreal[j] = (double)(cvRandReal(&rng)*(xreal_upper[j]-xreal_lower[j]) + xreal_lower[j]);
     }
      k++;
   }
  proyecta(pop_size_off, 2);
}

void opencv_abejas::preselect_tour()     //inicializa orden de individuos
{
    reset1();      //rompe el orden en que estan los individuos para que el torneo escoga de manera aleatoria
    tourneypos = 0;
}

/******************************************************
 * Función de selección por torneo
 * *************************************************/
int opencv_abejas::tour_select()
{
    int pick, winner, i;

    /* If remaining members not enough for a tournament, then reset list */
start_select :
    if((pop_size - tourneypos) < tourneysize)
    {
        reset1();
        tourneypos = 0;
    }

    /* Select tourneysize structures at random and conduct a tournament */
    winner=tourneylist[tourneypos];
    if(pop_size < 2)
      return (winner);
/* Added by RBA */
    if( winner < 0 || winner > pop_size-1)
	 {
   	cout<< "\n Warning !! ERROR1 tourpos = "<<tourneypos<<" winner = "<<winner<<endl;
      preselect_tour();
      goto start_select;
	 }
    for(i=1; i<tourneysize; i++)
    {
        pick=tourneylist[i+tourneypos];
/* Added by RBA */
        if (pick < 0 || pick > pop_size-1)
        {
            preselect_tour();
            printf("\n Warning !! ERROR in pick gen");
            cout<<pick<<","<<i<<","<<tourneypos<<","<<pop_size;
            getchar();
            goto start_select;
        }
      // case 1:
      if(MINM * oldpop[pick].obj < MINM * oldpop[winner].obj)
      	winner=pick;
    }

    /* Update tourneypos */
    tourneypos += tourneysize;
    return(winner);
}

void opencv_abejas::reset1()
/* Name changed from reset because of clash with lib. function - RBA */
/* Shuffles the tourneylist at random */
{
    int i, rand1, rand2, temp_site;

    for(i=0; i<pop_size; i++) tourneylist[i] = i;

    for(i=0; i < pop_size; i++)
    {

        rand1= cvRandInt(&rng) % pop_size;
        rand2= cvRandInt(&rng) % pop_size;
        temp_site = tourneylist[rand1];
        tourneylist[rand1]=tourneylist[rand2];
        tourneylist[rand2]=temp_site;
    }
}

/*====================================================================
CROSS - OVER  USING strategy of uniform 50% variables
  For one variable problem, it is crossed over as usual.
  For multivariables, each variable is crossed over with a probability
  of 50 % , each time generating a new random beta.
====================================================================*/
void opencv_abejas::cross_over(int first, int second, int childno1, int childno2)
{
    double difference,x_mean,beta;
    double u = 0.0;
    int site,k,x_s;

    x_s = 0;
    if (1)   /* Cross over has to be done */
    {
      no_xover++;
      if (nvar_real > 0)
      {
         for (site = 0; site < nvar_real; site++)
         {
             create_children(oldpop[first].xreal[site],oldpop[second].xreal[site],
									  &(newpop[childno1].xreal[site]),&(newpop[childno2].xreal[site]),
                             xreal_lower[site],xreal_upper[site],&u);
         }
       }
    }
    else /* Passing x-values straight */
    {
      for (site=0; site < nvar_real; site++)
      {
         newpop[childno1].xreal[site] = oldpop[first].xreal[site];
         newpop[childno2].xreal[site] = oldpop[second].xreal[site];
      }
    }
}

/*====================================================================
Creates two children from parents p1 and p2, stores them in addresses
pointed by c1 and c2.  low and high are the limits for x values and
rand_var is the random variable used to create children points.
====================================================================*/
void opencv_abejas::create_children(double p1, double p2, double *c1, double *c2, double low, double high, double *rand_var)
{
    double difference,x_mean,beta,v2,v1;
    double u, distance, umax, temp, alpha;
    int flag;

    if (c1 == NULL)  error_ptr_null("c1 in create_children");
    if (c2 == NULL)  error_ptr_null("c2 in create_children");
    if (rand_var == NULL)  error_ptr_null("rand_var in create_children");
    flag = 0;
    if ( p1 > p2)
    {
      temp = p1;
      p1 = p2;
      p2 = temp;
      flag = 1;
    }
    x_mean = (p1 + p2) * 0.5;
    difference = p2 - p1;
    if ( (p1-low) < (high-p2) )
      distance = p1- low;
    else
      distance = high - p2;
    if (distance < 0.0)
      distance = 0.0;
    if (RIGID && (difference > EPSILON))
    {
      alpha = 1.0 + (2.0 * distance / difference);
      umax = 1.0 - (0.5 / pow((double)alpha,(double)(n_distribution_c+1.0)));
      (*rand_var) = umax * cvRandReal(&rng);
    }
    else
      (*rand_var) = cvRandReal(&rng);
    beta = get_beta(*rand_var);
    if (fabs(difference * beta) > INFINITY)
      beta = INFINITY/difference;
    v2 = x_mean + (beta * 0.5 * difference);
    v1 = x_mean - (beta * 0.5 * difference);

    if (v2 < low) v2 = low;
    if (v2 > high) v2 = high;
    if (v1 < low) v2 = low;
    if (v1 > high) v2 = high;
    *c2 = v2; *c1 = v1;
    if (flag == 1) { temp = *c1; *c1 = *c2; *c2 = temp; }
}


/*===================================================================
Mutation Using polynomial probability distribution. Picks up a random
site and generates a random number u between -1 to 1, ( or between
minu to maxu in case of rigid boudaries) and calls the routine
get_delta() to calculate the actual shift of the value.
====================================================================*/
void opencv_abejas::mutation(INDIVIDUAL  *indiv)
{
   double distance1,x,delta_l,delta_u,delta,u;
   int k, site;

   if (indiv == NULL) error_ptr_null("indiv in mutation");

   if (nvar_real > 0)
     for (site = 0; site < nvar_real; site++)
       {
         no_mutation++;
         if(RIGID)
         {
            x = indiv->xreal[site];
            distance1 = xreal_lower[site] - x;
            delta_l = distance1/(xreal_upper[site] - xreal_lower[site]);
            if (delta_l < -1.0)
               delta_l = -1.0;
            distance1 = xreal_upper[site] - x;
            delta_u = distance1/(xreal_upper[site] - xreal_lower[site]);
            if (delta_u > 1.0)
               delta_u = 1.0;
            if (-1.0*delta_l < delta_u)
               delta_u = -1.0 * delta_l;
            else
               delta_l = -1.0 * delta_u;
         }
         else
         {
            delta_l = -1.0;
            delta_u = 1.0;
         }
			u = (double)cvRandReal(&rng);
         /* calculation of actual delta value */
         delta = get_delta(u, delta_l, delta_u) * (xreal_upper[site] - xreal_lower[site]);
         indiv->xreal[site] += delta;

         if (indiv->xreal[site] < xreal_lower[site])
            indiv->xreal[site] = xreal_lower[site];
         if (indiv->xreal[site] > xreal_upper[site])
            indiv->xreal[site] = xreal_upper[site];
       }
}

/*===================================================================
Calculates beta value for given random number u (from 0 to 1)
If input random numbers (u) are uniformly distributed for a set of
inputs, this results in uniform distribution of beta values in case
of BLX , and Binary Probability distribution simulation in case of
SBX.

====================================================================*/
double opencv_abejas::get_beta(double u)
{
   double beta;

   if (1.0-u < EPSILON ) u = 1.0 - EPSILON;
   if ( u < 0.0) u = 0.0;
   if (u < 0.5)
      beta = pow(2.0*u,(1.0/(n_distribution_c+1.0)));
   else
      beta = pow( (0.5/(1.0-u)),(1.0/(n_distribution_c+1.0)));
   return beta;
}

/*==================================================================
For given u value such that   -1 <= u <= 1, this routine returns a
value of delta from -1 to 1. Exact value of delta depends on specified
n_distribution. This is called by mutation().
====================================================================*/
double opencv_abejas::get_delta(double u, double delta_l, double delta_u)
{
  double delta, aa;

   if (u >= 1.0-1.0e-9)
      delta = delta_u;
   else if (u <= 0.0+1.0e-9)
      delta = delta_l;
   else
   {
      if (u < 0.5)
      {
         aa = 2.0*u + (1.0-2.0*u)*pow((1+delta_l),(n_distribution_m + 1.0));
         delta = pow(aa, (1.0 / (n_distribution_m + 1.0))) - 1.0;
      }
      else
      {
         aa = 2.0*(1-u) + 2.0*(u-0.5)*pow((1-delta_u),(n_distribution_m + 1.0));
         delta = 1.0 - pow(aa, (1.0 / (n_distribution_m + 1.0)));
      }
    }
  if (delta < -1.0 || delta > 1.0)
    {
      cout<<"Error in mutation!! delta = "<<delta<<endl;
      exit(-1);
    }
  return (delta);
}


//=======================================================================
//"Junta dos poblaciones: .";
//=======================================================================
void opencv_abejas::merge_pop()
{
   int k = 0;
   for(int i = 0; i < pop_size; i++, k++)
   {
      copy_individual(&oldpop[i],&tempop[k]);
   }

   for(int j = 0; j < pop_size_off; j++, k++)
   {
      copy_individual(&newpop[j],&tempop[k]);
   }
/*	for(int h = 0; h < pop_size+pop_size_off;h++)
		cout<<tempop[h].obj<<endl;
	getchar();*/
}


//=======================================================================
//Funciones para calcular el Sharing
//=======================================================================

void opencv_abejas::sharing2D()
{
	int Npix = buf_F1->width * buf_F1->height;
	double R = sqrt(Npix / pop_size_off) - 1;
	for(int v = 0; v < buf_F1->width; v+= (int)R)
	{
		for(int u = 0; u < buf_F1->height; u+= (int)R)
		{
			CvRect rect;
			rect. x = v; rect.y = u; rect.width = (int)R; rect.height = (int)R;
			sh(rect);
		}

	}
}

void opencv_abejas::sh(CvRect r)
{
	int sum_izq = 0, sum_der = 0;
	double Sum;
	//numero de proyecciones en la img izq - der

	for(int i = 0; i < pop_size+pop_size_off; i++)
	{
		if((tempop[i].v1 >= r.x && tempop[i].v1 < r.x + r.width) &&
			(tempop[i].u1 >= r.y && tempop[i].u1 < r.y + r.height))
			sum_izq ++;

		if((tempop[i].v2 >= r.x && tempop[i].v2 < r.x + r.width) &&
			(tempop[i].u2 >= r.y && tempop[i].u2 < r.y + r.height))
			sum_der ++;
	}

//penalizacion de cada individuo segun la celda en la que cae
	for(int i = 0; i < pop_size+pop_size_off; i++)
	{
		//if(/*tempop[i].obj - alpha_share*(sum_izq + abs(sum_izq - sum_der)) <= 0) // && */tempop[i].obj >= 0)
		{
			//cout<<"mayor"<<tempop[i].obj<<"= "<<sum_izq<<" , "<<sum_der<<" ALPHA"<<alpha_share;getchar();//tempop[i].obj = 0.0;
			tempop[i].obj -= alpha_share * (sum_izq + fabs(sum_izq - sum_der));
			//tempop[i].obj /= alpha_share * ((sum_izq + sum_der)/2.0);
			//cout<<"mayor"<<tempop[i].obj;getchar();
		}
	}
}


void opencv_abejas::sharing3D()
{
   bad_bees = 0; // modificacion
   for(int i = 0; i < pop_size + pop_size_off; i++)
   {
      if(tempop[i].obj != -100)
      {

      double sum = 0.0;
      for(int j = 0; j < pop_size + pop_size_off; j++)
      {
         sum += sh(i, j);
      }
      if(sum > 1) {
         bad_bees++; // modificacion
         tempop[i].obj  /= (sum);
      }
      else {
         tempop[i].obj  /= (sum);
      }
      }
   }
}

double opencv_abejas::sh(int uno, int dos)
{
   double dist = distanc(uno, dos, 3);
//   vcl_cout<<dist<<vcl_endl;
   if(dist <= sigma_share)
   {
      return (1 - dist/sigma_share);
   }
   else
      return 0.0;
}

double opencv_abejas::distanc(int one, int two, int op)
{
  int k;
  double sum, result;

  sum = 0.0;
  switch(op)
  {
     case 1:
      for (k=0; k<nvar_real; k++)
       sum += square((oldpop[one].xreal[k] - oldpop[two].xreal[k])/(xreal_upper[k]-xreal_lower[k]));
       result = sqrt(sum/nvar_real);
     break;
     case 2:
      for (k=0; k<nvar_real; k++)
       sum += square((newpop[one].xreal[k] - newpop[two].xreal[k])/(xreal_upper[k]-xreal_lower[k]));
       result = sqrt(sum/nvar_real);
     break;
     case 3:
      for (k=0; k<nvar_real; k++)
          sum += square(tempop[one].xreal[k] - tempop[two].xreal[k]);
      result = sqrt(sum);
     break;
     default:
      result = -1;
  }
   return result;
}

/***********************************************************************
 * Función que ordena la población de acuerdo al resultado en la función de aptitud
 * y regresa las "mu" mejores
 * ***********************************************************************/
void opencv_abejas::best_mu()
{

double increment = 0.01;
int bad_bee_percent = (bad_bees * 100) / (pop_size + pop_size_off);
int increment_pop = ((real_pop_size - pop_dif) * (bad_bee_percent / 2)) / 100;
// aumentar hijos aleatorios
if (bad_bee_percent > 50 && current_rate_cross > rate_dif + increment && current_rate_mut > rate_dif + increment)
   rate_dif += increment;
//disminuir hijos aleatorios
if (bad_bee_percent < 20 && rate_dif >= increment)
   rate_dif -= increment;
if (bad_bee_percent > 60 && real_pop_size - (pop_dif + increment_pop) > 0)
   pop_dif += increment_pop;
if (bad_bee_percent < 20 && pop_dif >= increment_pop)
   pop_dif -= increment_pop;
//cout << "-- rate_dif: " << rate_dif << endl;
//cout << "-- pop_dif: " << pop_dif << " real_pop_size: " << real_pop_size << endl;

	CvMat *M = cvCreateMat(pop_size+pop_size_off, 2, CV_32FC1);

   for(int i = 0; i < pop_size+pop_size_off; i ++)
   {
		CV_MAT_ELEM( *M, float, i, 0) = tempop[i].obj * MINM;
		CV_MAT_ELEM( *M, float, i, 1) = (float)i;
   }
	opencv_algos::quickSort(M, pop_size+pop_size_off);
   //for(int i = 0; i < pop_size; i++)
   for (int i = 0; i < real_pop_size - pop_dif; i++) // modificacion
   {
      copy_individual(&tempop[(int)CV_MAT_ELEM( *M, float, i, 1)], &oldpop[i]);
      if (stage == 0)
      {
         oldpop[i].state = m.giveFeedback(1, 1); // modificacion explora -> exp exitosa
         // hacer cambios al tamano de la poblacion dependiendo de los resultados de sharing3D
         // pasar algunas abejas (las peores) a inactividad: m.giveFeedback(4, 0);
         // cambiar el valor de pop_size segun corresponda
      } else if (stage == 2)
      {
         oldpop[i].state = m.giveFeedback(2, 1); // modificacion recolectora -> rec exitosa
         // hacer cambios al tamano de la poblacion dependiendo de los resultados de sharing3D
         // pasar algunas abejas (las peores) a inactividad: m.giveFeedback(5, 0);
         // cambiar el valor de pop_size segun corresponda
      }
   }
   for (int k = 0; k < pop_dif; k++) {
	if (stage == 0)
         {
            m.giveFeedback(1, 0); // modificacion explora -> inactividad
         } else if (stage == 2)
         {
            m.giveFeedback(2, 0); // modificacion recolectora -> inactividad
         }
   }
   pop_size = real_pop_size - pop_dif;
}




/*******************************************************************
 * función de la etapa de reclutamiento que define cuantas recolectoras
 * se le asignarán a cada sitio señalado por las exploradoras
 * ****************************************************************/
void opencv_abejas::asigna_recolectoras()
{
//asigna numero de recolectoras
   double val = 0.0, max, min, dimensiones, volumen = 1.0;
	CvMat *disp = cvCreateMat(pop_size, 1, CV_32FC1);

	for(int i = 0; i < pop_size; i++)
	{
		CV_MAT_ELEM(*disp, float, i, 0) = sqrt(square(oldpop[i].u1 - oldpop[i].u2) + square(oldpop[i].v1 - oldpop[i].v2));
		val += fabs(oldpop[i].obj);
	}
	cvMinMaxLoc(disp, &min, &max);
   for(int i = 0; i < pop_size; i++)
   {
		double num, res;
		num = fabs(oldpop[i].obj)/val * (double)pop_size_for;
		if(num <= 3.0)
		{
			num = 0;
		}
		if(modf(num, &res) >= 0.5)
			CV_MAT_ELEM(*sitios, float, i, 0) = (int)(num) + 1;
      else
         CV_MAT_ELEM(*sitios, float, i, 0) = (int)(num);

   }

//asigna espacio

   for(int i = 0; i < nvar_real; i++)
   {
      volumen *= fabs(xreal_lower[i]) + fabs(xreal_upper[i]);
   }
   dimensiones = pow(volumen/pop_size, (1.0/3.0));
   for(int i = 0; i < pop_size; i++)
   {
     double u = CV_MAT_ELEM(*disp, float, i, 0)/max, f, aux, res;
     int dim_penal = (int)dimensiones;

      f = 0.5 * (1-u) + 1 * u;
      aux = dimensiones * f;
      dim_penal = (int)aux;
      if(modf(aux, &res) >= 0.5)
         dim_penal ++;
      for(int j = 0; j < nvar_real; j++)
      {
			CV_MAT_ELEM(*espacios_lower, float, i, j) = oldpop[i].xreal[j] - dim_penal/2;
         if(CV_MAT_ELEM(*espacios_lower, float, i, j) < xreal_lower[j])
            CV_MAT_ELEM(*espacios_lower, float, i, j) = xreal_lower[j];

         CV_MAT_ELEM(*espacios_upper, float, i, j) = oldpop[i].xreal[j] + dim_penal/2;
         if(CV_MAT_ELEM(*espacios_upper, float, i, j) > xreal_upper[j])
            CV_MAT_ELEM(*espacios_upper, float, i, j) = xreal_upper[j];
      }
   }

	cvReleaseMat(&x1);
	cvReleaseMat(&x2);
}


/*====================================================================
Releases the memory for all news
====================================================================*/
void opencv_abejas::free_all()
{
   delete[] oldpop;
   delete[] newpop;
   delete[] tempop;

   oldpop = NULL;
   newpop = NULL;
   tempop = NULL;
}



/**************************************************************************
 * Función que dibuja cada abeja en la imagen izquierda y derecha
 * *********************************************************************/
void opencv_abejas::imprime_abejas()
{
	double Fe;
	CvPoint pt1, pt2;

        //result_file.open("resultados.txt", fstream::app);

	for(int k = 0; k < pop_size; k++)
	{
		pt1.y = (int)oldpop[k].u1;
  		if(modf(oldpop[k].u1, &Fe) >= 0.5)
     		pt1.y++;
	   pt1.x = (int)oldpop[k].v1;
  		if(modf(oldpop[k].v1, &Fe) >= 0.5)
     		pt1.x++;
		pt2.y = (int)oldpop[k].u2;
  		if(modf(oldpop[k].u2, &Fe) >= 0.5)
     		 pt2.y++;
	   pt2.x = (int)oldpop[k].v2;
  		if(modf(oldpop[k].v2, &Fe) >= 0.5)
     		 pt2.x++;

		cvCircle(buf_F1, pt1, 1, CV_RGB(247, 255, 12), -1);
		cvCircle(buf_F2, pt2, 1, CV_RGB(247, 255, 12), -1);

           //result_file << oldpop[k].u1;
           //result_file << ", " << oldpop[k].v1;
           //result_file << ", " << oldpop[k].u2;
           //result_file << ", " << oldpop[k].v2;
           //result_file << ", " << oldpop[k].xreal[0];
           //result_file << ", " << oldpop[k].xreal[1];
           //result_file << ", " << oldpop[k].xreal[2];
           //result_file << ", " << oldpop[k].obj << endl;
	}
        //result_file.close();

	cvShowImage("izquierda",buf_F1);
	cvShowImage("derecha",buf_F2);
	cvWaitKey(2);
}

void opencv_abejas::guarda_exploradoras()
{
  for(int i = 0; i < pop_size; i++)
  {
    copy_individual(&oldpop[i], &nextpop[i]);
  }
}

