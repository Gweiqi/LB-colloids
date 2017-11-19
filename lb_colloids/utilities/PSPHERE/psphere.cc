#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <glib.h>

#include "matrix.h"
#include "png.h"
#include "png_image.h"
#include "png_image_common.h"

#define N 32
#define NY N
#define NX N

#define DEFAULT_RADIUS     25
#define DEFAULT_POROSITY   0.5;
#define DEFAULT_DIMENSION  256


// these are not used yet
#define NORMAL_BOUNDARY 0
#define PERIODIC_BOUNDARY 1


gint iround(double val)
{
  if ((val-floor(val))>=0.5)
    return( (gint) ceil(val));
  else
    return( (gint) floor(val));   
}

int patch(MATRIXM *map, gint x, gint y,double r)
{
  
  gdouble circumference;
  gint nsteps;
  gint i;
  gdouble angle;
  gint ixpos,iypos,ix,newxpos,newypos;


  circumference=2*M_PI*r+1;
  nsteps=(unsigned int) ceil(circumference);
  //printf("%f %ld \n",circumference,nsteps);
  nsteps/=4;
  if (nsteps==1) nsteps=2;
  //if (r==1) nsteps=2;
  for (i=0;i<nsteps;i++)
    {
      angle=0.5*M_PI*((double)i/(double)nsteps);
      ixpos=iround(r*cos(angle));
      iypos=iround(r*sin(angle));
      for (ix=(-ixpos);ix<=ixpos;ix++) 
	{
	  newxpos=x+ix;
	  newypos=y+iypos;

	  // TODO: add NORMAL boundary condition
	  // this is a periodic boundary condition
	  if (newxpos>map->cols) newxpos-=map->cols;
	  if (newxpos<1) newxpos+=map->cols;
	  if (newypos>map->rows) newypos-=map->rows;
	  if (newypos<1) newypos+=map->rows;

	  map->matrix[newypos][newxpos]=PARTICLE_SPACE;

	  newypos=y-iypos;
	  if (newypos>map->cols) newypos-=map->cols;
	  if (newypos<1) newypos+=map->cols;

	  map->matrix[newypos][newxpos]=PARTICLE_SPACE;
	}
      //printf("%d %f %ld %ld \n",i,angle,ixpos,iypos); 
    }

  return(0);
}

int generate_plane(MATRIXM *map, gdouble porosity, gdouble radius, gint xyzdim)
{
  gdouble relative_radius,slice_radius,slice_distance;
  gint i;
  gint x,y,z;
  gint slice_location;
  gint lb,ub;
  gint number_of_spheres;
  gsl_rng *rng=gsl_rng_alloc(gsl_rng_taus);
  //time_t the_time;

  gsl_rng_set(rng,time((time_t)0));  // TODO set this to something random!, use systemtime

  slice_location=xyzdim/2;
  lb=slice_location-(int) radius;
  ub=slice_location+(int) radius;
  if ((lb<=0)||(ub>xyzdim))  // array offset 1
    { 
      printf("Radius to big or slice location wrong\n");
      exit(0);
    }
  relative_radius=radius/(double)xyzdim;
  number_of_spheres=iround(-3.0*log(porosity)/(4*M_PI*gsl_pow_3(relative_radius)));
  printf("Number of spheres %d\n",number_of_spheres);

  for (i=1;i<=number_of_spheres;i++)
    {
      z=1+gsl_rng_uniform_int(rng,xyzdim);
      if ((z>lb) && (z<ub))
	{
	  x=1+gsl_rng_uniform_int(rng,xyzdim);
	  y=1+gsl_rng_uniform_int(rng,xyzdim);
	  slice_distance=fabs((double) z- (double) slice_location);
	  slice_radius=sqrt(gsl_pow_2(radius)-gsl_pow_2(slice_distance))-0.5; // 0.5 because the radius is always 0.5 too large
	  if (slice_radius<0) slice_radius=0;
	  //	  printf("%d %d %d %f\n",x,y,z,slice_radius);
	  patch(map,y,x,slice_radius);
	}
    }

  gsl_rng_free(rng);
  return(0);
}

typedef struct XY_DATA
{
  guint x,y;
} XY_DATA;


void free_xy_data(gpointer listdata, gpointer userdata)
{
  free((XY_DATA *)listdata);
}

GList *burn_step_2d(MATRIXM *map,XY_DATA *xy_temp,GList *list,gint xstep, gint ystep, gint value1,gint value2)
{
  gint newxpos,newypos;
  XY_DATA *xy_data;
  GList *listnew=NULL;

  listnew=list;
  newxpos=xy_temp->x+xstep;
  newypos=xy_temp->y+ystep;
  
  // periodic boundary condition
  // TODO: add NORMAL boundary condition 
  // this is a periodic boundary condition
  if (newxpos>map->cols) newxpos-=map->cols;
  if (newxpos<1) newxpos+=map->cols;
  if (newypos>map->rows) newypos-=map->rows;
  if (newypos<1) newypos+=map->rows;
      
  if (map->matrix[newypos][newxpos]==value1)
    {
      xy_data=new (XY_DATA);
      xy_data->x=newxpos;
      xy_data->y=newypos;
      //      printf("%d %d\n",xy_data->x,xy_data->y);
      listnew=g_list_append(list,xy_data);
      map->matrix[newypos][newxpos]=value2;
    }

  return(listnew);
}


int burn_2d(MATRIXM *map,gboolean be_verbose)
{
  GList *listold=NULL;
  GList *listnew=NULL;
//    GList *listtemp=NULL;
  gint x,y;
  guint list_length;
  XY_DATA *xy_data,*xy_temp;
  guint iteration=0;
//    guint npix_percol=0;


  // initialize list with pixels of x=1,y=*
  x=1;
  for (y=1;y<=map->rows;y++)
    {
      if (map->matrix[y][x]==PORE_SPACE)
	{
	  xy_data=(XY_DATA *) malloc (sizeof(XY_DATA));
	  xy_data->x=x;
	  xy_data->y=y;
	  listnew=g_list_append(listnew,xy_data);
	  map->matrix[y][x]=CONNECTED_PORE_SPACE;
	  
	}
    }
  
  list_length=g_list_length(listnew);
  if (be_verbose) printf("Initialized: %d pixels\n",list_length);

  while (list_length)
    {
      iteration++;
      listold=g_list_copy(listnew);
//        printf("Num of pixels:%d\n",g_list_length(listold));
      listold=g_list_first(listold);
      g_list_free(listnew); // doesn't delete data;
      listnew=NULL;

      while (listold!= NULL)
	{
	  xy_temp=(XY_DATA *) listold->data;

	  listnew=burn_step_2d(map,xy_temp,listnew,1,0,PORE_SPACE,CONNECTED_PORE_SPACE);
	  listnew=burn_step_2d(map,xy_temp,listnew,1,1,PORE_SPACE,CONNECTED_PORE_SPACE);
	  listnew=burn_step_2d(map,xy_temp,listnew,1,-1,PORE_SPACE,CONNECTED_PORE_SPACE);

	  if (iteration>1)  // do not burn these pixels (want to go with positive x direction and check out whether the medium percolates)
	    {
	      listnew=burn_step_2d(map,xy_temp,listnew,0,1,PORE_SPACE,CONNECTED_PORE_SPACE);
	      listnew=burn_step_2d(map,xy_temp,listnew,0,-1,PORE_SPACE,CONNECTED_PORE_SPACE);
	      listnew=burn_step_2d(map,xy_temp,listnew,-1,1,PORE_SPACE,CONNECTED_PORE_SPACE);
	      listnew=burn_step_2d(map,xy_temp,listnew,-1,0,PORE_SPACE,CONNECTED_PORE_SPACE);
	      listnew=burn_step_2d(map,xy_temp,listnew,-1,-1,PORE_SPACE,CONNECTED_PORE_SPACE);
	    }
	  listold=g_list_next(listold);
	}
      list_length=g_list_length(listnew);
      if (be_verbose) printf("Iteration: %d number of pixels:%d\n",iteration,g_list_length(listnew));
      listold=g_list_first(listold);
      g_list_foreach(listold,(GFunc) free_xy_data,NULL); // delete data
      g_list_free(listold); 
    }

  // RIGHT TO LEFT, which pixels do really percolate?
  x=map->cols;
  for (y=1;y<=map->rows;y++)
    {
      if (map->matrix[y][x]==CONNECTED_PORE_SPACE)
	{
	  xy_data=(XY_DATA *) malloc (sizeof(XY_DATA));
	  xy_data->x=x;
	  xy_data->y=y;
	  listnew=g_list_append(listnew,xy_data);
	  map->matrix[y][x]=PERCOLATING_PORE_SPACE;
	  
	}
    }
  
  list_length=g_list_length(listnew);
  if (be_verbose) printf("Initialized: %d pixels\n",list_length);

  while (list_length)
    {
      iteration++;
      listold=g_list_copy(listnew);
//        printf("Num of pixels:%d\n",g_list_length(listold));
      listold=g_list_first(listold);
      g_list_free(listnew); // doesn't delete data;
      listnew=NULL;

      while (listold!= NULL)
	{
	  xy_temp=(XY_DATA *) listold->data;

	  listnew=burn_step_2d(map,xy_temp,listnew,-1,1,CONNECTED_PORE_SPACE,PERCOLATING_PORE_SPACE);
	  listnew=burn_step_2d(map,xy_temp,listnew,-1,0,CONNECTED_PORE_SPACE,PERCOLATING_PORE_SPACE);
	  listnew=burn_step_2d(map,xy_temp,listnew,-1,-1,CONNECTED_PORE_SPACE,PERCOLATING_PORE_SPACE);

	  if (iteration>1)  // do not burn these pixels (want to go with positive x direction and check out whether the medium percolates)
	    {
	      listnew=burn_step_2d(map,xy_temp,listnew,0,1,CONNECTED_PORE_SPACE,PERCOLATING_PORE_SPACE);
	      listnew=burn_step_2d(map,xy_temp,listnew,0,-1,CONNECTED_PORE_SPACE,PERCOLATING_PORE_SPACE);
	      listnew=burn_step_2d(map,xy_temp,listnew,1,0,CONNECTED_PORE_SPACE,PERCOLATING_PORE_SPACE);
	      listnew=burn_step_2d(map,xy_temp,listnew,1,1,CONNECTED_PORE_SPACE,PERCOLATING_PORE_SPACE);
	      listnew=burn_step_2d(map,xy_temp,listnew,1,-1,CONNECTED_PORE_SPACE,PERCOLATING_PORE_SPACE);
	    }
	  listold=g_list_next(listold);
	}
      list_length=g_list_length(listnew);
      if (be_verbose) printf("Iteration: %d number of pixels:%d\n",iteration,g_list_length(listnew));
      listold=g_list_first(listold);
      g_list_foreach(listold,(GFunc) free_xy_data,NULL); // delete data
      g_list_free(listold); 
    }
  


  g_list_foreach(listnew,(GFunc) free_xy_data,NULL); // delete data
  g_list_free(listnew); 

  return(0);
}

int invert_map(MATRIXM *map)
{
  gint ix,iy;

 for (ix=1;ix<=map->cols;ix++)
    { 
      for (iy=1;iy<=map->rows;iy++)
	{
	  if (map->matrix[iy][ix]==1) map->matrix[iy][ix]=0; 
	  if (map->matrix[iy][ix]==0) map->matrix[iy][ix]=1; 
	}
    }

 return(0);
}

int check_percolation_2d(MATRIXM *map)
{
  gint ix,iy;
//    gdouble connected_porosity=0.0;
//    gdouble real_porosity=0.0;
  MATRIXM *connectedness;
  gdouble min_conn,max_conn;

  connectedness=new MATRIXM(map->cols,2);

  min_conn=map->rows;
  max_conn=0.0;
  for (ix=1;ix<=map->cols;ix++)
    { 
      for (iy=1;iy<=map->rows;iy++)
	{
	  if (map->matrix[iy][ix]==PERCOLATING_PORE_SPACE) connectedness->matrix[ix][1]++;
	  if ((map->matrix[iy][ix]==PORE_SPACE) || (map->matrix[iy][ix]==CONNECTED_PORE_SPACE) || (map->matrix[iy][ix]==PERCOLATING_PORE_SPACE)) connectedness->matrix[ix][2]++;
	}
      min_conn=GSL_MIN(min_conn,connectedness->matrix[ix][1]);
      max_conn=GSL_MAX(max_conn,connectedness->matrix[ix][1]);
    }
  printf("Minimum percolation diameter %f \n",min_conn);
  printf("Maximum percolation diameter %f \n",max_conn);

  return(0);
}

int calc_porosity(MATRIXM *map)
{
  gint ix,iy;
  gdouble connected_porosity=0.0;
  gdouble real_porosity=0.0;

  for (iy=1;iy<=map->rows;iy++)
    { 
    for (ix=1;ix<=map->cols;ix++)
      {
	if (map->matrix[iy][ix]==PERCOLATING_PORE_SPACE) connected_porosity++;
	if ((map->matrix[iy][ix]==PORE_SPACE) || (map->matrix[iy][ix]==CONNECTED_PORE_SPACE)|| (map->matrix[iy][ix]==PERCOLATING_PORE_SPACE)) real_porosity++;
      }
    }

  real_porosity/=(double)(map->rows*map->cols);
  connected_porosity/=(double)(map->rows*map->cols);

  printf("Real porosity is %f\n",real_porosity);
  printf("Percolating porosity is %f\n",connected_porosity);

  return(0);
}

int main(int argc, char *argv[])
{
  MATRIXM *map;
  gdouble radius,porosity;
  gint xyzdim;
  char output_file_name[257];
  gint opt;
  gboolean be_verbose=FALSE;
  png_image *image=new png_image;
  png_bytep image_data;
  guint loc;
  png_uint_32 j,k;//, height, width;
  guint32 height,width;


  // default values (use define)
  radius=DEFAULT_RADIUS;
  porosity=DEFAULT_POROSITY;
  xyzdim=DEFAULT_DIMENSION;

  strcpy(output_file_name,"synthetic.png");

  while ((opt=getopt(argc,argv,"d:ho:p:r:v")) != -1)
    {
      switch (opt)
	{
	case 'h':
	  printf("psphere: generate an square of with d pixels and a porosity p of penetrable spheres of radius r\n");
	  printf("Usage:\n\n");
	  printf("psphere [-hv] [-d val] [-o file_name] [-p val] [-r val]\n");
	  printf("\n");
	  printf("-d val   : use size <val> as pixel size in x and y direction (default: %d)\n",xyzdim);
	  printf("-h       : this help\n");
	  printf("-o name  : write results to output file <name>\n");
	  printf("-p val   : use porosity <val>, default %5.2f \n",porosity);
	  printf("-r val   : use radius <val>, default %5.0f \n",radius);
	  printf("-v       : be verbose\n");
	  printf("\n");
	  exit(0);
	  break;

	case 'd':
	  xyzdim=atoi(optarg);
	  break;
	  
	case 'o':
	  strcpy(output_file_name,optarg);
	  break;

	case 'p':
	  porosity=atof(optarg);
	  if ((porosity<=0.0) || (porosity>=1.0))
	    {
	      printf("Porosity should be in the range <0.0 1.0>\n");
	      exit(1);
	    }
	  break;

	case 'r':
	  radius=atof(optarg);
	  break;

	case 'v':
	  be_verbose=TRUE;
	  break;

	case ':':
	  printf("Option %c needs a value\n",opt);
	  exit(2);
	  break;

	case '?':
	  //printf("Unknown option\n");
	  exit(3);
	  break;
	  
	}
    }

  if (be_verbose)
    {
      printf("Parsed options\n");
      printf("Dimension %d\n",xyzdim);
      printf("Porosity: %5.2f\n",porosity);
      printf("Radius %5.0f\n",radius);
      printf("Outputfile %s\n",output_file_name);
    }

  map=new MATRIXM(xyzdim,xyzdim);

  generate_plane(map, porosity,  radius, xyzdim);
  burn_2d(map,be_verbose);
  check_percolation_2d(map);
  calc_porosity(map);

  height=map->rows;
  width=map->cols;
  //bytes_per_pixel=BYTES_PER_PIXEL;
  image_data=new png_byte[height*width*BYTES_PER_PIXEL];
  for (k = 0; k < height; k++)
    {
      for (j = 0; j < width; j++)
	{
	  //printf("%ld %ld\n",k,j);
	  loc=k*width*BYTES_PER_PIXEL;
	  loc+=j*BYTES_PER_PIXEL;
	  if (map->matrix[k+1][j+1]==PERCOLATING_PORE_SPACE)
	    {	
	    image_data[loc+0]=COLOR_SCALE*(guint8)map->matrix[k+1][j+1];
	    image_data[loc+1]=COLOR_SCALE*(guint8)map->matrix[k+1][j+1];
	    image_data[loc+2]=COLOR_SCALE*(guint8)map->matrix[k+1][j+1];
	    }
	  else
	    {
	      image_data[loc]=0; //COLOR_SCALE*(guint8)map->matrix[k+1][j+1];
	      loc+=1;
	      image_data[loc]=0; //COLOR_SCALE*(guint8)map->matrix[k+1][j+1];
	      loc+=1;
	      image_data[loc]=0; //COLOR_SCALE*(guint8)map->matrix[k+1][j+1];
	    }
	}
    }
  printf("Writing image %s with dimension %dx%d\n",output_file_name,xyzdim,xyzdim);
  image->write_image(output_file_name,width,height,image_data,BYTES_PER_PIXEL);
  printf("Done\n");
  // delete image_data; // seg fault?

  return(0);

}
