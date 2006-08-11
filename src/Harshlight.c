/*
  compile it as a library, from the shell: R CMD SHLIB Harshlight.c it
  creates the file Harshlight.so From R then type:
  dyn.load("path_to_library/name_of_library")
  .C("name_of_the_function_called",
  as.type_arg1(arg1),as.type_arg2(arg2), DUP = FALSE)
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
//#include <malloc.h>
#include <math.h>
#include <unistd.h>
#include <R.h>
#define ERR3(x) { if(0) Rprintf("chk line %d [%s]\n",__LINE__, (x)); }

int ROW = 640;
int COL = 640;
int max_size = 0;
int c_id = -1;
int cluster_time = 1;
FILE *ps_out;
int ps_able;
int to_do_alloc = 16*640;
int *to_do_stack;
int to_do_items = 0;
int num_pages = 0;
//used in the recursive filling algorithm
//declared globally in order to avoid memory usage during the recursive step
int curr_x;
int curr_y;
int curr_i;
int curr_j;
int curr_spot;
int *curr_tmp;

time_t curr_time;

#define CUTOFF_1 1
#define NAN_SUB 100000
#define PSCANDO { if(!ps_able) return; }
#define PP fprintf(ps_out,
#define PS if(ps_able)fprintf(ps_out,
			      
/*
 *
 *
 *  Function definition
 *
 *
 *
*/

/* Functions used by the R program */

// Initialization of global variables

void init_harshlight(int *row_size, int *col_size, int *chips, int *status); //set the dimension of the analyzed chip

// Big defects

void extended_defects(double *list, double *med_obs, int *r, int* status); //calculate the median of a swliding window, used in the recognition of big defects, deals with boundary effects with matrix border duplication
double median(double *vector, int nelem, int valid_elem); //calculate the median
void sort_vector(double *vector, int nelem); //implementation of the quicksort algorithm to sort the values found in *vector array

// Small defects and clustering methods

void cluster_defects(int *array, int *obs_size, int *user_limit, int *connect, double *simulation_freq, double *user_alpha, int *type, int *status); //clusters the outliers in *array
int max(int a, int b); //max between two integers
int min(int a, int b); //min between two integers
void recursive_filling_four(int *curr_point, int *list, int *cluster_number, int *cluster_size, int *cut, int *status); //cluster points with a recursive algorithm, 4-neighbourhood
void recursive_filling_eight(int *curr_point, int *list, int *cluster_number, int *cluster_size, int *cut, int *status); //cluster points with a recursive algorithm, 8-neighbourhood
inline void add_stack(int *list, int *spot, int *cluster_number, int *cluster_size); //used by recursive_filling functions

// Diffuse defects

void diffuse_defects(double *list, double *quantile_bright, double *quantile_dark, int *r, double *p_bright, double *p_dark, double *q0, double *tresh_dark, double *tresh_bright, int *status); //diffuse_defects, duplicate borders
void trim_diagonal(int *list); //eliminates the points that have only one neighbour in their diagonal

// Image processing functions

void image_dilation(double *list, double *dil_list, int *r, int *status); //used to dilate images
void image_erosion(double *list, double *eros_list, int *r, int *status); //used to erase images

// Report file

void report_overall_header(char **fname, int *ext_rad, double *comp_q_br, double *comp_q_dr, int *comp_size, int *comp_conn, double *alpha, double *diff_br, double *diff_dr, double *diff_bin, int *diff_conn, int *diff_size, int *diff_rad, double *perc_empt, int *na_sub, int *interp, int *diff_close, char **chips, int *status); //overall header of the report
void chip_overall_header(); //general header for a chip
void openpsfile(char *fname); //open the postscript file
void grayimage_int(int wi, int he, int *gr); //prints a grayscale image
void chip_summary(double *var, int *num_c, int *num_d, double *perc_c, double *perc_d); //summary of chip's information
void chip_image(int *x, int *y, int *image, char **image_text, int *obs_size, int *obs_num, int *num_elem, int *type); //prints the original image, error image, compact defect image and diffuse defect image
void extended_stop(double *var); //called in case an extended defect is found

// General use functions

void free_memory(); //frees allocated memory
void simulations(int *phistogram, double *N, int *choice, int *status); //runs simulations of random matrices
int init_circular_mask(int **stack1, int radius, int col, int *status); // initializes arrays with the points inside a circular window of radius *r
int init_circ_mask_nocenter(int **stack1, int **stack2, int radius, int *status); //initializes queues with the points inside a circular window of radius *r (without its center)
double* duplicate_borders(double *list, int radius, int *status); //build a new matrix with duplicate borders for kernel issue at the boundaries

//Yupu's functions

void ErrorInt_row(double *list, int *size,int *status);
void norm(double *list, int *size,int *status);
int  handle_NA(double *list, int size);

void test(char **array, int *status);

void test(char **array, int *status){
  //Rprintf("%d\n", *(array + 1));
  *status = *status + 1;
  return;
}

int max(int a, int b){
  return (a > b ? a : b);
}

int min(int a, int b){
  return (a < b ? a : b);
}

void init_harshlight(int *row_size, int *col_size, int *chips, int *status){

  /* Set the dimension of the chip under analysis - default is 640X640 */

  ROW = *row_size;
  COL = *col_size;
  max_size = max(ROW, COL);
  to_do_alloc = 16*max_size;
  num_pages = *chips + 1;

  if ((to_do_stack = malloc(to_do_alloc * sizeof(int))) == NULL){
    Rprintf("init_harshlight: Cannot allocate memory!\n");
    Rprintf("to_do_alloc %d\n", to_do_alloc);
    fflush(stderr);
    *status = 1;
    return;
    //exit(EXIT_FAILURE);
  }

  ERR3( "opened harshlight" ); 

  return;

}

/* Functions */

inline void add_stack(int *list, int *spot, int *cluster_number, int *cluster_size){

  /* Add the spot to the cluster, increase the size of the cluster and add the spot to the queue */

  //Rprintf("b %d ", *(list + *spot));
  //fflush(stderr);
  *(list + *spot) = *cluster_number;
  //Rprintf("a %d\n", *(list + *spot));
  //fflush(stderr);
  *cluster_size += 1;
  to_do_stack[to_do_items++] = *spot;

  return;
  
}

void recursive_filling_four(int *curr_point, int *list, int *cluster_number, int *cluster_size, int *cut, int *status){

  /* Recursive flood-fill algorithm - 4-point neighbourhood */

  /*
   *curr_point          pixel in the array that is checked for neighbors
    *list               array of pixels
    *cluster_number      number of the cluster in which the curr_point is
    *cluster_size       size of the current cluster
    *cut                 if the cells around curr_point have value cut, they are considered part of the cluster
   */

  //int spot; //cell around curr_point that is evaluated
  //int curr_x; //coordinates of curr_point
  //int *tmp; //temporary pointer for realloc use

  ERR3(" Entering recursive FOUR, watch out!! "); 

  do{

  curr_x = (int)(*curr_point/COL);

  /* left of curr_point */
  curr_spot = *curr_point - 1;
  if((int)(curr_spot/COL) == curr_x && curr_spot >= 0 && *(list + curr_spot) == *cut){
    add_stack(list, &curr_spot, cluster_number, cluster_size);
  }

  /* right of curr_point */
  curr_spot = *curr_point + 1;
  if((int)(curr_spot/COL) == curr_x && curr_spot < ROW*COL && *(list + curr_spot) == *cut){
    add_stack(list, &curr_spot, cluster_number, cluster_size);
  }

  /* up of curr_point */
  curr_spot = *curr_point - COL;
  if(curr_spot >= 0 && *(list + curr_spot) == *cut){
    add_stack(list, &curr_spot, cluster_number, cluster_size);
  }

  /* down of curr_point */
  curr_spot = *curr_point + COL;
  if(curr_spot < ROW*COL && *(list + curr_spot) == *cut){
    add_stack(list, &curr_spot, cluster_number, cluster_size);
  }
  
  if(to_do_items > to_do_alloc - 50){
  ERR3( "Trying to reallocate memory" );
    if ((curr_tmp = realloc(to_do_stack, (to_do_alloc + 8*max_size) * sizeof(int))) == NULL){
      Rprintf("recursive_filling_four: Cannot reallocate memory!\n");
      //Rprintf("to_do_alloc %d max_size %d to_do_items %d\n", to_do_alloc, max_size, to_do_items);
      fflush(stderr);
      *status = 1;
      return;
    }
    to_do_stack = curr_tmp;
    to_do_alloc += 8*max_size;
  }
  
  //while(to_do_items){
  if(to_do_items){
    *curr_point = to_do_stack[--to_do_items];
    ERR3 ("decrementing to_do_items!"); 
    recursive_filling_four(curr_point, list, cluster_number, cluster_size, cut, status);
    if(*status){
      return;
    }
  }
  else{
    return;
  }
    //}
  }while(to_do_items >= 0);

  return;
}

void recursive_filling_eight(int *curr_point, int *list, int *cluster_number, int *cluster_size, int *cut, int *status){

  /* Recursive flood-fill algorithm - 8-point neighbourhood */

  /*
   *curr_point          pixel in the array that is checked for neighbors
    *list               array of pixels
    *cluster_number      number of the cluster in which the curr_point is
    *cluster_size       size of the current cluster
    *cut                 if the cells around curr_point have value cut, they are considered part of the cluster
   */

  //int curr_i, curr_j; //counters
  //curr_x, curr_y coordinates of curr_point, declared globally
  //int *curr_tmp; //temporary pointer for realloc use
  //int curr_spot;


  do{
  
    //Rprintf("coord ");
    //fflush(stderr);
    curr_x = (int)(*curr_point/COL);
    curr_y = *curr_point%COL;
    
    /* Check around curr_point */
    
    //Rprintf("x %d y %d\n", curr_x, curr_y);
    //fflush(stderr);
    
    for(curr_i = max(curr_x - 1, 0); curr_i <= min(curr_x + 1, ROW - 1); curr_i++){
      for(curr_j = max(curr_y - 1, 0); curr_j <= min(curr_y + 1, COL - 1); curr_j++){
	curr_spot = curr_i*COL + curr_j;
	if(*(list + curr_spot) == *cut){
	  add_stack(list, &curr_spot, cluster_number, cluster_size);
	}
      }
    }
    
    if(to_do_items > (to_do_alloc - 50)){
      ERR3( "Trying to reallocate memory\n");      
      //Rprintf("to_do_alloc %d max_size %d to_do_items %d\n", to_do_alloc, max_size, to_do_items);
      //fflush(stderr);
      
      curr_tmp = realloc(to_do_stack, (to_do_alloc + 8*max_size) * sizeof(int));
      if(curr_tmp == NULL){
	Rprintf("recursive_filling_eight: Cannot reallocate memory!\n");
	//Rprintf("to_do_alloc %d max_size %d to_do_items %d\n", to_do_alloc, max_size, to_do_items);
	fflush(stderr);
	*status = 1;
	//to_do_items = 0;
	return;
      }
      to_do_stack = curr_tmp;
      to_do_alloc += 8*max_size;
    }
    
    //while(to_do_items){
    if(to_do_items){
  
      if(0)if( to_do_items%1000 == 0) Rprintf("stack depth %d line %d\n", to_do_items,__LINE__); 
      //fflush(stderr);
      
      //Rprintf("to_do_items = %d\tto_do_stack = %d\tstatus = %d\n", to_do_items, to_do_stack[to_do_items - 1], *status);
      //fflush(stderr);
    
      *curr_point = to_do_stack[--to_do_items];
      
      //Rprintf("curr_point = %d cluster_number = %d cluster_size = %d cut = %d status = %d\n", *curr_point, *cluster_number, *cluster_size, *cut, *status);
      //fflush(stderr);

      //recursive_filling_eight(curr_point, list, cluster_number, cluster_size, cut, status);

      /* Handle recursive returns */
      
      //    if(*status){
      //return;
      //}

      //}
    }
    else{
      return;
    }
  }while(to_do_items >= 0);

  return;
}

			      /*
void recursive_filling_eight(int *curr_point, int *list, int *cluster_number, int *cluster_size, int *cut, int *status){
			      
  // Recursive flood-fill algorithm - 8-point neighbourhood

  
   *curr_point          pixel in the array that is checked for neighbors
    *list               array of pixels
    *cluster_number      number of the cluster in which the curr_point is
    *cluster_size       size of the current cluster
    *cut                 if the cells around curr_point have value cut, they are considered part of the cluster
   

  //int curr_i, curr_j; //counters
  //curr_x, curr_y coordinates of curr_point, declared globally
  //int *curr_tmp; //temporary pointer for realloc use
  //int curr_spot;

			        
  //Rprintf("coord ");
  //fflush(stderr);
  //curr_x = (int)(*curr_point/COL);
  //curr_y = *curr_point%COL;

  //Check around curr_point 
			      
  //Rprintf("x %d y %d\n", curr_x, curr_y);
  fflush(stderr);

  for(curr_i = max(curr_x - 1, 0); curr_i <= min(curr_x + 1, ROW - 1); curr_i++){
    for(curr_j = max(curr_y - 1, 0); curr_j <= min(curr_y + 1, COL - 1); curr_j++){
      curr_spot = curr_i*COL + curr_j;
      if(*(list + curr_spot) == *cut){
	add_stack(list, &curr_spot, cluster_number, cluster_size);
      }
    }
  }

  if(to_do_items > (to_do_alloc - 50)){
    ERR3( "Trying to reallocate memory\n");      
    //Rprintf("to_do_alloc %d max_size %d to_do_items %d\n", to_do_alloc, max_size, to_do_items);
    fflush(stderr);

    curr_tmp = realloc(to_do_stack, (to_do_alloc + 8*max_size) * sizeof(int));
    if(curr_tmp == NULL){
      Rprintf("recursive_filling_eight: Cannot reallocate memory!\n");
      Rprintf("to_do_alloc %d max_size %d to_do_items %d\n", to_do_alloc, max_size, to_do_items);
      fflush(stderr);
      *status = 1;
      //to_do_items = 0;
      return;
    }
    to_do_stack = curr_tmp;
    to_do_alloc += 8*max_size;
  }

  while(to_do_items){

    if(1)if( to_do_items%1000 == 0) Rprintf("stack depth %d line %d\n",to_do_items,__LINE__); 
    fflush(stderr);
    
    Rprintf("to_do_items = %d\tto_do_stack = %d\tstatus = %d\n", to_do_items, to_do_stack[to_do_items - 1], *status);
    fflush(stderr);
    
    *curr_point = to_do_stack[--to_do_items];
    
    Rprintf("curr_point = %d cluster_number = %d cluster_size = %d cut = %d status = %d\n", *curr_point, *cluster_number, *cluster_size, *cut, *status);
    fflush(stderr);

    recursive_filling_eight(curr_point, list, cluster_number, cluster_size, cut, status);
    
    // Handle recursive returns
			      
    if(*status){
      return;
    }
    
  }

  return;
}

*/

void sort_vector(double *vector, int nelem){

  /* Quicksort algorithm to sort a vector */
  /* modified from http://users.cs.cf.ac.uk/C.L.Mumford/tristan/QuickPage.html */

  int i, j;
  double start, tmp;

  if(nelem <= 1) return;

  start = vector[0];
  i = 0;
  j = nelem;
  do{
    do{
      i++;
    }while(vector[i] <start && i < nelem);

    do{
      j--;
    }while(vector[j] > start);

    if(i >= j) break;
    tmp = vector[i];
    vector[i] = vector[j];
    vector[j] = tmp;
  }while(i < j);

  tmp = vector[i-1];
  vector[i - 1] = vector[0];
  vector[0] = tmp;

  sort_vector(vector, i - 1);
  sort_vector(vector + i, nelem - i);
}

double median(double *vector, int nelem, int valid_elem){

  /*
    Calculates the median of the array *vector with nelem number of elements. The vector can include values that should not be considered
    (because, for example, they derive from missing values NA)
    valid_elem is the number of elements that must be considered in the array.
    Elements that do not have to be considered always have the highest number among all the elements, so that
    they are always at the end of the array when the array is sorted.
  */

  /* sort the array */
  sort_vector(vector, nelem);

  /* return the median */
  if(!(valid_elem%2)){
    return (double)((*(vector + valid_elem/2 - 1) + *(vector + valid_elem/2))/2);
  }
  else{
    return (double)(*(vector + valid_elem/2));
  }

}

void extended_defects(double *list, double *med_obs, int *r, int *status){

  /* This routine is used to detect big defects in the chip. *list is the error matrix, that is scanned with a sliding circular window, radius = *r.
     For each cell in the matrix, the median of the points included in the window centered in the cell is calculated and stored in *med_obs.
     To deal with border issues, the boundaries of the chip are duplicated.
  */

  /*
   *list      error matrix
   *med_obs   *list after convolution with a median kernel
   *r         radius used for the median kernel
  */

  int i, j, n; //counters
  int curr_point, curr_point_circle, area; //pixels that are analyzed inside the circle and number of pixels already analyzed inside the circle
  int num_points, *p_num_points; //number of points inside the median kernel
  int new_ROW, new_COL; //dimensions of the new image after border duplication
  double *circle_value; //value of the pixel inside *med_obs
  double *new_list; //*list with duplicated borders for window scanning
  int *mask; //the median kernel

  /* Build a new matrix with duplicate borders for window scanning */
  new_list = duplicate_borders(list, *r, status);
  if(*status){
    return;
  }

  new_ROW = ROW + 2*(*r); //do not initialise it as ROW + *r for the for loop, curr_point is dependent
  new_COL = COL + 2*(*r); //do not initialise it as COL + *r for the for loop, curr_point is dependent

  /* Initialize the two masks and inserts the values of the points present in the circular window */
  num_points = 0;
  p_num_points = &num_points;
  num_points = init_circular_mask(&mask, *r, new_COL, status);
  if(*status){
    return;
  }

  if ((circle_value = (double *)malloc(num_points*sizeof(double))) == NULL){
    Rprintf("extended_defects: Cannot allocate memory!\n");
    Rprintf("num_points %d\n", num_points);
    *status = 1;
    fflush(stderr);
    return;
    //exit(EXIT_FAILURE);
  }

  for(i = *r; i < new_ROW - *r; i++){
    for(j = *r; j < new_COL - *r; j++){
      
      area = 0; //stores the number of valid points inside the window (accounts for NaNs values in the matrix)
      curr_point = new_COL*i + j;

      for(n = 0; n < num_points; n++){

	area++;
	curr_point_circle = curr_point + mask[n];
	*(circle_value + n) = *(new_list + curr_point_circle);

	/* Account for NaN values in the matrix */
	if(isnan(*(circle_value + n))){
	  *(circle_value + n) = NAN_SUB;
	  area--;
	}
      }
      
      /* Calculate the median inside the window */

      *(med_obs + COL*(i - *r) + (j - *r)) = median(circle_value, num_points, area);
      
    }
  }
  
  /* Free the dynamically allocated memory */
  free(circle_value);
  free(mask);
  
  return;
  
}

void free_memory(){
  free(to_do_stack);
}

void cluster_defects(int *array, int *obs_size, int *user_limit, int *connect, double *simulation_freq, double *user_alpha, int *type, int *status)
{

  /* This function is used to flood-fill the *array image
     The number of cluster size is stored in *obs_size, but only if the cluster size is bigger than *user_limit, otherwise is not considered a cluster.
     *diagonal defines the neighbourhood: 0 = 4-neighbours, 1 = 8-neighbours
     *simulation_freq values of cluster size obtained through simulations
     *user_alpha statistical significance to be met in order to consider a cluster as a real cluster and not just occurring by chance only
     *type type of defects being analyzed: 
  */

  int i, j; //positions in the *array
  int cut; //cutoff to consider a pixel flagged
  int cluster_id; //cluster number
  int curr_point; //current point in the *array
  int cluster_size;//, *p_cluster_size; //size of the cluster found
  double alpha; //probability to find a cluster of the observed size by chance
  void (*p_recursive)(int*, int*, int*, int*, int*, int*) = NULL;

  to_do_items = 0;
  
  //p_cluster_size = &cluster_size;
  
  cluster_id = c_id;
  
  /* Eliminate the points that are lonely in the matrix, that is those points that are connected to only one other point in their diagonal */

  if(!*type){
    trim_diagonal(array);
  }
  if(!(*connect)){
    p_recursive = &recursive_filling_four;
  }
  else{
    p_recursive = &recursive_filling_eight;
  }

  /* Perform the scan throughout the whole image *array */
  for(i = 0; i < ROW; i++){
    for(j = 0; j < COL; j++){
      curr_point = COL*i + j;
      
      if(*(array + curr_point) == CUTOFF_1){
	cut = *(array + curr_point);
	*(array + curr_point) = cluster_id;
	cluster_size = 1;
	//cluster_size = 0;

	(*p_recursive)(&curr_point, array, &cluster_id, &cluster_size, &cut, status);
	if(*status){
	  return;
	}

	alpha = *(simulation_freq + cluster_size);
	  
	/* If the cluster is not big enough set the cluster_id to 0 */

	//*(obs_size + *p_cluster_size) += 1;
	if(cluster_size < *user_limit || alpha > *user_alpha){
	  //Rprintf("cluster_size %d\n", cluster_size);
	  //fflush(stderr);
	  *(array + curr_point) = 0;
	  (*p_recursive)(&curr_point, array, (array + curr_point), &cluster_size, &cluster_id, status);
	  if(*status){
	    return;
	  }
	    
	}
	else{
	  *(obs_size + cluster_size) += 1;
	  cluster_id--;
	}
      }
    }
  }
  
  if(cluster_time > 0){
    c_id = cluster_id;
  }
  else{
    c_id = -1;
  }
  cluster_time *= -1;
  
  return;
}

void trim_diagonal(int *list){

  /* Eliminates from the matrix *list the points that are connected only with one other pixel through the diagonals */

  int i, j, n; //positions of the pixle inside the *list and counters
  int nrow, ncol, nx, ny, circle_col, circle_points1, circle_points2, neighbours, curr_point, curr_point_circle; //parameters to define diagonals and number of neighbours for each pixel
  int num_diag; //parameters used while checking the diagonals
  int mask1[8], mask2[8], mask3[4], mask4[4]; //arrays in which to store the relative positions of diagonal pixels relative to the pixel under consideration
  /*
    mask1 contains the  position of the points in the mask, relative to a general point centered in (0,0)
    mask2 contains the values of the y positions of the points present in mask2, relative to its center (0,0)
    mask3 contains the  position of the points in the mask, relative to a general point centered in (0,0):
    the difference with mask1 is that the central point is not part of the mask, and only the diagonal points are considered
    mask4 contains the values of the y positions of the points present in mask4 relative to its center (0,0)
  */
  int mask_items1 = 0;
  int mask_items2 = 0;
  /*
    Nodes necessary to go through the values stored in the two masks
  */

  circle_col = 3; //points in the diameter = 2*r + 1, where r = 1

  /* Insert the values of the points present in the eight cells surrounding the cell under analysis */
  for(nrow = 0; nrow < circle_col; nrow++){
    for(ncol = 0; ncol < circle_col; ncol++){
      nx = nrow - 1;
      ny = ncol - 1;
      if(nx || ny){
	mask1[mask_items1] = COL*nx + ny;
	mask2[mask_items1++] = nx;
	if(nrow != 1 && ncol != 1){
	  mask3[mask_items2] = COL*nx + ny;
	  mask4[mask_items2++] = nx;
	}
      }
    }
  }

  circle_points1 = 8;
  circle_points2 = 4;

  /* Perform the pre-scan to eliminate the diagonal points in *list */

  /* Calculate the number of neighbours for each point in the matrix */
  for(i = 0; i < ROW; i++){
    for(j = 0; j < COL; j++){

      neighbours = 0;
      curr_point = COL*i + j;
      if(*(list + curr_point) == CUTOFF_1){
	
	for(n = 0; n < circle_points1; n++){
	  curr_point_circle = curr_point + mask1[n];
	  
	  if(curr_point_circle >= 0 && curr_point_circle < ROW*COL){
	    if((int)(curr_point_circle/COL) == (i + mask2[n])){
	      if(*(list + curr_point_circle) >= 1){
		neighbours++;
	      }
	    }
	  }
	}
	*(list + curr_point) = neighbours;
      }
    }
  }

  int counter = 0; //counter

  do{
    counter++;
    num_diag = 0;
    for(i = 0; i < ROW; i++){
      for(j = 0; j < COL; j++){
	curr_point = COL*i + j;
	
	/* Eliminate all points that have only 1 neighbour in the diagonal */
	if(*(list + curr_point) == 1){
	  
	  for(n = 0; n < circle_points2; n++){
	    curr_point_circle = curr_point + mask3[n];
	    
	    if(curr_point_circle >= 0 && curr_point_circle < ROW*COL){
	      if((int)(curr_point_circle/COL) == (i + mask4[n])){
		if(*(list + curr_point_circle) > 0){
		  *(list + curr_point) = 0;
		  *(list + curr_point_circle) = *(list + curr_point_circle) - 1;
		  num_diag = 1;
		  break;
		}
	      }
	    }
	  }
	}
      }
    }
  }while(num_diag);

  /* Old *list without the diagonal pixels */
  for(i = 0; i < ROW; i++){
    for(j = 0; j < COL; j++){
      curr_point = COL*j + i;
      if(*(list + curr_point) > 0){
	*(list + curr_point) = 1;
      }
    }
  }

}

void diffuse_defects(double *list, double *quantile_bright, double *quantile_dark, int *r, double *p_bright, double *p_dark, double *q0, double *tresh_dark, double *tresh_bright, int *status){

  /*
   *list                  image under analysis
   *r                     radius of the sliding window used to scan the chip
   *quantile_bright       value above which the pixel is considered an outlier
   *quantile_dark         value belowe which the pixel is considered an outlier
  */
  
  /*
    The function scans the *list image using a sliding window of radius *r.
    For each window, the number of outliers is calculated and a binomial test is performed. If the number of outliers is significant, the centre of the pixel is flagged as defected.
  */
  
  /* Duplicated borders for boundary effects */

  int i, j, curr_point, circle_points, n, curr_point_circle, num_quant_bright, num_quant_dark, curr_point_original, area;
  int new_ROW, new_COL; //dimensions of the new matrix after border duplication
  double cell_value, percent, num_bright, num_dark, limit_bright, limit_dark, limit_bright_tot, limit_dark_tot;
  double *new_list; //matrix with duplicated borders to deal with boundary issues

  int *mask1; //array with positions of points inside the window of radius *r

  curr_point = 0;

  /* Constructs another matrix with duplicated borders to deal with boundary effects */
  new_list = duplicate_borders(list, *r, status);
  if(*status){
    return;
  }

  new_ROW = ROW + 2*(*r); //do not initialise it as ROW + *r for the for loop, curr_point is dependent
  new_COL = COL + 2*(*r); //do not initialise it as COL + *r for the for loop, curr_point is dependent

  circle_points = init_circular_mask(&mask1, *r, new_COL, status);
  if(*status){
    return;
  }

  num_bright = *tresh_bright * (1 - *tresh_bright); //numerator for the p-null hypothesis (bright)
  num_dark = *tresh_dark * (1 - *tresh_dark); //numerator for the p_null hypothesis (dark)

  limit_bright_tot = *tresh_bright + (*q0 * sqrt(num_bright/circle_points)); //cutoff of the p-null hypothesis (bright) in case all the values inside the sliding window are valid
  limit_dark_tot = *tresh_dark + (*q0 * sqrt(num_dark/circle_points)); //cutoff of the p-null hypothesis (dark) in case all the values inside the sliding window are valid

  for(i = *r; i < new_ROW - *r; i++){
    for(j = *r; j < new_COL - *r; j++){

      area = 0;

      //Counts the number of outliers present in the window

      curr_point_original = COL*(i - *r) + (j - *r);
      curr_point = new_COL*i + j; //this is the reason why new_ROW and new_COL are not initialised like ROW + *r for the for loops

      num_quant_bright = 0;
      num_quant_dark = 0;

      for(n = 0; n < circle_points; n++){
	curr_point_circle = curr_point + mask1[n];
	cell_value = *(new_list + curr_point_circle);

	if(!isnan(cell_value)){
	  area++;// Account for NaN values in the matrix

	  if(cell_value >= *quantile_bright){
	    num_quant_bright++;
	  }
	  if(cell_value <= *quantile_dark){
	    num_quant_dark++;
	  }
	}

      }

      limit_bright = limit_bright_tot;
      limit_dark = limit_dark_tot;

      if(area != circle_points){
	limit_bright = *tresh_bright + (*q0 * sqrt(num_bright/area));
	limit_dark = *tresh_dark + (*q0 * sqrt(num_dark/area));
      }
      
      /* Fraction of outliers present in the window (bright) */
      percent = (double)((double)num_quant_bright/(double)area);
      *(p_bright + curr_point_original) = 0;

      /* If the fraction is statistically significant, the point is flagged */
      if(percent > limit_bright){
	*(p_bright + curr_point_original) = 1;
      }

      /* Fraction of outliers present in the window (dark) */
      percent = (double)((double)num_quant_dark/(double)area);
      *(p_dark + curr_point_original) = 0;

      /* If the fraction is statistically significant, the point is flagged */
      if(percent > limit_dark){
	*(p_dark + curr_point_original) = 1;
      }
    }
  }

  /* Free the dynamically allocated memory */  
  free(mask1);

  return;
}

void image_dilation(double *list, double *dil_list, int *r, int *status){

  /* Dilates an image passed though the array *list and stores the final result in the array *dil_list,
     using a circular kernel of radius *r */

  int *mask1, *mask2;
  int i, j, n, curr_point, curr_point_circle, circle_points;
  double cell_value;

  circle_points = init_circ_mask_nocenter(&mask1, &mask2, *r, status);
  if(*status){
    return;
  }

  for(i = 0; i < ROW; i++){
    for(j = 0; j < COL; j++){
      
      curr_point = COL*i + j;

      cell_value = 0;
      if(*(list + curr_point) >= 1){
	*(dil_list + curr_point) = 1;
	for(n = 0; n < circle_points; n++){
	  curr_point_circle = curr_point + mask1[n];

	  if(curr_point_circle >= 0 && curr_point_circle < ROW*COL){
	    if((int)(curr_point_circle/COL) == (i + mask2[n])){
	      *(dil_list + curr_point_circle) = 1;
	    }
	  }
	}
      }
    }
  }

  free(mask1);
  free(mask2);

}

void image_erosion(double *list, double *eros_list, int *r, int *status){

  /* Erodes an image passed though the array *list and stores the final result in the array *eros_list,
     using a circular kernel of radius *r */

  int *mask1, *mask2;
  int i, j, n, curr_point, curr_point_circle, circle_points;
  double cell_value;

  circle_points = init_circ_mask_nocenter(&mask1, &mask2, *r, status);
  if(*status){
    return;
  }

  for(i = 0; i < ROW; i++){
    for(j = 0; j < COL; j++){
      
      curr_point = COL*i + j;

      cell_value = 1;
      for(n = 0; n < circle_points; n++){
	curr_point_circle = curr_point + mask1[n];

	if(curr_point_circle >= 0 && curr_point_circle < ROW*COL){
	  if((int)(curr_point_circle/COL) == (i + mask2[n])){
	    if(*(list + curr_point_circle) == 0){
	      cell_value = 0;
	      break;
	    }
	  }
	}
      }
      *(eros_list + curr_point) = cell_value;
    }
  }
  
  free(mask1);
  free(mask2);

}

/* Initialise the circular masks with the coordinates of the points inside the circle, no problems with the boundaries */
int init_circular_mask(int **stack1, int radius, int col, int *status){

  int nrow, ncol, nx, ny, circle_col, r_square;
  int stack_items = 0;
  int area;
  

  circle_col = 2*radius + 1;
  r_square = radius*radius;
  area = 2*3.14*r_square + 50;

  if ((*stack1 = malloc(area * sizeof(int))) == NULL){
    Rprintf("init_circular_mask: Cannot allocate memory!\n");
    Rprintf("area %d\n", area);
    *status = 1;
    fflush(stderr);
    return(0);
    //exit(EXIT_FAILURE);
  }

  for(nrow = 0; nrow < circle_col; nrow++){
    for(ncol = 0; ncol < circle_col; ncol++){
      nx = nrow - radius;
      ny = ncol - radius;
      if((ny*ny) + (nx*nx) <= r_square){
	*(*stack1 + stack_items++) = col*nx + ny;
      }
    }
  }

  stack_items--;
  return stack_items;

}

/* Initialise the circular masks with the coordinates of the points inside the circle,
   without considering the center of the circle */
int init_circ_mask_nocenter(int **stack1, int **stack2, int radius, int *status){

  int nrow, ncol, nx, ny, circle_col, r_square;
  int stack_items = 0;
  int area;

  circle_col = 2*radius + 1;
  r_square = radius*radius;
  area = 2*3.14*r_square + 50;

  if ((*stack1 = malloc(area * sizeof(int))) == NULL){
    Rprintf("init_circ_mask_nocenter: Cannot allocate memory!\n");
    Rprintf("area stack1 %d\n", area);
    fflush(stderr);
    *status = 1;
    return(0);
    //exit(EXIT_FAILURE);
  }

  if ((*stack2 = malloc(area * sizeof(int))) == NULL){
    Rprintf("init_circ_mask_nocenter: Cannot allocate memory!\n");
    Rprintf("area stack2 %d\n", area);
    fflush(stderr);
    *status = 1;
    return(0);
    //exit(EXIT_FAILURE);
  }

  for(nrow = 0; nrow < circle_col; nrow++){
    for(ncol = 0; ncol < circle_col; ncol++){
      nx = nrow - radius;
      ny = ncol - radius;
      if(nx || ny){
	if((ny*ny) + (nx*nx) <= r_square){
	  *(*stack1 + stack_items) = COL*nx + ny;
	  *(*stack2 + stack_items++) = nx;
	}
      }
    }
  }
  stack_items--;

  return stack_items;

}

double* duplicate_borders(double *list, int radius, int *status){

  /* Builds a new matrix with the borders duplicated to deal with the boundary during kernel analysis */

  int area_new_list; //size of the new list
  int i, j; //counters
  int new_ROW, new_COL; //new size of the matrix
  int coord_row, coord_col; //position of the point in the new list in respect to the old list
  int curr_point = 0; //index of the cell in *list
  int newx_col, newy_col, add_val_row, add_val_col, change_row, change_col, cell_value;
  double *new_list; //new image with duplicated borders

  /* Calculates the dimensions for the new matrix */
  new_ROW = ROW + 2*radius;
  new_COL = COL + 2*radius;
  area_new_list = new_ROW * new_COL;
  
  /* Allocate the memory for the new matrix */
  if ((new_list = (double *)malloc(area_new_list*sizeof(double))) == NULL){
    Rprintf("duplicate_borders: Cannot allocate memory!\n");
    Rprintf("area_new_list %d\n", area_new_list);
    fflush(stderr);
    *status = 1;
    return(list);
    //exit(EXIT_FAILURE);
  }
  
  i = 0;
  j = 0;
  change_row = 0;

  newx_col = radius - 1;
  newy_col = newx_col;
  add_val_row = 0;

  for(i = 0; i < new_ROW; i++){

    newx_col = newx_col + add_val_row;
    cell_value = COL*newx_col + newy_col;
    *(new_list + new_COL*i) = *(list + cell_value); //elements in the first column of the new matrix
    coord_row = i - radius + 1;

    if(coord_row >= 0 && coord_row < ROW){
      if(!change_row){
	add_val_row = 0;
	change_row = 1;
      }
      else{
	add_val_row = 1;
      }
    }
    else{
      if(!change_row){
	add_val_row = -1;
      }
      else{
	add_val_row = 0;
	change_row = 0;
      }
    }

    change_col = 0;

    for(j = 1; j < new_COL; j++){
      coord_col = j - radius;

      /* Cells above or below the original matrix */
      if(coord_col >= 0 && coord_col < COL){
	if(!change_col){
	  add_val_col = 0;
	  change_col = 1;
	}
	else{
	  add_val_col = 1;
	}
      }

      /* Cells on the side of the original matrix */
      else{
	if(!change_col){
	  add_val_col = -1;
	}
	else{
	  add_val_col = 0;
	  change_col = 0;
	}
      }

      curr_point = new_COL*i + j;
      cell_value = cell_value + add_val_col;
      *(new_list + curr_point) = *(list + cell_value);
    }
  }

  return (new_list);

}

void simulations(int *phistogram, double *N, int *choice, int *status)
{
  
  /* The function is used to simulate during run-time the distribution of cluster size under the assumption that outliers are randomly distributed on the surface of the chip. 100.000 simulations are run. */
  /* *phistogram sotres the number of times at least one cluster of a certain size is found in each of the 100.000 randomized arrays with *N is the density of outliers and connected with a connectivity definition given by *choice (0 = 4-neighbourhood, 1 = 8-neighbourhood). */

  float x;
  int i, j, y, counter; //position of pixels in the array and counters
  int max_cluster_size;
  int *random_list; //simulated array
  int cluster_number, curr_point;
  int cluster;
  int cutoff = 1;
  void (*p_recursive)(int*, int*, int*, int*, int*, int*) = NULL;

  /* Try to allocate memory for the simulated array */
  if ((random_list = malloc(ROW*COL * sizeof(int))) == NULL){
    Rprintf("simulations: Cannot allocate memory!\n");
    Rprintf("ROW*COL %d\n", ROW*COL);
    fflush(stderr);
    *status = 1;
    return;
    //exit(EXIT_FAILURE);
  }
  
  //srandom( (unsigned)time( NULL ) );
  srand( (unsigned)time( NULL ) );

  if(!*choice){
    p_recursive = &recursive_filling_four;
  }
  else{
    p_recursive = &recursive_filling_eight;
  }

  Rprintf("Running simulations:\n");
  for(counter = 0; counter < 100000; counter++){
    max_cluster_size = 0;
    if(!(counter%1000)){
      Rprintf("%d of 100000\n", counter);
    }
    for (i = 0; i < ROW; i++){
      for(j = 0; j < COL; j++){
	//x = (float) random()/RAND_MAX;
	x = (float) rand()/RAND_MAX;
	if(x <= *N){
	  y =  1;
	}
	else{
	  y =  0;
	}
	random_list[COL*i + j] = y;
      }
    }    
    
    cluster_number = -1;

    for(i = 0; i <ROW; i++){
      for(j = 0; j < COL; j++){
	curr_point = COL*i + j;
	if(*(random_list + curr_point) == 1){
	  *(random_list + curr_point) = cluster_number;
	  cluster = 1;
	  (*p_recursive)(&curr_point, random_list, &cluster_number, &cluster, &cutoff, status);
	  if(*status){
	    free(random_list);
	    return;
	  }
	  if(cluster > max_cluster_size){
	    max_cluster_size = cluster;
	  }
	  cluster_number--;
	}
      }
    }
    
    for(i = 1; i <= max_cluster_size; i++){
      *(phistogram + i) = *(phistogram + i) + 1;
    }
  }
  
  //Rprintf("%d\n", (unsigned)time(NULL));
  free(random_list);

  return;
}

//Yupu's functions

int handle_NA(double *list, int size){
  int i;
  int valid = size;
  for(i = 0; i< size; i++){
    if(isnan(list[i])){
      list[i] = NAN_SUB;
      valid--;
    }
  }
  return valid;
}
      
void norm(double *list, int *size,int *status){
  
  int valid_count = handle_NA(list,*size);
  int i;
  double *sorted;

  if((sorted = (double *)malloc(*size*sizeof(double))) == NULL){
    Rprintf("norm: Cannot allocate memory for the sorted array!\n");
    Rprintf("size %d\n", *size);
    *status = 1;
    fflush(stderr);
    return;
  }

  for(i = 0; i < *size; i++){
    sorted[i] = list[i];
  }
  double med = median(sorted,*size,valid_count);
  
  for(i = 0; i < valid_count; i++){
    list[i] -= med;
  }
  free(sorted);
}

void ErrorInt_row(double *list, int *size, int *status){
  
  
  
  int valid_count = handle_NA(list,*size);
  int i;
  double *sorted;

  if((sorted = (double *)malloc(*size*sizeof(double))) == NULL){
    Rprintf("norm: Cannot allocate memory for the sorted array!\n");
    Rprintf("size %d\n", *size);
    *status = 1;
    fflush(stderr);
    return;
  }

  for(i = 0; i < *size; i++){
    sorted[i] = list[i];
  }
  double med = median(sorted,*size,valid_count);
  //Rprintf("median: %f\n",med);
  

  for(i = 0; i< *size; i++){
	list[i] = list[i] - med;
  }
  
  free(sorted);
  *status = 0;
  return;
}

//should worry about the NA point later
void ErrorInt_row2(double *list, int *size, int *status){
  double *sorted_row;
  int *rank_array;
  if((sorted_row = (double *)malloc(*size*sizeof(double))) == NULL){
    Rprintf("ErrorInt_row: Cannot allocate memory for the sorted array!\n");
    Rprintf("size %d\n", *size);
    *status = 1;
    fflush(stderr);
    return;
  }

  if((rank_array = (int *)malloc(*size*sizeof(int))) == NULL){
    Rprintf("ErrorInt_row: Cannot allocate memory for the index array!\n");
    Rprintf("size %d\n", *size);
    *status = 1;
    fflush(stderr);
    return;
  }
  int valid_count = handle_NA(list,*size);
 
  //copy the list
  int i;
  
  for(i = 0; i < *size; i++){
    sorted_row[i] = list[i];
  }
  
  //sort the list
  sort_vector(sorted_row,*size);
  
  
  //rank the list in running time O(n^2) 
  int j;
  for(i = 0; i < *size; i++){
    
    double ele = sorted_row[i];
    for(j =0; j < *size; j++){
      if(ele == list[j]){
	rank_array[j] = i;
      }
    }
  }
  
  //set the even boolean; get the 'middle'
  int magic;
  int even;
  if(valid_count %2){ 
    magic = (valid_count+1)/2-1;
    even = 0; 
  }
  else{
    magic = valid_count/2-1;
    even = 1;
  }
  
  /*test setting up value
    for(i = 0; i < *size;i++){
    list[i] = 10000;
    }
  */
  //get the error image
  for(i = 0; i< *size; i++){
    int index = rank_array[i]; //get rank
    if(index >= valid_count){
      
      //do nothing here for those NA elements
    }
    else if(even){
      //left wing
      if(index <= magic){
	list[i] = list[i] - sorted_row[magic+1];
      }
      //right wing
      else{
	list[i] = list[i] - sorted_row[magic];
      }
    }
    else{
      //middle point
      if(index == magic){
	list[i] =0;
      }
      //left wing
      else if(index < magic){
	list[i] = list[i]-((sorted_row[magic]+sorted_row[magic+1])/2);
      }
      //right wing
      else{
	list[i] = list[i] -((sorted_row[magic-1]+sorted_row[magic])/2);
      }
    }
  }
  
  free(sorted_row);
  free(rank_array);
  *status = 0;
  return;
}

/* REPORT */
/*********************
 *      612
 *
 * 
 * 7
 * 9
 * 2
 *
 *
 *
 *
 *********************/

void report_overall_header(char **fname, int *ext_rad, double *comp_q_br, double *comp_q_dr, int *comp_size, int *comp_conn, double *alpha, double *diff_br, double *diff_dr, double *diff_bin, int *diff_conn, int *diff_size, int *diff_rad, double *perc_empt, int *na_sub, int *interp, int *diff_close, char **chips, int *status){

  int x_text = 100;
  int y_text = 525;
  int x_text2 = x_text - 10;
  int x_text_val = x_text + 300;
  int c_conn = 8;
  int d_conn = 8;
  char interpol[] = "ON ";
  char closing_proc[] = "ON ";
  char sub_defect[] = "Median";

  if(!*comp_conn){
    c_conn = 4;
  }
  if(!*diff_conn){
    d_conn = 4;
  }
  if(!*interp){
    interpol[1] = 'F';
    interpol[2] = 'F';
  }
  if(!*diff_close){
    closing_proc[1] = 'F';
    closing_proc[2] = 'F';
  }
  if(*na_sub){
    sub_defect[0] = 'N';
    sub_defect[1] = 'A';
    sub_defect[2] = ' ';
    sub_defect[3] = ' ';
    sub_defect[4] = ' ';
    sub_defect[5] = ' ';
  }

  time(&curr_time);

  openpsfile(*fname);

  if(!ps_able)
    {
      Rprintf("pslib could not open, bailing out\n"); 
      *status = 1;
      fflush(stderr);
      return;
      //exit(1); 
    }

  PP "%%%%Pages: %d\n", num_pages);
  PP "%%%%EndComments\n");
  PP "%%%%Page: 1 1\n");
  PP "grestore gsave /Times-Roman findfont %d scalefont setfont\n %d %d moveto (Harshlight report: ", 12, 25, 755);
  PP "%s) show\n %d %d moveto (Version 1.1.2) show\n", asctime(localtime(&curr_time)), 490, 755);

  PP "%d %d moveto (EXTENDED defects:) show\n", x_text2, y_text - 25);
  PP "%d %d moveto (radius of the median kernel) show\n", x_text, y_text - 40);
  PP "%d %d moveto (%d pixels) show\n", x_text_val, y_text - 40, *ext_rad);

  PP "%d %d moveto (COMPACT defects::) show\n", x_text2, y_text - 65);
  PP "%d %d moveto (quantile for brigth outlier definition) show\n", x_text, y_text - 80);
  PP "%d %d moveto (%.3f) show\n", x_text_val, y_text - 80, *comp_q_br);

  PP "%d %d moveto (quantile for dark outlier definition) show\n", x_text, y_text - 95);
  PP "%d %d moveto (%.3f) show\n", x_text_val, y_text - 95, *comp_q_dr);

  PP "%d %d moveto (interpolation) show\n", x_text, y_text - 110);
  PP "%d %d moveto (%s) show\n", x_text_val, y_text - 110, interpol);

  PP "%d %d moveto (connectivity) show\n", x_text, y_text - 125);
  PP "%d %d moveto (%d-neighbourhood) show\n", x_text_val, y_text - 125, c_conn);

  PP "%d %d moveto (p_value for cluster size) show\n", x_text, y_text - 140);
  PP "%d %d moveto (%.3f) show\n", x_text_val, y_text - 140, *alpha);

  PP "%d %d moveto (minimum cluster size) show\n", x_text, y_text - 155);
  PP "%d %d moveto (%d pixels) show\n", x_text_val, y_text - 155, *comp_size);

  PP "%d %d moveto (minimum density) show\n", x_text, y_text - 170);
  PP "%d %d moveto (%.2f%%) show\n", x_text_val, y_text - 170, *perc_empt);

  PP "%d %d moveto (DIFFUSE defects) show\n", x_text2, y_text - 200);
  PP "%d %d moveto (percent of increase in intensity (bright outliers)) show\n", x_text, y_text - 215);
  PP "%d %d moveto (%.2f%%) show\n", x_text_val, y_text - 215, *diff_br);

  PP "%d %d moveto (percent of decrease in intensity (dark outliers)) show\n", x_text, y_text - 230);
  PP "%d %d moveto (%.2f%%) show\n", x_text_val, y_text - 230, *diff_dr);

  PP "%d %d moveto (p-value of the binomial test) show\n", x_text, y_text - 245);
  PP "%d %d moveto (%.3f) show\n", x_text_val, y_text - 245, *diff_bin);

  PP "%d %d moveto (radius of the circular mask in the binomial test) show\n", x_text, y_text - 260);
  PP "%d %d moveto (%d pixels) show\n", x_text_val, y_text - 260, *diff_rad);

  PP "%d %d moveto (connectivity) show\n", x_text, y_text - 275);
  PP "%d %d moveto (%d-neighbourhood) show\n", x_text_val, y_text - 275, d_conn);

  PP "%d %d moveto (minimium cluster size) show\n", x_text, y_text - 290);
  PP "%d %d moveto (%d pixels) show\n", x_text_val, y_text - 290, *diff_size);

  PP "%d %d moveto (closing procedure) show\n", x_text, y_text - 305);
  PP "%d %d moveto (%s) show\n", x_text_val, y_text - 305, closing_proc);

  PP "%d %d moveto (Defects substituted with) show\n", x_text2, y_text - 335);
  PP "%d %d moveto (%s) show\n", x_text_val, y_text - 335, sub_defect);

  PP "%d %d moveto (page 1 of %d) show\n", 25, 18, num_pages);
  PP "/Times-Roman findfont %d scalefont setfont\n %d %d moveto (PARAMETERS) show\n", 15, x_text + 150, y_text + 25);
  PP "%d %d moveto (Report for the chips %s) show\n", x_text - 50, y_text + 125, *chips);

  PP "newpath %d %d moveto %d %d lineto stroke\n", 20, 750, 550, 750);
  PP "newpath %d %d moveto %d %d lineto stroke\n", 20, 27, 550, 27);

  PP "showpage\n"); 
}

void chip_overall_header(int *chip_number, char **chip_name){

  PP "%%%%Page: %d %d\n", *chip_number + 1, *chip_number + 1);
  PP "grestore gsave /Times-Roman findfont %d scalefont setfont\n %d %d moveto (Harshlight report: ", 12, 25, 755);
  PP "%s) show\n %d %d moveto (Version 1.1.2) show\n", asctime(localtime(&curr_time)), 490, 755);
  PP "%d %d moveto (page %d of %d) show\n", 25, 18, *(chip_number) + 1, num_pages);
  PP "/Times-Roman findfont %d scalefont setfont\n %d %d moveto (Chip number %d, \"%s\") show\n ", 15, 200, 710, *chip_number, *chip_name);
  PP "newpath %d %d moveto %d %d lineto stroke\n", 20, 750, 550, 750);
  PP "newpath %d %d moveto %d %d lineto stroke\n", 20, 27, 550, 27);

}

void closepsfile(){
  if(ps_able){
    PP "%%EOF\n");
    fclose(ps_out);
    ps_able = 0;
  }
}

void chip_summary(double *var, int *num_c, int *num_d, double *perc_c, double *perc_d){

  int x_text = 50;
  int y_text = 150;
  int x_val_1 = x_text + 250;
  int x_val_2 = x_val_1 + 100;

  PP "grestore gsave /Times-Roman findfont %d scalefont setfont\n %d %d moveto (Chip summary:) show\n", 12, x_text, y_text);
  PP "%d %d moveto (Extended defects: the variance of the Error Image explained by the background is %.2f%%) show\n", x_text, y_text - 30, *var);

  PP "%d %d moveto (compact) show\n", x_text + 250, y_text - 60);
  PP "%d %d moveto (diffuse) show\n", x_text + 350, y_text - 60);

  PP "%d %d moveto (Number of clusters found:) show\n", x_text, y_text - 80);
  PP "%d %d moveto (%d) show\n", x_val_1, y_text - 80, *num_c);
  PP "%d %d moveto (%d) show\n", x_val_2, y_text - 80, *num_d);

  PP "%d %d moveto (Percent of the surface covered by the defects:) show\n", x_text, y_text - 100);
  PP "%d %d moveto (%.2f%%) show\n", x_val_1, y_text - 100, *perc_c);
  PP "%d %d moveto (%.2f%%) show\n", x_val_2, y_text - 100, *perc_d);
  PP "showpage\n");
}

void extended_stop(double *var){

  /* Called when an extended defect is found */

  int x_text = 50;
  int y_text = 150;

  PP "grestore gsave /Times-Roman findfont %d scalefont setfont\n %d %d moveto (Chip summary:) show\n", 12, x_text, y_text);
  PP "%d %d moveto (Extended defects: the variance of the Error Image explained by the background is %.2f%%) show\n", x_text, y_text - 30, *var);
  PP "%d %d moveto (The chip has been discarded due to the extended defect.) show\n", x_text, y_text - 50);
  PP "%d %d moveto (Please note: it is recommended to repeat the analysis after discarding the chip from the AffyBatch.) show\n", x_text, y_text - 70);
  PP "showpage\n");
}

void chip_image(int *x, int *y, int *image, char **image_text, int *obs_size, int *obs_num, int *num_elem, int *type){

  /* Prints the grayscale *image at position *x and *y on the page */

  int i;

  PS "grestore gsave %d %d translate 200 200 scale\n", *x, *y); 
  grayimage_int(ROW,COL,image);
  PP "grestore gsave /Times-Roman findfont %d scalefont setfont\n %d %d moveto (%s) show\n", 12, *x + 65, *y + 205, *image_text);
  if(*num_elem){
    for(i = 0; i < *num_elem; i++){
       PP "%% type %d %d %d\n", *type, *(obs_size + i), *(obs_num + i));
    }
  }
  else{
       PP "%% type %d 0 0\n", *type);
  }
}

void openpsfile(char *fname)
{
  ps_able=1;

  ps_out = fopen(fname,"w"); 

  if(ps_out == NULL){
    Rprintf("Could not open PS output file [%s] for writing\n", fname);
    fflush(stderr);
    ps_able=0; 
  }
  else{
    PP "%%!PS-Adobe-3.0\n");
  }
}

void grayimage_int(int wi, int he, int *gr)
{
/* uncompressed hexa bitmap of a grayscale 8bpp image */
int i,j;
                                                                                
PSCANDO;
                                                                                
PP "gsave\n/picstr %d string def\n%d %d 8\n",wi,wi,he);
PP "[ %d 0 0 %d 0 %d]\n",wi,-he,he);
PP "{ currentfile picstr readhexstring pop } image\n");
                                                                                
for(j=0;j<he;j++)
   {
   for(i=0;i<wi;i++)
      PP "%02x",*(gr++) & 0xff);
   PP "\n");
   }
PP "\n\ngrestore\n");
}
