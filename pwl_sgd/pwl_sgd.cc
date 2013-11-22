// pwl_sgd.cc alex berg march 2009

#include <string.h>
#include <math.h>
#include <mex.h>
#include <stdio.h>
#include <mclmcr.h>

int g_usesqrt=0;
int g_class_1_weight=0;
int g_verbosity = 1;

void vprint(const char *s) // wraps printing to control verbosity
{
  if (g_verbosity>0){
    mexPrintf(s);
  }
}

void encode(float *data, 
	    float *seg_data, 
	    int ndata, 
	    int ndims,
	    int nnz,
	    double *ind_i, 
	    double *ind_j, 
	    double *val)
{

  int count = 0;
  for (int i=0; i< ndata; i++){
    int dsofar = 0;

    for (int j=0; j< ndims; j++){
      float min_val = seg_data[0 + 4*j];
      float max_val = seg_data[1 + 4*j];
      float step = seg_data[2 + 4*j];
      float sqrt_step = sqrtf(step);
      float nsteps = seg_data[3 + 4*j];
      
      float v = data[i*ndims + j ];
      float msv = (v-min_val);
      float lval,uval,offset, omoffset;
      mxInt32 lind,uind;
      
      if (msv<=0){  // take case of lower end of range
	lind = 0;
	uind = 1;
	lval = 1;
	uval = 0;
      }else{
	lind = (mxInt32) (floor(msv/step));
	if (lind>=(nsteps-1)){
	  lind = nsteps-2;
	  uind = nsteps-1;
	  lval = 0;
	  uval = 1;
	}else{
	  offset = (msv-(lind*step))/step;
	  if (offset<0){
	    offset=0;
	  }else if (offset>1){
	    offset=1;
	  }
	  omoffset = 1-offset;
	  if (omoffset<0){
	    omoffset=0;
	  }else if (omoffset>1){
	    omoffset=1;
	  }
	  uind = lind+1;
	  lval = omoffset;
	  uval = offset;
	}
      }
      if (count>=nnz){
	mexPrintf("about to go past end of array in encoding data\n");
      }
      ind_i[count]=dsofar+lind+1;
      ind_j[count]=i+1;
      if (g_usesqrt>0){
	val[count]=lval*sqrt(step);
      }else{
	val[count]=lval;
      }
      count = count+1;

      if (count>=nnz){
	mexPrintf(" about to go past end of array in encoding data\n");
      }
      ind_i[count]=dsofar+uind+1;
      ind_j[count]=i+1;
      if (g_usesqrt>0){
	val[count]=uval*sqrt(step);
      }else{
	val[count]=uval;
      }
      count = count+1;
      dsofar = dsofar + nsteps;
    } 
  }
}


void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  srand48(0); //added for repeatability of code
  int IW_FLAG=0;
  int  mlen = 4097;
  char* message = (char *)malloc(sizeof(char)*mlen);
 
  if (!message){ mexErrMsgTxt("Memories... I have none.\n"); return; }
  


  int err=0;

  if (((nlhs!=1)&&(nlhs!=4)&&(nlhs!=6))|((nrhs!=10)&&(nrhs!=11))){
    mexPrintf("Usage : \n");
    mexPrintf("  [w [,ii,jj,vv [,delta,w_half ] ] ] = pwl_sgd( single(  data( ndims , nitems ) ),\n");
    mexPrintf("                                                single(   cat( 1, min_s( 1, ndims ) ,\n");
    mexPrintf("                                                                  max_s( 1, ndims ) ,\n");
    mexPrintf("                                                               seg_size( 1, ndims ) ,\n");
    mexPrintf("                                                                  nsegs( 1, ndims ) ) ),\n");
    mexPrintf("                                                single(labels( 1, ndims  ) ),\n");
    mexPrintf("                                                double(k), double(niters), double(lambda),\n");
    mexPrintf("                                                double(pwl_flag), double(init_iter), double(dbsqrt_flag),\n");
    mexPrintf("                                                double(class_1_weight)\n");
    mexPrintf("                                                [,double(i_w)] );\n");
    mexPrintf("\n");
    mexPrintf(" labels must be {-1, +1}\n");
    mexPrintf(" If you use an initial i_w and iters=0 then the output w will be\n");
    mexPrintf(" evaluation of the data against the input w.\n");
    mexPrintf(" your input had nlhs=%d  nrhs=%d\n",nlhs,nrhs);

    return;
  }

  

    // seg_vals ////////////////////////////////////////////////////////////
  if ((!mxIsSingle(prhs[0]))|(mxIsComplex(prhs[0]))){
    err = 1;
    vprint("data should be of type (real) single\n");
    return;
  }

  const int *sz = NULL;
  int ndims  = 0;
  int ndata  = 0;

  if (2!=mxGetNumberOfDimensions(prhs[0])){
    err = 1;
    vprint("data should be 2 dimensional\n");
    return;
  }else{
    sz = mxGetDimensions(prhs[0]);
    ndims  = sz[0];
    ndata  = sz[1];
  }

  float* data = (float *)mxGetData(prhs[0]);


  // seg_min_max_step_nsteps ///////////////////////////////////////////////////
  if ((!mxIsSingle(prhs[1]))|(mxIsComplex(prhs[1]))){
    err = 1;
    vprint("seg_min_max_step_nsteps should be of type (real) single\n");
    return;
  }

  const int *sz2 = NULL;
  int dummy_ndims  = 0;
  int dummy_4  = 0;

  if (2!=mxGetNumberOfDimensions(prhs[1])){
    err = 1;
    vprint("seg_data should be 2 dimensional\n");
    return;
  }else{
    sz2 = mxGetDimensions(prhs[1]);
    dummy_4  = sz2[0];
    dummy_ndims  = sz2[1];
  }
  
  if ((dummy_ndims!=ndims)|(dummy_4!=4)){
    err = 3;
    vprint("seg_data should be 4 x #data dimensions\n");
    return;
  }

  float* seg_data = (float *)mxGetData(prhs[1]);


  ///////////// LABELS ////////////////////////////////////////////////

  if ((!mxIsSingle(prhs[2]))|(mxIsComplex(prhs[2]))){
    err = 1;
    vprint("labels should be of type (real) single\n");
    return;
  }

  const int *sz3 = NULL;
  int dummy_1  = 0;
  int dummy_ndata  = 0;

  if (2!=mxGetNumberOfDimensions(prhs[2])){
    err = 1;
    vprint("labels should be 2 dimensional\n");
    return;
  }else{
    sz3 = mxGetDimensions(prhs[2]);
    dummy_1  = sz3[0];
    dummy_ndata  = sz3[1];
  }



  if ((dummy_ndata!=ndata)|(dummy_1!=1)){
    err = 3;
    mexPrintf("labes were %d by %d\n",dummy_1,dummy_ndata);
    mexPrintf("ndata is %d\n",ndata);
    vprint("labels should be 1 x #data\n");
    return;
  }

  float* labels = (float *)mxGetData(prhs[2]);

  //////////////// k /////////////
  if ((!mxIsDouble(prhs[3]))|(mxIsComplex(prhs[3]))){
    err = 1;
    vprint("k should be of type (real) double\n");
    return;
  }
  
  int k = (int)*(double *)mxGetData(prhs[3]);


  //////////////// niters /////////////
  if ((!mxIsDouble(prhs[4]))|(mxIsComplex(prhs[4]))){
    err = 1;
    vprint("niters should be of type (real) double\n");
    return;
  }
  
  int niters = (int)*(double *)mxGetData(prhs[4]);


  //////////////// lambda /////////////
  if ((!mxIsDouble(prhs[5]))|(mxIsComplex(prhs[5]))){
    err = 1;
    vprint("lambda should be of type (real) double\n");
    return;
  }
  
  double lambda = *(double *)mxGetData(prhs[5]);


  //////////////// pwl_flag /////////////
  if ((!mxIsDouble(prhs[6]))|(mxIsComplex(prhs[6]))){
    err = 1;
    vprint("pwl_flag should be of type (real) double\n");
    return;
  }
  
  int pwl_flag = ((int)*(double *)mxGetData(prhs[6]))>0;

  //////////////// init_iter /////////////
  if ((!mxIsDouble(prhs[7]))|(mxIsComplex(prhs[7]))){
    err = 1;
    vprint("init_iter should be of type (real) double\n");
    return;
  }
  
  int init_iter = (int)*((double *)mxGetData(prhs[7]));

  ////////////// use dbsqrt flag ////////////////////////
  if ((!mxIsDouble(prhs[8]))|(mxIsComplex(prhs[8]))){
    err = 1;
    vprint("dbsqrt_flag should be of type (real) double\n");
    return;
  }
  
  g_usesqrt = (int)(((int)*(double *)mxGetData(prhs[8]))>0);

  ////////////// class 1 weight  ////////////////////////
  if ((!mxIsDouble(prhs[9]))|(mxIsComplex(prhs[9]))){
    err = 1;
    vprint("class_1_weight should be of type (real) double\n");
    return;
  }
  
  g_class_1_weight = *((double *)mxGetData(prhs[9]));

  //////////////////////////////////////////////////////////////////////////////

  int tot_dims = 0;
  for (int i=0; i<ndims; i++){
    tot_dims = tot_dims + seg_data[3+i*4];
  }

  mexPrintf("====================================\n");
  mexPrintf("|       Training Parameters         |\n");
  mexPrintf("====================================\n");
  mexPrintf("      init_iter = %d\n",init_iter);
  mexPrintf("              k = %d\n",k);
  mexPrintf("           data = %d x %d\n", ndims, ndata);
  mexPrintf("       seg_data = %d x %d\n",dummy_4, dummy_ndims);
  mexPrintf("         labels = %d x %d\n",dummy_1, dummy_ndata);
  mexPrintf("         niters = %d\n",niters);
  mexPrintf("         lambda = %f\n",lambda);
  mexPrintf("       pwl_flag = %d\n",pwl_flag);
  mexPrintf("       use_sqrt = %d\n",g_usesqrt);
  mexPrintf(" class_1_weight = %d\n",g_class_1_weight);
  mexPrintf("     total dims = %d\n",tot_dims);
  mexPrintf("------------------------------------\n");
 



  double *i_w=NULL;
  if (11==nrhs){ // initialized w!!!

    IW_FLAG = 1;
    
    const int *sz_iw = NULL;
    dummy_1  = 0;
    int dummy_tot_dims  = 0;
    
    if ((2!=mxGetNumberOfDimensions(prhs[10]))|(!mxIsDouble(prhs[10]))){
      err = 1;
      vprint("input w should be 2 dimensional (tot_dims x 1) real double\n");
      return;
    }else{
      sz_iw = mxGetDimensions(prhs[10]);
      dummy_tot_dims  = sz_iw[0];
      dummy_1  = sz_iw[1];
    }
    

    
    if ((dummy_1!=1)|(dummy_tot_dims!=tot_dims)){
      err = 3;
      vprint("param 9==10 should be 2d tot_dims x 1\n");
      return;
    }
    printf("Using input w\n");
    i_w = (double *)mxGetData(prhs[10]);
    //    for (int i=0; i<tot_dims; i++){
    //      printf("%4.2f ",i_w[i]);
    //    }
    //    printf("\n");

  }else{
    IW_FLAG=0;
  }

  double *w = NULL;
  if ((0==niters)&&(IW_FLAG>0)){
    mexPrintf("niters = 0 so computing results of data against input i_w...\n");
    plhs[0]=mxCreateNumericMatrix( ndata, 1, mxDOUBLE_CLASS,  mxREAL); // w
    w = (double *)mxGetData(plhs[0]);
  }else{
    plhs[0]=mxCreateNumericMatrix( tot_dims, 1, mxDOUBLE_CLASS,  mxREAL); // w
    w = (double *)mxGetData(plhs[0]);
  }

  plhs[4]=mxCreateNumericMatrix( tot_dims, 1, mxDOUBLE_CLASS,  mxREAL); // w
  double *delta = (double *)mxGetData(plhs[4]);
//   double *delta = (double *)malloc(sizeof(double)*tot_dims);
//   if (NULL==delta){
//     printf("Out of memories delta\n");
//     return;
//   }
  plhs[5]=mxCreateNumericMatrix( tot_dims, 1, mxDOUBLE_CLASS,  mxREAL); // w
  double *w_half = (double *)mxGetData(plhs[5]);
//   double *w_half = (double *)malloc(sizeof(double)*tot_dims);
//   if (NULL==w_half){
//     printf("Out of memories w_half\n");
//     return;
//   }


  int nnz = ndims*ndata*2;
  plhs[1]=mxCreateNumericMatrix( nnz, 1, mxDOUBLE_CLASS,  mxREAL); // i,
  plhs[2]=mxCreateNumericMatrix( nnz, 1, mxDOUBLE_CLASS,  mxREAL); // j,
  plhs[3]=mxCreateNumericMatrix( nnz, 1, mxDOUBLE_CLASS, mxREAL); // value
  double *ind_i = (double *)mxGetData(plhs[1]);
  double *ind_j = (double *)mxGetData(plhs[2]);
  double *val = (double *)mxGetData(plhs[3]);


  printf(" encoding data...\n");

  encode(data, seg_data, ndata, ndims, nnz, ind_i, ind_j, val);



  if ((0==niters)&&(IW_FLAG>0)){
    printf("niters = 0 so computing results of data against input i_w...\n");
    for (int ind=0; ind<ndata; ind++){
      double res=0;
      int doff = ind*ndims*2;
      for (int i=0; i< ndims; i++){
	int ii = (int)(ind_i[2*i+doff]); 
	int ii1 = (int)(ind_i[2*i+doff+1]); 
	res = res + i_w[ii-1]*val[2*i+doff] + i_w[ii1-1]*val[2*i+1+doff];
      }
      w[ind]=res;
    }
    return;
  }



  printf(" starting sub gradient descent...\n");
  

  double nw=0;


  /// BEGIN -- NORMALIZATION BLOCK

  if (IW_FLAG>0){
    printf("using input w\n");
    for (int i=0; i<tot_dims; i++){
      w[i]=i_w[i];
    }
  }else{
    printf("initializing w\n");
    for (int i=0; i<tot_dims; i++){
      w[i]=2*drand48()-1;
    }
    nw=0;
    if (pwl_flag>0){
      int dc=0;
      for (int j=0; j<ndims; j++){
	dc=dc+1;
	for (int i=1; i<seg_data[3+j*4]; i++){
	  nw = nw+(w[dc]-w[dc-1])*(w[dc]-w[dc-1]);
	  dc=dc+1;
	}
      }
    }else{
      for (int i=0; i<tot_dims; i++){
      nw = nw + w[i]*w[i];
      }
    }

    nw = sqrt(nw)*sqrt(lambda);
    //  nw = sqrt(nw);
    for (int i=0; i<tot_dims; i++){
      w[i]=w[i]/nw;
    }

  }

  /// END  -- NORMALIZATION BLOCK




  printf("ndata = %d\n",ndata);

  //// BEGIN MAIN S.G.D. LOOP


  double sum_weights=0;
  int nwrong=0;
  double res=0;
  double step_size=0;
  for (int iter=0; iter<niters; iter++){
    step_size = 1/(lambda*((double)iter+(double)init_iter+1.0));
    for (int i=0; i< tot_dims; i++){
      delta[i]=0;
    }
    nwrong =0;
    sum_weights=0;
    for (int j=0; j<k; j++){
      // pick random data item
      int ind = (int)round((drand48()*ndata)-0.5);
      if (ind<0){ ind = 0; }
      //   if (ind>=(ndata-1)){ ind = ndata-1; }
   if (ind>=(ndata-1)){ ind = ndata-1; }
      // compute output
   //   printf("%d ",ind);
      res=0;
      int doff = ind*ndims*2;
      for (int i=0; i< ndims; i++){
	int ii = (int)(ind_i[2*i+doff]); 
	int ii1 = (int)(ind_i[2*i+doff+1]); 
	res = res + w[ii-1]*val[2*i+doff] + w[ii1-1]*val[2*i+1+doff];
      }
      if (labels[ind]>0){
	sum_weights = sum_weights+g_class_1_weight;
      }else{
	sum_weights = sum_weights+1;
      }

      if ((res*labels[ind])<1){  // was wrong
	nwrong = nwrong+1;
	if (labels[ind]>0){
	  for (int i=0; i<ndims; i++){
	    int ii = (int)(ind_i[2*i+doff]); 
	    int ii1 = (int)(ind_i[2*i+doff+1]); 
	    delta[ii-1]+= val[2*i+doff]*labels[ind]*g_class_1_weight;
	    delta[ii1-1]+= val[2*i+doff+1]*labels[ind]*g_class_1_weight;
	  }
	}else{
	  for (int i=0; i<ndims; i++){
	    int ii = (int)(ind_i[2*i+doff]); 
	    int ii1 = (int)(ind_i[2*i+doff+1]); 
	    delta[ii-1]+= val[2*i+doff]*labels[ind];
	    delta[ii1-1]+= val[2*i+doff+1]*labels[ind];
	  }
	}
	//	for (int i=0; i<tot_dims; i++){
	//	  w[i]=delta[i];
	//	}
	//	printf("ind = %d\n",ind);
	//	return;
      }
    }

// 		for (int i=0; i<tot_dims; i++){
// 		  w[i]=delta[i];
// 		}
//		return;

    printf("it = %10d       \t nwrong = %5d  \t step_size = %4.6f\r",iter,nwrong,step_size);
    if (nwrong>0){
      // compute w_half
      


      double nw_half=0;
      double om_step_size_lambda = 1-step_size*lambda;
      double step_size_lambda = step_size*lambda;
      // double step_size_k = step_size/(double)k;
 double step_size_k = step_size/(double)sum_weights;
 // printf("sum_weights = %4.2f k = %f \n",sum_weights,(double)k);
      if (0==pwl_flag){
	for (int i=0; i<tot_dims; i++){
	  w_half[i]=om_step_size_lambda*w[i]+ step_size_k*delta[i];
	  nw_half = nw_half+w_half[i]*w_half[i];
	}
	//printf("linear normalization\n");
      }else{
	int dc=0;
	for (int j=0; j<ndims; j++){
	  w_half[dc]=w[dc]-w[dc+1];
	  dc=dc+1;
	  for (int i=1; i<seg_data[3+j*4]-1; i++){
	    w_half[dc] = 2*w[dc]-w[dc-1]-w[dc+1];
	    dc=dc+1;
	  }
	  w_half[dc]=w[dc]-w[dc-1];
	  dc=dc+1;
	}
	//	if (7000==iter){
	//	  printf("debugging break\n");
	//	  return;
	//	}
	//	printf("last used dc was %d\n",dc-1);
// 	w_half[0] = w[0]-w[1];
// 	for (int i=1; i<(tot_dims-1); i++){
// 	  w_half[i] = 2*w[i]-w[i-1]-w[i+1];
// 	}
// 	w_half[tot_dims-1] = w[tot_dims-1]-w[tot_dims-2];

	for (int i=0; i<tot_dims; i++){
	  w_half[i]=w[i]-step_size_lambda*w_half[i]+ step_size_k*delta[i];
	}
	dc=0;
	for (int j=0; j<ndims; j++){
	  dc=dc+1;
	  for (int i=1; i<seg_data[3+j*4]; i++){
	    nw_half = nw_half+(w_half[dc]-w_half[dc-1])*(w_half[dc]-w_half[dc-1]);
	    dc=dc+1;
	  }
	}
	//	printf("PW linear normalization\n");
      }
      //      printf("nw_half=%4.2f\n",nw_half);
      nw_half = 1/(sqrt(nw_half)*sqrt(lambda));
      //      nw_half = 1/(sqrt(nw_half));
      //       printf("nw_half=%4.2f\n",nw_half);
      //nw_half = 1/(sqrt(nw_half));
      if (nw_half<1){
	//	printf("normalizing\n");
	for (int i=0; i<tot_dims; i++){
	  w[i]=w_half[i]*nw_half;
	}
	//	return;
      }else{
	for (int i=0; i<tot_dims; i++){
	  w[i]=w_half[i];
	}
      }
      /// END  -- NORMALIZATION BLOCK
      
    }
    if (((iter%5000)==0)|(niters==(iter+1))){
      int tot_wrong = 0;
      for (int mind=0; mind<ndata ; mind++){
	res=0;
	int doff = mind*ndims*2;
	for (int i=0; i< ndims; i++){
	  int ii = (int)(ind_i[2*i+doff]); 
	  int ii1 = (int)(ind_i[2*i+doff+1]); 
	  res = res + w[ii-1]*val[2*i+doff] + w[ii1-1]*val[2*i+1+doff];
	}
	if ((res*labels[mind])<1){
	  tot_wrong = tot_wrong+1;
	}
      }
      printf("\n\niter %d # over margin = %d\n",iter,tot_wrong);
      if (0==tot_wrong){
	return;
      }
    }
  }

  printf("\n");
  //  free(delta);
  //  free(w_half);
}


/*

Some useful functions from matrix.h
-----------------------------------

*** for 2d arrays ***

mxArray *mxCreateNumericMatrix(mwSize m, mwSize n, 
  mxClassID classid, mxComplexity ComplexFlag);

to get pointer to beginning of array (double *)
double *mxGetPr(const mxArray *pm);  

Get pointer to character array data
mxChar *mxGetChars(const mxArray *array_ptr);

Get pointer to array data
void *mxGetData(const mxArray *pm);

mxGetProperty returns the value at pa[index].propname
mxArray *mxGetProperty(const mxArray *pa, mwIndex index,
         const char *propname);

mxClassID mxGetClassID(const mxArray *pm);
const char *mxGetClassName(const mxArray *pm);

size_t mxGetM(const mxArray *pm);
size_t mxGetN(const mxArray *pm);


*** for nd array ***

mxArray *mxCreateNumericArray(mwSize ndim, const mwSize *dims, 
         mxClassID classid, mxComplexity ComplexFlag);

mwSize mxGetNumberOfDimensions(const mxArray *pm);

Use mxGetNumberOfDimensions to determine how many dimensions are in
the specified array. To determine how many elements are in each
dimension, call mxGetDimensions.


const mwSize *mxGetDimensions(const mxArray *pm);

The address of the first element in the dimensions array. Each integer
in the dimensions array represents the number of elements in a
particular dimension. The array is not NULL terminated.


*** for strings (in array single byte per char ) ***

int mxGetString(const mxArray *pm, char *str, mwSize strlen);


*** enums ***

typedef enum {
    mxREAL,
    mxCOMPLEX
} mxComplexity;


typedef enum {
	mxUNKNOWN_CLASS = 0,
	mxCELL_CLASS,
	mxSTRUCT_CLASS,
	mxLOGICAL_CLASS,
	mxCHAR_CLASS,
	mxVOID_CLASS,
	mxDOUBLE_CLASS,
	mxSINGLE_CLASS,
	mxINT8_CLASS,
	mxUINT8_CLASS,
	mxINT16_CLASS,
	mxUINT16_CLASS,
	mxINT32_CLASS,
	mxUINT32_CLASS,
	mxINT64_CLASS,
	mxUINT64_CLASS,
	mxFUNCTION_CLASS,
        mxOPAQUE_CLASS,
	mxOBJECT_CLASS
} mxClassID;


*/

