#include <time.h>
#include <math.h>
#include "feature.h"
#include "parameters.h"
#include "decisiontree.h"
#include "treeaccessory.h"
#include "mex.h"

// matlab entry point
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Currently, we do not have any input or output
  if (nrhs != 0)
    mexErrMsgTxt("Wrong number of inputs"); 
  if (nlhs != 3)
    mexErrMsgTxt("Wrong number of outputs");

  // set the seed for random numbers
  srand(time(NULL));

  // set parameters for learning
  PARAMETER param;
  SetParameters(param);
  mexPrintf("Finish setting up the parameters!\n");

  // get image feature descriptors from the file
  vector<SINGLEIM> imVec;
  int classNum;
 
   GetData(imVec, "test_0.txt", param);
   
  
  AdjustLabels(imVec, classNum);
  mexPrintf("Finish getting feature descriptors from the file!\n");

  // pre-pooling for the background information, if it is necessary
//  if (param.bgPrePooling == 1)
//    BgPrePooling(imVec, param);

  // the pre-allocated memory for all the trees
  int treeNum = GetTreeNum();
  int treeMemSize = int(pow(2, param.maxTreeDepth + 1)) * treeNum;
  TREENODE *treeMemory = new TREENODE[treeMemSize];
  InitTreeMem(treeMemory, treeMemSize);
  mexPrintf("Finish initializing tree memory.\n");
  vector<TREENODE*> trees;  // the set of decision trees
  ReadTrees(trees, treeMemory, treeNum, param);
  treeNum = trees.size();  // the final number of valid trees
  //treeNum = 84;
  mexPrintf("Finish reading trees from the memory!\n");

  int out[3];
  out[0] = imVec.size();  // the number of images
  out[1] = classNum;      // the number of classes
  out[2] = treeNum;       // the number of trees
  // create the return matrix
  plhs[0] = mxCreateNumericArray(3, out, mxSINGLE_CLASS, mxREAL);
  float *result = (float*)mxGetPr(plhs[0]);
  memset(result, 0, sizeof(float) * out[0] * out[1] * out[2]);
  float *resultPTR = result;
  for (int i = 0; i < treeNum; i++)
  {
    mexPrintf("Evaluating using the %d-th tree.\n", i);
    // initial the result with that from the previous trees
    if (i > 0)
      memcpy(resultPTR, resultPTR - out[0] * out[1], sizeof(float) * out[0] * out[1]);

    TestDecisionTree(imVec, trees[i], resultPTR, param);
    resultPTR += out[0] * out[1];
  }

  // the second output is the class label
  int classOut[2];
  classOut[0] = imVec.size();
  classOut[1] = 1;
  plhs[1] = mxCreateNumericArray(2, classOut, mxINT32_CLASS, mxREAL);
  int *classLabel1 = (int*)mxGetPr(plhs[1]);
  plhs[2] = mxCreateNumericArray(2, classOut, mxINT32_CLASS, mxREAL);
  int *classLabel2 = (int*)mxGetPr(plhs[2]);
  for (int i = 0; i < imVec.size(); i++)
  {
    classLabel1[i] = imVec[i].classLbl + 1;
    classLabel2[i] = imVec[i].classLbl2 + 1;
  }

  ReleaseData(imVec);
  for (int i = 0; i < treeNum; i++)
    ReleaseTreeMem(trees[i]);
  delete[] treeMemory;
  
}
