/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "denormalize.h"
#include "coverage.h"
#include <R.h>
#include <R_ext/Rdynload.h>

void denormalize_interface
(char **bedGraph_file_vec, char **out_file_vec){
  int status = denormalize(bedGraph_file_vec[0], out_file_vec[0]);
  if(status==ERROR_CANT_OPEN_INPUT_FILE){
    error("Unable to open input file");
  }
  if(status==ERROR_UNEXPECTED_INPUT_DATA){
    error("Unexpected input data");
  }
  if(status==ERROR_NON_NUMERIC_COVERAGE){
    error("Coverage data (fourth column) should be numeric");
  }
  if(status==ERROR_NEGATIVE_COVERAGE){
    error("Coverage data (fourth column) should be non-negative");
  }
  if(status==ERROR_CANT_OPEN_OUT_FILE){
    error("Unable to open out file");
  }
  if(status != 0){
    error("error code %d", status);
  }
}
  
void coverage_interface
(char **bedGraph_file_vec, double *coverage_vec){
  int status = coverage(bedGraph_file_vec[0], coverage_vec);
  if(status==ERROR_COVERAGE_CANT_OPEN_INPUT_FILE){
    error("Unable to open input file");
  }
  if(status==ERROR_COVERAGE_UNEXPECTED_INPUT){
    error("Unexpected input data");
  }
  if(status==ERROR_COVERAGE_NON_INTEGER){
    error("Coverage data (fourth column) should be integer");
  }
  if(status==ERROR_COVERAGE_NEGATIVE){
    error("Coverage data (fourth column) should be non-negative");
  }
  if(status != 0){
    error("error code %d", status);
  }
}
  
R_CMethodDef cMethods[] = {
  {"denormalize_interface",
   (DL_FUNC) &denormalize_interface, 2
   //,{CHARSXP, CHARSXP}
  },
  {"coverage_interface",
   (DL_FUNC) &coverage_interface, 2
   //,{CHARSXP, CHARSXP}
  },
  {NULL, NULL, 0}
};

extern "C" {
  void R_init_PeakSegPipeline(DllInfo *info) {
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    //R_useDynamicSymbols call says the DLL is not to be searched for
    //entry points specified by character strings so .C etc calls will
    //only find registered symbols.
    R_useDynamicSymbols(info, FALSE);
  }
}
