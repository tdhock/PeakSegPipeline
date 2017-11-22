/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "PeakSegFPOPLog.h"
#include "denormalize.h"
#include <R.h>
#include <R_ext/Rdynload.h>

void PeakSegFPOP_interface
(char **file_vec, char **pen_vec){
  char *bedGraph = file_vec[0];
  char *penalty = pen_vec[0];
  int status = PeakSegFPOP_disk(bedGraph, penalty);
  if(status==ERROR_PENALTY_NOT_FINITE){
    error("penalty=%s but must be finite", penalty);
  }
  if(status==ERROR_PENALTY_NEGATIVE){
    error("penalty=%s must be non-negative", penalty);
  }
  if(status==ERROR_UNABLE_TO_OPEN_BEDGRAPH){
    error("unable to open input file for reading %s", bedGraph);
  }
  if(status==ERROR_NOT_ENOUGH_COLUMNS){
    error("each line of input data file %s should have exactly four columns", bedGraph);
  }
  if(status==ERROR_NON_INTEGER_DATA){
    error("fourth column of input data file %s should be integer", bedGraph);
  }
  if(status==ERROR_INCONSISTENT_CHROMSTART_CHROMEND){
    error("there should be no gaps (columns 2-3) in input data file %s", bedGraph);
  }
  if(status != 0){
    error("error code %d", status);
  }
}
  
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
  
R_CMethodDef cMethods[] = {
  {"PeakSegFPOP_interface",
   (DL_FUNC) &PeakSegFPOP_interface, 2
   //,{REALSXP, REALSXP, INTSXP, INTSXP, REALSXP}
  },
  {"denormalize_interface",
   (DL_FUNC) &denormalize_interface, 2
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
