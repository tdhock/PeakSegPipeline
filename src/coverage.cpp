/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <fstream>
#include "coverage.h"
#include <R.h> //for Rprintf

int coverage(char *in_filename, double *coverage_vec){
  std::ifstream bedGraph_file(in_filename);
  if(!bedGraph_file.is_open()){
    return ERROR_COVERAGE_CANT_OPEN_INPUT_FILE;
  }
  std::string line;
  int items, line_i=0;
  char chrom[100];
  char extra[100] = "";
  double count, chromStart, chromEnd, count_from_int;
  int count_int;
  coverage_vec[0] = 0.0;
  while(std::getline(bedGraph_file, line)){
    line_i++;
    items = sscanf
      (line.c_str(),
       "%s %lf %lf %lf%s\n",
       chrom, &chromStart, &chromEnd, &count, extra);
    if(items < 4){
      Rprintf("error: expected '%%s %%d %%d %%f\\n' on line %d\n%s\n", line_i, line.c_str());
      return ERROR_COVERAGE_UNEXPECTED_INPUT;
    }
    count_int = count;
    count_from_int = count_int;
    if(0 < strlen(extra) || count_from_int != count){
      Rprintf("error: non-integer coverage data on line %d\n%s\n", line_i, line.c_str());
      return ERROR_COVERAGE_NON_INTEGER;
    }
    if(coverage < 0){
      Rprintf("error: negative coverage on line %d\n%s\n", line_i, line.c_str());
      return ERROR_COVERAGE_NEGATIVE;
    }else{
      coverage_vec[0] += count*(chromEnd-chromStart);
    }
  }
  return 0;
}

