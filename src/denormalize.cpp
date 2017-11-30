/* -*- compile-command: "R CMD INSTALL .." -*- */
 
#include <fstream>
#include <math.h> //for INFINITY, round
#include "denormalize.h"
#include <R.h> //for Rprintf

int denormalize(char *in_filename, char *out_filename){//data_count x 2
  std::ifstream bedGraph_file(in_filename);
  if(!bedGraph_file.is_open()){
    return ERROR_CANT_OPEN_INPUT_FILE;
  }
  std::string line;
  int chromStart, chromEnd, items, line_i=0, n_zeros=0;
  char chrom[100];
  char extra[100] = "";
  double coverage, min_coverage = INFINITY;
  while(std::getline(bedGraph_file, line)){
    line_i++;
    items = sscanf
      (line.c_str(),
       "%s %d %d %lf%s\n",
       chrom, &chromStart, &chromEnd, &coverage, extra);
    //printf("%s %d %d %d%s\n", chrom, chromStart, chromEnd, coverage, extra);
    if(items < 4){
      Rprintf("error: expected '%%s %%d %%d %%f\\n' on line %d\n%s\n", line_i, line.c_str());
      return ERROR_UNEXPECTED_INPUT_DATA;
    }
    if(0 < strlen(extra)){
      Rprintf("error: non-numeric coverage data on line %d\n%s\n", line_i, line.c_str());
      return ERROR_NON_NUMERIC_COVERAGE;
    }
    if(0 == coverage){
      n_zeros++;
    }else{
      if(coverage < 0){
	Rprintf("error: negative coverage on line %d\n%s\n", line_i, line.c_str());
	return ERROR_NEGATIVE_COVERAGE;
      }else{
	if(coverage < min_coverage){
	  min_coverage = coverage;
	}
      }
    }
  }
  bedGraph_file.clear();
  bedGraph_file.seekg(0, std::ios::beg);
  Rprintf("min_coverage=%f n_zeros=%d/%d\n", min_coverage, n_zeros, line_i);
  std::ofstream out_file;
  out_file.open(out_filename);
  if(!out_file.is_open()){
    return ERROR_CANT_OPEN_OUT_FILE;
  }
  line_i=0;
  while(std::getline(bedGraph_file, line)){
    line_i++;
    items = sscanf
      (line.c_str(),
       "%s %d %d %lf\n",
       chrom, &chromStart, &chromEnd, &coverage);
    int int_coverage = round(coverage/min_coverage);
    out_file << chrom <<
      "\t" << chromStart <<
      "\t" << chromEnd <<
      "\t" << int_coverage <<
      "\n";
  }
  return 0;
}

