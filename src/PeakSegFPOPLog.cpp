/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <iomanip> //for setprecision.
#include <fstream> //for ifstream etc.
// http://docs.oracle.com/cd/E17076_05/html/programmer_reference/arch_apis.html
#include <dbstl_vector.h>
#include <R.h>

#include "funPieceListLog.h"
#include "PeakSegFPOPLog.h"

u_int32_t PiecewiseFunSize(const PiecewisePoissonLossLog&fun){
  return sizeof(PoissonLossPieceLog)*fun.piece_list.size() +
    sizeof(int)*2;
}

void PiecewiseFunCopy(void *dest, const PiecewisePoissonLossLog&fun){
  char *p = (char*)dest;
  int n_pieces = fun.piece_list.size();
  memcpy(p, &n_pieces, sizeof(int));
  p += sizeof(int);
  memcpy(p, &(fun.chromEnd), sizeof(int));
  p += sizeof(int);
  for(PoissonLossPieceListLog::const_iterator it = fun.piece_list.begin();
      it != fun.piece_list.end(); it++){
    memcpy(p, &(*it), sizeof(PoissonLossPieceLog));
    p += sizeof(PoissonLossPieceLog);
  }
}

void PiecewiseFunRestore(PiecewisePoissonLossLog&fun, const void *src){
  int n_pieces;
  char *p = (char*)src;
  PoissonLossPieceLog piece;
  memcpy(&n_pieces, p, sizeof(int));
  p += sizeof(int);
  memcpy(&(fun.chromEnd), p, sizeof(int));
  p += sizeof(int);
  for(int piece_i=0; piece_i < n_pieces; piece_i++){
    memcpy(&piece, p, sizeof(PoissonLossPieceLog));
    p += sizeof(PoissonLossPieceLog);
    fun.piece_list.push_back(piece);
  }
}

int PeakSegFPOP_disk(char *bedGraph_file_name, char* penalty_str){
  double penalty = atof(penalty_str);
  if(penalty == INFINITY){
    //ok but maybe we should special case this, no need to run PDPA.
  }else if(!std::isfinite(penalty)){
    return ERROR_PENALTY_NOT_FINITE;
  }else if(penalty < 0){
    return ERROR_PENALTY_NEGATIVE;
  }
  std::ifstream bedGraph_file(bedGraph_file_name);
  if(!bedGraph_file.is_open()){
    return ERROR_UNABLE_TO_OPEN_BEDGRAPH;
  }
  std::string line;
  int chromStart, chromEnd, coverage, items, line_i=0;
  char chrom[100];
  char extra[100] = "";
  double cum_weight_i = 0.0, cum_weight_prev_i, cum_weighted_count=0.0;
  double min_log_mean=INFINITY, max_log_mean=-INFINITY, log_data;
  int data_i = 0;
  double weight;
  int first_chromStart, prev_chromEnd;
  while(std::getline(bedGraph_file, line)){
    line_i++;
    items = sscanf
      (line.c_str(),
       "%s %d %d %d%s\n",
       chrom, &chromStart, &chromEnd, &coverage, extra);
    //Rprintf("%s %d %d %d%s\n", chrom, chromStart, chromEnd, coverage, extra);
    if(items < 4){
      return ERROR_NOT_ENOUGH_COLUMNS;
    }
    if(0 < strlen(extra)){
      return ERROR_NON_INTEGER_DATA;
    }
    weight = chromEnd-chromStart;
    cum_weight_i += weight;
    cum_weighted_count += weight*coverage;
    if(line_i == 1){
      first_chromStart = chromStart;
    }else{
      if(chromStart != prev_chromEnd){
	return ERROR_INCONSISTENT_CHROMSTART_CHROMEND;
      }
    }      
    prev_chromEnd = chromEnd;
    log_data = log(coverage);
    if(log_data < min_log_mean){
      min_log_mean = log_data;
    }
    if(max_log_mean < log_data){
      max_log_mean = log_data;
    }
  }
  double best_cost, best_log_mean, prev_log_mean;
  // open segments file for writing.
  std::string penalty_prefix = bedGraph_file_name;
  penalty_prefix += "_penalty=";
  penalty_prefix += penalty_str;
  std::string segments_file_name = penalty_prefix + "_segments.bed";
  std::ofstream segments_file;
  segments_file.open(segments_file_name.c_str());
  // Also write loss file.
  std::string loss_file_name = penalty_prefix + "_loss.tsv";
  std::ofstream loss_file;
  // Also db file.
  std::string db_file_name = penalty_prefix + ".db";
  loss_file.open(loss_file_name.c_str());
  if(penalty == INFINITY){
    if(cum_weighted_count != 0){
      best_cost = cum_weighted_count *
	(1 - log(cum_weighted_count) + log(cum_weight_i)); 
    } else {
      best_cost = 0;
    }
    segments_file << chrom << "\t" << first_chromStart << "\t" << chromEnd << "\tbackground\t" << cum_weighted_count/cum_weight_i << "\n";
    segments_file.close();
    loss_file << std::setprecision(20) << penalty_str << //penalty constant
      "\t" << 1 << //segments
      "\t" << 0 << //peaks
      "\t" << (int)cum_weight_i << //total bases
      "\t" << best_cost/cum_weight_i << //mean penalized cost
      "\t" << best_cost << //total un-penalized cost
      "\t" << "feasible" <<
      "\t" << 0 <<
      "\t" << 0 <<
      "\n";
    loss_file.close();
    return 0;
  }
  int data_count = line_i;
  //Rprintf("data_count=%d min_log_mean=%f max_log_mean=%f\n", data_count, min_log_mean, max_log_mean);
  //return 0;
  bedGraph_file.clear();
  bedGraph_file.seekg(0, std::ios::beg);

  // Both Berkeley DB Backends need to know how to serialize
  // the FPOP solver classes:
  dbstl::DbstlElemTraits<PiecewisePoissonLossLog> *funTraits =
    dbstl::DbstlElemTraits<PiecewisePoissonLossLog>::instance();
  funTraits->set_size_function(PiecewiseFunSize);
  funTraits->set_copy_function(PiecewiseFunCopy);
  funTraits->set_restore_function(PiecewiseFunRestore);

  //Berkeley DB filesystem Backend:
  DbEnv *env = NULL;
  Db *db = dbstl::open_db(env, db_file_name.c_str(), DB_RECNO, DB_CREATE, 0);
  dbstl::db_vector<PiecewisePoissonLossLog> cost_model_mat(db, env);

  //Berkeley DB in-memory backend:
  //dbstl::db_vector<PiecewisePoissonLossLog> cost_model_mat;

  //STL in-memory backend:
  //std::vector<PiecewisePoissonLossLog> cost_model_mat;

  //Initialization of empty function pieces.
  PiecewisePoissonLossLog foo;
  for(int i=0; i<data_count*2; i++){
    cost_model_mat.push_back(foo);
  }
  
  PiecewisePoissonLossLog up_cost, down_cost, up_cost_prev, down_cost_prev;
  PiecewisePoissonLossLog min_prev_cost;
  int verbose=0;
  cum_weight_i = 0;
  double total_intervals = 0.0, max_intervals = 0.0;
  while(std::getline(bedGraph_file, line)){
    items = sscanf(line.c_str(), "%*s\t%d\t%d\t%d\n", &chromStart, &chromEnd, &coverage);
    weight = chromEnd-chromStart;
    cum_weight_i += weight;
    // if(data_i < 10 || data_i > 1192280){
    //   Rprintf("data_i=%d weight=%f cum=%f coverage=%d\n",
    // 	     data_i, weight, cum_weight_i, coverage);
    // }
    if(data_i==0){
      // initialization Cdown_1(m)=gamma_1(m)/w_1
      down_cost.piece_list.emplace_back
	(1.0, -coverage, 0.0,
	 min_log_mean, max_log_mean, -1, false);
    }else{
      // if data_i is up, it could have come from down_cost_prev.
      min_prev_cost.set_to_min_less_of(&down_cost_prev, verbose);
      int status = min_prev_cost.check_min_of(&down_cost_prev, &down_cost_prev);
      if(status){
	Rprintf("BAD MIN LESS CHECK data_i=%d status=%d\n", data_i, status);
	min_prev_cost.set_to_min_less_of(&down_cost_prev, true);
	Rprintf("=prev down cost\n");
	down_cost_prev.print();
	Rprintf("=min less(prev down cost)\n");
	min_prev_cost.print();
	throw status;
      }
      // C^up_t(m) = (gamma_t + w_{1:t-1} * M^up_t(m))/w_{1:t}, where
      // M^up_t(m) = min{
      //   C^up_{t-1}(m),
      //   C^{<=}_down_{t-1}(m) + lambda/w_{1:t-1}
      // in other words, we need to divide the penalty by the previous cumsum,
      // and add that to the min-less-ified function, before applying the min-env.
      min_prev_cost.set_prev_seg_end(data_i-1);
      // cost + lambda * model.complexity =
      // cost + penalty * peaks =>
      // penalty = lambda * model.complexity / peaks.
      // lambda is output by exactModelSelection,
      // penalty is input by PeakSegFPOP.
      min_prev_cost.add(0.0, 0.0, penalty/cum_weight_prev_i);
      if(data_i==1){
	up_cost = min_prev_cost;
      }else{
	up_cost.set_to_min_env_of(&min_prev_cost, &up_cost_prev, verbose);
	status = up_cost.check_min_of(&min_prev_cost, &up_cost_prev);
	if(status){
	  Rprintf("BAD MIN ENV CHECK data_i=%d status=%d\n", data_i, status);
	  up_cost.set_to_min_env_of(&min_prev_cost, &up_cost_prev, true);
	  Rprintf("=prev down cost\n");
	  down_cost_prev.print();
	  Rprintf("=min less(prev down cost) + %f\n", penalty);
	  min_prev_cost.print();
	  Rprintf("=prev up cost\n");
	  up_cost_prev.print();
	  Rprintf("=new up cost model\n");
	  up_cost.print();
	  throw status;
	}
      }
      up_cost.multiply(cum_weight_prev_i);
      up_cost.add
	(weight,
	 -coverage*weight,
	 0.0);
      up_cost.multiply(1/cum_weight_i);
      //Rprintf("computing down cost\n");
      // compute down_cost.
      if(data_i==1){
	//for second data point, the cost is only a function of the
	//previous down cost (there is no first up cost).
	down_cost = down_cost_prev;
      }else{
	// if data_i is down, it could have come from up_cost_prev.
	// if(data_i==2329683){
	//   Rprintf("computing cost data_i=%d\n", data_i);
	//   verbose=1;
	// }else{
	//   verbose=0;
	// }
	min_prev_cost.set_to_min_more_of(&up_cost_prev, verbose);
	status = min_prev_cost.check_min_of(&up_cost_prev, &up_cost_prev);
	if(status){
	  Rprintf("BAD MIN MORE CHECK data_i=%d status=%d\n", data_i, status);
	  min_prev_cost.set_to_min_more_of(&up_cost_prev, true);
	  Rprintf("=prev up cost\n");
	  up_cost_prev.print();
	  Rprintf("=min more(prev up cost)\n");
	  min_prev_cost.print();
	  throw status;
	}
	min_prev_cost.set_prev_seg_end(data_i-1);
	//NO PENALTY FOR DOWN CHANGE
	down_cost.set_to_min_env_of(&min_prev_cost, &down_cost_prev, verbose);
	status = down_cost.check_min_of(&min_prev_cost, &down_cost_prev);
	if(status){
	  Rprintf("BAD MIN ENV CHECK data_i=%d status=%d\n", data_i, status);
	  down_cost.set_to_min_env_of(&min_prev_cost, &down_cost_prev, true);
	  Rprintf("=prev up cost\n");
	  up_cost_prev.print();
	  Rprintf("=min more(prev up cost)\n");
	  min_prev_cost.print();
	  Rprintf("=prev down cost\n");
	  down_cost_prev.print();
	  Rprintf("=new down cost model\n");
	  down_cost.print();
	  throw status;
	}
      }
      down_cost.multiply(cum_weight_prev_i);
      down_cost.add
	(weight,
	 -coverage*weight,
	 0.0);
      down_cost.multiply(1/cum_weight_i);
    }//if(data_i initialization else update
    cum_weight_prev_i = cum_weight_i;
    total_intervals += up_cost.piece_list.size() + down_cost.piece_list.size();
    if(max_intervals < up_cost.piece_list.size()){
      max_intervals = up_cost.piece_list.size();
    }
    if(max_intervals < down_cost.piece_list.size()){
      max_intervals = down_cost.piece_list.size();
    }
    up_cost_prev = up_cost;
    down_cost_prev = down_cost;
    //Rprintf("data_i=%d data_i+data_count=%d\n", data_i, data_i+data_count);
    up_cost.chromEnd = chromEnd;
    cost_model_mat[data_i] = up_cost;
    down_cost.chromEnd = chromEnd;
    cost_model_mat[data_i + data_count] = down_cost;
    data_i++;
  }
  //Rprintf("AFTER\n");
  // Decoding the cost_model_vec, and writing to the output matrices.
  int prev_seg_end;
  int prev_seg_offset = 0;
  // last segment is down (offset N) so the second to last segment is
  // up (offset 0).
  down_cost = cost_model_mat[data_count*2-1];
  down_cost.Minimize
    (&best_cost, &best_log_mean,
     &prev_seg_end, &prev_log_mean);
  //Rprintf("mean=%f end_i=%d chromEnd=%d\n", exp(best_log_mean), prev_seg_end, down_cost.chromEnd);
  prev_chromEnd = down_cost.chromEnd;
  // mean_vec[0] = exp(best_log_mean);
  // end_vec[0] = prev_seg_end;
  bool feasible = true;
  line_i=1;
  while(0 <= prev_seg_end){
    line_i++;
    // up_cost is actually either an up or down cost.
    up_cost = cost_model_mat[prev_seg_offset + prev_seg_end];
    //Rprintf("decoding prev_seg_end=%d prev_seg_offset=%d\n", prev_seg_end, prev_seg_offset);
    segments_file << chrom << "\t" << up_cost.chromEnd << "\t" << prev_chromEnd << "\t";
    // change prev_seg_offset for next iteration.
    if(prev_seg_offset==0){
      //up_cost is actually up
      prev_seg_offset = data_count;
      segments_file << "background"; // prev segment is down.
    }else{
      //up_cost is actually down
      prev_seg_offset = 0;
      segments_file << "peak";
    }
    segments_file << "\t" << exp(best_log_mean) << "\n";
    prev_chromEnd = up_cost.chromEnd;
    if(prev_log_mean != INFINITY){
      //equality constraint inactive
      best_log_mean = prev_log_mean;
    }else{
      feasible = false;
    }
    up_cost.findMean
      (best_log_mean, &prev_seg_end, &prev_log_mean);
    //Rprintf("mean=%f end=%d chromEnd=%d\n", exp(best_log_mean), prev_seg_end, up_cost.chromEnd);
  }//for(data_i
  segments_file << chrom << "\t" << first_chromStart << "\t" << prev_chromEnd << "\tbackground\t" << exp(best_log_mean) << "\n";
  segments_file.close();
  //Rprintf("feasible=%d\n", feasible);
  int n_peaks = (line_i-1)/2;
  loss_file << std::setprecision(20) << penalty << //penalty constant
    "\t" << line_i << //segments
    "\t" << n_peaks << //peaks
    "\t" << (int)cum_weight_i << //total bases
    "\t" << best_cost << //mean penalized cost
    "\t" << best_cost*cum_weight_i-penalty*n_peaks << //total un-penalized cost
    "\t"; 
  if(feasible){
    loss_file << "feasible";
  }else{
    loss_file << "infeasible";
  }
  loss_file <<
    "\t" << total_intervals/(data_count*2) <<
    "\t" << max_intervals <<
    "\n";
  loss_file.close();
  // need to close the database file, otherwise it takes up disk space
  // until R exists.
  db->close(0);
  return 0;
}

