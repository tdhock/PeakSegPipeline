/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "funPieceListLog.h"
#include <list>
#include <math.h>
#include <stdio.h>
#include <R.h>

#define NEWTON_EPSILON 1e-12
#define NEWTON_STEPS 100
#define PREV_NOT_SET (-3)

#define ABS(x) ((x)<0 ? -(x) : (x))

PoissonLossPieceLog::PoissonLossPieceLog
(double li, double lo, double co, double m, double M, int i, double prev){
  Linear = li;
  Log = lo;
  Constant = co;
  min_log_mean = m;
  max_log_mean = M;
  data_i = i;
  prev_log_mean = prev;
}

PoissonLossPieceLog::PoissonLossPieceLog(){
}

bool PoissonLossPieceLog::has_two_roots(double equals){
  // are there two solutions to the equation Linear*e^x + Log*x +
  // Constant = equals ?
  if(Log == 0){//degenerate linear function.
    // f(mean) = Linear*mean + Constant - equals
    throw "problem";
    return false;
  }
  double optimal_mean = argmin_mean(); //min or max!
  double optimal_log_mean = log(optimal_mean); //min or max!
  double optimal_cost = getCost(optimal_log_mean);
  double optimal_cost2 = PoissonLoss(optimal_mean);
  // does g(x) = Linear*e^x + Log*x + Constant = equals? compare the
  // cost at optimum to equals.
  // g'(x)= Linear*e^x + Log,
  // g''(x)= Linear*e^x.
  if(0 < Linear){//convex.
    return optimal_cost + NEWTON_EPSILON < equals && optimal_cost2 + NEWTON_EPSILON < equals;
  }
  //concave.
  return equals + NEWTON_EPSILON < optimal_cost && equals + NEWTON_EPSILON < optimal_cost2;
}

double PoissonLossPieceLog::PoissonLoss(double mean){
  double loss_without_log_term = Linear*mean + Constant;
  if(Log==0){
    return loss_without_log_term;
  }
  double log_mean_only = log(mean);
  double log_coef_only = Log;
  double product = log_mean_only * log_coef_only;
  return loss_without_log_term + product;
}

double PoissonLossPieceLog::PoissonDeriv(double mean){
  return Linear + Log/mean;
}

// This function performs the root finding in the positive (mean)
// space, but needs to return a value in the log(mean) space.
double PoissonLossPieceLog::get_larger_root(double equals){
  double optimal_mean = argmin_mean(); //min or max!
  double optimal_cost = PoissonLoss(optimal_mean);
  double right_cost = getCost(max_log_mean);
  if(
     (optimal_cost < right_cost && right_cost < equals) ||
     (optimal_cost > right_cost && right_cost > equals)
     ){
    // intersection to the right of this interval, so just return some
    // value on the right side.
    return max_log_mean+1;
  }
  // Approximate the solution by the line through
  // (optimal_mean,optimal_cost) with the asymptotic slope. As m tends
  // to Inf, f'(m)=Linear+Log/m tends to Linear.
  //double candidate_root = optimal_mean + (equals-optimal_cost)/Linear;
  double candidate_root = optimal_mean + 1;
  // find the larger root of f(m) = Linear*m + Log*log(m) + Constant -
  // equals = 0.
  double candidate_cost, possibly_outside, deriv;
  double closest_positive_cost = INFINITY, closest_positive_mean=INFINITY;
  double closest_negative_cost = -INFINITY, closest_negative_mean=INFINITY;
  if(optimal_cost < 0){
    closest_negative_cost = optimal_cost;
    closest_negative_mean = optimal_mean;
  }else{
    closest_positive_cost = optimal_cost;
    closest_positive_mean = optimal_mean;
  }
  int step=0;
  do{
     candidate_cost = PoissonLoss(candidate_root) - equals;
     if(0 < candidate_cost && candidate_cost < closest_positive_cost){
       closest_positive_cost = candidate_cost;
       closest_positive_mean = candidate_root;
     }
     if(closest_negative_cost < candidate_cost && candidate_cost < 0){
       closest_negative_cost = candidate_cost;
       closest_negative_mean = candidate_root;
     }
     if(NEWTON_STEPS <= ++step){
       // Rprintf("larger root MAXSTEPS with equals=%e\n", equals);
       // print();
       // Rprintf("step=%d mean=%e cost=%e\n", step, candidate_root, candidate_cost);
       double between_closest = (closest_positive_mean + closest_negative_mean)/2;
       double between_cost = PoissonLoss(between_closest) - equals;
       if(ABS(between_cost) < ABS(candidate_cost)){
	 return log(between_closest);
       }else{
	 return log(candidate_root);
       }
     }
     deriv = PoissonDeriv(candidate_root);
     possibly_outside = candidate_root - candidate_cost/deriv;
     candidate_root = possibly_outside;
  }while(NEWTON_EPSILON < ABS(candidate_cost));
  //Rprintf("found root %e in %d steps!\n", candidate_root, step);
  return log(candidate_root);
}

double PoissonLossPieceLog::get_smaller_root(double equals){
  // root finding using the fact that g'(x)=Linear*e^x+Log is linear
  // as x tends to Inf.
  // find the smaller root of g(x) = Linear*e^x + Log*x + Constant -
  // equals = 0. x = log_mean
  double optimal_log_mean = argmin(); //min or max!
  double optimal_cost = getCost(optimal_log_mean);
  double left_cost = getCost(min_log_mean);
  if(
     (equals < left_cost && left_cost < optimal_cost) ||
     (equals > left_cost && left_cost > optimal_cost)
     ){
    // intersection to the left of this interval, so just return some
    // value on the left side-- it will be ignored later.
    return min_log_mean-1;
  }
  double candidate_root = optimal_log_mean - 1;
  double candidate_cost, possibly_outside, deriv;
  // as we search we will store bounds on the left and the right of
  // the zero point.
  double closest_positive_cost = INFINITY, closest_positive_log_mean=INFINITY;
  double closest_negative_cost = -INFINITY, closest_negative_log_mean=INFINITY;
  if(optimal_cost < 0){
    closest_negative_cost = optimal_cost;
    closest_negative_log_mean = optimal_log_mean;
  }else{
    closest_positive_cost = optimal_cost;
    closest_positive_log_mean = optimal_log_mean;
  }
  int step=0;
  double offset;
  do{
     candidate_cost = getCost(candidate_root) - equals;
     if(0 < candidate_cost && candidate_cost < closest_positive_cost){
       closest_positive_cost = candidate_cost;
       closest_positive_log_mean = candidate_root;
     }
     if(closest_negative_cost < candidate_cost && candidate_cost < 0){
       closest_negative_cost = candidate_cost;
       closest_negative_log_mean = candidate_root;
     }
     if(NEWTON_STEPS <= ++step){
       // Rprintf("smaller root MAXSTEPS equals=%e\n", equals);
       // print();
       // Rprintf("step=%d log_mean=%e cost=%e\n", step, candidate_root, candidate_cost);
       // Rprintf("neg_cost=%e neg_log_mean=%e pos_cost=%e pos_log_mean=%e\n", closest_negative_cost, closest_negative_
       //log_mean, closest_positive_cost, closest_positive_log_mean);
       double between_closest = (closest_positive_log_mean + closest_negative_log_mean)/2;
       double between_cost = getCost(between_closest) - equals;
       if(ABS(between_cost) < ABS(candidate_cost)){
	 return between_closest;
       }else{
	 return candidate_root;
       }
     }
     deriv = getDeriv(candidate_root);
     offset = candidate_cost/deriv;
     possibly_outside = candidate_root - offset;
     candidate_root = possibly_outside;
  }while(NEWTON_EPSILON < ABS(candidate_cost));
  return candidate_root;
}

double PoissonLossPieceLog::argmin_mean(){
  // f(m) = Linear*m + Log*log(m) + Constant,
  // f'(m)= Linear + Log/m = 0 means
  // m = -Log/Linear.
  return - Log / Linear;
}

double PoissonLossPieceLog::argmin(){
  // g(x) = Linear*e^x + Log*x + Constant,
  // g'(x)= Linear*e^x + Log = 0 means
  // x = log(-Log/Linear).
  return log(argmin_mean());
}

double PoissonLossPieceLog::getCost(double log_mean){
  // f(m) = Linear*m + Log*log(m) + Constant,
  // x = log(m),
  // g(x) = Linear*e^x + Log*x + Constant.
  double linear_term, log_term;
  if(log_mean == -INFINITY){
    linear_term = 0.0;
  }else{
    linear_term = Linear*exp(log_mean);
  }
  if(Log==0){
    log_term = 0.0;
  }else{
    log_term = Log*log_mean;
  }
  return linear_term + log_term + Constant;
}

double PoissonLossPieceLog::getDeriv(double log_mean){
  // g(x) = Linear*e^x + Log*x + Constant,
  // g'(x)= Linear*e^x + Log.
  double linear_term;
  if(log_mean == -INFINITY){
    linear_term = 0.0;
  }else{
    linear_term = Linear*exp(log_mean);
  }
  return linear_term + Log;
}

void PiecewisePoissonLossLog::set_to_min_less_of
(PiecewisePoissonLossLog *input, int verbose){
  piece_list.clear();
  PoissonLossPieceListLog::iterator it = input->piece_list.begin();
  PoissonLossPieceListLog::iterator next_it;
  double prev_min_cost = INFINITY;
  double prev_min_log_mean = it->min_log_mean;
  double prev_best_log_mean=INFINITY;
  while(it != input->piece_list.end()){
    double left_cost = it->getCost(it->min_log_mean);
    double right_cost = it->getCost(it->max_log_mean);
    if(prev_min_cost == INFINITY){
      // Look for min achieved in this interval.
      if(verbose){
	Rprintf("Searching for min in\n");
	it->print();
      }
      double next_left_cost=INFINITY;
      next_it = it;
      next_it++;
      if(it->Log==0){
	// degenerate linear function. since the Linear coef is never
	// negative, we know that this function must be increasing or
	// numerically constant on this interval.
	// g(x) = Linear*e^x + Constant,
	if(verbose)Rprintf("DEGENERATE LINEAR FUNCTION IN MIN LESS\n");
	// We used to check if(it->Linear==0) but there are some cases
	// when the function has a non-zero Linear coefficient, but is
	// numerically constant (e.g. Linear=156 between -inf and
	// -44). So now we check to see if the cost on the left and
	// right of the interval are equal.
	double right_left_diff = right_cost - left_cost;
	if(verbose)Rprintf("right_cost-left_cost=%e\n", right_left_diff);
	bool right_left_equal = right_left_diff < NEWTON_EPSILON;
	bool next_cost_more_than_left;
	if(next_it == input->piece_list.end()){
	  next_cost_more_than_left = true;
	}else{
	  next_left_cost = next_it->getCost(next_it->min_log_mean);
	  double next_left_diff = next_left_cost-left_cost;
	  if(verbose)Rprintf("next_left_cost-left_cost=%e\n", next_left_diff);
	  next_cost_more_than_left = NEWTON_EPSILON < next_left_diff;
	}
	// next_cost_more_than_left is true if the cost on the left of
	// the next piece is numerically greater than the cost on the
	// left of this interval => this is a NOT sufficient condition
	// for finding a minimum on the left and starting a constant
	// interval. We also need to make sure the left and right
	// sides are not equal -- there are some cases where this
	// interval is numerically constant, and the next interval's
	// left side is slightly larger in cost, but we can't find a
	// root in that next interval.
	if(next_cost_more_than_left && !right_left_equal){
	  // don't store this interval, but store its min cost as a
	  // constant.
	  prev_min_cost = left_cost;
	  prev_best_log_mean = it->min_log_mean;
	  if(verbose){
	    Rprintf("Increasing interval left_cost=%e(stored) right_cost=%e diff=%e\n", left_cost, right_cost, right_cost-left_cost);
	    it->print();
	  }
	}else{
	  if(verbose){
	    Rprintf("Numerically constant convex piece\n");
	    it->print();
	  }
	  // store this numerically constant interval.
	  piece_list.emplace_back
	    (it->Linear, it->Log, it->Constant,
	     prev_min_log_mean, it->max_log_mean,
	     PREV_NOT_SET, INFINITY); // equality constraint active on convex piece.
	  prev_min_log_mean = it->max_log_mean;
	}
      }else{//not degenerate linear
	double mu = it->argmin();
	double mu_cost = it->getCost(mu);
	bool next_ok;
	if(next_it == input->piece_list.end()){
	  next_ok = true;
	}else{
	  next_left_cost = next_it->getCost(next_it->min_log_mean);
	  next_ok = NEWTON_EPSILON < next_left_cost-mu_cost;
	}
	// Compute the cost at the next interval (interval to the
	// left), to check if the cost at the minimum is less than the
	// cost on the edge of the next function piece. This is
	// necessary because sometimes there are numerical issues.
 	if(verbose){
	  Rprintf("min cost=%f at log_mean=%f\n", mu_cost, mu);
	  Rprintf("next-mu=%e right-mu=%e\n", next_left_cost-mu, right_cost-mu);
	}
	bool cost_ok = NEWTON_EPSILON < right_cost-mu_cost && next_ok;
	if(mu <= it->min_log_mean && cost_ok){
	  /* The minimum is achieved on the left or before this
	     interval, so this function is always increasing in this
	     interval. We don't need to store it, but we do need to keep
	     track of the minimum cost, which occurs at the min mean
	     value in this interval. */
	  if(verbose)Rprintf("min before interval\n");
	  prev_min_cost = it->getCost(it->min_log_mean);
	  prev_best_log_mean = it->min_log_mean;
	}else if(mu < it->max_log_mean && cost_ok){
	  // Minimum in this interval, so add a convex piece up to the
	  // min, and keep track of the min cost to create a constant
	  // piece later. NB it is possible that prev_min_log_mean==mu in
	  // which case we do not need to store the convex piece.
	  if(verbose){
	    Rprintf("min in this interval at log_mean=%f cost=%f\n", mu, mu_cost);
	    Rprintf("right_cost=%f right-constant=%e\n", right_cost, right_cost-mu_cost);
	    Rprintf("next_left_cost=%f next-constant=%e\n", next_left_cost, next_left_cost-mu_cost);
	  }
	  if(prev_min_log_mean < mu){
	    piece_list.emplace_back
	      (it->Linear, it->Log, it->Constant, prev_min_log_mean, mu,
	       PREV_NOT_SET, INFINITY); // equality constraint active on convex piece.
	  }
	  prev_min_log_mean = mu;
	  prev_best_log_mean = mu;
	  prev_min_cost = mu_cost;
	  if(verbose)Rprintf("prev_min_cost=%f\n", prev_min_cost);
	}else{
	  // Minimum after this interval, so this function is
	  // decreasing on this entire interval, and so we can just
	  // store it as is.
	  if(verbose)Rprintf("min after interval\n");
	  piece_list.emplace_back
	    (it->Linear, it->Log, it->Constant, prev_min_log_mean, it->max_log_mean,
	     PREV_NOT_SET, INFINITY); // equality constraint active on convex piece.
	  prev_min_log_mean = it->max_log_mean;
	}//if(non-degenerate mu in interval
      }//if(degenerate linear cost.
    }else{//prev_min_cost is finite
      // Look for a function with prev_min_cost in its interval.
      if(verbose){
	Rprintf("Searching for intersection with %f\n", prev_min_cost);
	Rprintf("cost at limits=[%f,%f] cost-constant=[%e,%e]\n",
	       left_cost, right_cost,
	       left_cost-prev_min_cost, right_cost-prev_min_cost);
	it->print();
      }
      if(it->Log==0){
	//degenerate Linear case
	if(it->Linear < 0){
	  //decreasing linear function.
	  throw(500);//this should never happen.
	}else{
	  //increasing linear function, so will not intersect the
	  //constant below.
	}
      }else{// Log is not zero.
	// Theoretically there can be zero, one, or two intersection
	// points between the constant function prev_min_log_mean and a
	// non-degenerate Poisson loss function piece. 
	if(it->has_two_roots(prev_min_cost)){
	  // There are two mean values where the Poisson loss piece
	  // intersects the constant function prev_min_log_mean, but we
	  // are only concerned with the first mean value (the
	  // lesser of the two).
	  double mu = it->get_smaller_root(prev_min_cost);
	  double cost_mu = it->getCost(mu);
	  if(verbose)Rprintf("smaller_root log_mean=%f root_cost=%f root_cost-side_cost=[%f,%f]\n", mu, cost_mu, cost_mu-left_cost, cost_mu-right_cost);
	  if(it->min_log_mean < mu && mu < it->max_log_mean){
	    // The smaller intersection point occurs within the
	    // interval, so the constant interval ends here, and we
	    // can store it immediately.
	    piece_list.emplace_back
	      (0, 0, prev_min_cost,
	       prev_min_log_mean, mu, PREV_NOT_SET,
	       prev_best_log_mean);// equality constraint inactive.
	    prev_min_cost = INFINITY;
	    prev_min_log_mean = mu;
	    it--;
	  }//if(mu in interval)
	}//if(has two roots
	if(right_cost <= prev_min_cost+NEWTON_EPSILON && prev_min_cost < INFINITY){
	  //ends exactly/numerically on the right.
	  if(verbose)Rprintf("constant numerically equal on right\n");
	  piece_list.emplace_back
	    (0, 0, prev_min_cost,
	     prev_min_log_mean, it->max_log_mean, 
	     PREV_NOT_SET,
	     prev_best_log_mean);
	  prev_min_cost = INFINITY;
	  prev_min_log_mean = it->max_log_mean;
	}
      }//if(Log is zero
    }//if(prev_min_cost is finite
    it++;
    if(verbose){
      Rprintf("current min-less-------------------\n");
      print();
    }
  }//while(it
  if(prev_min_cost < INFINITY){
    // ending on a constant piece.
    it--;
    piece_list.emplace_back
      (0, 0, prev_min_cost,
       prev_min_log_mean, it->max_log_mean, PREV_NOT_SET,
       prev_best_log_mean);//equality constraint inactive on constant piece.
  }
}

void PiecewisePoissonLossLog::set_to_min_more_of
(PiecewisePoissonLossLog *input, int verbose){
  piece_list.clear();
  PoissonLossPieceListLog::iterator it = input->piece_list.end();
  PoissonLossPieceListLog::iterator prev_it;
  it--;
  double prev_min_cost = INFINITY;
  double prev_max_log_mean = it->max_log_mean;
  double prev_best_log_mean=INFINITY;
  it++;
  if(verbose)print();
  while(it != input->piece_list.begin()){
    it--;
    if(prev_min_cost == INFINITY){
      // Look for min achieved in this interval.
      if(verbose){
	Rprintf("Searching for min in\n");
	it->print();
      }
      if(it->Log==0){
	if(verbose)Rprintf("DEGENERATE LINEAR FUNCTION IN MIN MORE\n");
	//degenerate Linear function. since the Linear coef is never
	//negative, we know that this function must be increasing or
	//numerically constant on this interval. In both cases we
	//should just store this interval.
	piece_list.emplace_front
	  (it->Linear, it->Log, it->Constant, it->min_log_mean, prev_max_log_mean,
	   PREV_NOT_SET, INFINITY); // equality constraint active on convex piece.
	prev_max_log_mean = it->min_log_mean;
      }else{
	double mu = it->argmin();
	double mu_cost = it->getCost(mu);
	// Compute the cost at the prev interval (interval to the
	// left), to check if the minimum is really a minimum. This is
	// necessary because sometimes there are numerical issues.
	bool prev_ok;
	if(it == input->piece_list.begin()){
	  prev_ok = true;
	}else{
	  prev_it = it;
	  prev_it--;
	  double prev_cost_right = prev_it->getCost(prev_it->max_log_mean);
	  prev_ok = NEWTON_EPSILON < prev_cost_right-mu_cost;
	}
	double this_cost_left = it->getCost(it->min_log_mean);
	if(it->max_log_mean <= mu){
	  double this_cost_right = it->getCost(it->max_log_mean);
	  // the cost diff below is for the case when the min is to the
	  // right of the interval, in which case the coefficients say
	  // that the function is decreasing on this entire
	  // interval. however we need to test if it is numerically
	  // decreasing or constant, so we compare the cost on the left
	  // of the interval to the cost on the right. this_cost_diff
	  // should always be positive, but if it is numerically very
	  // close to zero we treat the interval as constant.
	  double this_cost_diff = this_cost_left-this_cost_right;
	  if(verbose){
	    Rprintf("min after this interval\n");
	    Rprintf("cost_left_right=[%f,%f] diff=%f\n",
		    this_cost_left, this_cost_right, this_cost_diff);
	  }
	  if(NEWTON_EPSILON < this_cost_diff){
	  /* The minimum is achieved after this interval, so this
	     function is always decreasing in this interval. This
	     happens when the minimum occurs exactly on the right of
	     this interval (the previous step -- the next interval --
	     must have concluded that there was a min before the
	     interval). We don't need to store this interval, because
	     we start creating a constant piece here. */
	    if(verbose)Rprintf("decreasing interval, starting constant piece from right.\n");
	    prev_min_cost = this_cost_right;
	    prev_best_log_mean = it->max_log_mean;
	  }else{
	  /* NOTE: if we test only it->max_log_mean <= mu, then that
	     is bad for intervals that are numerically constant. For
	     example sometimes we have three intervals: increasing,
	     decreasing, increasing according to the coefficients, but
	     we should test the cost values at the edges and treat it
	     as increasing, constant, increasing. */
	    if(verbose)Rprintf("constant interval, storing numerically constant convex piece.\n");
	    piece_list.emplace_front
	      (it->Linear, it->Log, it->Constant, it->min_log_mean, prev_max_log_mean,
	       PREV_NOT_SET, INFINITY); // equality constraint active on convex piece.
	    prev_max_log_mean = it->min_log_mean;
	  }
	}else if(it->min_log_mean < mu && NEWTON_EPSILON < this_cost_left-mu_cost && prev_ok){
	  // Minimum in this interval, so add a convex piece up to the
	  // min, and keep track of the min cost to create a constant
	  // piece later. NB it is possible that mu==prev_max_log_mean, in
	  // which case we do not need to save the convex piece.
	  if(verbose)Rprintf("min in this interval at mu=%f\n", mu);
	  if(mu < prev_max_log_mean){ 
	    piece_list.emplace_front
	      (it->Linear, it->Log, it->Constant, mu, prev_max_log_mean, 
	       PREV_NOT_SET, INFINITY); // equality constraint active on convex piece.
	  }
	  prev_max_log_mean = mu;
	  prev_best_log_mean = mu;
	  prev_min_cost = mu_cost;
	}else{
	  // Minimum before this interval, so this function is
	  // increasing on this entire interval, and so we can just
	  // store it as is.
	  if(verbose)Rprintf("min before this interval\n");
	  piece_list.emplace_front
	    (it->Linear, it->Log, it->Constant, it->min_log_mean, prev_max_log_mean,
	     PREV_NOT_SET, INFINITY); // equality constraint active on convex piece.
	  prev_max_log_mean = it->min_log_mean;
	}
      }//if(degenerate linear)else
    }else{//prev_min_cost is finite
      // Look for a function with prev_min_cost in its interval.
      double left_cost = it->getCost(it->min_log_mean);
      double right_cost = it->getCost(it->max_log_mean);
      if(verbose){
	Rprintf("Searching for intersection with %f\n", prev_min_cost);
	Rprintf("cost at limits=[%f,%f] cost-constant=[%e,%e]\n",
	       left_cost, right_cost,
	       left_cost-prev_min_cost, right_cost-prev_min_cost);
	it->print();
      }
      double mu = INFINITY;
      if(it->Log==0){
	//degenerate Linear case, there is one intersection point.
	mu = log((prev_min_cost - it->Constant)/it->Linear);
	if(verbose)Rprintf("degenerate linear intersection at log_mean=%f\n", mu);
      }else{// Log is not zero.
	// Theoretically there can be zero, one, or two intersection
	// points between the constant function prev_min_log_mean and a
	// non-degenerate Poisson loss function piece. 
	if(it->has_two_roots(prev_min_cost)){
	  // There are two mean values where the Poisson loss piece
	  // intersects the constant function prev_min_log_mean, but we
	  // are only concerned with the second mean value (the
	  // greater of the two). 
	  mu = it->get_larger_root(prev_min_cost);
	  if(verbose)Rprintf("large root log_mean=%f\n", mu);
	}//if(there are two roots
      }//if(Log is zero
      if(it->min_log_mean < mu && mu < it->max_log_mean){
	// The intersection point occurs within the interval, so the
	// constant interval ends here, and we can store it
	// immediately.
	if(verbose)Rprintf("%f in interval\n", mu);
	piece_list.emplace_front
	  (0, 0, prev_min_cost,
	   mu, prev_max_log_mean,
	   PREV_NOT_SET,
	   prev_best_log_mean);// equality constraint inactive on constant piece.
	prev_min_cost = INFINITY;
	prev_max_log_mean = mu;
	it++;
      }else if(left_cost <= prev_min_cost+NEWTON_EPSILON){
	//ends exactly/numerically on the left.
	if(verbose)Rprintf("constant numerically equal on left\n");
	piece_list.emplace_front
	  (0, 0, prev_min_cost,
	   it->min_log_mean, prev_max_log_mean,
	   PREV_NOT_SET,
	   prev_best_log_mean);// equality constraint inactive on constant piece.
	prev_min_cost = INFINITY;
	prev_max_log_mean = it->min_log_mean;
      }
    }//if(prev_min_cost is finite
    if(verbose){
      Rprintf("current min-more-------------------\n");
      print();
    }
  }//while(it
  if(prev_min_cost < INFINITY){
    // ending on a constant piece.
    piece_list.emplace_front
      (0, 0, prev_min_cost,
       it->min_log_mean, prev_max_log_mean,
       PREV_NOT_SET,
       prev_best_log_mean);//equality constraint inactive on constant piece.
  }
}

void PiecewisePoissonLossLog::add(double Linear, double Log, double Constant){
  PoissonLossPieceListLog::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->Linear += Linear;
    it->Log += Log;
    it->Constant += Constant;
  }
}

void PiecewisePoissonLossLog::multiply(double x){
  PoissonLossPieceListLog::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->Linear *= x;
    it->Log *= x;
    it->Constant *= x;
  }
}

void PiecewisePoissonLossLog::set_prev_seg_end(int prev_seg_end){
  PoissonLossPieceListLog::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->data_i = prev_seg_end;
  }
}

void PiecewisePoissonLossLog::findMean
(double log_mean, int *seg_end, double *prev_log_mean){
  PoissonLossPieceListLog::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    if(it->min_log_mean <= log_mean && log_mean <= it->max_log_mean){
      *seg_end = it->data_i;
      *prev_log_mean = it->prev_log_mean;
      return;
    }
  }
}

double PiecewisePoissonLossLog::findCost(double log_mean){
  PoissonLossPieceListLog::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    if(it->min_log_mean <= log_mean && log_mean <= it->max_log_mean){
      return it->getCost(log_mean);
    }
  }
  // we should never get here, but we include this return in order to
  // avoid clang++ warning control may reach end of non-void function.
  return INFINITY;
}

void PiecewisePoissonLossLog::print(){
  PoissonLossPieceListLog::iterator it;
  Rprintf("%10s %10s %15s %15s %15s %15s %s\n",
	 "Linear", "Log", "Constant",
	 "min_log_mean", "max_log_mean",
	 "prev_log_mean", "data_i");
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->print();
  }
}

void PoissonLossPieceLog::print(){
  Rprintf("%.20e %.20e %.20e %15f %15f %15f %d\n",
	 Linear, Log, Constant,
	 min_log_mean, max_log_mean,
	 prev_log_mean, data_i);
}

bool PoissonLossPieceLog::equality_constraint_active(){
  return prev_log_mean == INFINITY;
}

void PiecewisePoissonLossLog::Minimize
(double *best_cost,
 double *best_log_mean,
 int *data_i,
 double *prev_log_mean){
  double candidate_cost, candidate_log_mean;
  PoissonLossPieceListLog::iterator it;
  *best_cost = INFINITY;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    candidate_log_mean = it->argmin();
    if(candidate_log_mean < it->min_log_mean){
      candidate_log_mean = it->min_log_mean;
    }else if(it->max_log_mean < candidate_log_mean){
      candidate_log_mean = it->max_log_mean;
    }
    candidate_cost = it->getCost(candidate_log_mean);
    if(candidate_cost < *best_cost){
      *best_cost = candidate_cost;
      *best_log_mean = candidate_log_mean;
      *data_i = it->data_i;
      *prev_log_mean = it->prev_log_mean;
    }
  }
}

// check that this function is the minimum on all pieces.
int PiecewisePoissonLossLog::check_min_of
(PiecewisePoissonLossLog *prev, PiecewisePoissonLossLog *model){
  PoissonLossPieceListLog::iterator it;
  return 0;// disables checks!
  for(it = piece_list.begin(); it != piece_list.end(); it++){
    if(it != piece_list.begin()){
      PoissonLossPieceListLog::iterator pit = it;
      pit--;
      if(pit->max_log_mean != it->min_log_mean){
	Rprintf("prev->max_log_mean != it->min_log_mean min\n");
	return 3;
      }
      // double cost_prev = pit->getCost(pit->max_log_mean);
      // double cost_here = it->getCost(it->min_log_mean);
      // if(0.1 < ABS(cost_prev - cost_here)){
      // 	Rprintf("discontinuity detected at %f, %f != %f\n", pit->max_log_mean, cost_prev, cost_here);
      // 	pit->print();
      // 	it->print();
      // 	return 4;
      // }
    }
    if(it->max_log_mean <= it->min_log_mean){
      Rprintf("max_log_mean<=min_log_mean=%15.10f min\n", it->min_log_mean);
      return 2;
    }
    double mid_log_mean = (it->min_log_mean + it->max_log_mean)/2;
    if(-INFINITY == mid_log_mean){
      mid_log_mean = it->max_log_mean - 1;
    }
    if(INFINITY == mid_log_mean){
      mid_log_mean = it->min_log_mean + 1;
    }
    double cost_min = it->getCost(mid_log_mean);
    double cost_prev = prev->findCost(mid_log_mean);
    if(cost_prev+1e-6 < cost_min){
      Rprintf("prev(%f)=%f\n", mid_log_mean, cost_prev);
      prev->print();
      Rprintf("min(%f)=%f\n", mid_log_mean, cost_min);
      print();
      return 1;
    }
    double cost_model = model->findCost(mid_log_mean);
    double cost_model_diff = cost_min - cost_model;
    if(NEWTON_EPSILON < cost_model_diff){
      Rprintf("min-model=%e\n", cost_model_diff);
      Rprintf("model(%f)=%f\n", mid_log_mean, cost_model);
      model->print();
      Rprintf("min(%f)=%f\n", mid_log_mean, cost_min);
      print();
      return 1;
    }
  }
  for(it = prev->piece_list.begin(); it != prev->piece_list.end(); it++){
    if(it != prev->piece_list.begin()){
      PoissonLossPieceListLog::iterator pit = it;
      pit--;
      if(pit->max_log_mean != it->min_log_mean){
	Rprintf("prev->max_log_mean != it->min_log_mean prev\n");
	return 3;
      }
    }
    if(it->max_log_mean <= it->min_log_mean){
      Rprintf("max_log_mean<=min_log_mean=%15.10f prev\n", it->min_log_mean);
      return 2;
    }
    double mid_log_mean = (it->min_log_mean + it->max_log_mean)/2;
    if(-INFINITY == mid_log_mean){
      mid_log_mean = it->max_log_mean - 1;
    }
    if(INFINITY == mid_log_mean){
      mid_log_mean = it->min_log_mean + 1;
    }
    double cost_prev = it->getCost(mid_log_mean);
    double cost_min = findCost(mid_log_mean);
    if(cost_prev+1e-6 < cost_min){
      Rprintf("prev(%f)=%f\n", mid_log_mean, cost_prev);
      prev->print();
      Rprintf("min(%f)=%f\n", mid_log_mean, cost_min);
      print();
      return 1;
    }
  }
  for(it = model->piece_list.begin(); it != model->piece_list.end(); it++){
    if(it != model->piece_list.begin()){
      PoissonLossPieceListLog::iterator pit = it;
      pit--;
      if(pit->max_log_mean != it->min_log_mean){
	Rprintf("prev->max_log_mean != it->min_log_mean model\n");
	return 3;
      }
    }
    if(it->max_log_mean <= it->min_log_mean){
      Rprintf("max_log_mean<=min_log_mean=%15.10f model\n", it->min_log_mean);
      return 2;
    }
    double mid_log_mean = (it->min_log_mean + it->max_log_mean)/2;
    if(-INFINITY == mid_log_mean){
      mid_log_mean = it->max_log_mean - 1;
    }
    if(INFINITY == mid_log_mean){
      mid_log_mean = it->min_log_mean + 1;
    }
    double cost_model = it->getCost(mid_log_mean);
    double cost_min = findCost(mid_log_mean);
    double cost_model_diff = cost_min - cost_model;
    if(NEWTON_EPSILON < cost_model_diff){
      Rprintf("min-model=%e\n", cost_model_diff);
      Rprintf("model(%f)=%f\n", mid_log_mean, cost_model);
      model->print();
      Rprintf("min(%f)=%f\n", mid_log_mean, cost_min);
      print();
      return 1;
    }
  }
  return 0;
}

void PiecewisePoissonLossLog::set_to_min_env_of
(PiecewisePoissonLossLog *fun1, PiecewisePoissonLossLog *fun2, int verbose){
  PoissonLossPieceListLog::iterator
    it1 = fun1->piece_list.begin(),
    it2 = fun2->piece_list.begin();
  if(verbose){
    Rprintf("computing min env of:\n");
    Rprintf("=min-less/more\n");
    fun1->print();
    Rprintf("=cost model\n");
    fun2->print();
  }
  piece_list.clear();
  while(it1 != fun1->piece_list.end() &&
	it2 != fun2->piece_list.end()){
    push_min_pieces(fun1, fun2, it1, it2, verbose);
    if(verbose){
      print();
      Rprintf("------\n");
    }
    double last_max_log_mean = piece_list.back().max_log_mean;
    if(it1->max_log_mean == last_max_log_mean){
      it1++;
    }
    if(it2->max_log_mean == last_max_log_mean){
      it2++;
    }
  }
}

bool sameFuns
(PoissonLossPieceListLog::iterator it1,
 PoissonLossPieceListLog::iterator it2){
  return it1->Linear == it2->Linear &&
    it1->Log == it2->Log &&
    ABS(it1->Constant - it2->Constant) < NEWTON_EPSILON;
}

void PiecewisePoissonLossLog::push_min_pieces
(PiecewisePoissonLossLog *fun1,
 PiecewisePoissonLossLog *fun2,
 PoissonLossPieceListLog::iterator it1,
 PoissonLossPieceListLog::iterator it2,
 int verbose){
  bool same_at_left;
  double last_min_log_mean;
  PoissonLossPieceListLog::iterator prev2 = it2;
  prev2--;
  PoissonLossPieceListLog::iterator prev1 = it1;
  prev1--;
  if(it1->min_log_mean < it2->min_log_mean){
    //it1 function piece starts to the left of it2.
    same_at_left = sameFuns(prev2, it1);
    last_min_log_mean = it2->min_log_mean;
  }else{
    //it1 function piece DOES NOT start to the left of it2.
    last_min_log_mean = it1->min_log_mean;
    if(it2->min_log_mean < it1->min_log_mean){
      //it2 function piece starts to the left of it1.
      same_at_left = sameFuns(prev1, it2);
    }else{
      //it1 and it2 start at the same min_log_mean value.
      if(it1==fun1->piece_list.begin() &&
	 it2==fun2->piece_list.begin()){
	same_at_left = false;
      }else{
	same_at_left = sameFuns(prev1, prev2);
      }
    }
  }
  PoissonLossPieceListLog::iterator next2 = it2;
  next2++;
  PoissonLossPieceListLog::iterator next1 = it1;
  next1++;
  bool same_at_right;
  double first_max_log_mean;
  if(it1->max_log_mean < it2->max_log_mean){
    if(verbose)Rprintf("it2 function piece continues to the right of it1.\n");
    same_at_right = sameFuns(next1, it2);
    first_max_log_mean = it1->max_log_mean;
  }else{
    first_max_log_mean = it2->max_log_mean;
    if(it2->max_log_mean < it1->max_log_mean){
      if(verbose)Rprintf("it2 function piece ends before it1.\n");
      same_at_right = sameFuns(it1, next2);
    }else{
      if(verbose)Rprintf("it2 and it1 end at same max_log_mean.\n");
      if(next1==fun1->piece_list.end() &&
	 next2==fun2->piece_list.end()){
	if(verbose)Rprintf("at the end so next can't be the same.\n");
	same_at_right = false;
      }else{
	if(verbose){
	  Rprintf("comparing next function pieces.\n");
	  next1->print();
	  next2->print();
	}
	same_at_right = sameFuns(next1, next2);
      }
    }
  }
  if(last_min_log_mean == first_max_log_mean){
    // we should probably never get here, but if we do, no need to
    // store this interval.
    if(verbose){
      Rprintf("prev\n");
      fun1->print();
      Rprintf("model\n");
      fun2->print();
      Rprintf("interval size 0!-----------------\n");
    }
    return;
  }
  if(sameFuns(it1, it2)){
    // The functions are exactly equal over the entire interval so we
    // can push either of them.
    push_piece(it1, last_min_log_mean, first_max_log_mean);
    if(verbose)Rprintf("exactly equal over entire interval\n");
    return;
  }
  PoissonLossPieceLog diff_piece
    (it1->Linear - it2->Linear,
     it1->Log - it2->Log,
     it1->Constant - it2->Constant,
     last_min_log_mean, first_max_log_mean,
     -5, false);
  // Evaluate the middle in the original space, to avoid problems when
  // first_max_log_mean is -Inf.
  double mid_mean = (exp(first_max_log_mean) + exp(last_min_log_mean))/2;
  double cost_diff_mid = diff_piece.getCost(log(mid_mean));
  // Easy case of equality on both left and right.
  if(same_at_left && same_at_right){
    if(verbose)Rprintf("Same on both the left and the right\n");
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_log_mean, first_max_log_mean);
    }else{
      push_piece(it2, last_min_log_mean, first_max_log_mean);
    }
    return;
  }
  // Easy degenerate cases that do not require root finding.
  if(diff_piece.Log == 0){
    // g(x) = Linear*e^x + Constant = 0,
    // x = log(-Constant/Linear).
    if(diff_piece.Linear == 0){
      // They are offset by a Constant.
      if(diff_piece.Constant < 0){
	push_piece(it1, last_min_log_mean, first_max_log_mean);
      }else{
	push_piece(it2, last_min_log_mean, first_max_log_mean);
      }
      if(verbose)Rprintf("offset by a constant=%e\n", diff_piece.Constant);
      return;
    }
    if(diff_piece.Constant == 0){
      // The only difference is the Linear coef.
      if(diff_piece.Linear < 0){
	push_piece(it1, last_min_log_mean, first_max_log_mean);
      }else{
	push_piece(it2, last_min_log_mean, first_max_log_mean);
      }
      if(verbose)Rprintf("only diff is linear coef\n");
      return;
    }
    double log_mean_at_equal_cost = log(-diff_piece.Constant/diff_piece.Linear);
    if(last_min_log_mean < log_mean_at_equal_cost &&
       log_mean_at_equal_cost < first_max_log_mean){
      // the root is in the interval, so we need to add two intervals.
      if(0 < diff_piece.Linear){
	push_piece(it1, last_min_log_mean, log_mean_at_equal_cost);
	push_piece(it2, log_mean_at_equal_cost, first_max_log_mean);
      }else{
	push_piece(it2, last_min_log_mean, log_mean_at_equal_cost);
	push_piece(it1, log_mean_at_equal_cost, first_max_log_mean);
      }
      if(verbose)Rprintf("Log zero with one root in interval\n");
      return;
    }
    // the root is outside the interval, so one is completely above
    // the other over this entire interval.
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_log_mean, first_max_log_mean);
    }else{
      push_piece(it2, last_min_log_mean, first_max_log_mean);
    }
    if(verbose)Rprintf("Log zero with no roots in interval\n");
    return;
  }//if(diff->Log == 0
  double cost_diff_left = diff_piece.getCost(last_min_log_mean);
  double cost_diff_right = diff_piece.getCost(first_max_log_mean);
  bool two_roots = diff_piece.has_two_roots(0.0);
  double smaller_log_mean=INFINITY, larger_log_mean=INFINITY;
  if(two_roots){
    smaller_log_mean = diff_piece.get_smaller_root(0.0);
    larger_log_mean = diff_piece.get_larger_root(0.0);
    if(verbose)Rprintf("Computed crossing points: %e %e\n", smaller_log_mean, larger_log_mean);
  }
  if(same_at_right){
    // they are equal on the right, but we don't know if there is
    // another crossing point somewhere to the left.
    if(two_roots){
      // there could be a crossing point to the left.
      double log_mean_at_crossing = smaller_log_mean;
      double log_mean_between_zeros = (log_mean_at_crossing + first_max_log_mean)/2;
      double cost_between_zeros = diff_piece.getCost(log_mean_between_zeros);
      double log_mean_at_optimum = diff_piece.argmin();
      if(verbose){
	Rprintf("Determining if there is a crossing point in interval...\n", last_min_log_mean, cost_diff_left);
	diff_piece.print();
	Rprintf("cost_diff(left:%e)=%e\n", last_min_log_mean, cost_diff_left);
	Rprintf("cost_diff(cross:%e)=%e\n", log_mean_at_crossing, diff_piece.getCost(log_mean_at_crossing));
	Rprintf("cost_diff(between:%e)=%e\n", log_mean_between_zeros, cost_between_zeros);
	Rprintf("cost_diff(optimum:%e)=%e\n", log_mean_at_optimum, diff_piece.getCost(log_mean_at_optimum));
	Rprintf("cost_diff(right:%e)=%e\n", first_max_log_mean, cost_diff_right);
      }
      if(last_min_log_mean < log_mean_at_crossing &&
	 log_mean_at_crossing < log_mean_at_optimum && log_mean_at_optimum < first_max_log_mean){
	//the cross point is in the interval.
	if(cost_diff_left < 0){
	  push_piece(it1, last_min_log_mean, log_mean_at_crossing);
	  push_piece(it2, log_mean_at_crossing, first_max_log_mean);
	}else{
	  push_piece(it2, last_min_log_mean, log_mean_at_crossing);
	  push_piece(it1, log_mean_at_crossing, first_max_log_mean);
	}
	if(verbose)Rprintf("equal on the right with one crossing in interval\n");
	return;
      }
      // at log_mean=-Inf (mean=0) we can just test the Log
      // coefficient to determine which function is greater. If
      // cost(log_mean)=Linear*exp(log_mean)+Log*log_mean+Constant
      // then cost(-Inf)=Linear*0+Log*-Inf+Constant so
      // cost(-Inf)=f1-f2<0 implies Log*-Inf<0 implies 0<Log implies
      // it1<it2.
      bool it1_smaller_at_mean0 = 0 < diff_piece.Log;
      if(log_mean_at_crossing < last_min_log_mean){
	if(verbose)Rprintf("equal on the right, cross before interval\n");
	if(it1_smaller_at_mean0){
	  push_piece(it2, last_min_log_mean, first_max_log_mean);
	}else{
	  push_piece(it1, last_min_log_mean, first_max_log_mean);
	}
      }else{
	if(verbose)Rprintf("equal on the right, no cross before interval\n");
	if(it1_smaller_at_mean0){
	  push_piece(it1, last_min_log_mean, first_max_log_mean);
	}else{
	  push_piece(it2, last_min_log_mean, first_max_log_mean);
	}
      }
      return;
    }//if(two_roots
    if(verbose){
      Rprintf("equal on the right, but no roots anywhere\n");
    }
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_log_mean, first_max_log_mean);
    }else{
      push_piece(it2, last_min_log_mean, first_max_log_mean);
    }
    return;
  }
  if(same_at_left){
    // equal on the left.
    if(two_roots){
      // There could be a crossing point to the right.
      double log_mean_at_crossing = larger_log_mean;
      double log_mean_at_optimum = diff_piece.argmin();
      if(verbose)Rprintf("larger_log_mean=%f\n", log_mean_at_crossing);
      if(last_min_log_mean < log_mean_at_optimum &&
	 log_mean_at_optimum < log_mean_at_crossing &&
      	 log_mean_at_crossing < first_max_log_mean){
	// the crossing point is in this interval.
	if(cost_diff_right < 0){
	  push_piece(it2, last_min_log_mean, log_mean_at_crossing);
	  push_piece(it1, log_mean_at_crossing, first_max_log_mean);
	}else{
	  push_piece(it1, last_min_log_mean, log_mean_at_crossing);
	  push_piece(it2, log_mean_at_crossing, first_max_log_mean);
	}
	if(verbose)Rprintf("equal on the left with crossing in interval\n");
	return;
      }
    }//if(there may be crossing
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_log_mean, first_max_log_mean);
    }else{
      push_piece(it2, last_min_log_mean, first_max_log_mean);
    }
    if(verbose)Rprintf("equal on the left with no crossing in interval\n");
    return;
  }
  // The only remaining case is that the curves are equal neither on
  // the left nor on the right of the interval. However they may be
  // equal inside the interval, so let's check for that.
  double first_log_mean = INFINITY, second_log_mean = INFINITY;
  if(two_roots){
    bool larger_inside =
      last_min_log_mean < larger_log_mean && larger_log_mean < first_max_log_mean;
    if(verbose)Rprintf("smaller_log_mean=%f\nlarger_log_mean=%f\n",
		      smaller_log_mean,
		      larger_log_mean);
    bool smaller_inside =
      last_min_log_mean < smaller_log_mean &&
      0 < exp(smaller_log_mean) &&
      smaller_log_mean < first_max_log_mean;
    if(larger_inside){
      if(smaller_inside && smaller_log_mean < larger_log_mean){
	// both are in the interval.
	first_log_mean = smaller_log_mean;
	second_log_mean = larger_log_mean;
	if(verbose){
	  diff_piece.print();
	  Rprintf("%f and %f in [%f,%f]\n",
		 smaller_log_mean, larger_log_mean,
		 last_min_log_mean, first_max_log_mean);
	}
      }else{
	// smaller mean is not in the interval, but the larger is.
	first_log_mean = larger_log_mean;
	if(verbose){
	  Rprintf("%f in [%f,%f]\n",
		 first_log_mean,
		 last_min_log_mean, first_max_log_mean);
	}
      }
    }else{
      // larger mean is not in the interval
      if(smaller_inside){
	// smaller mean is in the interval, but not the larger.
	first_log_mean = smaller_log_mean;
	if(verbose){
	  Rprintf("%f in [%f,%f]\n",
		 first_log_mean,
		 last_min_log_mean, first_max_log_mean);
	}
      }
    }
  }//if(two_roots
  if(second_log_mean != INFINITY){
    // two crossing points.
    bool it1_larger_before;
    if(second_log_mean-first_log_mean < first_log_mean-last_min_log_mean){
      //bigger distance between start of last interval and first root,
      //so test which function is bigger at a point before the first
      //root.
      double before_mean = (exp(last_min_log_mean) + exp(first_log_mean))/2;
      double cost_diff_before = diff_piece.getCost(log(before_mean));
      if(cost_diff_before < 0){
	it1_larger_before = true;
      }else{
	it1_larger_before = false;
      }
    }else{
      // more space between the two roots, so test which function is
      // bigger at a midpoint between the two.
      double log_mean_between = (first_log_mean+second_log_mean)/2;
      double cost_diff_between = diff_piece.getCost(log_mean_between);
      if(cost_diff_between < 0){
	it1_larger_before = false;
      }else{
	it1_larger_before = true;
      }
    }
    if(it1_larger_before){
      push_piece(it1, last_min_log_mean, first_log_mean);
      push_piece(it2, first_log_mean, second_log_mean);
      push_piece(it1, second_log_mean, first_max_log_mean);
    }else{
      push_piece(it2, last_min_log_mean, first_log_mean);
      push_piece(it1, first_log_mean, second_log_mean);
      push_piece(it2, second_log_mean, first_max_log_mean);
    }
    if(verbose)Rprintf("not equal on the sides, 2 crossing points\n");
  }else if(first_log_mean != INFINITY){
    // "one" crossing point. actually sometimes we have last_min_log_mean
    // < first_log_mean < first_max_log_mean but cost_diff_before and
    // cost_diff_after have the same sign! In that case we need to
    // just push one piece.
    double before_mean = (exp(last_min_log_mean) + exp(first_log_mean))/2;
    double cost_diff_before = diff_piece.getCost(log(before_mean));
    if(verbose){
      Rprintf("cost_diff_before(%f)=%f\n", log(before_mean), cost_diff_before);
    }
    double after_mean = (first_max_log_mean + first_log_mean)/2;
    double cost_diff_after = diff_piece.getCost(after_mean);
    if(verbose)Rprintf("cost_diff_after(%f)=%f\n", after_mean, cost_diff_after);
    if(cost_diff_before < 0){
      if(cost_diff_after < 0){
	// f1-f2<0 meaning f1<f2 on the entire interval, so just push it1.
	push_piece(it1, last_min_log_mean, first_max_log_mean);
      }else{
	push_piece(it1, last_min_log_mean, first_log_mean);
	push_piece(it2, first_log_mean, first_max_log_mean);
      }
    }else{//f1(before)-f2(before)>=0 meaning f1(before)>=f2(before)
      if(cost_diff_after < 0){
	//f1(after)-f2(after)<0 meaning f1(after)<f2(after)
	push_piece(it2, last_min_log_mean, first_log_mean);
	push_piece(it1, first_log_mean, first_max_log_mean);
      }else{
	//f1(after)-f2(after)>=0 meaning f1(after)>=f2(after)
	push_piece(it2, last_min_log_mean, first_max_log_mean);
      }
    }
    if(verbose)Rprintf("not equal on the sides, 1 crossing point\n");
  }else{
    // "zero" crossing points. actually there may be a crossing point
    // in the interval that is numerically so close as to be identical
    // with last_min_log_mean or first_max_log_mean.
    if(verbose){
      Rprintf("not equal on the sides, zero crossing points\n");
      Rprintf("cost_diff left=%e mid=%e right=%e\n",
	     cost_diff_left, cost_diff_mid, cost_diff_right);
    }
    double cost_diff;
    if(ABS(cost_diff_mid) < NEWTON_EPSILON){
      cost_diff = cost_diff_right;
    }else{
      cost_diff = cost_diff_mid;
    }
    if(cost_diff < 0){
      push_piece(it1, last_min_log_mean, first_max_log_mean);
    }else{
      push_piece(it2, last_min_log_mean, first_max_log_mean);
    }
  }
}

void PiecewisePoissonLossLog::push_piece
(PoissonLossPieceListLog::iterator it, double min_log_mean, double max_log_mean){
  if(max_log_mean <= min_log_mean){
    // why do we need this? TODO: figure out where this is called, and
    // just do not call push_piece at all.
    return;
  }
  PoissonLossPieceListLog::iterator last=piece_list.end();
  --last;
  if(piece_list.size() &&
     sameFuns(last, it) &&
     it->prev_log_mean == last->prev_log_mean &&
     it->data_i == last->data_i){
    //it is the same function as last, so just make last extend
    //further to the right.
    last->max_log_mean = max_log_mean;
  }else{
    //it is a different function than last, so push it on to the end
    //of the list.
    piece_list.emplace_back
      (it->Linear, it->Log, it->Constant,
       min_log_mean, max_log_mean,
       it->data_i, it->prev_log_mean);
  }
}
