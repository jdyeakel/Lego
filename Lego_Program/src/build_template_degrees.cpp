#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
List build_template_degrees(
  double num_play,
  NumericVector probs;
  double tr_avoid
){
  double p_a = probs(0);
  double p_n = probs(1);
  double p_i = probs(2);
  double p_m = probs(3);
  double p_e = probs(4);

  //Calculate pairwise probabilities

  double pr_ne = p_n*(p_e/(p_e+p_n+p_i+p_m));
  double pr_nn = p_n*(p_n/(p_e+p_n+p_i+p_m));
  double pr_ni = p_n*(p_i/(p_e+p_n+p_i+p_m));
  double pr_nm = p_n*(p_m/(p_e+p_n+p_i+p_m));
  double pr_ia = p_i*(p_a/(p_e+p_a+p_n+p_i));
  double pr_ie = p_i*(p_e/(p_e+p_a+p_n+p_i));
  double pr_ii = p_i*(p_i/(p_e+p_a+p_n+p_i));
  double pr_aa = p_a*(p_a/(p_a+p_i));
  double pr_ee = p_e*(p_e/(p_i+p_n+p_e));

  NumericVector pw_prob(9);
  pw_prob(0) = pr_ne;
  pw_prob(1) = pr_nn;
  pw_prob(2) = pr_ni;
  pw_prob(3) = pr_nm;
  pw_prob(4) = pr_ia;
  pw_prob(5) = pr_ie;
  pw_prob(6) = pr_ii;
  pw_prob(7) = pr_aa;
  pw_prob(8) = pr_ee;

  NumericVector pw_prob_sort = sort(pw_prob);

  NumericVector prob_line += pw_prob_sort;

  //  A = 1
  //  N = 2
  //  I = 3
  //  M = 4
  //  E = 5

  double N_e = p_e * num_play.L * (num_play.L - 1);
  double mean_k = N_e / (num_play - 1);
  int aux = 1;
  while (aux == 1) {
    NumericVector degrees_draw = round(rexp((num_play-1),1/mean_k)),0);
    NumericVector degrees(num_play);
    degrees(0) = 0;
    for (int i=1;i<num_play;i++) {
      degrees(i) = degrees_draw(i);
    }
    if(max(degrees)<(num_play-1)){
      aux=0;
    }
  }s


}
