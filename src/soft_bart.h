#ifndef SOFT_BART_H
#define SOFT_BART_H

#include <RcppArmadillo.h>
#include "functions.h"

struct Hypers;
struct Node;
struct Forest;
struct Opts;

typedef std::vector<Node*> Nodevec;


class Hypers {

  public:

  // Quantities
  double alpha;
  double beta;
  double gamma;
  double sigma;
  double sigma_hat;
  double sigma_mu;
  double sigma_mu_hat;
  double shape;
  double tau_rate;
  double temperature;
  int num_tree;
  int num_groups;
  arma::vec s;
  arma::vec logs;
  arma::vec logZ;
  arma::uvec group;
  std::vector<std::vector<unsigned int> > group_to_vars;

  // Basis Function Hypers
  double length_scale;
  double shape_length_scale;
  double rate_length_scale;

  // Updates
  void UpdateSigma(const arma::vec& residuals);
  void UpdateSigmaMu(const arma::vec& means);
  void UpdateAlpha();

  // Construct/Destructure
  Hypers(Rcpp::List hypers);

  // Utilities
  int SampleVar() const;
  
  // Stuff for update alpha
  double a_hat;
  double b_hat;
  double alpha_shape_1;
  double alpha_scale;

};

class Node {

  public:

  bool is_leaf;
  bool is_root;
  bool is_left;
  Node* left;
  Node* right;
  Node* parent;

  // Branch params
  int depth;
  int var;
  double val;
  arma::sp_vec lowers;
  arma::sp_vec uppers;
  arma::sp_uvec parent_split_indicators;
  double tau;

  // Leaf parameters
  double mu;
  // Data for computing weights
  double current_weight;

  // Basis function parameters for targeted smoothing
  double frequency;
  double offset;
  // bool is_sine;

  // Constructor / Destructor
  Node(const Hypers& hypers);
  Node(Node* parent, const Hypers& hypers, bool is_left);
  ~Node();

  // Updates and such
  void BirthLeaves(const Hypers& hypers);
  void GenBelow(const Hypers& hypers);
  void get_limits_below();
  void GetW(const arma::mat& X, int i);
  void DeleteLeaves();
  void UpdateMu(const arma::vec& Y, const arma::vec& t, const arma::mat& X, const Hypers& hypers);
  void UpdateTau(const arma::vec& Y, const arma::vec& t, const arma::mat& X, const Hypers& hypers);
  void SetTau(double tau_new);
  double loglik_tau(double tau_new, const arma::mat& X, const arma::vec& Y, const arma::vec& t, const Hypers& hypers);
  
};

struct Opts {

  int num_burn;
  int num_thin;
  int num_save;
  int num_print;

  bool update_sigma_mu;
  bool update_sigma;
  bool update_s;
  bool update_alpha;
  bool update_tau;

  Opts(Rcpp::List opts);

};

class Forest {

  private:

  Nodevec trees;
  Hypers hypers;
  Opts opts;

  public:

  Forest(Rcpp::List hypers, Rcpp::List opts);
  ~Forest();
  void IterateGibbs(arma::vec& Y_hat, const arma::mat&X, const arma::vec& Y, const arma::vec& t);
  arma::mat do_gibbs(const arma::mat& X,
                     const arma::vec& Y,
                     const arma::vec& t,
                     const arma::mat& X_test,
                     const arma::vec& t_test,
                     int num_iter);

  // Getters and Setters, from R
  void set_s(const arma::vec& s);
  void set_sigma(double sigma);
  double get_sigma();
  arma::vec get_s();
  arma::uvec get_counts();
  Rcpp::List get_params();
  
  // Prediction interfact
  arma::vec predict_vec(const arma::mat& X_test, const arma::vec& t_test);

};

// Various tree probabilities
double SplitProb(Node* node, const Hypers& hypers);
double growth_prior(int leaf_depth, const Hypers& hypers);


// Functions for dealing with trees
void leaves(Node* x, std::vector<Node*>& leafs);
std::vector<Node*> leaves(Node* x);
std::vector<Node*> not_grand_branches(Node* tree);
void not_grand_branches(std::vector<Node*>& ngb, Node* node);
void branches(Node* n, Nodevec& branch_vec);
std::vector<Node*> branches(Node* root);
Node* rand(std::vector<Node*> ngb);

// Sampling and selecting nodes
Node* birth_node(Node* tree, double* leaf_node_probability);
double probability_node_birth(Node* tree);
Node* death_node(Node* tree, double* p_not_grand);

// Functions for collecting quantities from the forest
arma::vec get_means(std::vector<Node*>& forest);
void get_means(Node* node, std::vector<double>& means);
arma::vec predict(const std::vector<Node*>& forest,
                  const arma::mat& X,
                  const arma::vec& t,
                  const Hypers& hypers);
arma::vec predict(Node* node,
                  const arma::mat& X,
                  const arma::vec& t,
                  const Hypers& hypers);
arma::uvec get_var_counts(std::vector<Node*>& forest, const Hypers& hypers);
void get_var_counts(arma::uvec& counts, Node* node, const Hypers& hypers);

// Functions for computing with the forest
void GetSuffStats(Node* n, const arma::vec& y,
                  const arma::vec& t,
                  const arma::mat& X, const Hypers& hypers,
                  arma::vec& mu_hat_out, arma::mat& Omega_inv_out);
double LogLT(Node* n, const arma::vec& Y,
             const arma::vec& t,
             const arma::mat& X, const Hypers& hypers);


// MCMC on the trees
void birth_death(Node* tree, const arma::mat& X, const arma::vec& Y, const arma::vec& t,
                 const Hypers& hypers);
void node_birth(Node* tree, const arma::mat& X, const arma::vec& Y, const arma::vec& t,
                const Hypers& hypers);
void node_death(Node* tree, const arma::mat& X, const arma::vec& Y, const arma::vec& t,
                const Hypers& hypers);
Node* draw_prior(Node* tree, const arma::mat& X, const arma::vec& Y, const arma::vec& t,
                 const Hypers& hypers);

// Functions for the perturb algorithm
double calc_cutpoint_likelihood(Node* node);
std::vector<double> get_perturb_limits(Node* branch);
void perturb_decision_rule(Node* tree,
                           const arma::mat& X,
                           const arma::vec& Y,
                           const arma::vec& t,
                           const Hypers& hypers);

// MCMC Functions for tree parameters
double logprior_tau(double tau, double tau_rate);
double tau_proposal(double tau);
double log_tau_trans(double tau_new);
arma::vec get_tau_vec(const std::vector<Node*>& forest);

// The backfitting algorithm
void TreeBackfit(std::vector<Node*>& forest, arma::vec& Y_hat,
                 const Hypers& hypers, const arma::mat& X, const arma::vec& Y, const arma::vec& t,
                 const Opts& opts);

/* arma::vec loglik_data(const arma::vec& Y, const arma::vec& Y_hat, const Hypers& hypers); */
/* void IterateGibbsWithS(std::vector<Node*>& forest, arma::vec& Y_hat, */
/*                        Hypers& hypers, const arma::mat& X, const arma::vec& Y, */
/*                        const Opts& opts); */
/* void IterateGibbsNoS(std::vector<Node*>& forest, arma::vec& Y_hat, */
/*                      Hypers& hypers, const arma::mat& X, const arma::vec& Y, */
/*                      const Opts& opts); */
/* double alpha_to_rho(double alpha, double scale); */
/* double rho_to_alpha(double rho, double scale); */
/* double forest_loglik(std::vector<Node*>& forest, double gamma, double beta); */
/* double tree_loglik(Node* node, int node_depth, double gamma, double beta); */
void UpdateS(std::vector<Node*>& forest, Hypers& hypers);

/* // PERTURB STUFF */

#endif
