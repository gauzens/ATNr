#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// #include <Rcpp.h>
using namespace arma;

// generate documentation with:
// Rcpp::compileAttributes()           # this updates the Rcpp layer from C++ to R
// roxygen2::roxygenize(roclets="rd")  # this updates the documentation based on roxygen comments

//' @name Scaled
//' @title Store parameters and functions associated to the scaled version of ATN
//' @description Type the name of the class to see its methods
//' @field nb_s Total number of species
//' @field nb_b Number of basal species
//' @field c double: interference competition
//' @field X Vector of metabolic rates (length = number of species)
//' @field max_feed Vector of maximum feeding rates (length = number of consumers)
//' @field e Vector of assimilation efficiencies (length = number of species)
//' @field r Vector of producers maximum growth rates (length = number of basal species)
//' @field BM Vector of body masses (length = number of species)
//' @field dB Vector of local derivatives (length = number of species)
//' @field B0 Vector of half saturation densities (length = number of consumers)
//' @field fw Adjacency matrix of the food-web (dim = number of species * number of species)
//' @field w Matrix of relative consumption rates (dim = number of species * number of consumers)
//' @field F Matrix of per-capita feeding rates (dim = number of species * number of consumers)
//' @field q hill exponent for the type of functional response
//' @field K Carrying capacity of basal species
//' @field ext Extinction threshold for species
//' @field alpha Plant resource competition
//' @field ODE Calculate the derivatives for the scaled version of the ATN model \itemize{
//' \item Parameter: bioms -  Local species biomasses
//' \item Parameter: t - Integration time point
//' \item Returns a vector of growth rate for each species at time t
//' }


class Scaled{
public:
  // number of species
  int nb_s;
  // number of basal species
  int nb_b;
  // hill coefficient
  double q;

  // extinction threshold
  double ext;
  // carrying capacity of plants
  double K;

  
  // metabolic rates
  vec X;
  vec total_X; // population metabolism
  vec c; // interference competition
  // maximum feeding rate. Could be a vetor of length = nb of animals (maybe more intuitive for users?)
  vec max_feed;
  // multiplication of feeding rate and metabolism
  vec xy;
  // half stauration density
  vec B0;
  // assimilation efficiencies
  vec e;
  // mass specific growth rates of plants
  vec r;
  // net growth rate of plant
  vec G;

  // plant competition for resources
  mat alpha;


  // body masses:
  vec BM;
  // vec log_BM;
    
  // vector of derivatives
  vec dB;

    
  Mat<int> fw; 

  // functional response
  mat F;
  // consumption rates
  mat w;
  mat wt;

  vec pow_bioms;
  vec pow_B0;


  // contains lower part of the feeding rate ratio
  vec low;

  vec out_fluxes;

    // plant competition
  vec s;

  // slicers 
  uvec animals;
  uvec plants;
  uvec extinct;

  double out = 0.0;
  int i = 0;
  
  Scaled(int s, int b):
    nb_s(s), nb_b(b) {
    int n_cons = nb_s - nb_b;
    ext = 1e-6;

    X.zeros(nb_s);
    e.zeros(nb_s);
    r.zeros(nb_b);
    G.zeros(nb_b);
    pow_bioms.zeros(nb_s);
    low.zeros(n_cons);
    dB.zeros(nb_s);
    out_fluxes.zeros(nb_s);
    BM.ones(nb_s);
    max_feed.zeros(n_cons);
    c.zeros(n_cons);
    xy.zeros(n_cons);
    total_X.zeros(n_cons);
    B0.zeros(n_cons);
    pow_B0.zeros(n_cons);
    // matrices
    F.zeros(nb_s, n_cons);
    w.zeros(nb_s, n_cons);
    alpha.zeros(nb_b, nb_b);
    wt = w.t();

    animals = linspace<uvec>(nb_b, nb_s-1, n_cons);
    plants = linspace<uvec>(0, nb_b-1, nb_b);
  }


  void print(){
    Rcpp::Rcout << "plants: "  << plants << std::endl;
    Rcpp::Rcout << "animals: "  <<  animals << std::endl;
    // Rcpp::Rcout << "F: "  <<  F << std::endl;
  }
  
  
  void initialisations(){
    // intermediate matrices for the functional response:
    wt = w.t(); // is that truly needed? I'm even not sure that transpose are calculated
    xy = X(animals)%max_feed;
    pow_B0 = pow(B0, q);
  }



  // NumericVector ODE(double t, NumericVector bioms, NumericVector p){  // for sundials
  vec ODE(vec bioms, double t){ // for odeintr

    // checkUserInterrupt();
    extinct = find(bioms < ext);
    bioms.elem(extinct).fill(0.0);

    pow_bioms = pow(bioms, q);
    // Rcpp::Rcout << "aa " << std::endl;
    // calculate values for feeding rates
    F = w.each_col() % pow_bioms;
    // Rcpp::Rcout << "bb " << std::endl;
    // at this point: Fij = wij * Bi
    // now compute the lower part of the ratio:
    low = w.t()*pow_bioms + c%bioms(animals) + pow_B0;

    // and make the division
    F.each_row() /=low.t();

    // out_fluxes: sum of out flux for each resource species, col vector
    out_fluxes = F*(bioms(animals)%xy);

    // realised met. rate
    total_X = X%bioms;
    
    // realised growth rate for plants:
    s = alpha * bioms(plants);
    G = 1 - s / K;
    // Rcpp::Rcout << "dd " << std::endl;
    dB(plants) = r%bioms(plants)%G - out_fluxes(plants) - total_X(plants);
    dB(animals) = -out_fluxes(animals) - total_X(animals) + 
                          (F.t()*e)%(bioms(animals)%xy);

    dB.elem(extinct).fill(0.0);

    return dB;
  }

};



RCPP_MODULE(ScaledModule){
using namespace Rcpp;
  class_<Scaled>("Scaled")
  .constructor<int, int>("constructor")
  .method("print", &Scaled::print)
  .method("ODE", &Scaled::ODE)
  .method("initialisations", &Scaled::initialisations)
  .field("nb_s", &Scaled::nb_s)
  .field("nb_b", &Scaled::nb_b)
  .field("BM", &Scaled::BM)
  // .field("log_BM", &Scaled::log_BM)
  .field("r", &Scaled::r)
  .field("X", &Scaled::X)
  .field("e", &Scaled::e)
  .field("w", &Scaled::w)
  .field("B0", &Scaled::B0)
  .field("c", &Scaled::c)
  .field("q", &Scaled::q)
  .field("dB", &Scaled::dB)
  .field("F", &Scaled::F)
  .field("fw", &Scaled::fw)
  .field("max_feed", &Scaled::max_feed)
  .field("K", &Scaled::K)
  .field("ext", &Scaled::ext)
  .field("alpha", &Scaled::alpha)
  ;  
}
