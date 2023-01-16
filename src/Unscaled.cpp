#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// #include <Rcpp.h>
using namespace arma;

// generate documentation with:
// Rcpp::compileAttributes()           # this updates the Rcpp layer from C++ to R
// roxygen2::roxygenize(roclets="rd")  # this updates the documentation based on roxygen comments

//' @name Unscaled
//' @title Store parameters and functions associated to the unscaled version of ATN
//' @description Type the name of the class to see its methods
//' @field nb_s Total number of species
//' @field nb_b Number of basal species
//' @field c double: interference competition
//' @field X Vector of metabolic rates (length = number of species)
//' @field a Matrix of attack rates (dim = number of species * number of consumers)
//' @field h Matrix of handling times (dim = number of species * number of consumers)
//' @field e Vector of assimilation efficiencies (length = number of species)
//' @field r Vector of producers maximum growth rates (length = number of basal species)
//' @field BM Vector of body masses (length = number of species)
//' @field dB Vector of local derivatives (length = number of species)
//' @field fw Adjacency matrix of the food-web (dim = number of species * number of species)
//' @field F Matrix of per-capita feeding rates (dim = number of species * number of consumers)
//' @field q hill exponent for the type of functional response
//' @field K Carrying capacity of basal species
//' @field alpha Plant resource competition
//' @field ext Extinction threshold for species
//' @field ODE Calculate the derivatives for the scaled version of the ATN model \itemize{
//' \item Parameter: bioms -  Local species biomasses
//' \item Parameter: t - Integration time point
//' \item Returns a vector of growth rate for each species at time t
//' }


class Unscaled{
public:
  int nb_s; // number of species
  int nb_b; // number of basal species

  double ext;


  vec q;
  vec X; // metabolic rates
  vec total_X; // population metabolism
  vec e; // assimilation efficiencies
  vec r; // growth rates of plants
  // net growth rate of plant
  vec G;
  vec c; // interference competition

  // NumericVector p = NumericVector(2);
  // vector of derivatives
  vec dB;
  
  Mat<int> fw; 
  
  mat a;
  // handling times
  mat h;
  // functional response
  mat F;
  // plants carying capacity
  vec K; 
  // plant competition for resources
  mat alpha;

  // internal variables for optimisation
  // index of plants (optimisation purpose)
  mat ah;
  mat pow_bioms;
    // contains lower part of the feeding rate ratio
  vec low;

  vec out_fluxes;

  // plant competition
  vec s;

  // body masses
  vec BM;

  // slicers 
  uvec animals;
  uvec plants;
  uvec extinct;

  double out = 0;
  int i = 0;
  int n_cons;
  
  Unscaled(int ns, int nb):
    nb_s(ns), nb_b(nb) {
    n_cons = nb_s - nb_b;
    X.zeros(nb_s);
    e.zeros(nb_s);
    r.zeros(nb_b);
    G.zeros(nb_b);
    c.zeros(n_cons);
    low.zeros(n_cons);
    dB.zeros(nb_s);
    out_fluxes.zeros(nb_s);
    total_X.zeros(n_cons);
    BM.zeros(nb_s);
    q.zeros(n_cons);

    // matrices
    pow_bioms.zeros(nb_s, n_cons);
    a.zeros(nb_s, n_cons);
    h.zeros(nb_s, n_cons);
    F.zeros(nb_s, n_cons);
    alpha.zeros(nb_b, nb_b);
	  // t_ah.zeros(n_cons, nb_s);
    ah.zeros(nb_s, n_cons);

    animals = linspace<uvec>(nb_b, nb_s-1, n_cons);
    plants = linspace<uvec>(0, nb_b-1, nb_b);

    // q = 0;

    ext = 0;

    }


  void print(){
    Rcpp::Rcout << "plants: "  << plants << std::endl;
    Rcpp::Rcout << "animals: "  <<  animals << std::endl;
  }

  void initialisations(){
    // intermediate matrices for the functional response:

    // t_ah = (a%h).t();
    ah = (a%h);
  }

  
  vec ODE(vec bioms, double t){ 
    
    extinct = find(bioms < ext);
    bioms.elem(extinct).fill(0.0);
    pow_bioms.each_col() = bioms;
    pow_bioms = pow(pow_bioms.each_row(), q.t());
    // calculate the upper part of the feeding rate function
    F = a % pow_bioms;
    // and the lower part
  	low = sum(ah%pow_bioms,0).t() + c%bioms(animals) + 1;
    // divide both to obtained the matrix of feeding rates
   	F.each_row() /=low.t();

    // now, compute outfluxes for every species:
  	out_fluxes = F*bioms(animals);

  	// realised met. rate
  	total_X = X%bioms;

  	// realised growth rate for plants:
    s = alpha * bioms(plants);
    // s = bioms(plants);
  	G = 1 - s / K;

    // Rcpp::Rcout << G.t();

  	dB(plants) = r%bioms(plants)%G - out_fluxes(plants) - total_X(plants);
  	dB(animals) = -out_fluxes(animals) - total_X(animals) + 
                          (F.t()*e)%bioms(animals);
    
    // derivates for plants
    dB.elem(extinct).fill(0.0);
    
    return dB;
  }
  
};



RCPP_MODULE(UnscaledModule){
  using namespace Rcpp;
  class_<Unscaled>("Unscaled")
    .constructor<int, int>("constructor") //constructor
    .method("print", &Unscaled::print)
    .method("ODE", &Unscaled::ODE)
    .method("initialisations", &Unscaled::initialisations)
    .field("nb_s", &Unscaled::nb_s)
    .field("nb_b", &Unscaled::nb_b)
    .field("K", &Unscaled::K)
    .field("r", &Unscaled::r)
    .field("X", &Unscaled::X)
    .field("e", &Unscaled::e)
    .field("a", &Unscaled::a)
    .field("c", &Unscaled::c)
    .field("h", &Unscaled::h)
    .field("q", &Unscaled::q)
    .field("dB", &Unscaled::dB)
    .field("BM", &Unscaled::BM)
    .field("F", &Unscaled::F)
    .field("fw", &Unscaled::fw)
    .field("ext", &Unscaled::ext)
    .field("alpha", &Unscaled::alpha)
    ;  
}
