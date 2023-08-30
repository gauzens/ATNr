#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// #include <Rcpp.h>
using namespace arma;

// generate documentation with:
// Rcpp::compileAttributes()           # this updates the Rcpp layer from C++ to R
// roxygen2::roxygenize(roclets="rd")  # this updates the documentation based on roxygen comments

//' @name Unscaled_nuts
//' @title Store parameters and functions associated to the unscaled version of ATN including nutrient dynamics
//' @description Type the name of the class to see its methods
//' @field nb_s Total number of species
//' @field nb_b Number of basal species
//' @field nb_n Number of nutrient pool
//' @field c double: interference competition
//' @field b Matrix of attack rates (dim = number of species * number of consumers)
//' @field h Matrix of handling times (dim = number of species * number of consumers)
//' @field X vector of metabolic rates (length = number of species)
//' @field K matrix of plant nutrient efficiencies (dim = number of nutrients * number of plants)
//' @field V matrix of plant relative nutrient content (dim = number of nutrients * number of plants)
//' @field S Vector of maximum nutrient concentration (length = number of plants)
//' @field r Vector of maximum growth rate of plant species (length = number of plant species)
//' @field e Vector of assimilation efficiencies (length = number of species)
//' @field BM Vector of body masses (length = number of species)
//' @field dB Vector of local derivatives (length = number of species)
//' @field fw Adjacency matrix of the food-web (dim = number of species * number of species)
//' @field w Matrix of relative consumption rates (dim = number of species * number of consumers)
//' @field F Matrix of per-capita feeding rates (dim = number of species * number of consumers)
//' @field q hill exponent for the type of functional response
//' @field ext Extinction threshold for species
//' @field ODE Calculate the derivatives for the scaled version of the ATN model \itemize{
//' \item Parameter: bioms -  Local species biomasses
//' \item Parameter: t - Integration time point
//' \item Returns a vector of growth rate for each species at time t
//' }

class Unscaled_nuts_eco{
public:
  int nb_s; // number of species
  int nb_b; // number of basal species
  int nb_n; // number of nutrients
  int nb_f; // number of fisheries
  int n_tot; // bool prefs MORE PRECISE?
  int n_cons;
  // global nutrient turn over rate (rate of replenishment)
  // used in calculating change in nutrient concentration
  double D;

  double ext;

  vec X; // metabolic rates
  vec total_X; // population metabolism
  vec e; // assimilation efficiencies
  vec r; // growth rates of plants
  vec S; // maximal nutrient level
  vec c; // interference competition
  vec out_fluxes; // out fluxes for all species
  // body masses
  vec BM ;
  vec q;
  
  // biomasses
  vec bioms;
  // r*G for plants, as there is no need to compute it at each ODE call
  
  // NumericVector test;
  vec p;
  // Coltor of derivatives
  vec dB;
  // 
  
  Mat<int> fw; 
  
  mat b;
  // handling times
  mat h;
  // functional response
  mat F;
  // consumption rates
  mat w;
  
  //harvested by boats
  vec Harvest;
  


  // plants nutrient uptake efficiency (K(i,j): plant i on nutrient j)
  mat K; //!!!!!!!! change to same row*col in inits in R
  // relative content in the plant species' biomass
  mat V; //rows = nuts, cols = plants

  // fishery part:

  double boat_eff; // slope of the linear scaling between boat harvest and fish density
  double price_scaling; 
  double price_exp; 
  double boat_cost; // cost of a boat

  uvec targets; // indices of fisheries
  uvec internal_targets;

  // internal variables for optimisation
  vec G; // species specific growth factor Gi
  mat pow_bioms;
  vec low;
  vec uptake;
  vec bioms_non_nut;
  vec bioms_comunity;

  mat wbh_mat;
  mat wb_mat;

  // iterators
  int res;
  vec::iterator res_end;

  // slicers
  uvec animals;
  uvec plants;
  uvec nut;
  uvec non_nut;
  uvec extinct;
  uvec biology;
  uvec fisheries;
  
  Unscaled_nuts_eco(int ns, int nb, int nn, int nf):
    nb_s(ns), nb_b(nb), nb_n(nn), nb_f(nf) {
      n_tot = nb_s + nb_n;
      n_cons = nb_s - nb_b;
      // initialise vectors to 0
      X.zeros(nb_s);
      e.zeros(nb_s);
      r.zeros(nb_b);
      S.zeros(nb_n);
      G.zeros(nb_b);
      c.zeros(n_cons);
      low.zeros(n_cons);
      dB.zeros(n_tot + nb_f);
      uptake.zeros(nb_b);
      out_fluxes.zeros(nb_s);
      bioms_non_nut.zeros(nb_s);
      BM.ones(nb_s);
      total_X.zeros(n_cons);
      q.zeros(n_cons);

      // matrices
      pow_bioms.zeros(nb_s, n_cons);
      b.zeros(nb_s, n_cons);
      h.zeros(nb_s, n_cons);
      F.zeros(nb_s, n_cons);
      w.zeros(nb_s, n_cons);
      K.zeros(nb_n, nb_b);
      V.zeros(nb_n, nb_b);
      wb_mat.zeros(nb_s, n_cons);
      wbh_mat.zeros(n_cons, nb_s);

      // scalars
      D = 0.0;

      // iterator
      res_end = G.end();

      // slicers TODO check optional argument N, =100 by default
      // but inconsistant with this: https://stackoverflow.com/questions/25790634/how-to-create-a-vector-from-1n-in-c-armadillo
      animals = linspace<uvec>(nb_n + nb_b, n_tot-1, n_cons);
      plants = linspace<uvec>(nb_n, nb_n + nb_b-1, nb_b);
      nut = linspace<uvec>(0, nb_n-1, nb_n);
      non_nut = linspace<uvec>(nb_n, n_tot-1, nb_s);
      biology = linspace<uvec>(0, n_tot-1, nb_n + nb_s);
      fisheries = linspace<uvec>(n_tot, n_tot + nb_f -1, nb_f);
    }

  void initialisations(){
    // intermediate matrices for the functional response:
    // matrices used for the computation of IN feeding rates
    // wb_mat; w * b (upper part of F)
    wb_mat = (w%b);
    // adding efficiencies
    // wb_mat.each_col() *= e; !!! check that option later
    // wbh_mat: used in the below part of F, rows = cons, cols = res
    wbh_mat = (w%b%h);
    internal_targets = targets + nb_n - 1;

  }

  void print(){
    Rcpp::Rcout << "nb_s:"  << std::endl << nb_s << std::endl;
    Rcpp::Rcout << "nb_b:"  << std::endl << nb_b << std::endl; 
    Rcpp::Rcout << "bioms: " << bioms << std::endl; 
    Rcpp::Rcout << "G: " << G << std::endl; 
    // Rcout << " prey" << prey << std::endl;
  }
  
  
  

  // NumericVector ODE(double t, NumericVector bioms, NumericVector p){  // for sundials
  vec ODE(vec bioms, double t){ // for odeintr
    // here, the matricies of attack rates, handling times, feeding rates ...
    // are of diemnsions: nb species * nb consumers
    // so non square matrices
    // this is because I removed all the elements relative to a plant on a plant
    // it implies less claculations for the feeding rates
    // but more subsetting afterwards. 
    // not sure of the best option though 
    // (would it not be faster to use square matrices everywhere?)
    
    // Rcpp::Rcout << "aaa " << bioms.n_elem << std::endl;

    extinct = find(bioms < ext);
    uvec not_nut_extinct = find(extinct >= nb_n);
    extinct = extinct(not_nut_extinct);
    bioms.elem(extinct).zeros();
    // Rcpp::Rcout << bioms.t()  << std::endl;

    bioms_comunity = bioms(biology);
    bioms_non_nut = bioms.elem(non_nut);
    pow_bioms.each_col() = bioms_non_nut;
    pow_bioms = pow(pow_bioms.each_row(), q.t());

    // calculate values for feeding rates
    // F contains first the upper part of the feeding rates
    F = wb_mat % pow_bioms;

    // wbh_mat*bioms: gives for each consumer i the sum over prey j of
    // wij*hij*bij*Bj
    low = sum(wbh_mat%pow_bioms,0).t() + c%bioms_comunity(animals) + 1;
    low = low%BM(animals-nb_n);

    F.each_row() /=low.t();

    // out_fluxes: sum of out flux for each resource species, col vector
    out_fluxes = F*bioms_comunity(animals);

    // realised met. rate
    total_X = X%bioms_comunity(non_nut);

    // species specific growth factor for each plant
    // iterators are pointers, so not good here as I have to access to 
    // the ith element of different vectors
    for (res = 0; res != nb_b; ++res){
      G(res) = min(bioms_comunity(nut) / (K.col(res) + bioms_comunity(nut)));
    }
    // G could be calculated using matrix inversion, 
    // not sure there is much to win here though
    // Cholesky algs are in general in O(n^3)
    // KandBioms = k.each_col() + bioms(nut);
    // G = min(inv(KandBioms.each_col() / bioms(nut)), DIMENSION)

    // plant uptake:
    uptake = r%bioms_comunity(plants)%G;
    // derivatives for non nutrients
    dB(plants) = uptake - out_fluxes(plants-nb_n) - total_X(plants-nb_n);
 
    dB(animals) = -out_fluxes(animals-nb_n) - total_X(animals-nb_n) + 
                          (F.t()*e)%bioms_comunity(animals);
     
    // note: e may be preintegrated as constant ober time

    dB.elem(extinct).zeros();

    // derivate for nutrients
    dB(nut) = D * (S - bioms_comunity(nut)) - V*uptake;

    //fisheries
    
    Harvest = boat_eff * bioms(fisheries) % bioms(internal_targets);
    dB(fisheries) = price_scaling * pow(boat_eff * bioms(fisheries) % bioms(internal_targets), price_exp + 1) - boat_cost*boat_eff*bioms(fisheries);
    dB(internal_targets) = dB(internal_targets) - Harvest;

    // Rcpp::Rcout << "tagets id " << internal_targets.t() << std::endl;
    // Rcpp::Rcout << "db fisheries " << dB(fisheries).t()<< std::endl;
    // Rcpp::Rcout << "fisheries id" << fisheries.t()<< std::endl;
    // Rcpp::Rcout << "biology id" << biology.t()<< std::endl;
    // Rcpp::Rcout << "biology id" << biology.t()<< std::endl;

    // vec xx;
    // xx = price_scaling * pow(boat_eff*bioms(fisheries) % bioms(internal_targets), price_exp + 1);
    // Rcpp::Rcout << "gains" << xx.t()<< std::endl;
    // Rcpp::Rcout << "harvest" << Harvest.t()<< std::endl;
    // Rcpp::Rcout << "biom fisheries" << bioms(fisheries).t()<< std::endl;
    // Rcpp::Rcout << "bioms target" << bioms(internal_targets).t()<< std::endl;

    return(dB);
  }
  
};


RCPP_MODULE(Unscaled_nuts_ecoModule){
  using namespace Rcpp;
  class_<Unscaled_nuts_eco>("Unscaled_nuts_eco")
    .constructor<int, int, int, int>("constructor") //constructor
    .method("print", &Unscaled_nuts_eco::print)
    .method("ODE", &Unscaled_nuts_eco::ODE)
    .method("initialisations", &Unscaled_nuts_eco::initialisations)
    .field("nb_s", &Unscaled_nuts_eco::nb_s)
    .field("nb_b", &Unscaled_nuts_eco::nb_b)
    .field("nb_n", &Unscaled_nuts_eco::nb_n)
    .field("nb_f", &Unscaled_nuts_eco::nb_f)
    .field("BM", &Unscaled_nuts_eco::BM)
    .field("K", &Unscaled_nuts_eco::K)
    .field("D", &Unscaled_nuts_eco::D)
    .field("S", &Unscaled_nuts_eco::S)
    .field("r", &Unscaled_nuts_eco::r)
    .field("X", &Unscaled_nuts_eco::X)
    .field("e", &Unscaled_nuts_eco::e)
    .field("w", &Unscaled_nuts_eco::w)
    .field("b", &Unscaled_nuts_eco::b)
    .field("c", &Unscaled_nuts_eco::c)
    .field("q", &Unscaled_nuts_eco::q)
    .field("V", &Unscaled_nuts_eco::V)
    .field("dB", &Unscaled_nuts_eco::dB)
    .field("F", &Unscaled_nuts_eco::F)
    .field("h", &Unscaled_nuts_eco::h)
    .field("fw", &Unscaled_nuts_eco::fw)
    .field("ext", &Unscaled_nuts_eco::ext)
    .field("boat_eff", &Unscaled_nuts_eco::boat_eff)
    .field("price_scaling", &Unscaled_nuts_eco::price_scaling)
    .field("price_exp", &Unscaled_nuts_eco::price_exp)
    .field("boat_cost", &Unscaled_nuts_eco::boat_cost)
    .field("fishery_targets", &Unscaled_nuts_eco::targets)
    ;  
}
