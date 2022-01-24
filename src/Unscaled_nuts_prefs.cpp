#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// #include <Rcpp.h>
using namespace arma;
using namespace Rcpp;

// generate documentation with:
// Rcpp::compileAttributes()           # this updates the Rcpp layer from C++ to R
// roxygen2::roxygenize(roclets="rd")  # this updates the documentation based on roxygen comments

//' @name Unscaled_nuts_prefs
//' @title Store parameters and functions associated to the unscaled version of ATN including nutrient dynamics
//' @description Type the name of the class to see its methods
//' @field nb_s Total number of species
//' @field nb_b Number of basal species
//' @field nb_n Number of nutrient pool
//' @field X Coltor of metabolic rates (length = number of species)
//' @field K1 Vector of maximum feeding rates (length = number of consumers)
//' @field K2 Vector of producers maximum growth rates (length = number of basal species)
//' @field e Vector of assimilation efficiencies (length = number of species)
//' @field BM Vector of body masses (length = number of species)
//' @field dB Vector of local derivatives (length = number of species)
//' @field B0 Vector of half saturation densities (length = number of consumers)
//' @field fw Adjacency matrix of the food-web (dim = number of species * number of species)
//' @field w Matrix of relative consumption rates (dim = number of species * number of consumers)
//' @field F Matrix of per-capita feeding rates (dim = number of species * number of consumers)
//' @field q parameter for the type of functional response (hill exponent = 1 + q)
//' @field K Carrying capacity of basal species
//' @field ext extinction threshold for species
//' @field ODE Calculate the derivatives for the scaled version of the ATN model \itemize{
//' \item Parameter: bioms -  Local species biomasses
//' \item Parameter: t - Integration time point
//' \item Returns a Coltor of growth rate for each species at time t
//' }

class Unscaled_nuts_prefs{
public:
  int nb_s; // number of species
  int nb_b; // number of basal species
  int nb_n = 2; // number of nutrients
  int n_tot = nb_s + nb_n; // bool prefs MORE PRECISE?
  // global nutrient turn over rate (rate of replenishment)
  // used in calculating change in nutrient concentration
  double D = 0.25;
  double q;
  double c;
  double ext;
  double temperature;

  vec X; // metabolic rates
  vec total_X; // population metabolism
  vec e; // assimilation efficiencies
  vec r; // growth rates of plants
  vec S; // maximal nutrient level
  vec out_fluxes; // out fluxes for all species
  // body masses
  vec BM ;
  vec log_BM;
  
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

  // plants nutrient uptake efficiency (K(i,j): plant i on nutrient j)
  mat K; //!!!!!!!! change to same row*col in inits in R
  // relative content in the plant species' biomass
  mat V; //rows = nuts, cols = plants

  // internal variables for optimisation
  vec G; // species specific growth factor Gi
  vec pow_bioms;
  vec low;
  vec uptake;
  vec bioms_non_nut;
  uvec prey;

  mat wbh_mat;
  mat wb_mat;

  // iterators
  int res;
  int cons;
  vec::iterator res_end;

  // slicers
  uvec animals;
  uvec plants;
  uvec nut;
  uvec non_nut;
  uvec extinct;
  
  Unscaled_nuts_prefs(int ns, int nb, int nn):
    nb_s(ns), nb_b(nb), nb_n(nn) {
      int n_tot = nb_s + nb_n;
      int n_cons = nb_s - nb_b;
      // initialise vectors to 0
      X.zeros(nb_s);
      e.zeros(nb_s);
      r.zeros(nb_b);
      S.zeros(nb_n);
      G.zeros(nb_b);
      pow_bioms.zeros(nb_s);
      low.zeros(n_cons);
      dB.zeros(n_tot);
      uptake.zeros(nb_b);
      out_fluxes.zeros(nb_s);
      low.zeros(n_cons);
      bioms_non_nut.zeros(nb_s);
      BM.ones(nb_s);

      // matrices
      b.zeros(nb_s, n_cons);
      h.zeros(nb_s, n_cons);
      F.zeros(nb_s, n_cons);
      w.zeros(nb_s, n_cons);
      K.zeros(nb_n, nb_b);
      V.zeros(nb_n, nb_b);
      wb_mat.zeros(nb_s, n_cons);
      wbh_mat.zeros(n_cons, nb_s);

      // iterator
      res_end = G.end();


      // slicers TODO check optional argument N, =100 by default
      // but inconsistant with this: https://stackoverflow.com/questions/25790634/how-to-create-a-vector-from-1n-in-c-armadillo
      animals = linspace<uvec>(nb_n + nb_b, n_tot-1, n_cons);
      plants = linspace<uvec>(nb_n, nb_n + nb_b-1, nb_b);
      nut = linspace<uvec>(0, nb_n-1, nb_n);
      non_nut = linspace<uvec>(nb_n, n_tot-1, nb_s);
    }

  void initialisations(){
    // intermediate matrices for the functional response:
    // matrices used for the computation of IN feeding rates
    // wb_mat; w * b (upper part of F)
    wb_mat = (w%b);
    // adding efficiencies
    // wb_mat.each_col() *= e; !!! check that option later
    // wbh_mat: used in the below part of F, rows = cons, cols = res
    wbh_mat = (w%b%h).t();
  }

  void print(){
    Rcpp::Rcout << "nb_s:"  << std::endl << nb_s << std::endl;
    Rcpp::Rcout << "nb_b:"  << std::endl << nb_b << std::endl; 
    Rcpp::Rcout << "bioms: " << bioms << std::endl; 
    Rcpp::Rcout << "G: " << G << std::endl; 
    // Rcout << " prey" << prey << std::endl;
  }
  
    mat compute_w(vec bioms, vec log_BM, double Temperature){
      // I loop over cons first for simplicity
    // could be optimised by looping over resources first to match C storage of matricies
    double prod, alpha, omega, eps, total;
    
    // I rescale the biomasses to nb ind. that fall into the values that I used to predict sdn paramters.
    // note: scale to hve species with the max generalsim to reach 4
    vec scaled_bioms = (bioms / sum(bioms))*4*nb_s/6;
    vec bioms_noNut = scaled_bioms(non_nut); //the biomass vector contains also nutrients, but they are not in fw... to change?
    // Rcout << "scale done " << std::endl;
    for (cons = nb_b; cons < nb_s; cons++){
      prey = find(fw.col(cons)>0);
      vec bioms_prey = bioms_noNut(prey);
      prod = sum(bioms_prey);
      total = 0.0;
      for (res = 0; res < nb_s; res++){
        // Rcout << *res << " " << *cons << "  ||  ";
        alpha = 0.286013574 + 0.001660838*prod + 0.102041062*Temperature;
        omega = 0.5200681304 + 0.0884532552*log_BM(cons) - 0.0005778026*prod + 0.0002432819*prod*Temperature;
        eps = -1.2291493541 + 0.7314166498*log_BM(cons) + 0.0005495436*prod + 0.0833382603*Temperature - 0.0384199792*log_BM(cons)*Temperature - 0.0000650646*prod*Temperature;
        if (fw(res, cons) > 0){
          w(res, cons-nb_b) = skew_norm(log_BM(res), eps, omega, alpha);
          // Rcout << "resBM " << log_BM(*res) << "  alpha" << alpha << "  omega: " << omega << "  eps: " << eps << std::endl;
          // Rcout << "w: " << w(*res, *cons - nb_b) << "  cons: " << *cons - nb_b << "  res: " << *res << std::endl;
          total = total +  w(res, cons-nb_b);
        }
      }
      
      // Rcout << std::endl;
      
      // Rcout << "w_tot: " << total << "   cons: " << *cons << std::endl;
      if (total > 0.0) {
        w.col(cons-nb_b) = w.col(cons-nb_b) / total;
      }else{
        w.col(cons-nb_b) = 0.0;
      }
    }
    return(w);
  }
  
    double skew_norm(double x, double xi, double omega, double alpha){
    // no skew normal function in Rcpp
    // have to call the R function from C directly
    // the R function return a pointer (SEXP) that has to be converted to double with REAL
    //  !!! library sn needs to be loaded in the R session
    Rcpp::Environment base("package:sn");
    Rcpp::Function dsn_r = base["dsn"];
    SEXP a =  dsn_r(Rcpp::_["x"] = x, Rcpp::_["xi"] = xi, Rcpp::_["omega"] = omega, Rcpp::_["alpha"] = alpha);
    return *REAL(a);
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
    bioms.elem(extinct).fill(0.0);
    // Rcpp::Rcout << bioms.t()  << std::endl;
    bioms_non_nut = bioms.elem(non_nut);
    pow_bioms = pow(bioms_non_nut, 1+q);

    w = compute_w(bioms, log_BM, temperature);
    wb_mat = (w%b);
    wbh_mat = (w%b%h).t();

    F = wb_mat.each_col() % pow_bioms;
    // Rcpp::Rcout << "aaaa"  << std::endl;
    // wbh_mat*bioms: gives for each consumer i the sum over prey j of
    // wij*hij*bij*Bj
    low = wbh_mat*pow_bioms + c*bioms(animals) + 1;
    // Rcpp::Rcout << "bbbb"  << std::endl;
    low = low%BM(animals-nb_n);

    // Rcpp::Rcout << "low calculated: " << low.n_elem << std::endl;
    F.each_row() /=low.t();
    // Rcpp::Rcout << "F calculated"  << std::endl;
    // out_fluxes: sum of out flux for each resource species, col vector
    out_fluxes = F*bioms(animals);
    // Rcpp::Rcout << "out_fluxes"  << std::endl;
    // realised met. rate
    total_X = X%bioms(non_nut);
    // Rcpp::Rcout << "X"  << std::endl;

    // species specific growth factor for each plant
    // iterators are pointers, so not good here as I have to access to 
    // the ith element of different vectors
    for (res = 0; res != nb_b; ++res){
      // Rcpp::Rcout << res << "  ";
      G(res) = min(bioms(nut) / (K.col(res) + bioms(nut)));
      // Rcpp::Rcout << "!!" << res << "  " << G(res) << "  " << bioms(res+nb_n) << "-------";
    }
    // G could be calculated like that, not sure there is uch to win here though
    // Cholesky algs are in general in O(n^3)
    // KandBioms = k.each_col() + bioms(nut);
    // G = min(inv(KandBioms.each_col() / bioms(nut)), DIMENSION)

    // Rcpp::Rcout << "G: " <<  G.t() << std::endl;
    // plant uptake:
    uptake = r%bioms(plants)%G;
    // Rcpp::Rcout << std::endl << uptake.t() << std::endl;
    // derivatives for non nutrients
    dB(plants) = uptake - out_fluxes(plants-nb_n) - total_X(plants-nb_n);
     // Rcpp::Rcout << "db_plant"  << std::endl;  
    // dB(animals) = -out_fluxes(animals-nb_n) - total_X(animals-nb_n) + 
    //                       F.t()*(e%bioms(non_nut));
    dB(animals) = -out_fluxes(animals-nb_n) - total_X(animals-nb_n) + 
                      (F.t()*e)%bioms(animals);

    // Rcpp::Rcout << "V: " << V << std::endl;
    // Rcpp::Rcout << "db"  << std::endl;       
    // note: e may be preintegrated as constant ober time

    dB.elem(extinct).fill(0.0);
    // derivate for nutrients
    // dB(nut) = V*uptake;
     // Rcpp::Rcout << "V*uptake"  << std::endl;
    dB(nut) = D * (S - bioms(nut)) - V*uptake;
    // Rcpp::Rcout << dB(nut).t()  << std::endl; 
    // Rcpp::Rcout << "db_nuts"  << std::endl; 
    
    return(dB);
  }
  
};



RCPP_MODULE(Unscaled_nuts_prefsModule){
  using namespace Rcpp;
  class_<Unscaled_nuts_prefs>("Unscaled_nuts_prefs")
    .constructor<int, int, int>("constructor") //constructor
    .method("print", &Unscaled_nuts_prefs::print)
    .method("ODE", &Unscaled_nuts_prefs::ODE)
    .method("initialisations", &Unscaled_nuts_prefs::initialisations)
    .field("nb_s", &Unscaled_nuts_prefs::nb_s)
    .field("nb_b", &Unscaled_nuts_prefs::nb_b)
    .field("nb_n", &Unscaled_nuts_prefs::nb_n)
    .field("BM", &Unscaled_nuts_prefs::BM)
    .field("log_BM", &Unscaled_nuts_prefs::log_BM)
    .field("K", &Unscaled_nuts_prefs::K)
    .field("D", &Unscaled_nuts_prefs::D)
    .field("S", &Unscaled_nuts_prefs::S)
    .field("r", &Unscaled_nuts_prefs::r)
    .field("X", &Unscaled_nuts_prefs::X)
    .field("e", &Unscaled_nuts_prefs::e)
    .field("w", &Unscaled_nuts_prefs::w)
    .field("b", &Unscaled_nuts_prefs::b)
    .field("c", &Unscaled_nuts_prefs::c)
    .field("h", &Unscaled_nuts_prefs::h)
    .field("q", &Unscaled_nuts_prefs::q)
    .field("V", &Unscaled_nuts_prefs::V)
    .field("temperature", &Unscaled_nuts_prefs::temperature)
    .field("dB", &Unscaled_nuts_prefs::dB)
    .field("D", &Unscaled_nuts_prefs::D)
    .field("F", &Unscaled_nuts_prefs::F)
    .field("uptake", &Unscaled_nuts_prefs::uptake)
    .field("h", &Unscaled_nuts_prefs::h)
    .field("fw", &Unscaled_nuts_prefs::fw)
    .field("ext", &Unscaled_nuts_prefs::ext)
    ;  
}
