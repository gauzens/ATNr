#include <Rcpp.h>
// #include <RcppArmadillo.h>
using namespace Rcpp;

// generate documentation with:
// Rcpp::compileAttributes()           # this updates the Rcpp layer from C++ to R
// roxygen2::roxygenize(roclets="rd")  # this updates the documentation based on roxygen comments

//' @name Unscaled_loops
//' @title Store parameters and functions associated to the unscaled version of ATN
//' @description To not use. For testing purpose only. Please use Rcpp_Unscaled instead. 

class Unscaled_loops{
public:
  int nb_s; // number of species
  int nb_b; // number of basal species

  double q;
  double test_double;
  double ext;
  double s;


  NumericVector X; // metabolic rates
  NumericVector e; // assimilation efficiencies
  NumericVector r; // growth rates of plants
  NumericVector c;  // interference competitio  
  // body masses
  NumericVector BM;
  // NumericVector log_BM;
  
  // NumericVector test;
  // NumericVector p;
  NumericVector dB;
  
  LogicalMatrix fw; 
  
  NumericMatrix a;
  // handling times
  NumericMatrix h;
  // functional response
  NumericMatrix F;
  // plants carying capacity
  NumericVector K;
  NumericMatrix alpha;
  // internal variables for optimisation
  // index of plants (optimisation purpose)

  IntegerVector plants;
  IntegerVector animals;
  IntegerVector all;

  NumericVector pow_bioms;
  
  // LogicalVector prey = fw[_,1] = 1;
  IntegerVector::iterator cons;
  IntegerVector::iterator cons2;
  IntegerVector::iterator res;
  IntegerVector::iterator nut; //not properly needed, but more readable
  NumericVector uptake;
  double out;
  int i;
  
  Unscaled_loops(int ns, int nb):
    nb_s(ns), nb_b(nb) {


    X = NumericVector(nb_s); // metabolic rates
    e = NumericVector(nb_s); // assimilation efficiencies
    r = NumericVector(nb_b); // growth rates of plants
    c = NumericVector(nb_s - nb_b);  // interference competition
    
    // body masses
    BM = NumericVector(nb_s);
    // log_BM = NumericVector(nb_s);
    
    // p = NumericVector(2);
    // vector of derivatives
    dB = NumericVector(nb_s);
    
    fw = LogicalMatrix(nb_s, nb_s); 
    
    a = NumericMatrix(nb_s,nb_s - nb_b);
    // handling times
    h = NumericMatrix(nb_s,nb_s - nb_b);
    // functional response
    F = NumericMatrix(nb_s,nb_s - nb_b);
    // plants carying capacity
    K = NumericVector(nb_b); 

    alpha = NumericMatrix(nb_b, nb_b);
    // internal variables for optimisation
    // index of plants (optimisation purpose)

    plants = Range(0, nb_b - 1);
    animals = Range(nb_b, nb_s - 1);
    all = Range(0, nb_s - 1);

    pow_bioms = NumericVector(nb_s);
    
    // LogicalVector prey = fw[_,1] = 1;
    IntegerVector::iterator cons;
    IntegerVector::iterator cons2;
    IntegerVector::iterator res;
    IntegerVector::iterator nut; //not properly needed, but more readable
    uptake = NumericVector(nb_b);
    out = 0;
    i = 0;
    q = 0.0;
    ext = 0.0;
    // s = 0.0;

    }
  
  
  void print(){
    Rcout << "nb_s:"  << std::endl << nb_s << std::endl;
    Rcout << "nb_b:"  << std::endl << nb_b << std::endl;
    Rcout << "plants: " << plants << std::endl; 
    Rcout << "dbplant " << dB[plants] << std::endl;
    Rcout << "r[plants]" << r[plants] << std::endl;
    // Rcout << " prey" << prey << std::endl;
  }
  
  
  
  
  double F_rate(int prey, int pred, NumericVector bioms){
    double tot = 0;
    int i;

    for (i=0; i<nb_s; i++){
      tot += h(i,pred)*a(i,pred) * pow_bioms[i];
    }
    return ((a(prey,pred)*pow_bioms(prey)) / 
            (1 + c(pred)*bioms(pred + nb_b) + tot));
  }
  
  // NumericVector ODE(double t, NumericVector bioms, NumericVector p){  // for sundials
  NumericVector ODE(NumericVector bioms, double t){ // for odeintr
    
    bioms[bioms < ext] = 0.0;
    pow_bioms = pow(bioms, q);
    
    for (res = all.begin(); res != all.end(); res++){
      for (cons = animals.begin(); cons != animals.end(); cons++){
        if ((a(*res, *cons - nb_b) > 0) & (bioms(*res) > 0.0) & (bioms(*cons)>0.0)){
          F(*res, *cons - nb_b) = F_rate(*res, *cons-nb_b, bioms);
        }
        else{
          F(*res, *cons-nb_b) = 0.0;
        }
      }
    }
    // Rcout << "F done " << std::endl; 

    
    // derivates for plants
    for (res = plants.begin(); res != plants.end(); res++){
      out = 0;
      for (cons = animals.begin(); cons != animals.end(); cons++){
        out += bioms[*cons] * F(*res, *cons - nb_b);
      }
      // plant resource competition 
      // s = 0;
      // for (i=0; i<nb_b; i++){
      //   s += alpha(*res, i)*bioms[i];
      // }

      // Rcout << 1-s/K[*res] << "  ";
      // s = bioms[*res];
      dB[*res] = r[*res]*bioms[*res]*(1-bioms[*res]/K[*res]) - out - X[*res]*bioms[*res];
    }
    // Rcout << " plants done " << std::endl;


    // derivative for animals
    for (cons = animals.begin(); cons != animals.end(); cons++){
      out = 0;
      // int in = 0;
      // in = sum(e * F(_,*cons - nb_b)) * bioms[*cons + nb_n];
      for (cons2 = animals.begin(); cons2 != animals.end(); cons2++){
        out += bioms[*cons2]*F(*cons,*cons2-nb_b);
      }

      dB[*cons] = sum(e * F(_,*cons - nb_b)) * bioms[*cons] - out - X[*cons]*bioms[*cons];
    }
    // Rcout << "db: " << dB << std::endl;
    dB[bioms < ext] = 0.0;
    
    return dB;
  }
  
};



RCPP_MODULE(Unscaled_loopsModule){
  using namespace Rcpp;
  class_<Unscaled_loops>("Unscaled_loops")
    .constructor<int, int>("constructor") //constructor
    .method("print", &Unscaled_loops::print)
    .method("ODE", &Unscaled_loops::ODE)
    .field("nb_s", &Unscaled_loops::nb_s)
    .field("nb_b", &Unscaled_loops::nb_b)
    .field("BM", &Unscaled_loops::BM)
    // .field("log_BM", &Unscaled_loops::log_BM)
    .field("K", &Unscaled_loops::K)
    .field("r", &Unscaled_loops::r)
    .field("X", &Unscaled_loops::X)
    .field("e", &Unscaled_loops::e)
    .field("a", &Unscaled_loops::a)
    .field("c", &Unscaled_loops::c)
    .field("h", &Unscaled_loops::h)
    .field("q", &Unscaled_loops::q)
    .field("dB", &Unscaled_loops::dB)
    .field("F", &Unscaled_loops::F)
    .field("fw", &Unscaled_loops::fw)
    .field("ext", &Unscaled_loops::ext)
    .field("alpha", &Unscaled_loops::alpha)
    ;  
}
