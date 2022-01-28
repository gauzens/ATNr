#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// generate documentation with:
// Rcpp::compileAttributes()           # this updates the Rcpp layer from C++ to R
// roxygen2::roxygenize(roclets="rd")  # this updates the documentation based on roxygen comments

//' @name Unscaled_nuts_loops
//' @title Store parameters and functions associated to the unscaled version of ATN
//' @description To not use. For testing purpose only. Please use Rcpp_Unscaled_nuts instead. 


class Unscaled_nuts_loops{
public:
  int nb_s; // number of species
  int nb_b; // number of basal species
  int nb_n = 2; // number of nutrients
  int n_tot = nb_s + nb_n; // bool prefs MORE PRECISE?
  double temperature = 20;
  // global nutrient turn over rate (rate of replenishment)
  // used in calculating change in nutrient concentration
  double D = 0.25;

  // half saturation density of nutrient, or nutrient uptake efficiency
  double q;
  double test_double;
  double ext;


  NumericVector X = NumericVector(nb_s); // metabolic rates
  NumericVector e = NumericVector(nb_s); // assimilation efficiencies
  NumericVector r = NumericVector(nb_b); // growth rates of plants
  NumericVector S = NumericVector(nb_n); // maximal nutrient level
  // interference competition
  NumericVector c = NumericVector(nb_s - nb_b); 
  
  // body masses
  NumericVector BM = NumericVector(nb_s);
  NumericVector log_BM = NumericVector(nb_s);
  
  // biomasses
  NumericVector bioms = NumericVector(nb_s + nb_n);
  // r*G for plants, as there is no need to compute it at each ODE call
  
  NumericVector test;
  NumericVector p = NumericVector(2);
  // vector of derivatives
  NumericVector dB = NumericVector(nb_n + nb_s);
  // 
  
  LogicalMatrix fw = LogicalMatrix(); 
  
  NumericMatrix b = NumericMatrix(nb_s,nb_s - nb_b);
  // handling times
  NumericMatrix h = NumericMatrix(nb_s,nb_s - nb_b);
  // functional response
  NumericMatrix F = NumericMatrix(nb_s,nb_s - nb_b);
  // consumption rates
  NumericMatrix w = NumericMatrix(nb_s,nb_s - nb_b);

  // relative content in the plant species' biomass
  NumericMatrix V = NumericMatrix(nb_n, nb_b); 
  // plants nutrient uptake efficiency (K(i,j): plant i on nutrient j)
  NumericMatrix K = NumericMatrix(nb_n, nb_b); 

  // internal variables for optimisation
  // index of plants (optimisation purpose)

  IntegerVector plants = Range(0, nb_b - 1);
  IntegerVector animals = Range(nb_b, nb_s - 1);
  IntegerVector plants_bioms = Range(nb_n, nb_b - 1 + nb_n);
  IntegerVector animals_bioms = Range( nb_b + nb_n, nb_s - 1 + nb_n);
  IntegerVector non_nut = Range(0, nb_s - 1);
  IntegerVector nuts = Range(0, nb_n - 1);
  NumericVector G = NumericVector(nb_b);
  NumericVector g_temp = NumericVector(nb_b);
  NumericVector zeros = NumericVector(nb_s); // accordingly to documentation that should work
  NumericVector pow_bioms = NumericVector(nb_s);
  NumericVector temp_nut = NumericVector(nb_n);
  NumericVector temp_b = NumericVector(nb_n);
  
  // LogicalVector prey = fw[_,1] = 1;
  IntegerVector::iterator cons;
  IntegerVector::iterator cons2;
  IntegerVector::iterator res;
  IntegerVector::iterator nut; //not properly needed, but more readable
  NumericVector uptake = NumericVector(nb_b);
  double out = 0;
  int i = 0;
  
  Unscaled_nuts_loops(int ns, int nb, int nn):
    nb_s(ns), nb_b(nb), nb_n(nn) {}
  
  // Rcpp::XPtr<parameters_prefs> make_model(int s, int b, int n){
  //   parameters_prefs* m = new parameters_prefs(s, b, n);
  //   Rcpp::XPtr<parameters_prefs> ptr(m);
  //   return(ptr);
  // }
  
  void print(){
    Rcout << "nb_s:"  << std::endl << nb_s << std::endl;
    Rcout << "nb_b:"  << std::endl << nb_b << std::endl;
    Rcout << "plants: " << plants << std::endl; 
    Rcout << "bioms: " << bioms << std::endl; 
    Rcout << "bioms plants: " << bioms[plants] << std::endl; 
    Rcout << "G: " << G << std::endl; 
    Rcout << "Gplant: " << G[plants] << std::endl;
    Rcout << "dbplant " << dB[plants] << std::endl;
    Rcout << "r[plants]" << r[plants] << std::endl;
    // Rcout << " prey" << prey << std::endl;
    test_double = pow(0,1.5);
    test = F(_, 0);
    Rcout << "test " << test << std::endl;
  }
  
  
  
  
  double F_rate(int prey, int pred, NumericVector bioms){
    double tot = 0;
    int i;

    for (i=0; i<nb_s; i++){
      tot += w(i, pred)*h(i,pred)*b(i,pred) * pow_bioms[i+nb_n];
    }
    return ((w(prey, pred)*b(prey,pred)*pow_bioms(prey + nb_n)) / 
            ((1 + c(pred)*bioms(pred + nb_n + nb_b) + tot)*BM(pred+nb_b)));
  }
  
  // NumericVector ODE(double t, NumericVector bioms, NumericVector p){  // for sundials
  NumericVector ODE(NumericVector bioms, double t){ // for odeintr
    
    bioms[bioms < ext] = 0.0;
    pow_bioms = pow(bioms, q);
    
    for (res = non_nut.begin(); res != non_nut.end(); res++){
      for (cons = animals.begin(); cons != animals.end(); cons++){
        if ((b(*res, *cons-nb_b) > 0) & (bioms(*res + nb_n) > 0.0) & (bioms(*cons + nb_n)>0.0)){
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
        out += bioms[*cons + nb_n] * F(*res, *cons - nb_b);
      }

      for (nut=nuts.begin(); nut!=nuts.end(); nut++){
        temp_nut[*nut] = bioms[*nut] / (K(*nut,*res)+ bioms[*nut]);
      }

      G[*res] = min(temp_nut);
      uptake[*res] = r[*res]*G[*res]*bioms[*res + nb_n];
      dB[*res + nb_n] = uptake[*res] - out - X[*res]*bioms[*res + nb_n];
    }
    // Rcout << " plants done " << std::endl;
    // derivative for animals
    for (cons = animals.begin(); cons != animals.end(); cons++){
      out = 0;
      // int in = 0;
      // in = sum(e * F(_,*cons - nb_b)) * bioms[*cons + nb_n];
      for (cons2 = animals.begin(); cons2 != animals.end(); cons2++){
        out += bioms[*cons2 + nb_n]*F(*cons,*cons2-nb_b);
      }
      dB[*cons + nb_n] = sum(e * F(_,*cons - nb_b)) * bioms[*cons + nb_n] - out - X[*cons]*bioms[*cons + nb_n];
    }
    // Rcout << "db: " << dB << std::endl;
    dB[bioms < ext] = 0.0;

    // derivate for nutrients
    for (nut=nuts.begin(); nut!=nuts.end(); nut++){
      out = 0;
      for (res = plants.begin(); res != plants.end(); res++){
        out += V(*nut,*res)*uptake[*res];
      }
      dB[*nut] = D * (S[*nut] - bioms[*nut]) - out;  
    }
    
    
    return dB;
  }
  
};



RCPP_MODULE(Unscaled_nuts_loopsModule){
  using namespace Rcpp;
  class_<Unscaled_nuts_loops>("Unscaled_nuts_loops")
    .constructor<int, int, int>("constructor") //constructor
    .method("print", &Unscaled_nuts_loops::print)
    .method("ODE", &Unscaled_nuts_loops::ODE)
    .field("nb_s", &Unscaled_nuts_loops::nb_s)
    .field("nb_b", &Unscaled_nuts_loops::nb_b)
    .field("nb_n", &Unscaled_nuts_loops::nb_n)
    .field("BM", &Unscaled_nuts_loops::BM)
    .field("log_BM", &Unscaled_nuts_loops::log_BM)
    .field("K", &Unscaled_nuts_loops::K)
    .field("D", &Unscaled_nuts_loops::D)
    .field("S", &Unscaled_nuts_loops::S)
    .field("r", &Unscaled_nuts_loops::r)
    .field("X", &Unscaled_nuts_loops::X)
    .field("e", &Unscaled_nuts_loops::e)
    .field("w", &Unscaled_nuts_loops::w)
    .field("b", &Unscaled_nuts_loops::b)
    .field("c", &Unscaled_nuts_loops::c)
    .field("h", &Unscaled_nuts_loops::h)
    .field("q", &Unscaled_nuts_loops::q)
    .field("V", &Unscaled_nuts_loops::V)
    // .field("bioms", &Unscaled_nuts_loops::bioms)
    .field("dB", &Unscaled_nuts_loops::dB)
    .field("D", &Unscaled_nuts_loops::D)
    .field("F", &Unscaled_nuts_loops::F)
    .field("uptake", &Unscaled_nuts_loops::uptake)
    .field("fw", &Unscaled_nuts_loops::fw)
    .field("ext", &Unscaled_nuts_loops::ext)
    ;  
    ;  
}
