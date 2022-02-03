#include <Rcpp.h>
// #include <RcppArmadillo.h>
using namespace Rcpp;

// generate documentation with:
// Rcpp::compileAttributes()           # this updates the Rcpp layer from C++ to R
// roxygen2::roxygenize(roclets="rd")  # this updates the documentation based on roxygen comments

//' @name Scaled_loops
//' @title Store parameters and functions associated to the scaled version of ATN
//' @description To not use. For testing purpose only. Please use Rcpp_Scaled instead. 


class Scaled_loops{
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
  // net growth rate of plant
  double G;
  // plant competition
  double s;
  
  // metabolic rates
  NumericVector X;
  // maximum feeding rate. Could be a vetor of length = nb of animals (maybe more intuitive for users?)
  NumericVector max_feed;

  // assimilation efficiencies
  NumericVector e;

  // mass specific growth rates of plants
  NumericVector r;

  // interference competition
  NumericVector c;  
  // body masses:
  NumericVector BM;
  // NumericVector log_BM;
    
  // vector of derivatives
  NumericVector dB;
  // 
    
  LogicalMatrix fw;
  
  NumericVector B0;

  // plant competition
  NumericMatrix alpha;

  // functional response
  NumericMatrix F;
  // consumption rates
  NumericMatrix w;
  
  // internal variables for optimisation
  // index of plants (optimisation purpose). Is it truly needed?
  IntegerVector plants;
  IntegerVector animals;
  IntegerVector all_sp;

  NumericVector pow_bioms;
  NumericVector pow_B0;
  
  // NumericVector zeros = NumericVector(nb_s); // accordingly to documentation that should work
  
  // LogicalVector prey = fw[_,1] = 1;
  IntegerVector::iterator cons;
  IntegerVector::iterator cons2;
  IntegerVector::iterator res;

  double out;
  int i;

  // temporary:
  NumericVector low;
  NumericVector tot2;
  NumericVector cs;
  
  Scaled_loops(int s, int b):
    nb_s(s), nb_b(b) {

       // metabolic rates
  X = NumericVector(nb_s);
  // maximum feeding rate. Could be a vetor of length = nb of animals (maybe more intuitive for users?)
  max_feed = NumericVector(nb_s-nb_b);

  // assimilation efficiencies
  e = NumericVector(nb_s);

  // mass specific growth rates of plants
  r = NumericVector(nb_b);

  // interference competition
  c = NumericVector(nb_s - nb_b);  
  // body masses:
  BM = NumericVector(nb_s);
  // log_BM = NumericVector(nb_s);
    
  // vector of derivatives
  dB = NumericVector(nb_s);
  // 
    
  fw = LogicalMatrix(nb_s,nb_s); 
  
  B0 = NumericVector(nb_s-nb_b);

  // plant competition
  alpha = NumericMatrix(nb_b, nb_b);

  // functional response
  F = NumericMatrix(nb_s,nb_s-nb_b);
  // consumption rates
  w = NumericMatrix(nb_s,nb_s-nb_b);
  
  // internal variables for optimisation
  // index of plants (optimisation purpose). Is it truly needed?
  plants = Range(0, nb_b - 1);
  animals = Range(nb_b, nb_s - 1);
  all_sp = Range(0, nb_s - 1);

  pow_bioms = NumericVector(nb_s);
  

  low = NumericVector(nb_s-nb_b);
  tot2 = NumericVector(nb_s-nb_b);
  cs = NumericVector(nb_s-nb_b);

  K = 0.0;
  out = 0.0;
  i = 0;
  ext = 0.0;
  q = 0.0;
  s = 0.0;
  
  // NumericVector zeros = NumericVector(nb_s); // accordingly to documentation that should work
  
  // LogicalVector prey = fw[_,1] = 1;
  // IntegerVector::iterator cons;
  // IntegerVector::iterator cons2;
  // IntegerVector::iterator res;
    }
  
  // Rcpp::XPtr<parameters_prefs> make_model(int s, int b, int n){
  //   parameters_prefs* m = new parameters_prefs(s, b, n);
  //   Rcpp::XPtr<parameters_prefs> ptr(m);
  //   return(ptr);
  // }

  void print(){
    Rcout << "pow_bioms: "  << pow_bioms << std::endl;
    // Rcout << "out_fluxes: "  <<  out << std::endl;
    Rcout << "F: "  <<  F << std::endl;
  }
  
    

  double F_rate(int prey, int pred, NumericVector bioms){
    double tot = 0.0;
    double res;

    // is the result of w(_,pred) * pow_bioms reallocated at every call ad freed each time?
    // if yes, then I should make the intermediate vector a class attribute maybe?
    
    tot = sum(w(_,pred) * pow_bioms);
    res = (w(prey, pred) * pow_bioms[prey]) / 
           (pow_B0[pred] + c[pred]*bioms[pred+nb_b] + tot);
    // low[pred] = pow_B0[pred] + c[pred]*bioms[pred+nb_b] + tot;
    // tot2[pred] = tot;
    // cs[pred] = c*bioms[pred+nb_b];
    // Rcout << "    pow_bioms prey: " <<  pow_bioms[prey] << "   res: " << res << std::endl;;
    // res = 0;
    return (res);
  }



  // NumericVector ODE(double t, NumericVector bioms, NumericVector p){  // for sundials
  NumericVector ODE(NumericVector bioms, double t){ // for odeintr

    checkUserInterrupt();
    bioms[bioms < ext] = 0.0;
    pow_bioms = pow(bioms, q);

    // note: no need to be recalculated at each call of ODE
    pow_B0 = pow(B0, q);

    // calculate the matrix of feeding rates
    // the if statement might prevent all loop optimisations
    // better to calculate it for everything?
    // calculate everything before hand helps to vectorise some calculations later on
    for (res = all_sp.begin(); res != all_sp.end(); res++){
      for (cons = animals.begin(); cons != animals.end(); cons++){
        // Rcout << "res: " << *res << "  cons: " << *cons << std::endl ;
        if ((fw(*res, *cons) > 0) & (bioms(*res) > 0.0) & (bioms(*cons)>0.0)){
            F(*res, *cons-nb_b) = F_rate(*res, *cons-nb_b , bioms);
        }
        else{
          F(*res, *cons-nb_b) = 0.0;
        }
      }
    }
    // Rcout << "low: " << low << std::endl;
    // Rcpp::Rcout << "tot2 " << tot2 << std::endl;
    // Rcpp::Rcout << "c " << cs << std::endl;
    // Rcpp::Rcout << "B0 " << pow_B0 << std::endl;
    // derivates for plants

    // G = 1 - bioms[plants]*(double(nb_b)/K); multiplying a sclice of ea vector by something is not permitted? 
    for (res = plants.begin(); res != plants.end(); res++){
      // first, sum up what is eaten by consumers
      out = 0;
      for (cons = animals.begin(); cons != animals.end(); cons++){
        out += X[*cons]*max_feed[*cons-nb_b]*bioms[*cons] * F(*res, *cons - nb_b);
      }
            // plant resource competition 
      s = 0;
      for (i=0; i<nb_b; i++){
        s += alpha(*res, i)*bioms[i];
      }
      G = 1 - s/K;
      // dB[*res] = r[*res]*G[*res]*bioms[*res] - out - X[*res]*bioms[*res];
      dB[*res] = r[*res]*G*bioms[*res] - out - X[*res]*bioms[*res];

    }

    // derivative for animals
    for (cons = animals.begin(); cons != animals.end(); cons++){
      // first, sum up what is eaten by consumers (cons2 here)
      out = 0;
      for (cons2 = animals.begin(); cons2 != animals.end(); cons2++){
        // Rcout << "\tconsumer2: " << *cons2 << std::endl ;
        out += X[*cons2]*max_feed[*cons2-nb_b]*bioms[*cons2]*F(*cons,*cons2-nb_b);
      }
      // Rcout << out << " ";
      dB[*cons] = X[*cons]*max_feed[*cons-nb_b]*sum(e * F(_,*cons-nb_b)) * bioms[*cons] - out - X[*cons]*bioms[*cons];
    }
    // Rcout << "db: " << dB << std::endl;
    dB[bioms < ext] = 0.0; // should be useless
    // Rcout << "animals done " << std::endl;
    // derivate for nutrients

    return dB;
  }

};



RCPP_MODULE(Scaled_loopsModule){
using namespace Rcpp;
  class_<Scaled_loops>("Scaled_loops")
  .constructor<int, int>("constructor")
  .method("print", &Scaled_loops::print)
  .method("ODE", &Scaled_loops::ODE)
  .field("nb_s", &Scaled_loops::nb_s)
  .field("nb_b", &Scaled_loops::nb_b)
  .field("BM", &Scaled_loops::BM)
  // .field("log_BM", &Scaled_loops::log_BM)
  .field("r", &Scaled_loops::r)
  .field("X", &Scaled_loops::X)
  .field("e", &Scaled_loops::e)
  .field("w", &Scaled_loops::w)
  .field("B0", &Scaled_loops::B0)
  .field("c", &Scaled_loops::c)
  .field("q", &Scaled_loops::q)
  .field("dB", &Scaled_loops::dB)
  .field("F", &Scaled_loops::F)
  .field("fw", &Scaled_loops::fw)
  .field("max_feed", &Scaled_loops::max_feed)
  .field("K", &Scaled_loops::K)
  .field("ext", &Scaled_loops::ext)
  .field("alpha", &Scaled_loops::alpha)
  ;  
}
