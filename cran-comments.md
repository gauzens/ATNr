## R CMD check results

There were 2 warnings:

* checking for code/documentation mismatches ... WARNING
Variables with usage in documentation object â€˜schneiderâ€™ but not in code:
  â€˜schneiderâ€™

 ==> Solved

 * checking for missing documentation entries ... WARNING
Undocumented S4 classes:

==> Here, the S4 classes are C++ classes from Rcpp. I followed Rcpp recomendations to document classes and it is possible to acess to a documentation for these classes. Could that be a false positive? 