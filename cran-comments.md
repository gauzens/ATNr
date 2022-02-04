# comment from last submission:

Thanks, we still see with valgrind enabled checks:


==3227138== Conditional jump or move depends on uninitialised value(s)
==3227138==    at 0x4BCD5B5: __ieee754_pow_sse2 (in /usr/lib64/libm-2.31.so)
==3227138==    by 0x4BE17D7: pow@@GLIBC_2.29 (in /usr/lib64/libm-2.31.so)
==3227138==    by 0x17C6DFDE: operator[] (R-devel/site-library/Rcpp/include/Rcpp/sugar/functions/pow.h:36)
==3227138==    by 0x17C6DFDE: void Rcpp::Vector<14, Rcpp::PreserveStorage>::import_expression<Rcpp::sugar::Pow<14, true, Rcpp::Vector<14, Rcpp::PreserveStorage>, double> >(Rcpp::sugar::Pow<14, true, Rcpp::Vector<14, Rcpp::PreserveStorage>, double> const&, long) (R-devel/site-library/Rcpp/include/Rcpp/vector/Vector.h:1085)
==3227138==    by 0x17C6E496: import_sugar_expression<true, Rcpp::sugar::Pow<14, true, Rcpp::Vector<14, Rcpp::PreserveStorage>, double> > (R-devel/site-library/Rcpp/include/Rcpp/vector/Vector.h:1071)
==3227138==    by 0x17C6E496: Vector<true, Rcpp::sugar::Pow<14, true, Rcpp::Vector<14, Rcpp::PreserveStorage>, double> > (R-devel/site-library/Rcpp/include/Rcpp/vector/Vector.h:165)
==3227138==    by 0x17C6E496: Scaled_loops::Scaled_loops(int, int) (test/ATNr.Rcheck/00_pkg_src/ATNr/src/Scaled_loops.cpp:89)
...
==3227138==  Uninitialised value was created by a heap allocation
==3227138==    at 0x483AE7D: operator new(unsigned long) (builddir/build/BUILD
valgrind-3.16.1/coregrind/m_replacemalloc/vg_replace_malloc.c:342)
==3227138==    by 0x17C713F5: get_new (R-devel/site-library/Rcpp/include/Rcpp/module/Module_generated_Constructor.h:57)
==3227138==    by 0x17C713F5: Rcpp::class_<Scaled_loops>::newInstance(SEXPREC**,
 int) (R-devel/site-library/Rcpp/include/Rcpp/module/class.h:131)
==3227138==    by 0x17BEF969: class__newInstance(SEXPREC*) (/tmp/RtmprdPLXC/R.INSTALL21fb3a55b5e55/Rcpp/src/module.cpp:143)

and lots more like that.

Please fix and resubmit. 



# response:

Thanks for raising this issues, and my apologise for the time you invest on my package. 
I check the code and ran valgrim on the different files to correct issues. 

The final check was made by creating a file (memcheck.R) at the root of the package containing the code of the examples, a source() to the R file associated to the vignette in /doc and a source() on all test files. 

results was :

R -d "valgrind -s --tool=memcheck --leak-check=full --track-origins=yes" --vanilla < memcheck.R 

==12545== Memcheck, a memory error detector
==12545== Copyright (C) 2002-2017, and GNU GPL'd, by Julian Seward et al.
==12545== Using Valgrind-3.15.0 and LibVEX; rerun with -h for copyright info
==12545== Command: /usr/lib/R/bin/exec/R --vanilla
==12545== 

...

==23967== 
==23967== HEAP SUMMARY:
==23967==     in use at exit: 98,642,028 bytes in 22,963 blocks
==23967==   total heap usage: 2,031,079 allocs, 2,008,116 frees, 2,245,895,650 bytes allocated
==23967== 
==23967== LEAK SUMMARY:
==23967==    definitely lost: 0 bytes in 0 blocks
==23967==    indirectly lost: 0 bytes in 0 blocks
==23967==      possibly lost: 0 bytes in 0 blocks
==23967==    still reachable: 98,642,028 bytes in 22,963 blocks
==23967==                       of which reachable via heuristic:
==23967==                         newarray           : 4,264 bytes in 1 blocks
==23967==         suppressed: 0 bytes in 0 blocks
==23967== Reachable blocks (those to which a pointer was found) are not shown.
==23967== To see them, rerun with: --leak-check=full --show-leak-kinds=all
==23967== 
==23967== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)


I hope this check is enough to guaranty that all problems that can be detected by valgrind are properly solved. 