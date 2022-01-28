## R CMD check results

2 notes:

N  checking sizes of PDF files under 'inst/doc' (1.2s)
   Unable to find GhostScript executable to run checks on size reduction
N  checking installed files from 'inst/doc'
   The following files should probably not be installed:
     'model_descriptions.tex'

## r_hub

Installation of Rcpp and Rcpp Armadillo packages failed in that local environment:

ATNr 1.0: PREPERROR
Build ID:   ATNr_1.0.tar.gz-68f2a5c7816d4483b2fa1dad08b6a1b9
Platform:   Debian Linux, R-devel, GCC ASAN/UBSAN
Submitted:  19 minutes 4 seconds ago
Build time:   18 minutes 19.4 seconds

#> attributes.cpp: In function ‘bool Rcpp::attributes::checkRSignature(const Rcpp::attributes::Function&, std::string)’:
1781#> attributes.cpp:2813:45: error: ‘>>’ should be ‘> >’ within a nested template argument list
1782#> 2813 | Rcpp::as>(pargs_cv);
1783#> | ^~
1784#> | > >
1785#> make: *** [/usr/local/lib/R/etc/Makeconf:177: attributes.o] Error 1
1786#> ERROR: compilation failed for package ‘Rcpp’

#########################################################################################

Build ID:   ATNr_1.0.tar.gz-7362a6c0e8184568b9f666a80cbbf4ef
Platform:   Ubuntu Linux 20.04.1 LTS, R-release, GCC
Submitted:  14 minutes 26.9 seconds ago
Build time:   14 minutes 5 seconds

* checking installed package size ... NOTE
  installed size is 24.1Mb
  sub-directories of 1Mb or more:
    libs  23.3Mb

* checking installed files from ‘inst/doc’ ... NOTE
The following files should probably not be installed:
  ‘model_descriptions.tex’

Consider the use of a .Rinstignore file: see ‘Writing R Extensions’,
or move the vignette sources from ‘inst/doc’ to ‘vignettes’.

#########################################################################################

Build ID:   ATNr_1.0.tar.gz-88551b23f3cb4df892a568bd8d478c4a
Platform:   Windows Server 2022, R-devel, 64 bit
Submitted:  39 minutes 27.4 seconds ago
Build time:   39 minutes 15.2 seconds

* checking sizes of PDF files under 'inst/doc' ... NOTE
Unable to find GhostScript executable to run checks on size reduction

* checking installed files from 'inst/doc' ... NOTE
The following files should probably not be installed:
  'model_descriptions.tex'

or move the vignette sources from 'inst/doc' to 'vignettes'.
Consider the use of a .Rinstignore file: see 'Writing R Extensions',

#########################################################################################

Build ID:   ATNr_1.0.tar.gz-e500421680884ec0a4dccfbae6076a0a
Platform:   Fedora Linux, R-devel, clang, gfortran
Submitted:  15 minutes 15.1 seconds ago
Build time:   14 minutes 41.7 seconds

* checking re-building of vignette outputs ... WARNING
Error(s) in re-building vignettes:
--- re-building ‘ATNr.Rmd’ using rmarkdown
--- finished re-building ‘ATNr.Rmd’

--- re-building ‘model_descriptions.tex’ using tex
Error: processing vignette 'model_descriptions.tex' failed with diagnostics:
Running 'texi2dvi' on 'model_descriptions.tex' failed.
LaTeX errors:
! LaTeX Error: File `authblk.sty' not found.

Type X to quit or <RETURN> to proceed,
or enter new name. (Default extension: sty)

! Emergency stop.
<read *> 
         
l.13 \usepackage
                {hyperref}^^M
!  ==> Fatal error occurred, no output PDF file produced!
--- failed re-building ‘model_descriptions.tex’

SUMMARY: processing the following file failed:
  ‘model_descriptions.tex’

Error: Vignette re-building failed.
Execution halted