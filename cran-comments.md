Dear Julia Haider (or new reviewer),

thanks a lot for taking the time to review the package. You can find below the comments made on the last submission and how the issues were soved.

Best.

#Please omit the redundant "in R" from the title.
==> new title: run Allometric Trophic Networks models

# Please do not start the description with "This package", package name, title or similar.
==> new description: Implements the differential equations associated to different versions of Allometric Trophic Models (ATN) to estimate the temporal dynamics of species biomasses in food webs. It offers several features to generate synthetic food webs and to parametrise models as well as a wrapper to the ODE solver deSolve.

# Please add some more details about the package functionality and implemented methods in your Description text.
==> see above. I hope it is fine know. 


# If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

==> No references for now. 

# Please rather use the Authors@R field and declare Maintainer, Authors and Contributors with their appropriate roles with person() calls.
e.g. something like:
Authors@R: c(person("Alice", "Developer", role = c("aut", "cre","cph"),
                     email = "alice.developer@some.domain.net"),
              person("Bob", "Dev", role = "aut") )
              
==> changed authors to: 
authors@R: 
    c(person("Benoit", "Gauzens", role = c("cre", "aut"), email = "benoit.gauzens@gmail.com"), 
      person("Emilio", "Berti", role = "aut"))


# Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. \value{No return value, called for side effects} or similar)
Missing Rd-tags e.g.:
    create_matrix_parameter.Rd: \value
    plot_odeweb.Rd: \value
    run_checks.Rd: \value

==> added return value for `create_matrix_parameter`.
for plot_odeweb.Rd: added "#' @return No return value, called for side effects"
for run_checks.Rd: added "#' @return No return value, only throw an error if parameters are inconsistent. "

# \dontrun{} should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ('# Not run:') as a warning for the user.
Does not seem necessary.

==> No more \dontrun{} in the descriptions. all examples are executed. 

# Please unwrap the examples if they are executable in < 5 sec, or replace \dontrun{} with \donttest{}.

==> same as above, No more \dontrun{} in the descriptions.

# Please always make sure to reset to user's options(), working directory or par() after you changed it in examples and vignettes and demos.
e.g.:
oldpar <- par(mfrow = c(1,2))
...
par(oldpar)
e.g. inst/doc/ATNr.R

==> Changed as suggested

# You are using installed.packages() in your code. As mentioned in the notes of installed.packages() help page, this can be very slow. Therefore please do not use installed.packages(). 

==> No more use of `installed.packages()` in the code. All packages that were called by `installed.packages()` are now in declared in Imports

