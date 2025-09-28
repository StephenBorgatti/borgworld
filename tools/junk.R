devtools::document()
devtools::load_all()

devtools::install()

usethis::use_package_doc()

devtools::document()
devtools::build()   # creates ezborgstats_0.0.0.9000.tar.gz


# clean the masking object and unload the pkg
rm(ez_vars, envir = .GlobalEnv)
devtools::unload("ezborgstats")       # ignore message if not loaded

# rebuild docs and load from source
devtools::document()
devtools::load_all()                  # now ez_vars comes from the package
