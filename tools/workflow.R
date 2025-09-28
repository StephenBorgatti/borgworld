# 1. Document the package (generates man files from roxygen)
devtools::document()

# 2. Check if everything is working (optional but recommended)
devtools::check()  # Full check - takes longer
# OR for a quicker check:
#devtools::test()  # If you have tests
devtools::load_all()  # Just load to test functions work

# 3. Install the package locally
devtools::install()

# 4. Build the package (if you want to share it)
devtools::build()  # Creates a .tar.gz file
devtools::build(binary = TRUE)  # Creates a binary package (optional)

