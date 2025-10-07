# 1. Load your changes immediately to test them
devtools::load_all()  # Ctrl+Shift+L

# 2. Test your function informally in console
#your_function(test_data)

# 3. If roxygen comments changed, regenerate docs
devtools::document()  # Ctrl+Shift+D

# 4. Iterate on 1-3 until happy, then...

# 5. Run checks (can be slow, so not every edit)
devtools::check()     # Ctrl+Shift+E

# 6. If releasing/milestone, update version in DESCRIPTION
usethis::use_version("minor")  # or "minor", "major", "dev"
usethis::use_version("patch")  # or "minor", "major", "dev"

# 7. Commit and push via Git tab

##occasionally
# Periodically, for full testing:
devtools::install()   # Install the package
# Restart R (Session > Restart R in RStudio)
library(borgworld)    # Test as users would use it
help(package = "borgworld")  # Check documentation
