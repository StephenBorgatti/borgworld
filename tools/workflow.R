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
usethis::use_version("patch", push = TRUE)


# 7. Commit and push via Git tab

##occasionally
# Periodically, for full testing:
devtools::install()   # Install the package
# or from different folder
devtools::install("C:/Users/sborg2/R/borgworld")
# Restart R (Session > Restart R in RStudio)
library(borgworld)    # Test as users would use it
help(package = "borgworld")  # Check documentation


#when rdb is corrupt

# Clean the package completely
devtools::clean_dll()  # Remove compiled code
devtools::unload()     # Unload if loaded

# Delete the corrupt files manually
unlink("C:/Users/sborg2/AppData/Local/R/win-library/4.5/borgworld", recursive = TRUE)

# Rebuild
devtools::document()   # Regenerate documentation
devtools::install()    # Reinstall the package
devtools::load_all()   # Load it again

#when claude code pushes something to github

# fetch changes
system("git fetch origin")
#merge changes from git branch learned from system()
system("git merge origin/claude/document-bwordcount-function-01Jr8e32wzpVEVBg5aAG14Sy")

#push back after edits
system("git push origin main")

# after doing the pull in github
remotes::install_github("StephenBorgatti/borgworld")

