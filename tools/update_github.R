#update github project


# WORKFLOW FOR MAKING CHANGES TO BORGWORLD

# 1. Make changes to your package (like adding bpca)

# 2. Document and test
devtools::document()
devtools::check()  # Or devtools::load_all() for quick testing

# 3. Commit your changes FIRST
gert::git_add(".")
gert::git_commit("Add bpca function for social science PCA")

# 4. THEN bump version (creates its own commit)
usethis::use_version()  # Choose "minor" for new functions

# 5. Push to GitHub (pushes both commits)
gert::git_push()

# 6. Now others can install the new version
# Others would run: devtools::install_github("StephenBorgatti/borgworld")


#------------------------------------------------------
# Or use the RStudio Git pane (easier!)
# - Check boxes next to changed files to stage
# - Click "Commit"
# - Write message and commit
# - Click "Push" (green arrow)
