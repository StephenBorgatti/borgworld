Here's how to get your GitHub Personal Access Token (PAT):
Step 1: Create the token on GitHub
Run this in R - it will open the correct GitHub page:
r
usethis::create_github_token()

OR manually go to: https://github.com/settings/tokens/new
Step 2: On the GitHub page that opens:

Note/Name: Give it a name like "R Package Development" or "borgworld"
Expiration: Choose how long the token lasts (30 days, 90 days, or custom)
Select scopes - Check these boxes at minimum:

✅ repo (this checks all repo sub-boxes)
✅ workflow (if you plan to use GitHub Actions)
✅ user (for your user info)


Scroll to bottom and click "Generate token" (green button)

Step 3: CRITICAL - Copy the token immediately!
⚠️ The token will look like: ghp_xxxxxxxxxxxxxxxxxxxx
Copy it NOW - GitHub will only show it once! If you lose it, you'll need to create a new one.
Step 4: Store it in R
r# Run this and paste your token when prompted
gitcreds::gitcreds_set()

# It will show:
# ? Enter password or token:
# Paste your token here (it will be hidden) and press Enter
Step 5: Verify it worked
r# Check that your token is stored
gitcreds::gitcreds_get()
# Should show your GitHub username and a masked token

# Now try using GitHub
usethis::use_github()  # This should work now!
Tips:

Save your token somewhere secure (like a password manager) as backup
The token is like a password - don't share it or commit it to code
If you get authentication errors later, create a new token
On Mac/Linux, the token is stored in your system keychain
On Windows, it's stored in the Windows Credential Manager

Once this is done, you won't need to enter the token again for future GitHub operations (until it expires)!RetryClaude can make mistakes. Please double-check responses.Research Opus 4.1
