# 2  Download handbook and data

# reference: https://www.epirhandbook.com/en/new_pages/data_used.html#download-offline-handbook

#-------------------------------------------------------------------------------

# temporary settings adjustments, 0ne time
# should not be run next time

# Unset the GitHub PAT for this session
Sys.unsetenv("GITHUB_PAT")

# managing git credentials

# Install the remotes package if you don't have it
if(!require(remotes)) install.packages("remotes")

install.packages("gitcreds")
library(gitcreds)

# Adding your git credentials
gitcreds::gitcreds_set()

# deleting git credentials
gitcreds::gitcreds_delete()

#-------------------------------------------------------------------------------
# install the latest version of the Epi R Handbook package
# 1. Install the pak package
install.packages("pak")
library(pak)

# 2. Use pak to install the handbook
pak::pkg_install("appliedepi/epirhandbook")

# checking if epirhandbook package has been installed
# this can also be used for other packages
"epirhandbook" %in% installed.packages()

# loading the epirhandbook library
library(epirhandbook)

# download the offline handbook to your computer
download_book()



# download all the example data into a folder on your computer
get_data("all")






