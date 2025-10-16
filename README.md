# BRRR
Brain Regions Results Repository


Open R and run the following code:

# You will first need to install these packages. This should take about 10 minutes (but you'll never have to do it again). 
# Once you've run this you don't have to do it again, only the second part. 
pkgs <- c("shiny","ape","ratematrix","ggplot2","dplyr","tidyr","purrr","DT","colourpicker","readr", "devtools")
install.packages(setdiff(pkgs, rownames(installed.packages())), repos = "https://cloud.r-project.org")
devtools::install_github("joannabaker/BayesTraitR")

# Then, open the app. This will take about 5-8 minutes. Don't accidentally close it...
library(shiny)
runGitHub("BRRR", "joannabaker")
