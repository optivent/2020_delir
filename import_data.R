#### prep steps ####
library(usethis)
?use_github
edit_r_environ()
use_github(protocol = 'https', auth_token = Sys.getenv("GITHUB_PAT"))

#### import the data ####

