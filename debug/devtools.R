WD <- '~/Desktop/'
setwd(WD)

if(F){
	install.packages("formatR")	
}

devtools::create("BayesMPSimulation") 
WD2 <- '~/Desktop/BayesMPSimulation'

## licenses

## copy your code into /R folder


setwd(WD2)
## make the code neat
formatR::tidy_dir("R")


## dependency


## check the package
devtools::check()

## build the package
devtools::build()


devtools::install()


