library(plspm)

# Benchmark the performance of plspm and matrixpls with 1000 bootstrap samples using the customer satisfaction example

# load dataset satisfaction
data(satisfaction)

# inner model matrix
IMAG = c(0,0,0,0,0,0)
EXPE = c(1,0,0,0,0,0)
QUAL = c(0,1,0,0,0,0)
VAL = c(0,1,1,0,0,0)
SAT = c(1,1,1,1,0,0)
LOY = c(1,0,0,0,1,0)
sat_inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
colnames(sat_inner) <- rownames(sat_inner)

# outer model list
sat_outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
# vector of modes (reflective indicators)
sat_mod = rep(c("A"), 6)

# Set up the native matrixpls model
reflective <- matrix(0, max(unlist(sat_outer)),nrow(sat_inner))
rownames(reflective) <- colnames(satisfaction)[1:nrow(reflective)]
colnames(reflective) <- colnames(sat_inner)

col <- 1		

for(outerSpec in sat_outer){
	reflective[outerSpec,col] <- 1
	col <- col + 1 
}

formative <- t(reflective)
formative[] <- 0

nativeModel <- list(inner=sat_inner, 
										reflective=reflective, 
										formative=formative) 

weightRelations <- t(reflective)
	
#
# Single thread comparisons
#

# Store the old option values

opt_boot.ncpus <- getOption("boot.ncpus")
opt_boot.parallel <- getOption("boot.parallel")
opt_cores <- getOption("mc.cores")

options(boot.ncpus = 1)
options(boot.parallel = "no")
options(mc.cores = 1)

plspm_boot <- system.time(plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val = TRUE, br=1000))
plspm_no_boot <- system.time(replicate(100,plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE)))

matrixpls.plspm_boot <- system.time(matrixpls.plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val = TRUE, br=1000))
matrixpls.plspm_no_boot <- system.time(replicate(100,matrixpls.plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE)))

matrixpls_boot <- system.time(matrixpls.boot(satisfaction[,1:27], R = 1000, model = nativeModel, weightRelations = weightRelations,
																				innerEstimator = matrixpls.inner.centroid,
																				tol = 0.00001, iter = 100, validateInput = FALSE))
S <- cov(satisfaction[,1:27])

matrixpls_no_boot <- system.time(replicate(100,matrixpls(S, model = nativeModel, weightRelations = weightRelations,
																	 						innerEstimator = matrixpls.inner.centroid,
																	 						tol = 0.00001, iter = 100, validateInput = FALSE)))
#
# Multithreaded bootstrapping
#

cores = detectCores()

options(boot.ncpus = cores)

if(.Platform$OS.type == "windows") {
	options(boot.parallel = "snow")
} else {
	options(boot.parallel = "multicore")
}
options(mc.cores = cores)

matrixpls.plspm_boot_multicore <- system.time(matrixpls.plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val = TRUE, br=1000))

matrixpls_boot_multicore <- system.time(matrixpls.boot(satisfaction[,1:27], R = 1000, model = nativeModel, weightRelations = weightRelations,
																 																						 innerEstimator = matrixpls.inner.centroid,
																 																						 tol = 0.00001, iter = 100, validateInput = FALSE))
																 
#
# Set the options back to default values
#


options(boot.ncpus = opt_boot.ncpus)
options(boot.parallel = opt_boot.parallel)
options(mc.cores = opt_cores)

# Print the times

rbind(plspm_no_boot,matrixpls.plspm_no_boot, matrixpls_no_boot)

rbind(plspm_boot,matrixpls.plspm_boot, matrixpls_boot,matrixpls.plspm_boot_multicore, matrixpls_boot_multicore)
																 