remove.packages("rstan")
if (file.exists(".RData")) file.remove(".RData")
install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
pkgbuild::has_build_tools(debug = TRUE)
dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, ifelse(.Platform$OS.type == "windows", "Makevars.win", "Makevars"))
if (!file.exists(M)) file.create(M)
cat("\nCXX14FLAGS=-O3 -march=native -mtune=native",
    if( grepl("^darwin", R.version$os)) "CXX14FLAGS += -arch x86_64 -ftemplate-depth-256" else 
      if (.Platform$OS.type == "windows") "CXX11FLAGS=-O3 -march=native -mtune=native" else
        "CXX14FLAGS += -fPIC",
    file = M, sep = "\n", append = TRUE)
library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
model_string<-"
// saved as 8schools.stan
data {
int<lower=0> J;         // number of schools 
real y[J];              // estimated treatment effects
real<lower=0> sigma[J]; // standard error of effect estimates 
}
parameters {
real mu;                // population treatment effect
real<lower=0> tau;      // standard deviation in treatment effects
vector[J] eta;          // unscaled deviation from mu by school
}
transformed parameters {
vector[J] theta = mu + tau * eta;        // school treatment effects
}
model {
target += normal_lpdf(eta | 0, 1);       // prior log-density
target += normal_lpdf(y | theta, sigma); // log-likelihood
}
"
schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
fit <- stan(model_code=model_string, data = schools_dat)
print(fit)
plot(fit)
pairs(fit, pars = c("mu", "tau", "lp__"))

la <- extract(fit, permuted = TRUE) # return a list of arrays 
mu <- la$mu 

### return an array of three dimensions: iterations, chains, parameters 
a <- extract(fit, permuted = FALSE) 

### use S3 functions on stanfit objects
a2 <- as.array(fit)
m <- as.matrix(fit)
d <- as.data.frame(fit)
