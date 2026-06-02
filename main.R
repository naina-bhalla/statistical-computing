# MTH210 Course Project
# Name: Naina Bhalla
# Roll Number: 240674
# Data File: fort.22

# Setting seed for reproducibility of the results written in the report
set.seed(0)

# 1. Loading the Data

# Ensure the file 'fort.22' is in your working directory.
# x <- scan("fort.22")

# OR directly putting the data from the file
x <- c(1.080, 1.245, 0.111, 1.970, 2.137, 2.643, 
        1.440, 4.519, 1.914, 0.203, 3.369, 1.210, 
        0.162, 0.138, 3.060, 0.678, 1.879, 0.221, 
        1.702, 8.334, 0.878, 1.513, 1.290, 0.610, 
        2.343)
n <- length(x)

# 2. Calculations

# Log-Likelihood Gradient and Hessian for f(x)
grad <- function(alpha, lambda, x) {
  n <- length(x)
  d_alpha <- n/alpha + sum(log(1 - exp(-lambda * x)))
  d_lambda <- n/lambda - sum(x) +
    (alpha - 1) * sum((x * exp(-lambda * x)) / (1 - exp(-lambda * x)))       
  return(c(d_alpha, d_lambda))
}

hessian <- function(alpha, lambda, x) {
  n <- length(x)
  d2_alpha <- -n / (alpha^2)
  d2_lambda <- -n / (lambda^2) - (alpha - 1) * sum((x^2 * exp(-lambda*x)) / (1 - exp(-lambda*x))^2)
  d_alpha_lambda <- sum((x * exp(-lambda*x)) / (1 - exp(-lambda*x)))
  
  H <- matrix(c(d2_alpha, d_alpha_lambda,
                d_alpha_lambda, d2_lambda), 2, 2)
  return(H)
}

# Newton-Raphson function
newton_raphson <- function(x, tol=1e-6, max_iter=100, quiet=FALSE) {
  
  # Initial guesses
  n <- length(x)
  
  # Since the PDF is a generalized exponential, 
  # for exp(lambda) (with alpha = 1),
  # the MLE of lambda is 1/mean
  
  lambda_0 <- 1 / mean(x)
  
  # Equating the equation for alpha to 0 using the initial value of lambda
  # That will give us an initial value of alpha
  
  alpha_0 <- -n / sum(log(1 - exp(-lambda_0 * x))) 
  
  # putting both the initial values into a theta vector
  theta <- c(alpha_0, lambda_0)
  
  for (i in 1:max_iter) {
    g <- grad(theta[1], theta[2], x)
    H <- hessian(theta[1], theta[2], x)
    
    # Update step
    step <- solve(H, g)
    theta_new <- theta - step
    
    # Enforce positivity since both alpha>0 and lambda>0
    # according to the PDF and its parametric space
    
    theta_new[1] <- max(theta_new[1], 1e-6)
    theta_new[2] <- max(theta_new[2], 1e-6)
    
    # Convergence check using tol=1e-6
    if (max(abs(theta_new - theta)) < tol) {
      if(!quiet) cat("Converged in", i, "iterations\n")
      return(list(est=theta_new, iter=i))
    }
    
    theta <- theta_new
  }
  stop("Newton-Raphson did not converge")
}

# 3. MLE Estimation
mle_result <- newton_raphson(x)
alpha_hat <- mle_result$est[1]
lambda_hat <- mle_result$est[2]


# 4. Generation of x using Inverse Transform from the PDF
generate_from_f <- function(n, alpha, lambda) {
  u <- runif(n)
  return(-log(1 - u^(1/alpha)) / lambda)
}

# 5. Bootstrap Initialization
B <- 1000   # number of replications

boot_alpha_np <- numeric(B)
boot_lambda_np <- numeric(B)

boot_alpha_p <- numeric(B)
boot_lambda_p <- numeric(B)

# 6. Non-parametric Bootstrap

# Sampling from the initial set of values provided in the data

for (b in 1:B) {
  x_boot <- sample(x, replace=TRUE)
  
  # quiet=TRUE prevents the 1000 "Converged in..." lines from printing
  res <- try(newton_raphson(x_boot, quiet=TRUE), silent=TRUE)
  
  if (!inherits(res, "try-error")) {
    boot_alpha_np[b] <- res$est[1]
    boot_lambda_np[b] <- res$est[2]
  }
}

# Remove zeros 
boot_alpha_np <- boot_alpha_np[boot_alpha_np > 0]
boot_lambda_np <- boot_lambda_np[boot_lambda_np > 0]


# 7. Parametric Bootstrap
# Sampling from the original form of the PDF
for (b in 1:B) {
  x_boot <- generate_from_f(n, alpha_hat, lambda_hat)
  
  # quiet=TRUE prevents the 1000 "Converged in..." lines from printing
  res <- try(newton_raphson(x_boot, quiet=TRUE), silent=TRUE)
  
  if (!inherits(res, "try-error")) {
    boot_alpha_p[b] <- res$est[1]
    boot_lambda_p[b] <- res$est[2]
  }
}

boot_alpha_p <- boot_alpha_p[boot_alpha_p > 0]
boot_lambda_p <- boot_lambda_p[boot_lambda_p > 0]

# 8. FINAL RESULTS
cat("MLE Estimates:\n")
cat("alpha =", alpha_hat, "\n")
cat("lambda =", lambda_hat, "\n")

cat("\n95% Confidence Intervals\n")

cat("\nNon-Parametric Bootstrap:\n")
cat("alpha CI:", quantile(boot_alpha_np, c(0.025, 0.975)), "\n")
cat("lambda CI:", quantile(boot_lambda_np, c(0.025, 0.975)), "\n")

cat("\nParametric Bootstrap:\n")
cat("alpha CI:", quantile(boot_alpha_p, c(0.025, 0.975)), "\n")
cat("lambda CI:", quantile(boot_lambda_p, c(0.025, 0.975)), "\n")
