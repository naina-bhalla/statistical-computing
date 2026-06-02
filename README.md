# MTH210: Statistical Computing

## Project Overview
This repository contains an R-based statistical computing pipeline designed to estimate the parameters ($\alpha$ and $\lambda$) of a Generalized Exponential distribution. The project implements a custom 2D **Newton-Raphson optimization algorithm** from scratch to find the Maximum Likelihood Estimates (MLEs). Furthermore, it evaluates the variance of these estimates by constructing 95% confidence intervals using both **parametric** and **non-parametric bootstrap** techniques.

---

## Algorithmic Features & Methodology

### 1. Newton-Raphson Optimizer
* **Analytical Derivatives:** The script computes the exact gradient vector and 2x2 Hessian matrix of the log-likelihood function.
* **Data-Driven Initialization:** To guarantee stable convergence, the algorithm initializes $\lambda_0$ as the inverse of the sample mean (mirroring a standard exponential distribution where $\alpha=1$). The initial $\alpha_0$ is then dynamically computed by setting its partial derivative to zero.
* **Strict Convergence Criteria:** The optimization loop terminates when the $L_\infty$ norm (maximum absolute difference) of the parameter update step falls below a strict $10^{-6}$ tolerance threshold.

### 2. Computational Safeguards
* **Domain Positivity:** The support of the Generalized Exponential distribution requires strictly positive parameters. The algorithm uses a bounding safeguard (`max(theta_new, 1e-6)`) to prevent parameter steps from jumping into negative space.
* **Singularity Protection:** During bootstrap resampling, heavily skewed datasets can cause the Hessian matrix to become singular. Optimization calls inside the loop are wrapped in `try()` blocks to catch errors silently, discard failed iterations, and prevent the 1,000-replication loop from crashing.

### 3. Bootstrap Confidence Intervals
* **Non-Parametric Bootstrap:** Resamples the original dataset with replacement to build an empirical distribution of the MLEs.
* **Parametric Bootstrap:** Utilizes an **Inverse Transform Sampling** function (`generate_from_f`) to simulate new datasets directly from the fitted Generalized Exponential probability density function.
* Both methods execute $B=1000$ replications to construct robust 95% confidence intervals via the 2.5th and 97.5th percentiles.

---

## How to Run

### Prerequisites
This script requires **Base R**. No external packages or dependencies are needed.

### Execution
1. Clone the repository and ensure your working directory contains the script.
2. The dataset (`fort.22`) is currently hardcoded into the `main.R` script as a numeric vector for ease of execution. If you wish to read directly from a file, uncomment the `x <- scan("fort.22")` line.
3. Run the script from the command line or RStudio:
   ```bash
   Rscript main.R
