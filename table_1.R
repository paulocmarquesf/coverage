# Universal coverage tolerance table for split conformal prediction

# Table 1 of https://arxiv.org/abs/2303.02770

begin <- Sys.time()

calibration_size <- function(epsilon, tau, alpha) {
    # minimum n_0 such that the pair (n_0, alpha) is feasible
    n_0 <- ceiling((1 - alpha) / alpha) # <=> ( ceiling((1 - alpha)*(n_0 + 1)) <= n_0 )
    repeat { # simple linear search; R is fast enough
        b <- ceiling((1 - alpha)*(n_0 + 1))
        g <- floor(alpha*(n_0 + 1))
        pr <- pbeta(1 - alpha + epsilon, b, g) - pbeta(1 - alpha - epsilon, b, g)
        if (pr >= tau) break
        n_0 <- n_0 + 1
    }
    cat(sprintf("epsilon = %.3f, tau = %.2f, 1 - alpha = %.2f, n_0 = %i\n", epsilon, tau, 1 - alpha, n_0))
}

for (epsilon in c(0.1, 0.05, 0.01, 0.005, 0.001)) {
    cat("\n")
    for (tau in c(0.9, 0.95, 0.99)) {
        cat("\n")
        for (alpha in c(0.2, 0.15, 0.10, 0.05)) {
            calibration_size(epsilon, tau, alpha)
        }
    }
}

print(Sys.time() - begin)
