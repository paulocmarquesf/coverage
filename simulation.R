# Universal distribution of the empirical coverage in split conformal prediction
#   Paulo C. Marques F. <PauloCMF1@insper.edu.br>
#   https://arxiv.org/abs/2303.02770

begin <- Sys.time()

library(ranger)

friedman <- function(n, p = 10) {
    w <- rexp(1)
    X <- matrix(runif(n*p), nrow = n, ncol = p, dimnames = list(1:n, paste0("x_", 1:p)))
    y <- 10*sin(pi*X[, 1]*X[, 2]) + 20*(X[, 3] - 0.5)^2 + 10*X[, 4] + 5*X[, 5] + rnorm(n) + w
    data.frame(cbind(y, X))
}

alpha <- 0.2

N <- 10^4 # replications

t <- 10^2 # training size
n <- 10   # calibration size
m <- 10^3 # future batch size

set.seed(42)

coverage <- numeric(N)
pb <- txtProgressBar(min = 1, max = N, style = 3)
for (i in 1:N) {
    db <- friedman(t + n + m)
    trn <- db[1:t, ]
    cal <- db[(t + 1):(t + n), ]
    tst <- db[(t + n + 1):(t + n + m), ]
    rf <- ranger(y ~ ., data = trn, num.trees = 100)
    S <- abs(cal$y - predict(rf, data = cal)$predictions)
    s_hat <- sort(S)[ceiling((1 - alpha)*(n + 1))]
    y_hat <- predict(rf, data = tst)$predictions
    coverage[i] <- mean((y_hat - s_hat < tst$y) & (tst$y < y_hat + s_hat))
    setTxtProgressBar(pb, i)
}
close(pb)

print(Sys.time() - begin)

round(c(1 - alpha, 1 - alpha + 1 / (n  + 1)), 3) # MVP

min(coverage)
which.min(coverage)

hist(coverage, prob = TRUE, breaks = "FD",
     col = "steelblue", xlim = c(0, 1), ylim = c(0, 4),
     xlab = "Empirical Coverage", main = "Split conformal prediction")

g <- function(x) dbeta(x, ceiling((1 - alpha)*(n + 1)), floor(alpha*(n + 1)))
curve(g, col = "red", lwd = 2, add = TRUE)
