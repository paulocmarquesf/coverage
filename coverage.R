# On the Universal distribution of the coverage in split conformal prediction

# Paulo C. Marques F. <PauloCMF1@insper.edu.br>

# Paper: https://bit.ly/3muCZJC

library(tidyverse)
library(ggExtra)
library(ranger)

friedman <- function(n, p = 10) {
    w <- rexp(1)
    X <- matrix(runif(n*p), nrow = n, ncol = p, dimnames = list(1:n, paste0("x_", 1:p)))
    y <- 10*sin(pi*X[, 1]*X[, 2]) + 20*(X[, 3] - 0.5)^2 + 10*X[, 4] + 5*X[, 5] + rnorm(n) + w
    data.frame(cbind(y, X))
}

alpha <- 1 - 0.8

N <- 10^4 # replications

r <- 10^2 # training
n <- 10   # calibration
m <- 10^3 # future horizon

set.seed(42)

coverage <- numeric(N)
pb <- txtProgressBar(min = 1, max = N, style = 3)
for (t in 1:N) {
    db <- friedman(r + n + m)
    training <- db[1:r, ]
    calibration <- db[(r + 1):(r + n), ]
    test <- db[(r + n + 1):(r + n + m), ]
    model <- ranger(y ~ ., data = training)
    R <- abs(calibration$y - predict(model, data = calibration)$predictions)
    s_hat <- sort(R)[ceiling((1 - alpha)*(n + 1))]
    y_hat <- predict(model, data = test)$predictions
    coverage[t] <- mean((y_hat - s_hat < test$y) & (test$y < y_hat + s_hat))
    setTxtProgressBar(pb, t)
}
close(pb)

round(c(1 - alpha, 1 - alpha + 1 / (n  + 1)), 3)

which.min(coverage)
min(coverage)

###

theme_set(theme_bw())

plt <- tibble(x = seq_along(coverage), y = coverage) %>%
    ggplot(aes(x, y)) +
        geom_point(size = 0.5, alpha = 0.25) +
        geom_hline(yintercept = 1 - alpha, size = 0.95, linetype = "dashed", alpha = 0.85) +
        geom_hline(yintercept = 1 - alpha + 1 / (n + 1), size = 0.95, linetype = "dashed", alpha = 0.85) +
        scale_x_continuous(limits = c(1, N), breaks = c(1, 2500, 5000, 7500, 10000), expand = c(0.015, 0.015),
                           labels = scales::comma_format(big.mark = " ")) +
        scale_y_continuous(limits = c(0, 1), expand = c(0.025, 0.025)) +
        labs(x = "Replication", y = "Future coverage") +
        theme(axis.title.y = element_text(vjust = 2.5),
              axis.title.x = element_text(vjust = -0.75))

fd <- 2 * IQR(coverage) / length(coverage)^(1/3) # Freedman-Diaconis rule

ggMarginal(plt, type = "histogram", margins = "y", size = 15,
           yparams = list(binwidth = fd, fill = "light gray"))

###

nominal_alpha <- 1 - 0.9
eps <- 0.02
gamma <- 0.95

g <- function(n, alpha, eps, gamma) {
    b <- ceiling((1 - alpha)*(n + 1))
    g <- floor(alpha*(n + 1))
    abs(pbeta(1 - alpha + eps, b, g) - pbeta(1 - alpha - eps, b, g) - gamma)
}

(calibration_size <- ceiling(optimize(function(n) g(n, nominal_alpha, eps, gamma), lower = 10, upper = 10^4)$minimum))
