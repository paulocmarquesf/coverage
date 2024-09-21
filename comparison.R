# Table 1 of https://arxiv.org/abs/2303.02770
# compared to page 15 of https://arxiv.org/abs/2107.07511v6

tau <- 0.9 # delta = 0.1
alpha <- 0.1

eps_AB <- c(0.1, 0.05, 0.01, 0.005, 0.001)
n_0_AB <- c(22, 102, 2491, 9812, 244390)

n_0 <- c(11, 90, 2429, 9733, 243492) # from table_1.R

b_AB <- ceiling((1 - alpha)*(n_0_AB + 1))
g_AB <- floor(alpha*(n_0_AB + 1))
pr_AB <- pbeta(1 - alpha + eps_AB, b_AB, g_AB) - pbeta(1 - alpha - eps_AB, b_AB, g_AB)

b <- ceiling((1 - alpha)*(n_0 + 1))
g <- floor(alpha*(n_0 + 1))
pr <- pbeta(1 - alpha + eps_AB, b, g) - pbeta(1 - alpha - eps_AB, b, g)

cbind(pr, pr_AB)

pr < pr_AB
