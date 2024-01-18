library(torch)

`%*%.torch_tensor` <- function(e1, e2) {
  torch_matmul(e1, e2)
}

t.torch_tensor <- function(e) {
  torch_transpose(e, -1, -2)
}

diagonal <- function(e) {
  torch_diagonal(e, dim1 = -1, dim2 = -2)
}

sym <- function(e) {
  (e + t(e)) / 2
}

prox <- function(x, lambda) {
  sign(x) * (abs(x) - lambda)$maximum(0)
}

distance <- function(x, y) {
  mapply(function(x, y) {
    max(abs(x - y))
  }, x, y)
}

EM.gvemm <- function(Y, D, X, lambda, iter, eps) {
  init <- function() {
    with(parent.frame(), {
      N <- nrow(Y)
      J <- ncol(Y)
      K <- ncol(D)
      G <- max(X)

      Y <- torch_tensor(Y)
      eta <- torch_tensor(matrix(0.125, N, J))
      Sigma <- torch_stack(replicate(G, diag(K), simplify = F))
      Mu <- torch_zeros(c(G, K))
      a.mask <- torch_tensor(D * 1.0)
      a.mask.diag <- a.mask$diag_embed()
      a <- a.mask$clone()
      b <- torch_zeros(J)
      gamma.mask <- torch_stack(c(torch_zeros_like(a.mask), replicate(G - 1, a.mask)))
      gamma <- torch_zeros_like(gamma.mask)
      beta.mask <- torch_cat(list(torch_zeros(c(1, J)), torch_ones(c(G - 1, J))))
      beta <- torch_zeros_like(beta.mask)
    })
  }

  parameters <- function() {
    lapply(list(SIGMA = SIGMA, MU = MU, Sigma = Sigma, Mu = Mu, a = a, b = b, gamma = gamma, beta = beta), as.array)
  }

  update <- function() {
    with(parent.frame(), {
      aG <- (a + gamma)$unsqueeze(4)
      aG.t <- t(aG)
      AG <- aG[X]
      AG.t <- t(AG)
      BB <- (b - beta)[X]

      Sigma.inv <- Sigma$inverse()
      SIGMA.inv <- Sigma.inv[X] + 2 * (eta$view(c(N, -1, 1, 1)) * (aG %*% aG.t)[X])$sum(2)
      SIGMA <- sym(SIGMA.inv$inverse())
      MU <- (SIGMA %*% (((Y - 0.5 + 2 * eta * BB)$view(c(N, -1, 1, 1)) * AG)$sum(2) + (Sigma.inv %*% Mu$unsqueeze(3))[X]))$squeeze(3)
      Mu <- torch_stack(tapply(1:N, X, function(n) {
        MU[n]$mean(1)
      }))
      mu <- MU - Mu[X]
      sigma.mu <- SIGMA + mu$unsqueeze(3) %*% mu$unsqueeze(2)
      Sigma <- sym(torch_stack(tapply(1:N, X, function(n) {
        sigma.mu[n]$mean(1)
      })))

      mu <- a.mask.diag %*% MU$view(c(N, 1, -1, 1))
      sigma.mu <- a.mask.diag %*% SIGMA$unsqueeze(2) %*% a.mask.diag + mu %*% t(mu)
      xi <- sqrt(BB$square() - 2 * BB * (AG.t %*% mu)$view(c(N, -1)) + (AG.t %*% sigma.mu %*% AG)$view(c(N, -1)))
      eta <- torch_where(abs(xi) < 1e-3, 0.125, (1 / (1 + exp(-xi)) - 0.5) / (2 * xi))
      a <- ((2 * eta$view(c(N, -1, 1, 1)) * sigma.mu)$sum(1)$pinverse() %*% ((Y - 0.5)$view(c(N, -1, 1, 1)) * mu + 2 * eta$view(c(N, -1, 1, 1)) * (BB$view(c(N, -1, 1, 1)) * mu - sigma.mu %*% gamma[X]$unsqueeze(4)))$sum(1))$squeeze(3) * a.mask
      b <- (0.5 - Y + 2 * eta * (beta[X] + (AG.t %*% mu)$view(c(N, -1))))$sum(1) / (2 * eta$sum(1))

      gamma.beta  <- torch_stack(tapply(1:N, X, function(n) {
        N <- length(n)
        torch_cat(list(prox(((Y[n] - 0.5)$unsqueeze(3) * mu[n]$squeeze(4) + 2 * eta[n]$unsqueeze(3) * (BB[n]$unsqueeze(3) * mu[n]$squeeze(4) - (sigma.mu[n] %*% a$unsqueeze(3))$squeeze(4)))$sum(1), lambda) / diagonal((2 * eta[n]$view(c(N, -1, 1, 1)) * sigma.mu[n])$sum(1)),
                       (prox(((Y[n] - 0.5) + 2 * eta[n] * (b - (AG.t[n] %*% mu[n])$view(c(N, -1))))$sum(1), lambda) / (2 * eta[n]$sum(1)))$unsqueeze(2)), 2)
      }))
      gamma$set_data(gamma.beta[, , 1:K]$masked_fill(gamma.mask == 0, 0))
      beta$set_data(gamma.beta[, , (K + 1)] * beta.mask)

      mu <- Mu[1]$clone()
      MU$sub_(mu)
      Mu$sub_(mu)
      b$sub_(a %*% mu)
      beta$add_(gamma %*% mu)
      sigma <- Sigma[1]$diag()$sqrt()
      a$mul_(sigma)
      gamma$mul_(sigma)
      sigma.inv <- (1 / sigma)$diag()
      SIGMA$set_data(sym(sigma.inv %*% SIGMA %*% sigma.inv))
      Sigma$set_data(sym(sigma.inv %*% Sigma %*% sigma.inv))
    })
  }

  init()
  params.old <- NULL
  for (i in 1:iter) {
    update()
    params <- parameters()
    if (!is.null(params.old) && all(distance(params, params.old) < eps))
      break
    params.old <- params
  }
  lambda <- 0
  gamma.mask <- gamma != 0
  beta.mask <- beta != 0
  params.old <- NULL
  for (i in 1:iter) {
    update()
    params <- parameters()
    if (!is.null(params.old) && all(distance(params, params.old) < eps))
      break
    params.old <- params
  }
  params
}

IC.gvemm <- function(Y, X, SIGMA, MU, Sigma, Mu, a, b, gamma, beta, c) {
  N <- nrow(Y)
  K <- max(X)
  Y <- torch_tensor(Y)
  SIGMA <- torch_tensor(SIGMA)
  MU <- torch_tensor(MU)$unsqueeze(2)
  Sigma <- torch_tensor(Sigma)
  Mu <- torch_tensor(Mu)
  a <- torch_tensor(a)
  b <- torch_tensor(b)
  gamma <- torch_tensor(gamma)
  beta <- torch_tensor(beta)

  AG <- (a + gamma)[X]$unsqueeze(4)
  AG.t <- t(AG)
  BB <- (b - beta)[X]
  xi <- sqrt(BB$square() - 2 * BB * (AG.t %*% MU$view(c(N, 1, -1, 1)))$view(c(N, -1)) + (AG.t %*% (SIGMA + t(MU) %*% MU)$unsqueeze(2) %*% AG)$view(c(N, -1)))
  MU$unsqueeze_(4)
  mu <- MU$squeeze(2) - Mu[X]$unsqueeze(3)
  Q <- as.array((nnf_logsigmoid(xi) + (0.5 - Y) * (BB - (AG.t %*% MU)$view(c(N, -1))) - xi / 2)$sum() - (Sigma$logdet()[X]$sum() + diagonal(linalg_solve(Sigma[X], SIGMA + mu %*% t(mu)))$sum()) / 2)
  l0 <- as.array(sum(gamma != 0) + sum(beta != 0))
  c(ll = Q, l0 = l0, AIC = -2 * Q + l0 * 2, BIC = -2 * Q + l0 * log(N), GIC = -2 * Q + c * l0 * log(N) * log(log(N)))
}

gvemm <- function(Y, D, X, c, progress, Lambda0 = seq(0.2, 0.8, by = 0.1), iter = 1000, eps = 1e-3) {
  if (is.character(X))
    X <- as.factor(X)
  if (!is.integer(X))
    X <- as.integer(X)
  N <- nrow(Y)
  
  get_ic <- function(result) {
    as.data.frame(t(sapply(result, function(r) {
      c(lambda0 = r$lambda0, ic = r$GIC)
    })))
  }
  
  run <- function(Lambda0) {
    lapply(Lambda0, function(lambda0) {
      progress$set(detail = paste0(lambda0))
      lambda <- sqrt(N) * lambda0
      em <- EM.gvemm(Y, D, X, lambda, iter, eps)
      list2env(em, environment())
      list2env(as.list(IC.gvemm(Y, X, SIGMA, MU, Sigma, Mu, a, b, gamma, beta, c)), environment())
      progress$inc(1 / length(Lambda0))
      list(lambda0 = lambda0, lambda = lambda, SIGMA = SIGMA, MU = MU, Sigma = Sigma, Mu = Mu, a = a, b = b, gamma = gamma, beta = beta, ll = ll, l0 = l0, AIC = AIC, BIC = BIC, GIC = GIC)
    })
  }
  
  result <- run(Lambda0)
  ic <- get_ic(result)
  ic_min <- which.min(ic$ic)
  if (ic_min == 1) {
    progress$set(0, 'Trying additional lambdas: ')
    result <- c(run(seq(0, Lambda0[1] - 0.1, by = 0.1)), result)
    ic <- get_ic(result)
    ic_min <- which.min(ic$ic)
  }
  lambda <- Lambda0[length(Lambda0)] + 0.1
  while (ic_min == nrow(ic)) {
    progress$set(0, 'Trying additional lambdas: ')
    result <- c(result, run(seq(lambda, lambda + 0.5, by = 0.1)))
    ic <- get_ic(result)
    ic_min <- which.min(ic$ic)
    lambda <- lambda + 0.4
  }
  list(result[[ic_min]], ic)
}
