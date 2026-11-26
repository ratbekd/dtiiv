test_that("moment derivatives match finite differences", {
  M <- dtiiv:::.moment_registry()
  set.seed(3); s <- rnorm(50); u <- rnorm(50); h <- 1e-6
  for (m in M) {
    expect_equal(m$gs(s, u), (m$g(s + h, u) - m$g(s - h, u)) / (2 * h), tolerance = 1e-6, info = m$id)
    expect_equal(m$gu(s, u), (m$g(s, u + h) - m$g(s, u - h)) / (2 * h), tolerance = 1e-6, info = m$id)
  }
})

test_that("the sets cover the true slope and exclude OLS in a skewed design", {
  set.seed(11)
  d <- dtiiv_simulate(4000, beta = 1, lambda = 0.5, clusters = 100)
  fit <- dtiiv(y ~ x, d, controls = ~ w, cluster = ~ cluster, sets = c("VT", "four", "mean", "error"))
  expect_s3_class(fit, "dtiiv")
  expect_gt(fit$ols[["estimate"]], 1.2)                      # OLS is biased upwards
  for (nm in names(fit$sets)) {
    cs <- fit$sets[[nm]]$confset
    expect_true(any(cs$lower <= 1 & cs$upper >= 1), info = nm)
    expect_false(any(cs$lower <= fit$ols[["estimate"]] & cs$upper >= fit$ols[["estimate"]]), info = nm)
    expect_equal(unname(coef(fit, nm)), 1, tolerance = 0.1, info = nm)
  }
})

test_that("the statistic at theta = 0 is the OLS slope and the set is invariant to rescaling", {
  set.seed(5)
  d <- dtiiv_simulate(1500)
  f1 <- dtiiv(y ~ x, d, controls = ~ w, sets = "four")
  d2 <- transform(d, y = 3 * y + 7, w = 10 * w)
  f2 <- dtiiv(y ~ x, d2, controls = ~ w, sets = "four")
  expect_equal(3 * f1$sets$four$estimate, f2$sets$four$estimate, tolerance = 1e-5)
  expect_equal(3 * as.matrix(f1$sets$four$confset[, 1:2]), as.matrix(f2$sets$four$confset[, 1:2]),
               tolerance = 1e-5)
  expect_equal(f1$sets$four$profile$beta[(f1$ngrid + 1) / 2], f1$ols[["estimate"]])
})

test_that("dtiiv_test agrees with the confidence set", {
  set.seed(7)
  d <- dtiiv_simulate(2000, clusters = 50)
  fit <- dtiiv(y ~ x, d, controls = ~ w, cluster = ~ cluster, sets = "two-sided")
  cs <- fit$sets$`two-sided`$confset
  inside <- (cs$lower[1] + cs$upper[1]) / 2
  tt <- dtiiv_test(fit, c(inside, fit$ols[["estimate"]]))
  expect_true(tt$p_value[1] > 0.05)
  expect_true(tt$p_value[2] < 0.05)
  ends <- dtiiv_test(fit, c(cs$lower[1], cs$upper[1]))
  if (!cs$lower_at_edge[1]) expect_equal(ends$stat[1], ends$crit[1], tolerance = 1e-6)
  if (!cs$upper_at_edge[1]) expect_equal(ends$stat[2], ends$crit[2], tolerance = 1e-6)
})

test_that("user-defined moments and custom sets work", {
  set.seed(2)
  d <- dtiiv_simulate(2000)
  m5 <- dtiiv_moment("us5", function(s, u) u * s^5, function(s, u) 5 * u * s^4, function(s, u) s^5)
  fit <- dtiiv(y ~ x, d, controls = ~ w, sets = list(custom = list("us2", m5), poly = c("us2", "us3")))
  expect_named(fit$sets, c("custom", "poly"))
  expect_equal(fit$sets$custom$moments, c("us2", "us5"))
  expect_error(dtiiv(y ~ x, d, sets = "nonsense"), "unknown set")
  expect_error(dtiiv(y ~ x, d, sets = list(a = c("us2", "zzz"))), "unknown moment")
})

test_that("input checks and warnings", {
  set.seed(1)
  d <- dtiiv_simulate(500)
  expect_error(dtiiv(y ~ x + w, d), "exactly one regressor")
  expect_error(dtiiv(y ~ x, d, controls = ~ nothere), "not found")
  d$xb <- as.numeric(d$x > 1); d$yb <- as.numeric(d$y > 3)
  expect_error(dtiiv(yb ~ xb, d, sets = "mean"), "binary")
  expect_warning(dtiiv(y ~ xb, d, sets = "VT"), "binary")
  d$x[1:5] <- NA
  fit <- dtiiv(y ~ x, d, sets = "VT")
  expect_equal(fit$n_dropped, 5)
  expect_equal(fit$n, 495)
})

test_that("methods run", {
  set.seed(4)
  d <- dtiiv_simulate(1000, clusters = 40)
  fit <- dtiiv(y ~ x, d, controls = ~ w, cluster = ~ cluster, sets = c("VT", "mean"), instruments = ~ w)
  expect_output(print(fit), "DT-IIV")
  expect_output(print(summary(fit)), "variance tendency")
  expect_true(is.matrix(confint(fit, set = "mean")))
  expect_error(confint(fit, level = 0.9), "refit")
  expect_length(coef(fit, "all"), 2)
  pdf(NULL); expect_silent(plot(fit)); dev.off()
  expect_s3_class(dtiiv_moments()$moments, "data.frame")
})
