# pam.rsph's methods are now registered with S3method(), so UseMethod finds
# them from any environment. FINDING 17 fixed. The direct-access helper is kept
# so the method-level tests stay independent of dispatch.
rsph_method <- function(cls) {
  get(paste0("pam.rsph.", cls), envir = asNamespace("TimeMetric"))
}

test_that("pam.rsph methods are registered with S3method (FINDING 17 fixed)", {
  # All three methods are now declared with @exportS3Method, so NAMESPACE
  # carries S3method(pam.rsph, <class>) and UseMethod resolves them anywhere.
  root <- skip_without_source_tree()
  ns <- readLines(file.path(root, "NAMESPACE"))

  expect_true(any(grepl("S3method(pam.rsph,coxph)", ns, fixed = TRUE)))
  expect_true(any(grepl("S3method(pam.rsph,survreg)", ns, fixed = TRUE)))
  expect_true(any(grepl("S3method(pam.rsph,aareg)", ns, fixed = TRUE)))
})

test_that("pam.rsph dispatches in-process for coxph and survreg", {
  d <- fx_surv()

  expect_no_error(TimeMetric:::pam.rsph(fx_cox(), test_data = d))
  expect_no_error(TimeMetric:::pam.rsph(fx_survreg(), test_data = d))
  # dispatch and direct method call must agree
  expect_equal(
    TimeMetric:::pam.rsph(fx_cox(), test_data = d)$Re,
    rsph_method("coxph")(fx_cox(), test_data = d)$Re,
    tolerance = 1e-8
  )
})

test_that("a clean session can now dispatch pam.rsph (FINDING 17 fixed)", {
  # The user-facing check. Run in a fresh subprocess because testthat's own test
  # environment inherits the package namespace and would mask a regression.
  skip_on_cran()
  skip_if_covr()
  pkg_root <- skip_without_source_tree()

  script <- sprintf('
    suppressWarnings(pkgload::load_all(%s, quiet = TRUE, attach_testthat = FALSE))
    d <- sim_cox_weibull_censored(n = 50, pi_c = 0.3, v = 2,
                                  beta = c(0.5, -0.5), seed = 1001)
    d <- d[, c("time", "status", "x1", "x2")]
    m <- survival::coxph(survival::Surv(time, status) ~ x1 + x2,
                         data = d, x = TRUE, y = TRUE)
    r <- try(TimeMetric:::pam.rsph(m, test_data = d), silent = TRUE)
    cat(if (inherits(r, "try-error")) as.character(r) else "DISPATCH WORKED")
  ', shQuote(pkg_root))

  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    args = c("--vanilla", "-e", shQuote(script)),
    stdout = TRUE, stderr = TRUE
  ))
  out <- paste(out, collapse = "\n")

  expect_no_match(out, "no applicable method")
  expect_match(out, "DISPATCH WORKED")
})

test_that("pam.rsph.coxph works when called directly", {
  d <- fx_surv()

  res <- rsph_method("coxph")(fx_cox(), test_data = d)

  expect_type(res, "list")
  # components pam.summary.rsph consumes
  expect_true(all(c("meanr", "ranks", "perfr", "Re", "times") %in% names(res)))
  expect_true(all(is.finite(res$ranks)))
  expect_snapshot_value(sort(names(res)), style = "serialize")
  expect_snapshot_value(snap_num(head(res$ranks, 10)), style = "serialize")
  expect_snapshot_value(snap_num(res$Re), style = "serialize")
})

test_that("pam.rsph.survreg works when called directly", {
  d <- fx_surv()

  res <- rsph_method("survreg")(fx_survreg(), test_data = d)

  expect_type(res, "list")
  expect_true(all(c("meanr", "ranks", "perfr") %in% names(res)))
  expect_snapshot_value(sort(names(res)), style = "serialize")
})

test_that("the coxph and survreg methods are genuinely different code paths", {
  d <- fx_surv()

  cox_res <- rsph_method("coxph")(fx_cox(), test_data = d)
  reg_res <- rsph_method("survreg")(fx_survreg(), test_data = d)

  expect_false(isTRUE(all.equal(cox_res$ranks, reg_res$ranks, tolerance = 1e-6)))
})

test_that("pam.summary.rsph converts an rsph object into R_E over time", {
  d <- fx_surv()
  obj <- rsph_method("coxph")(fx_cox(), test_data = d)

  res <- TimeMetric:::pam.summary.rsph(obj, times = stats::median(d$time))

  expect_s3_class(res, "data.frame")
  expect_identical(names(res), c("times", "Rti", "dRti"))
  expect_identical(nrow(res), 1L)
  expect_true(all(is.finite(res$Rti)))
  expect_snapshot_value(snap_num(c(res$times, res$Rti, res$dRti)),
                        style = "serialize")
})

test_that("pam.rsph.aareg exists as a third method", {
  ns <- asNamespace("TimeMetric")

  expect_true(exists("pam.rsph.aareg", envir = ns, inherits = FALSE))
  expect_true(is.function(get("pam.rsph.aareg", envir = ns)))
})

test_that("no pam.print generic exists, so pam.print.rsph is unreachable", {
  # Documents FINDING: pam.print.rsph can never be dispatched.
  ns <- asNamespace("TimeMetric")

  expect_false(exists("pam.print", envir = ns, inherits = FALSE))
  expect_true(exists("pam.print.rsph", envir = ns, inherits = FALSE))
})
