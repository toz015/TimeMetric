# pam.rsph's methods are defined in the namespace but never registered with
# S3method(), so UseMethod cannot find them from outside the package. Tests that
# need a method call it directly out of the namespace, which is what these
# helpers do. See FINDING 17.
rsph_method <- function(cls) {
  get(paste0("pam.rsph.", cls), envir = asNamespace("TimeMetric"))
}

test_that("pam.rsph methods are never registered with S3method (FINDING 17)", {
  # NAMESPACE contains no S3method() directives at all, so pam.rsph.coxph and
  # friends are invisible to UseMethod outside TimeMetric's own namespace.
  # Inside the package -- and inside testthat, whose test environment inherits
  # the namespace -- dispatch resolves, which is why R_E still computes in
  # pam.predicted_survial_eval. A user calling from the global environment gets
  # "no applicable method". Reproduced in a clean subprocess below.
  root <- skip_without_source_tree()

  expect_false(
    any(grepl("S3method", readLines(file.path(root, "NAMESPACE"))))
  )
})

test_that("a clean session cannot dispatch pam.rsph (FINDING 17)", {
  # The user-facing consequence of the missing S3method registration. Run in a
  # fresh subprocess because testthat's own test environment inherits the
  # package namespace and would mask the failure.
  #
  # When the methods are registered (spec step 3), invert this test.
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

  expect_match(out, "no applicable method")
  expect_match(out, "pam.rsph")
  expect_no_match(out, "DISPATCH WORKED")
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
