plot_input <- function() {
  pam.coxph_restricted(model = fx_cox(), covs = fx_covs(),
                       new_data = fx_surv(), tau = 10e10)
}

test_that("plot_pred requires linear.pred, pred, times and status", {
  pred <- plot_input()

  expect_true(all(c("linear.pred", "pred", "times", "status") %in% names(pred)))
  expect_error(
    plot_pred(pred[c("pred", "times", "status")]),
    "Input data must contain columns: linear.pred, pred, times, status"
  )
})

test_that("plot_pred with the default sample_index draws NOTHING (FINDING 18)", {
  # sample_index defaults to NULL, and the function then subsets with
  # x_var[sample_index]. x[NULL] is a zero-length vector, so the data frame
  # behind the plot is empty and the chart is blank -- silently, with no
  # warning, on the documented default call.
  #
  # When this is fixed, invert the test to expect 200 rows.
  p <- plot_pred(plot_input())
  built <- ggplot2::ggplot_build(p)

  expect_s3_class(p, "ggplot")
  expect_identical(nrow(p$data), 0L)
  expect_identical(nrow(built$data[[1]]), 0L)
  expect_identical(nrow(built$data[[2]]), 0L)
})

test_that("plot_pred renders every subject when sample_index is supplied", {
  p <- plot_pred(plot_input(), sample_index = seq_len(200))
  built <- ggplot2::ggplot_build(p)

  expect_identical(length(p$layers), 2L)
  expect_identical(nrow(built$data[[1]]), 200L)
  expect_identical(nrow(built$data[[2]]), 200L)
  expect_snapshot_value(
    vapply(built$data, nrow, integer(1)), style = "serialize"
  )
})

test_that("plot_pred uses the supplied labels", {
  p <- plot_pred(plot_input(), sample_index = seq_len(200),
                 title = "characterization", xlab = "RS", ylab = "Days")

  expect_identical(p$labels$title, "characterization")
  expect_identical(p$labels$x, "RS")
  expect_identical(p$labels$y, "Days")
})

test_that("plot_pred defaults xlab to Risk Score and ylab to Days", {
  p <- plot_pred(plot_input(), sample_index = seq_len(200))

  expect_identical(p$labels$x, "Risk Score")
  expect_identical(p$labels$y, "Days")
})

test_that("plot_pred restrict_time caps the plotted times", {
  cap <- 1
  p <- plot_pred(plot_input(), sample_index = seq_len(200), restrict_time = cap)

  expect_true(all(p$data$times <= cap))
  expect_gt(sum(p$data$times == cap), 0)
})

test_that("summary_pred_plot combines one panel per input", {
  pred <- plot_input()

  p <- summary_pred_plot(list(a = pred, b = pred), ncol = 2)

  expect_s3_class(p, "patchwork")
  expect_identical(length(p$patches$plots) + 1L, 2L)
  expect_no_error(suppressWarnings(print(p)))
})

test_that("summary_pred_plot accepts explicit panel titles", {
  pred <- plot_input()

  p <- summary_pred_plot(list(a = pred, b = pred),
                         titles = c("first", "second"), ncol = 2)

  expect_s3_class(p, "patchwork")
  expect_no_error(suppressWarnings(print(p)))
})
