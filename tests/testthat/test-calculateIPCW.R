#test_that("calculateIPCW produced expected IP weights in diabetes.complications", {
#  library(mets)
#  data(diabetes.complications)
#  output <- calculateIPCW_new(formula=Event(t,epsilon)~+1, data=diabetes.complications, code.censoring=0, strata_name='strata', specific.time=8)
#  tested <- round(output[1:6], digit=3)
#  expected <- c(1.656, 1.656, 0, 1.656, 1.656, 1.004)
#  expect_equal(tested, expected)
#})

#test_that("calculateIPCW produced expected IP weights in diabetes.complications", {
#  library(mets)
#  data(diabetes.complications)
#  output <- calculateIPCW_new(formula=Event(t,epsilon)~+1, data=diabetes.complications, code.censoring=0, strata_name='strata', specific.time=8)
#  tested <- round(output[1:6], digit=3)
#  names(tested) <- NULL
#  output <- calculateIPCW(formula=Event(t,epsilon)~+1, data=diabetes.complications, code.censoring=0, strata_name='strata', specific.time=8)
#  expected <- round(output[1:6], digit=3)
#  names(expected) <- NULL
#  expect_equal(tested, expected)
#})

