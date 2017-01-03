
context("test systemdf")

if(Sys.info()["sysname"] %in% c("Linux", "Darwin")) {

	test_that("test systemdf", {

		df = data.frame(x = sample(1:10, 10), y = sample(11:20, 10))
		df2 = systemdf("sort -k1,1n `df`")
		expect_that(nrow(df2), equals(nrow(df)))
		expect_that(sort(df[[1]]), equals(df2[[1]]))

		expect_error(systemdf("ls -$"))
	})
}
