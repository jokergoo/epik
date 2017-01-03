context("test reduce2")

test_that("test reduce2", {
	gr = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 4, 8 ,12), end = c(2, 5, 10, 20)))
	gr2 = reduce2(gr, gap = bp(2))
	expect_that(ranges(gr2), equals(IRanges(start = c(1, 8), end = c(5, 20))))

	gr2 = reduce2(gr, gap = 0.1)
	expect_that(ranges(gr2), equals(IRanges(start = c(1, 4, 8), end = c(2, 5, 20))))
})
