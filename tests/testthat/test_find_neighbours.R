context("test find_neighbours")

gr1 = GRanges(seqname = "chr1", ranges = IRanges(start = c(4, 10), end = c(6, 16)))
gr2 = GRanges(seqname = "chr1", ranges = IRanges(start = c(7, 13), end = c(8, 20)))


test_that("test find_neighbours", {
	expect_that(length(find_neighbours(gr1, gr2, upstream = 3, downstream = 3)), equals(3))
	expect_that(length(find_neighbours(gr1, gr2, upstream = 10, downstream = 10)), equals(4))
})
