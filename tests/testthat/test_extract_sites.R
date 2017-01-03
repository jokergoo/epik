context("test extract_sites")

site = c(2, 5, 9, 10, 15, 20)

test_that("test binary_search()", {
	expect_that(binary_search(site, 1, TRUE), equals(NA+0))
	expect_that(binary_search(site, 1, FALSE), equals(1))
	expect_that(binary_search(site, 5, TRUE), equals(2))
	expect_that(binary_search(site, 5, FALSE), equals(2))
	expect_that(binary_search(site, 12, TRUE), equals(4))
	expect_that(binary_search(site, 12, FALSE), equals(5))
	expect_that(binary_search(site, 30, TRUE), equals(6))
	expect_that(binary_search(site, 30, FALSE), equals(NA+0))

	expect_that(binary_search(site, c(1, 5, 12, 30), FALSE), equals(c(1, 2, 5, NA+0)))
	expect_that(binary_search(site, c(1, 5, 12, 30), TRUE), equals(c(NA+0, 2, 4, 6)))
})

test_that("test extract_sites", {

	expect_that(extract_sites(0, 1, site, FALSE, 0), equals(integer(0)))
	expect_that(extract_sites(0, 2, site, FALSE, 0), equals(2))
	expect_that(extract_sites(1, 8, site, FALSE, 0), equals(c(2, 5)))
	expect_that(extract_sites(2, 5, site, FALSE, 0), equals(c(2, 5)))
	expect_that(extract_sites(3, 4, site, FALSE, 0), equals(integer(0)))
	expect_that(extract_sites(1, 21, site, FALSE, 0), equals(c(2, 5, 9, 10, 15, 20)))
	expect_that(extract_sites(9, 16, site, FALSE, 0), equals(c(9, 10, 15)))
	expect_that(extract_sites(15, 20, site, FALSE, 0), equals(c(15, 20)))
	expect_that(extract_sites(20, 25, site, FALSE, 0), equals(c(20)))
	expect_that(extract_sites(25, 26, site, FALSE, 0), equals(integer(0)))
	expect_that(extract_sites(c(1, 9, 15), c(3, 12, 20), site, FALSE, 0), equals(c(2, 9, 10, 15, 20)))
	expect_that(extract_sites(c(1, 3), c(9, 12), site, FALSE, 0), equals(extract_sites(c(1, 7), c(6, 12), site, FALSE, 0)))
	for(i in 1:10) {
		site = sort(sample(1000, 100))
		pos = do.call("rbind", lapply(1:10, function(i) sort(sample(max(site), 2))))
		ir = IRanges(pos[, 1], pos[, 2])
		ir_site = IRanges(site, site)
		mtch = as.matrix(findOverlaps(ir_site, ir))
		s1 = start(ir_site[unique(mtch[, 1])])
		s2 = extract_sites(pos[, 1], pos[, 2], site, FALSE, 0)
		expect_that(s1, equals(s2))

		site = sort(sample(1000, 100))
		ir_site = IRanges(site, site)
		mtch = as.matrix(findOverlaps(ir_site, ir))
		s1 = start(ir_site[unique(mtch[, 1])])
		s2 = extract_sites(pos[, 1], pos[, 2], site, FALSE, 0)
		expect_that(s1, equals(s2))
	}
})

if(Sys.getenv("IS_PBS") != "") {

site = sort(sample(10000000, 1000000))
pos = do.call("rbind", lapply(1:1000, function(i) sort(sample(max(site), 2))))

system.time(for(i in 1:1000) {
	site[site >= pos[i, 1] & site <= pos[i, 2]]
})

system.time(for(i in 1:1000) {
	extract_sites(pos[i, 1], pos[i, 2], site, FALSE, 0)
})

system.time(extract_sites(pos[, 1], pos[, 2], site, FALSE, 0))

ir = IRanges(pos[, 1], pos[, 2])
ir_site = IRanges(site, site)

system.time({
	mtch <- as.matrix(findOverlaps(ir_site, ir))
	s1 <- start(ir_site[unique(mtch[, 1])])
})		

}
