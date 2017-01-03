context("Test `logical_segment`")

test_that("Test logical_segment", {
    
    l = c(TRUE, FALSE)
    expect_that(logical_segment(l), equals(data.frame(start_index = 1, end_index = 1)))
    
    l = c(TRUE, FALSE, TRUE)
    expect_that(logical_segment(l), equals(data.frame(start_index = c(1, 3), end_index = c(1, 3))))

    l = c(TRUE, FALSE, TRUE, TRUE, FALSE)
    expect_that(logical_segment(l), equals(data.frame(start_index = c(1, 3), end_index = c(1, 4))))

    l = c(FALSE, TRUE, FALSE, TRUE, TRUE, FALSE)
    expect_that(logical_segment(l), equals(data.frame(start_index = c(2, 4), end_index = c(2, 5))))

})

