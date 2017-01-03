gr = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 4, 8 ,16), 
    end = c(2, 5, 10, 30)), value = 1:4)
reduce2(gr, gap = bp(2))
reduce2(gr, gap = 0.6)
reduce2(gr, gap = 0.6, max_gap = 4)
