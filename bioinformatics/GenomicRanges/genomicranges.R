
library(GenomicRanges)
library(IRanges)

# creating IRanges
ir1 <- IRanges(c(1, 5, 10), c(3, 7, 15))
ir2 <- IRanges(c(2,4,6), c(7, 9, 15))


ir1
plotGR(ir1)
plotGR(ir2)

start(ir1)
width(ir1)

# inter-range operations
# reduce and disjoin
reduce(ir2)
disjoin(ir2)
ir22 <- disjoin(ir2)
plotGR(ir22)

gaps(ir2)

# intra-range operations
shift(ir2, shift = 4)


# advanced operations
?findOverlaps
findOverlaps(ir1, ir2)
