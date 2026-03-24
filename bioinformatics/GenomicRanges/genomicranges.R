
#===============================================================================
# GenomicRanges
# IRanges
# GenomicFeatures

# Lawrence, M., Huber, W., Pagès, H., Aboyoun, P., Carlson, M., 
# Gentleman, R., Morgan, M. T., & Carey, V. J. (2013). Software for computing 
# and annotating genomic ranges. PLoS computational biology, 9(8), e1003118. 
# https://doi.org/10.1371/journal.pcbi.1003118
#===============================================================================

library(GenomicRanges)
library(IRanges)
library(GenomicFeatures)

#-------------------------------------------------------------------------------
# 1. IRanges

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
# comparisons and overlaps
?findOverlaps
findOverlaps(ir1, ir2)


#===============================================================================
# GenomicRanges

library(BiocManager)
BiocManager::install('karyoploteR')

library(karyotype)# required to plot GRange objects

# creating GRange object
gr1 <- GRanges(seqnames = 'chr1', strand = c('+', '-', '+'), 
               ranges = IRanges(start = c(1,3, 5), width = 3 ) )

gr1

plotGR(gr1)

?GRanges

flank(gr1) # function creates new ranges that are located to the side of the 
# original ranges. It is commonly used to isolate upstream regulatory regions or downstream termination regions.
promoters(gr1) # extract the promoter region around the Transcription Start Site (TSS)


#-------------------------------------------------------------------------------
# Example using karyoploteR
# 1. Define the genomic context (e.g., human genome hg19)
# We limit it to chr1 since that is where gr1 is located
kp <- plotKaryotype(genome = "hg19", chromosomes = "chr1")

# 2. Plot your gr1 object
plotGR(kp, data = gr1)
