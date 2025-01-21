#! /usr/bin/env julia -t 15
using PopGen, PooledArrays, CSV

# outlier
yft = read_from("YFT.kinrm.adjustedpop.outlier.gen")

# workaround from bug that doesn't compress the indices #
yft.loci.name = PooledArray(yft.loci.name, compress = true)
yft.loci.population = PooledArray(yft.loci.population, compress = true)
yft.loci.locus = PooledArray(yft.loci.locus, compress = true)

newpops = ["ATL", "GOM","IVC","SEN","VZ"]
populations!(yft, newpops)

# quickly precompile the fst stuff using the outlier data#
fst_outlier = pairwise_fst(yft, method = "Hudson92", iterations = 10000)
CSV.write("hudsonfst_outlier.csv", fst_outlier.results)
# end of precompile stuff #

# purge memory of the PopData
yft = nothing
GC.gc()

#  neutral data
yft = read_from("YFT.kinrm.adjustedpop.neutral.gen")

# workaround from bug that doesn't compress the indices #
yft.loci.name = PooledArray(yft.loci.name, compress = true)
yft.loci.population = PooledArray(yft.loci.population, compress = true)
yft.loci.locus = PooledArray(yft.loci.locus, compress = true)

fst_neutral = pairwise_fst(yft, method = "Hudson92", iterations = 10000)
CSV.write("hudsonfst_neutral.csv", fst_neutral.results)

