using PopGen
using DataFrames
using CSV
using PooledArrays

data = read_from("/home/pdimens/Omega/USM PhD/Projects/Active/Yellowfin Tuna/Popgen/Analyses/nonkin/data/YFT.snp.kinrm.pcrelate.gen")
populations!(data, ["wATL","GOM","IVC", "GOM", "SEN", "GOM", "VZ"])
data.loci.population = PooledArray(data.loci.population, compress = true)

global_fst = summary_stats(data, by = "locus")
CSV.write("locbyloc.global.nei.fst", global_fst)

locbyloc = _pairwise_Hudson(data)
CSV.write("locbyloc.pairwise.hudson.fst", locbyloc)
function _pairwise_Hudson(data::PopData)
    !isbiallelic(data) && throw(error("Data must be biallelic to use the Hudson estimator"))
    idx_pdata = groupby(data.loci, :population)
    pops = getindex.(keys(idx_pdata), :population)
    npops = length(idx_pdata)
    n_loci = size(data)[2]
    locnames = loci(data)
    results = zeros(npops, npops)
    res = Float64[]
    p1 = String[]
    p2 = String[]
    locs = String[]
    for i in 2:npops
        for j in 1:(i-1)
            pop1 = reshape(idx_pdata[i].genotype, :, n_loci)
            pop2 = reshape(idx_pdata[j].genotype, :, n_loci)
            #results[i,j] = hudson_fst(pop1,pop2)
            append!(res, hudson_fst(pop1, pop2))
            append!(p1, fill(pops[i], n_loci))
            append!(p2, fill(pops[j], n_loci))
            append!(locs, locnames)
        end
    end
    return DataFrame(:pop1 => p1, :pop2 => p2,:locus => locs ,:fst => res)
end

function hudson_fst(population_1::T, population_2::T) where T<:AbstractMatrix
    fst_perloc = @views [_hudson_fst(population_1[:,i], population_2[:,i]) for i in 1:size(population_1,2)]
end


# helper function to do the math for Hudson FST on a locus
function _hudson_fst(pop1::T, pop2::T) where T<:GenoArray
    p1_frq = PopGen.allele_freq(pop1)
    p2_frq = PopGen.allele_freq(pop2)
    # find the shared allele(s) and choose one of them to be "P"
    # this is a safeguard if one population is completely homozygous for an allele
    p_allele = intersect(keys(p1_frq), keys(p2_frq)) |> first
    p1 = p1_frq[p_allele]
    q1 = 1.0 - p1
    p2 = p2_frq[p_allele]
    q2 = 1.0 - p2 
    # degrees of freedom is the number of alleles - 1
    df1 = (count(!ismissing, pop1) * 2) - 1
    df2 = (count(!ismissing, pop2) * 2) - 1
    numerator = (p1 - p2)^2 - (p1*q1/df1) - (p2*q2/df2)
    denominator = (p1*q2) + (p2*q1)
    return numerator/denominator
end