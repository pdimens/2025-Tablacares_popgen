using Base: project_file_name_uuid
using PopGen, DataFrames, CSV, NaturalSort, GeneticVariation

# the genepop of haplotypes (probably need to remake this)
x = genepop("../../all_samples_haplotyped/YFT_maf01haplo.gen")

miss = missing_data(x)

kin = CSV.File("kin_pcrelate.txt") |> DataFrame

offenders = Vector{String}(undef, 3)
for i in 1:3
    tmp = sort(filter(j -> j.name in [kin.ID1[i], kin.ID2[i]], miss), :missing)
    offenders[i] = tmp.name[end]
end

newdata = omit(x, name = offenders)

genepop(newdata, filename = "YFT.pcrelate.kinrm.haplo.gen")


# the vcf of SNPs
offenders_manual = ["IVC_2708", "VZ_2606", "LA_2348"]

y = vcf("../../all_samples/data/out.14.recode,maf01.vcf", rename_loci = true)

# generate population names to overwrite missings fields
popnames = getindex.(split.(samples(y), "_"),1)
replace!(popnames, "MSAL" => "GOM", "LA" => "GOM", "TXL" => "TX")
populations!(y, String.(y.meta.name), String.(popnames))

# naturally sort loci names (convenience)
sort!(y.loci, [:name, :locus], lt = natural)

# remove the offenders
y_rm = omit(y, name = offenders_manual)

# write to file
genepop(y_rm, filename = "../data/YFT.snp.kinrm.pcrelate.gen")