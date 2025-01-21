using PopGen: keep
using PopGen, CSV, DataFrames

data = genepop("YFT.snp.kinrm.pcrelate.gen")

outlier_table = CSV.File("outFLANK_kinless_putatives.csv") |> DataFrame
outliers = getindex.(split.(outlier_table.LocusName, "."),1)

neutral_data = omit(data, locus = outliers)
genepop(filtered_data, filename = "YFT.snp.kinrm.pcrelate.neutral.gen")
structure(neutral_data, filename = "YFT.snp.kinrm.pcrelate.neutral.str", faststructure = true)

outlier_data = PopGen.keep(data, locus = outliers)
genepop(outlier_data, filename = "YFT.snp.kinrm.pcrelate.outlier.gen")
structure(outlier_data, filename = "YFT.snp.kinrm.pcrelate.outlier.str", faststructure = true)
