suppressPackageStartupMessages({
  library("SingleR")
  library("Seurat")
  library("data.table")
  library("doFuture")
})

registerDoFuture();
options('future.globals.maxSize' = 10 * 1024^3)

getdata <- function(...)
{
  e <- new.env()
  name <- data(..., envir = e)[1]
  e[[name]]
}

load("test/pbmc68k_subset_seuratV3.rda")

singleRNumCpu = 8
singleRefSpecies = "";
refdata =  getdata(list = "blueprint_encode", package = "SingleR")
sc.data = GetAssayData(object = slimobj, slot = "counts", assay = "RNA");



# "fineTure_0.05_doFuture_1CPU"
print("result_1 => fine tune thres:0.05, doFuture: 1CPU")
print(system.time({
  plan(multiprocess, workers = 1)
  result_1 = SingleR::SingleR(
    "single",sc.data,refdata$data,types=refdata$types,
    sd.thres = refdata$sd.thres,genes = refdata$de.genes,
    fine.tune.thres = 0.05,numCores = 4, step = 10000
  )
}))

print("result_2 => fine tune thres:0.01, doFuture: 1CPU")
print(system.time({
  plan(multiprocess, workers = 1)
  result_2 = SingleR::SingleR(
    "single",sc.data,refdata$data,types=refdata$types,
    sd.thres = refdata$sd.thres,genes = refdata$de.genes,
    fine.tune.thres = 0.01,numCores = 4, step = 10000
  )
}))

print("result_3 => fine tune thres:0.01, doFuture: 4CPU")
print(system.time({
  plan(multiprocess, workers = 4)
  result_3 = SingleR::SingleR(
    "single",sc.data,refdata$data,types=refdata$types,
    sd.thres = refdata$sd.thres,genes = refdata$de.genes,
    fine.tune.thres = 0.01,numCores = 4, step = 1000
  )
}))
print("Adjusted Rand Index between result_1 and result2:")
print(mclust::adjustedRandIndex(result_1$labels, result_2$labels))
print("Adjusted Rand Index between result_2 and result3:")
print(mclust::adjustedRandIndex(result_2$labels, result_3$labels))