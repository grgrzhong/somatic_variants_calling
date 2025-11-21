
datafile <- "/lustre1/g/path_my/pipeline/somatic_variants_calling/data/DFSP/FACETS/DFSP-001-T/DFSP-001-T.snp_pileup.gz"

datafile=system.file("extdata", "stomach.csv.gz", package = "facets")
head(read.csv(datafile))[, c(1:2, 5:12)]

library(facets)
library(facetsSuite)

rcmat = readSnpMatrix(datafile)
xx = preProcSample(
    rcmat,
    ndepth = 15,
    het.thresh = 0.25,
    snp.nbhd = 250,
    cval = 25,
    gbuild = "hg38",
    ugcpct = NULL,
    hetscale = TRUE,
    unmatched = FALSE
)

oo = procSample(
    xx,
    cval = 150,
    min.nhet = 15,
    dipLogR = NULL
)
oo$dipLogR

fit= emcncf(
    oo,
    unif=FALSE, 
    min.nhet=15, 
    maxiter=20, 
    eps=1e-3
)

head(fit$cncf)
fit$purity
fit$ploidy
fit$loglik
fit$seglen
fit$cncf
fit$emflags


plotSample(x=oo, emfit=fit)
logRlogORspider(oo$out, oo$dipLogR)

?facetsSuite::gene_level_changes()
