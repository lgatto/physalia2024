library(msdata)
library(rpx)
library(Spectra)
library(PSMatch)
library(tidyverse)


idf <- ident(full.names = TRUE)
basename(idf)

id <- PSM(idf)
id
dim(id)

names(id)

head(id$spectrumID)
head(id$spectrumFile)
head(id$sequence)
head(id$DatabaseAccess)

## Verify that this table contains 5802 matches for 5343 scans and
## 4938 peptides sequences.

length(unique(id$spectrumID))
length(unique(id$sequence))
length(unique(id$DatabaseAccess))

as.data.frame(id) |>
    as_tibble() |>
    group_by(isDecoy) |>
    summarise(ms = mean(MS.GF.RawScore))


## Target: protein.fasta -> peptide
## Decoy:  reverse protein.fasta -> peptide

table(id$isDecoy)

## Compare the distribution of raw identification scores of the decoy
## and non-decoy hits. Interpret the figure.

as.data.frame(id) |>
    as_tibble() |>
    ggplot(aes(x = MS.GF.RawScore,
               colour = isDecoy)) +
    geom_density()


table(table(id$spectrumID))

i <- which(id$spectrumID == "controllerType=0 controllerNumber=1 scan=1774")

id[i, ] |>
    as.data.frame() |>
    DT::datatable()

id2 <- reducePSMs(id, id$spectrumID)

id2

(j <- which(id2$spectrumID == "controllerType=0 controllerNumber=1 scan=1774"))

id2[j, "DatabaseAccess"]

## Filtering

idtbl <- as_tibble(id)

## - Remove decoy hits

idtbl <- idtbl |>
    filter(!isDecoy)

## - Keep first rank matches

idtbl <- idtbl |>
    filter(rank == 1)

## - Remove shared peptides. Start by identifying scans that match
##   different proteins. For example scan 4884 matches proteins
##   XXX_ECA3406 and ECA3415. Scan 4099 match XXX_ECA4416_1,
##   XXX_ECA4416_2 and XXX_ECA4416_3. Then remove the scans that match
##   any of these proteins.

mltm <- group_by(spectrumID) |>
    mutate(nProts = length(unique(DatabaseAccess))) |>
    filter(nProts > 1) |>
    pull(spectrumID)

mltm

idtbl <- idtbl |>
    filter(!spectrumID %in% mltm)

idtbl

idf <- filterPSMs(id)

idf <- id |>
    filterPsmDecoy() |>
    filterPsmRank()

## PSMatch vignette: Understanding protein groups with adjacency
##                   matrices

data.frame(idf[1:10, c("sequence", "DatabaseAccess")])

data.frame(idf) |>
    as_tibble() |>
    filter(DatabaseAccess == 'ECA2006')

data.frame(idf) |>
    as_tibble() |>
    filter(sequence == 'RQCRTDFLNYLR')


adj <- makeAdjacencyMatrix(idf)
adj[1:15, 1:5]

dim(adj)

describePeptides(idf)

describeProteins(idf)

cc <- ConnectedComponents(adj)

connectedComponents(cc, 1)

connectedComponents(cc, 527)

connectedComponents(cc, 38)

connectedComponents(cc, 920)

i <- which(nrows(cc) > 2 & ncols(cc) > 2)
dims(cc)[i, ]

cx <- connectedComponents(cc, 1082)


plotAdjacencyMatrix(cx)

## https://lgatto.github.io/2023_06_15_CSAMA_Brixen


## Combining raw and id data


sp <- Spectra(f)
spectraVariables(sp)

head(sp$spectrumId)

idf <- filterPSMs(id)
names(idf)

head(idf$spectrumID)

table(table(idf$spectrumID))


which(table(idf$spectrumID) == 4)

idf[idf$spectrumID == "controllerType=0 controllerNumber=1 scan=5490", ] |>
    as.data.frame() |>
    DT::datatable()

idf <- reducePSMs(idf, idf$spectrumID)



spid <- joinSpectraData(sp, idf,
                        by.x = "spectrumId",
                        by.y = "spectrumID")

spectraVariables(spid)

all(is.na(filterMsLevel(spid, 1L)$sequence))


table(is.na(filterMsLevel(spid, 2L)$sequence))

## Visualise MS2 scans

i <- which(spid$MS.GF.RawScore > 100)[1]

plotSpectra(spid[i])

spid[i]$sequence

calculateFragments("THSQEEMQHMQR")

mz(spid[i])

mz(spid[i])[[1]]

pdi <- data.frame(peaksData(spid[i])[[1]])

pdi$label <- addFragments(spid[i])

addFragments(spid[i])

plotSpectra(spid[i], labels = addFragments,
            labelCol = "steelblue",
            labelPos = 3)

filter(pdi,
       !is.na(label)) |>
    arrange(intensity)

spid[i] |>
    filterIntensity(200) |>
    plotSpectra(labels = addFragments,
                labelCol = "steelblue",
                labelPos = 3)


spid <- countIdentifications(spid)

table(msLevel(spid),
      spid$countIdentifications)

spid |>
    filterMsLevel(1) |>
    spectraData() |>
    as_tibble() |>
    ggplot(aes(x = rtime,
               y = totIonCurrent)) +
    geom_line() +
    geom_point(aes(colour =
                       ifelse(countIdentifications == 0,
                              NA, countIdentifications)),
               size = 4) +
    labs(colour = "Number of ids")

## Comparing spectra - distances


## - Create a new Spectra object containing the MS2 spectra with
##   sequences "SQILQQAGTSVLSQANQVPQTVLSLLR" and
##   "TKGLNVMQNLLTAHPDVQAVFAQNDEMALGALR".


k <- which(spid$sequence %in% c("SQILQQAGTSVLSQANQVPQTVLSLLR", "TKGLNVMQNLLTAHPDVQAVFAQNDEMALGALR"))

spk <- spid[k]

plotSpectra(spk)

spk$sequence

## - Calculate the 5 by 5 similarity matrix between all spectra using
##   compareSpectra(). See the ?Spectra man page for details. Draw a
##   heatmap of that matrix, for example pheatmap::pheatmap()

mat <- compareSpectra(spk)
colnames(mat) <- rownames(mat) <- strtrim(spk$sequence, 2)
mat

library(pheatmap)
pheatmap(mat)


plotSpectraMirror(spk[1], spk[2])

plotSpectraOverlay(spk[3:5],
                   col = c("red", "steelblue", "green"))

plotSpectraMirror(spk[3], spk[4])

## Recap exercise

## Download the 3 first mzML and mzID files from the PXD022816
## project (Morgenstern, Barzilay, and Levin 2021).

library(rpx)

px <- PXDataset("PXD022816")
pxfiles(px)
pxtax(px)
pxref(px)

fmzml <- pxget(px, grep("mzML", pxfiles(px), value = TRUE)[1:3])
basename(fmzml)

fid <- pxget(px, grep("mzID", pxfiles(px), value = TRUE)[1:3])
basename(fid)

## Generate a Spectra object and a table of filtered PSMs. Visualise
## the total ion chromatograms and check the quality of the
## identification data by comparing the density of the decoy and
## target PSMs id scores for each file.

sp <- Spectra(fmzml)


sp

sp$file <- basename(dataOrigin(sp))

table(basename(dataOrigin(sp)), msLevel(sp))

filterMsLevel(sp, 1) |>
    spectraData() |>
    as_tibble() |>
    ggplot(aes(x = rtime,
               y = totIonCurrent,
               colour = file)) +
    geom_line() +
    facet_wrap(~ file)


id <- PSM(fid)
id$file <- sub("^.+QEP2", "QEP2", id$spectrumFile)

data.frame(id) |>
    as_tibble() |>
    count(file)


data.frame(id) |>
    as_tibble() |>
    ggplot(aes(x = MetaMorpheus.score,
               colour = isDecoy)) +
    geom_density() +
    facet_wrap(~ file)


summary(id$PSM.level.q.value)

data.frame(id) |>
    as_tibble() |>
    count(file, isDecoy)

idf <- filterPSMs(id)


## Join the raw and identification data. Beware though that the
## joining must now be performed by spectrum ids and by files.

## primary keys

sp$pkey <- paste0(sub("^.+QEP", "QEP", sp$file),
                  sub("^.+scan=", "::", sp$spectrumId))


idf$pkey <- paste0(idf$file,
                   sub("^.+scan=", "::", idf$spectrumID))

idf[which(idf$pkey == idf$pkey[7]), ] |>
    data.frame() |>
    DT::datatable()

table(table(idf$pkey))

sp <- joinSpectraData(sp, idf, by.x = "pkey")

table(filterMsLevel(sp, 2)$file,
      !is.na(filterMsLevel(sp, 2)$sequence))

length(unique(sp$sequence))
