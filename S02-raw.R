## Raw MS data

## binary - vendor-specific formats
## open formats: mzML, mzXML

## proteowizard msconvert: https://proteowizard.sourceforge.io/
##
## ThermoRawFileParser: https://github.com/compomics/ThermoRawFileParser

library(Spectra)

## data.frame
## tibble

spd <- DataFrame(msLevel = c(1L, 2L),
                 rtime = c(1.1, 1.2))

spd$mz <- list(
    c(100, 103.3, 132, 210),
    c(45, 100, 200)
)

spd$intensity <- list(
    c(45, 12, 345, 20),
    c(45, 122, 12)
)

sp <- Spectra(spd)

sp

## - spectraVariables() and spectraData()
## - peaksData()
## - sp[]

spectraVariables(sp)

spectraData(sp)

peaksData(sp)[[1]]

peaksData(sp)[[2]]

sp[c(1, 2, 1, 1)]


sp <- Spectra(f)

length(sp)
spectraVariables(sp)
pd <- peaksData(sp)

spectraVariables(sp)

msLevel(sp)

sp$msLevel

msLevel(sp)[[1234]]

plot(pd[[1234]], type = "h")


## How many MS level are there, and how many scans of each level?
table(msLevel(sp))


filterMsLevel(sp, 2L)

sp[msLevel(sp) == 2L]

## Extract the index of the MS2 spectrum with the highest base peak
## intensity.

sp2 <- filterMsLevel(sp, 2L)

sp2[which.max(sp2$basePeakIntensity)]


plotSpectra(sp2[4192])
plotSpectra(sp2[1234])
plotSpectra(sp2[1230])
plotSpectra(sp2[1230:1233])
plotSpectra(sp[1])
plotSpectra(sp[1:4])


## The chromatogram can be created by extracting the totIonCurrent and
## rtime variables for all MS1 spectra. Annotate the spectrum of
## interest.

spectraVariables(sp)

## plot(..., type = "l") ## line plot
## plot(..., type = "h") ## 'histogram' plot


plot(rtime(sp), tic(sp), type = "l")
plot(sp$rtime, sp$totIonCurrent, type = "l")


sp1 <- filterMsLevel(sp, 1L)
plot(rtime(sp1), tic(sp1), type = "l")
abline(v = rtime(sp)[2807], col = "red")

MsCoreUtils::formatRt(rtime(sp)[2800:2820])

sp[2807]

library(tidyverse)

spectraData(sp) |>
    as.data.frame() |>
    as_tibble() |>
    filter(msLevel == 1) |>
    ggplot(aes(x = rtime,
               y = totIonCurrent)) +
    geom_line()


## The filterPrecursorScan() function can be used to retain a set
## parent (MS1) and children scans (MS2), as defined by an acquisition
## number. Use it to extract the MS1 scan of interest and all its MS2
## children.

sp2 <- filterPrecursorScan(sp, 2807)

## Plot the MS1 spectrum of interest and highlight all the peaks that
## will be selected for MS2 analysis.

plotSpectra(sp2[1], xlim = c(400, 1000))

abline(v = precursorMz(sp2)[-1], col = "grey")

## Use plotSpectra() function to plot all 10 MS2 spectra in one call.

plotSpectra(sp2[-1])

plotSpectra(sp2[2:11])

plotSpectra(filterMsLevel(sp2, 2L))


## Focus of mz range

plotSpectra(sp[2807], xlim = c(521.2, 522.5))

plotSpectra(sp[2807], xlim = c(521.25, 521.4))

par(mfrow = c(2, 1))

## Processing

plotSpectra(sp[2807], xlim = c(521.2, 522.5))
Spectra::pickPeaks(sp[2807]) |>
    filterIntensity(1e7) |>
    plotSpectra(xlim = c(521.25, 522.5))

table(msLevel(sp), centroided(sp))


## More visualisation

plotSpectra(sp2[7],
            xlim = c(126, 132))

mzLabel <- function(z) {
    ## z is an instance of class Spectra
    z <- peaksData(z)[[1L]]
    lab <- format(z[, "mz"], digits = 4)
    lab[z[, "intensity"] < 1e5] <- ""
    lab
}

plotSpectra(sp2[7],
            labels = mzLabel,
            xlim = c(126, 132))


sp2 <- filterMsLevel(sp, 2L)

anyDuplicated(precursorMz(sp2))

i <- which(precursorMz(sp2) == precursorMz(sp2)[37])

plotSpectra(sp2[i])


plotSpectraMirror(sp2[31], sp2[37])

plotSpectraOverlay(sp2[i], col = c("red", "steelblue"))

## BiocManager::install("RforMassSpectrometry/SpectraVis")

library(SpectraVis)

plotlySpectra(sp2[31])

browseSpectra(sp)


BiocManager::install("MsBackendMgf")


(fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE))

basename(fls)

sciex <- Spectra(fls)

dataOrigin(sciex)

table(dataOrigin(sciex))



####################################################3
library(mzR)

Spectra(f)

x <- openMSfile(f)

hd <- header(x) ## like spectraData from Spectra

pk <- mzR::peaks(x) ## like peaksData from Spectra

Spectra(DataFrame(hd))
