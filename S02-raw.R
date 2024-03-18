## Raw MS data

## binary - vendor-specific formats
## open formats: mzML

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

## The filterPrecursorScan() function can be used to retain a set
## parent (MS1) and children scans (MS2), as defined by an acquisition
## number. Use it to extract the MS1 scan of interest and all its MS2
## children.

## ?Spectra
