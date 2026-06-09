library(XML)
read.xrdml <- function(filename) {
    #' read an RDXML file from Malvern Panalytical
    #' @details
    #' read an RDXML file from Malvern Panalytical and turn it
    #' into a list of metadata and data
    #' To be consistent between devices, make a data dataframe
    #' containing: I, theta, q, d for
    #' intensity, 2Theta, q in 1/nm, and characteristic thickness
    #' d in nm (2pi/q)
    #' (theta is actually 2theta, double the scattering angle)
    #'
    #' @param filename  char    the name / path to the file
    #'
    #' @return list with metadata and data in the $data subfield
    #'
    #' @export

    if (!file.exists(filename)) {
        cat('File does not exist!', filename, '\n')
        return(NULL)
    }

    a <- xmlToList(filename)
    cat('loaded file\n')

    meas <- a$xrdMeasurement
    scans <- meas$scan
    # measured positions and intensities are in here:
    data <- scans$dataPoints

    # get the wavelengths of the device:
    tmp <- unlist(meas$usedWavelength)
    lambda <- list(values = as.numeric(tmp[seq(1,8,2)]),
                   units = tmp[seq(2,8,2)])

    indx <- grepl('Angstrom', lambda$units)
    if (sum(indx) > 0) {
        lambda$values[indx] <- lambda$values[indx]/10
        lambda$units[indx] <- 'nm'
        names(lambda$units) <- gsub('\\.attrs\\.', '', names(lambda$units))
    }

    # sample information
    sample = a$sample
    sample$.attrs <- NULL
    sample$type <- a$sample$.attrs[1]


    cat('defining results\n')
    results = list(sample= sample)
    results$lambda <- lambda

    if (!is.null(sample$id)){
        results$SampleName <- sample$id
    }

    for (i in grep('position', names(data))) {
        tmp <- data[[i]]
        n.pos <- tmp$.attr['axis']
        n.unit <- tmp$.attr['unit']

        if(grepl('startPosition', names(tmp))[1] == TRUE) {
            axis.start <- as.numeric(tmp[[1]])
            axis.end <- as.numeric(tmp[[2]])
        } else {
            axis.start <- as.numeric(tmp[[1]])
            axis.end <- axis.start
        }
        results[[n.pos]] <- list(start=axis.start, end=axis.end, unit = n.unit)
    }

    results$comment <- unlist(a$comment)
    results$program.file <- meas$comment
    results$time.start <- scans$header$startTimeStamp
    results$time.end <- scans$header$endTimeStamp
    results$author <- scans$header$author
    results$scan.axis <- scans$.attrs['scanAxis']
    results$scan.mode <- scans$.attrs['mode']
    results$scan.status <- scans$.attrs['status']
    results$scan.append <- as.numeric(scans$.attrs['appendNumber'])
    cat('scan appended\n')
    results$software <- with(scans$header$source,
                             list(
                             name = applicationSoftware$text,
                             version = applicationSoftware$.attrs['version'],
                             control.name = instrumentControlSoftware$text,
                             control.version= instrumentControlSoftware$.attrs['version'],
                             instrument.ID = instrumentID
                                )
                             )

    results$incident.beam.path <- meas$incidentBeamPath
    results$diffracted.beam.path <- meas$diffractedBeamPath
    results$measurement.type <- meas$.attrs['measurementType']
    results$status <- meas$.attrs['status']
    results$mode <- meas$.attr['sampleMode']
    results$countingTime <- list(time=as.numeric(data$commonCountingTime$text),
                                 unit= data$commonCountingTime$.attrs)

    results$data <- data.frame(
                               I = as.numeric(
                                              unlist(
                                                     strsplit(data$intensities$text, ' ')
                                                     )
                               )
                            )

    # results$I <- as.numeric(unlist(strsplit(data$intensities$text, ' ')))
    results$I.unit <- data$intensities$.attr['unit']


    # provide an average wavelength
    results$wavelength <- mean(lambda$values[1:2])

    # some cleaning up
    if (!is.null(results[['2Theta']])) {
            results$theta.2 <- results[['2Theta']]
            results[['2Theta']] <- NULL

            results$data$theta <- with(results,
                                       seq(theta.2$start,
                                         theta.2$end,
                                         length = length(data$I)
                                          )
                                       )
            results$data$theta.rad <- with(results$data, pi*theta/180)

            stheta <- with(results$data, sin(theta.rad/2))

            # l.mean <- mean(results$lambda$values[1:2])
            # results$data$q <- 4*pi/l.mean*stheta
            results$data$q <- 4*pi/results$wavelength*stheta

            # 2 pi / q simplifies to
            results$data$d <- (results$wavelength/2)/stheta
            results$d.unit <- lambda$units[1]
            results$q.unit <- paste('1', lambda$units[1], sep='/')
    }

    results$type= 'XRD'
    return(results)
}
