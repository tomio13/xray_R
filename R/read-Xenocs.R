read.xenocs <- function(filename) {
    #' read a xenocs scattering data file
    #'
    #' @details
    #' the Xenocs system exports an angle-averaged scattering
    #' intensity profile into a tabular text file with a large
    #' header. This function is built to import this data for
    #' further analysis.
    #'
    #' ## reading the parameters
    #' the header is a set of commented lines.
    #' These commented lines contain parameters in the form of
    #' # key             value
    #' the value can be text, numeric or boolean
    #'
    #' This script tries to be non-deterministic pulling all
    #' up, and testing conversions.
    #'
    #' Then pull in the intensity table into a dataframe
    #' under name: 'data'. The q values are converted to 1/nm
    #' unit, intensity and error remains as is.
    #' Also generate a column for theta and d equivalent values.
    #' The column theta.rad is twice the scattering angle in radians.
    #' (that is 'theta' is actually 2 theta)
    #'
    #' @param filename  string, the file to be read
    #'
    #' @return a list with all parameters and data
    #'          data is in the $data subfield as a data frame
    #'
    #' @export

    if (!file.exists(filename)) {
        cat('File:', filename, 'not found\n')
        return(NULL)
    }

    res = list()

    fp <- file(filename, 'rt')
    is.head <- TRUE
    # the first line is comments only
    txt.lines <- readLines(fp, -1, ok= TRUE, warn= FALSE)
    head.indx <- grepl('^#', txt.lines)
    i <- 1
    res <- list()

    while(head.indx[i]) {
        # a line is like:
        # text_label                   value
        # so simplify multiple spaces to single ones
        this.line <- gsub('[[:space:]]{2,}', ' ', gsub('#\\s', '', txt.lines[i]))
        # cat(this.line,'\n')

        if (grepl(' ', this.line)[1]) {
            split.line <- strsplit(this.line, ' ')[[1]]
            label <- split.line[1]
            if (label == 'Wavelength') {
                label = 'wavelength'
            }
            # then we may have a text, so merge it back in a nice way
            value <- paste(split.line[-1], collapse=' ')
            # cat(label, ':', value,'dim:', dim(value), 'type:', typeof(value),'\n')

            if (!is.na(value) && (value == 'True' || value == 'False')) {
                value = ifelse('True', TRUE, FALSE)
            } else {
                try(num.value <- suppressWarnings(as.numeric(value)), TRUE)
                # as.numeric throws an error and returns NA
                # for things it cannot interpret
                if(!is.na(num.value)){
                    value <- num.value
                }
            }
            res[[label]] <- value
        }
        i <- i+1
    }

    # Now, pull up the data table in the simple way
    seek(fp, 0L, 'start')
    res$data  <- read.table(fp, header= TRUE, skip= i-1)
    # it tends to have some ending dots in the column name
    # names(res$data) <- gsub('\\.$', '', names(res$data))
    # since sometimes the second and third column has different name,
    # depending how the data file was created...
    names(res$data) <- c('q.A.1', 'I', 'sigma')
    close(fp)

    res$wavelength <- res$wavelength*1E9 # turn to nanometers
    # theta is actually 2 Theta that is 2 * 180/pi asin(q lambda/(4*pi))
    # res$wavelength is now in nm, and q in 1/Angstrom
    res$data$q <- with(res$data, q.A.1*10.0)
    res$data$theta.rad <- with(res$data, 2*asin(q*res$wavelength/(4*pi)))
    res$data$theta <- with(res$data, 180/pi*theta.rad)
    # create a more SI q in nm from the one in Angstroms:
    res$data$d <- with(res$data, 2*pi/q)
    res$data$q.A.1 <- NULL

    # determine type:
    if (res$data$q[1] > 1.0) {
        res$type <- 'WAXS'
    } else {
        res$type <- 'SAXS'
    }

    return(res)
}
