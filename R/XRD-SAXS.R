library(tomio.signals)

subtract.background <- function(signal,
                                bg,
                                f= NULL,
                                normalize= FALSE,
                                I.N= 3,
                                new.q = TRUE) {
    #' subtract a (measured) background signal from the signal
    #' @details
    #' Subtract two signals from one another: the actual signal
    #' and one measured on the background (empty sample holder)
    #' Normally do a 1:1 subtraction, or scale the signals to one-another
    #' If normalize is set, then use the zero intensities to scale tha background
    #' to the signal. (I - I.0/I.bg.0 * bg)
    #' If 'f' as for factor is set then use I - f*bg
    #' If normalize is set, f is ignored.
    #'
    #' @param signal    an X-ray list, where data has I for intensities
    #' @param bg    an X-ray list, where data has I for background
    #' @param f     float, scalar factor for scaling the background to the signal
    #' @param normalize Boolean if TRUE, normalize to the zero intensities
    #' @param I.N   integer, how many points of the starting part of
    #'                          intensities to use in calculating the zero
    #'                          intensity? I.0 = mean(I[1:I.N])
    #' @param new.q Boolean if we should identify the valid q-range, as where
    #'                          the corrected signal is first above the background
    #'
    #' @return  an X-ray list data$I replaced with the subtracted signal
    #'
    #' @export

    if (is.null(signal$data) || is.null(bg$data) ||
        is.null(signal$data$I) || is.null(bg$data$I)) {
        cat('Incorrect list type for signal or background!\n')
        return();
    }
    q <- signal$data$q
    sig <- signal$data$I
    backg <- bg$data$I

    if (normalize) {
        I.0 <- mean(sig[1:I.N])
        b.0 <- mean(backg[1:I.N])
        if (b.0 > 0) {
            backg <- I.0/b.0 * backg
        } else {
            cat('Invalid background zero intensity!', bg, '\n')
            return()
        }

    } else if (!is.null(f) && f > 0) {
        backg <- f*backg
    }

    # we will store the data in a copy of sig,
    res <- signal
    sig.bg <- sig - backg
    res$data$I <- sig.bg
    N <- length(sig)

    if (new.q) {
        # now, estimate the validity of the data
        #
        # the minimum q is where the signal is above the background
        indx <- which(sig.bg > backg)
        if (length(indx) > 0) {
        } else {
            cat('subtracted signal remains below bg\n')
            cat('Try signal above background\n')
            indx <- which(sig > backg)
            print(indx)
            if (length(indx) < 1) {
                # most probably trying to subtract background from itself
                # or swapped curve to bg.
                cat('really bad, set i.0 to 1\n')
                indx <- 1
            }
        }
        i.0 <- min(indx)
        res$data <- res$data[i.0:N,]
    }
    with(signal$data,
         plot(q, I,
              xlab='q, 1/nm',
              ylab='intensity, AU',
              type='l',
              log='xy',
              lwd=2)
         )
    with(bg$data, lines(q, I, col='red', lwd=2))
    with(res$data, lines(q, I, col='blue', lwd=2))
    legend('topright', legend=c('signal', 'background', 'subtracted'),
           fill=c('black', 'red', 'blue'))

    return(res)
}


plot.xray <- function(sp.1, sp.2= NULL,
                      theta= FALSE,
                      drop= 3,
                      add= FALSE,
                      type= 'l',
                      xlim= NULL,
                      log= NULL,
                      lwd= 2,
                      ...) {
    #' Plot a SAXS and/or WAXS (XRD) data set
    #' @details
    #' Plot together the spectra of two halves, ideally one is WAXS
    #' the other SAXS.
    #' If the second set is NULL, plot only a single curve, like
    #' the XRD data from the XRDML reader.
    #'
    #' In both cases take the theta or q.A.1 signals for the X-axis.
    #'
    #' @param sp.1 X-ray list objects
    #' @param sp.2 X-ray list objects, optional. If present the
    #'              plot is the merge of the two curves together
    #' @param theta     Boolean, whether plot 2 theta axis for x
    #'                  if TRUE, log= 'y' else log='xy' by default
    #'
    #' @param drop  integer, how many points to drop from the end
    #'                       of the first curve and the beginning
    #'                       of the second one
    #'                       ignored for single curves
    #' @param add       Boolean, if TRUE, do not make a new plot
    #' @param type string, the type parameter of plot
    #' @param xlim array, the xlim parameter of plot
    #' @param log string, the log parameter of plot
    #' @param lwd integer, the line width parameter of plot
    #' @param ... other parameters are passed to plot
    #'
    #' @return Nothing
    #'
    #' @export

    data.1 <- sp.1$data
    N <- length(data.1$I)

    if (!is.null(sp.2)) {
        data.2 <- sp.2$data
    } else {
        data.2 <- list(theta= NULL, q= NULL, I= NULL)
        drop <- 0
    }

    if (drop > 0){
        N <- N - drop
    }

    if (theta) {
        if (drop > 0) {
            x <- c(data.1$theta[1:N], data.2$theta[-c(1:drop)])
        } else {
            x <- c(data.1$theta, data.2$theta)
        }
        xlabel <- expression(paste('2', Theta, sep=''))
        logset = 'y'
    } else {

        # if drop is 0, then -c(0:1) still drops the first point
        if (drop > 0) {
            x <- c(data.1$q[1:N], data.2$q[-c(1:drop)])
        } else {
            x <- c(data.1$q, data.2$q)
        }

        xlabel <- 'q, 1/nm'
        if (!is.null(sp.1$type)){
            # we do log-log for SAXS
            if (sp.1$type == 'SAXS') {
                logset= 'xy'
            } else {
                # else it is WAXS or XRD
                logset= ''
            }
        }
    }
    # depending on how the dat file was created,
    # intensity may have a key as I.q or I.abs.units

    if (drop > 0){
        y <- c(data.1$I[1:N], data.2$I[-c(1:drop)])
    } else {
        y <- c(data.1$I, data.2$I)
    }
    if (!is.null(sp.2)){
        # we also have to sort this data
        # for clear plotting
        s.indx <- order(x)
        x <- x[s.indx]
        y <- y[s.indx]
    }

    if (!is.null(log)) {
        logset= log
    }

    if (is.null(xlim)){
        xrange <- range(x)
    }

    if (add) {
        cat('adding', length(x), length(y), 'in range', min(x), max(x), '\n')
        lines(x, y, type= type, lwd= lwd, ...)
    } else {
        cat('plotting', length(x), length(y), 'points as', as.character(xlabel), '\n')
        plot(x, y,
            type= type, xlim= xlim, log= logset,
            xlab= xlabel, ylab= 'intensity, counts',
            lwd= lwd, ...)
    }
    invisible(list(x= x, y=y))
}


plot.xray.group <- function(pattern,
                            complement= TRUE,
                            cols= NULL,
                            alpha= 1,
                            legends= NULL,
                            envir= .GlobalEnv,
                            ...) {
    #' take a variable name pattern and plot them as
    #' X-ray data sets using plot.xray().
    #'
    #' @details
    #' Plot based on a search pattern, where the found data are all X-ray
    #' scattering curves in the format provided by the reader files.
    #' Especially, have a $data subset which contains q/theta and intensity
    #' columns, suitable for plot.xray().
    #'
    #' @param pattern   string, a pattern to get the variables via ls()
    #'                  if an array of strings is provided, it is considered
    #'                  the names of variables to load (or a hit list)
    #' @param complement Boolean, find the complementer data set
    #' @param cols  string array, name of colors to be used
    #' @param alpha numeric 0 - 1, a transparency value for generating colors
    #' @param legends    string array, legends to be used
    #' @param envir     environment, where to find the data
    #' @param ... all further arguments are passed to plot.xray
    #'
    #' @export

    if (length(pattern) > 1) {
        lst <- pattern
    } else {
        lst <- ls(pattern= pattern, envir= envir)
    }

    if (length(lst) < 1) {
        cat('pattern', pattern, 'not found\n')
        return()
    }

    cat('found', length(lst), 'data sets\n')

    if (is.null(cols)) {
        N <- length(lst)
        if (N > 5) {
            cols<- hcl.colors(N, 'RdYlBu', alpha= alpha)
        } else {
            # add simple colors and alpha channel:
            # only rgb can add the alpha channels,
            # so we convert first to rgb then back
            cols <- col2rgb(c('black', 'blue', 'green', 'red', 'orange')[1:N])/255
            # rgb needs red, green, blue parameters separately
            cols =  rgb(cols[1,], cols[2,], cols[3,], alpha= alpha)
        }
    }
    # cat('color values:', cols, '\n')

    not.first <- FALSE
    for (i in seq_along(lst)) {
        n <- lst[i]
        a <- get(n, envir= envir)

        n.2 <- NULL
        if (complement) {
            if (grepl('.*\\.0\\.', n)) {
                n.2 <- gsub('\\.0\\.', '.1.', n)
            } else if (grepl('.*\\.1\\.', n)){
                n.2 <- gsub('\\.1\\.', '.0.', n)
            }
        }
        # cat('current names', n, '2:', n.2, '\n')
        if (!is.null(n.2) && exists(n.2, envir= envir)) {
            b <- get(n.2, envir= envir)
            # cat('plotting', n, 'c:', cols[i], '\n')
            plot.xray(a, b, add= not.first, col=cols[i], ...)
        } else {
            # cat('plotting', n, 'c:', cols[i], '\n')
            plot.xray(a, add= not.first, col=cols[i], ...)
        }
        not.first <- TRUE
    }
    if (!is.null(legends)) {
        legend('topright', fill= cols, legend= legends, text.width= NA, cex= 0.8)
    }
}


find.xray.peak <- function(data,
                           sigma= 20,
                           q.min = NULL,
                           q.max = NULL,
                           max.I = NULL,
                           log.I = FALSE,
                           bg.I  = FALSE,
                           peak.threshold= 0.05,
                           peak.window= 2,
                           bg.window= 80,
                           verbose= FALSE
                           ) {
    #' find peaks using the Gaussian find.peaks function from tomio.signals
    #'
    #' @details
    #' Find the peaks using the sign change of the first derivative, employing
    #' the find.peaks() from tomio.signals package.
    #' From the X-ray data set, use the $data$q and $data$I arrays.
    #'
    #' @param data  an X-ray data list
    #' @param sigma standard deviation of the Gaussian filter
    #'              used for differentiation and smoothing
    #'              (this is in index values)
    #' @param q.min if set, keep only q > q.min
    #' @param q.max if set, keep only q < q.max
    #' @param max.I cut the intensities to be below this value
    #'              useful for SAXS data, where the starting
    #'              part can span 6 orders of magnitude
    #'
    #' @param log.I Boolean if set, take log(I) instead of I
    #'
    #' @param bg.I  Boolean if set, subtract a background from I
    #'              using the generic background function from tomio.signals
    #'              with a window parameter of bg.window
    #'
    #' @param peak.threshold    relative threshold for peaks
    #'              passed to find.peaks
    #'
    #' @param peak.window   an index value; how many points
    #'              to be used around the peak to refine the
    #'              position (passed to find.peaks)
    #'
    #' @param bg.window     the window parameter for background
    #'              subtraction, if bg.I is set.
    #'
    #' @param verbose   Boolean, if TRUE, plot the results
    #'
    #' @return  data frame     containing the q,
    #'                                        theta,
    #'                                        2theta,
    #'                                        d=2pi/w and peak intensities
    #'
    #' @export

    this.data <- data$data
    wavelength <- data$wavelength # in nm hopefully
    if (is.null(this.data$I) || is.null(this.data$q)){
        cat('one of the arrays is missing!\n')
        return();
    }

    if (!is.null(q.min)) {
        this.data <- this.data[this.data$q > q.min,]
    }

    if (!is.null(q.max)) {
        this.data <- this.data[this.data$q < q.max,]
    }

    if (!is.null(max.I)) {
        indx <- this.data$I < max.I
        this.data <- this.data[indx,]

    }

    if (log.I){
        this.data$I <- log(this.data$I)
        this.data <- this.data[!is.na(this.data$I),]
    }

    if(nrow(this.data) < 2) {
        cat('No data remained after filtering!\n')
        return();
    }

    if (verbose == TRUE) {
        plot(I ~ q, data= this.data,
             type='o',
             xlab='q, 1/nm',
             ylab='intensity, au')
    }
    if (bg.I) {
        bg <- with(this.data, background(I, bg.window))
        this.data$I <- this.data$I - bg
        if (verbose) {
            lines(this.data$q, bg, col='green', lwd= 2)
        }
    }

    q.peaks <- find.peaks(this.data[,c('q','I')],
                          sigma= sigma,
                          peak.threshold= peak.threshold,
                          peak.window= peak.window)
    d.peaks <- 2*pi/q.peaks
    I.peaks <- rep(0, length(q.peaks))
    # q <- 4 pi / lambda sin(theta)
    # theta.peaks <- asin(q.peaks/(4*pi)*wavelength)*180.0/pi
    theta.peaks <- asin(wavelength/(2*d.peaks))*180.0/pi

    for (i in seq_along(q.peaks)) {
        q.i <- q.peaks[i]
        indx <- which(this.data$q == q.i)
        if (length(indx) > 0) {
            I.peaks[i] <- this.data$I[indx][1]
        } else {
            I.min <- tail(this.data$I[this.data$q < q.i], 1)
            I.max <- this.data$I[this.data$q > q.i][1]
            I.peaks[i] <- mean(c(I.min, I.max))
        }
    }

#    cat('w:', wavelength,'\n')
#    cat('q:', q.peaks, '\n')
#    cat('d:', d.peaks, '\n')
#    cat('I:', I.peaks, '\n')
#    cat('t:', theta.peaks, '\n')

    if (verbose) {
        abline(v = q.peaks, col='red')
        legend('topright', fill='black', legend= data$SampleName)
    }

    return(data.frame(
                q= q.peaks,
                theta= theta.peaks,
                theta.2= 2*theta.peaks,
                d= d.peaks,
                I = I.peaks
                )
    )
}


# for backwards compatibility
#plot.xenocs <- plot.xray
# plot.xenocs.group <- plot.xray.group

read.xray.files <- function(folder='./',
                           filetype= c('dat', 'xrdml'),
                           prefix= '',
                           remove= '',
                           recursive= FALSE,
                           envir= .GlobalEnv
                           ) {
    #' read all data files from a folder and put them into the
    #' main memory.
    #'
    #' @details
    #' Run in a similar manner as read.all.files in tomio.signals
    #' Turn the file names to variable names such:
    #' - remove dates as [0-9]\{4\}-\[0-9\]\{2\}-\[0-9\]\{2\}
    #' - convert - and _ characters to dots
    #' - add a prefix with a '.' before
    #'
    #' @param folder    text, where to search for the files
    #' @param filetype  text, one of dat or xrdml
    #' @param prefix    text, if non-empty, prepend to variable
    #'                      names with a dot
    #' @param remove    text, remove this pattern from the file name
    #'                      during conversion to variable name
    #' @param recursive Boolean, whether the searhc should be recursive
    #' @param envir     environment where to dump the results
    #'
    #' @return nothing
    #'
    #' @export

    # ensure file type can be only 'dat' or 'xrdml'
    filetype= match.arg(filetype)

    lst <- dir(folder,
               pattern= paste('\\', filetype, sep='.'),
               recursive= recursive,
               full.names= TRUE)

    if (length(lst) == 0) {
        cat('No files found!\n')
        return()
    }

    # do the clean up in steps
    # first the extension,
    # then custom content
    namelst <- gsub(paste('\\', filetype, sep='.'), '',
                              basename(lst)
                    )

    if (remove != '') {
        namelst <- gsub(remove, '', namelst)
        cat('after removing', remove,'we have', namelst, '\n')
    }

    # then standard also converting '-' and '_' to '.'.
    namelst <- gsub('[_-]+', '.',
                    # kill the date:
                    gsub('[0-9]{4}-[0-9]{2}-[0-9]{2}[_-]+', '',
                         namelst
                         )
                    )
    namelst <- gsub('\\.{2,}', '.', namelst)

    cat('Starting namelist:', namelst, '\n')
    if (prefix !='' && !is.null(prefix)){
        namelst <- paste(prefix, namelst, sep='.')
    }

    cat('Loading files to\n')
    cat(namelst, sep='\n')

    if (filetype == 'xrdml') {
        reader = read.xrdml
    } else {
        reader = read.xenocs
    }

    # load the files
    for (i in seq_along(lst)){
        a <- reader(lst[i])
        assign(namelst[i], a, envir= envir)
    }
    return()
}


write.xray <- function(x, file, sep=',') {
    #' use write.table to export the intensity table of an x-ray record
    #' @details
    #' this function just dumps I, dI, q, theta, theta.rad and d into a text file,
    #' not exporting the metadata.
    #' It is envisioned as an easy export for further plotting or processing of
    #' scattering curves.
    #'
    #' @param x the variable to be exported
    #' @param file the file name or pointer to be exported to
    #' @param sep   separator character passed to write .table
    #'
    #' @return Nothing
    #'
    #' @export

    if (is.null(x$data)) {
        cat('Invalud data, not data subset is found\n')
        return()
    }

    with(x, write.table(data, file= file, sep=sep, row.names= FALSE, col.names= TRUE))
    return()
}
