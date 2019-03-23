#' loint function
#'
#' This function generates a random raster file via repeated LOESS interpolation
#'
#' @author Marco Zanon , \email{marco[at]zanon.xyz}
#'
#' @param size Raster size
#' @param min_z,mid_z,max_z These parameters give you a bit of control on the initial population of the grid. Initial values are randomly picked among these parameters (NOTE: the final ratser will aways be rescaled between 1 and 100).
#' @param secondary_passes_num Number of secondary passes (rowise and columnwise interpolation across a fully populated grid.
#' @param secondary_pass_type Decides whether performing only rowise passes or columnwise too. Values: 1 only one rowise, 2: rowise and columnwise, 0: no secondary pass (equivalent to secondary_passes_num = 0).
#' @param diagonal_pass_type Decides how and if performing diagonal passes. Values: 1 only in one direction, 2: diagonal and inverted diagonal passes. Diagonal interpolation is quite time consuming, so currelty there is no option to perform more than one.
#' @param primary_span,secondary_span,diagonal_span LOESS smoothing spans for the first grid population, for the secondary pass and for the diagonal pass respectively.
#' @param k Number of rows and columns used for the first interpolation step.
#' @param noise_level Noisiness of the final product.
#'
#' @return Returns a raster layer
#'
#' @import raster
#'
#' @examples
#'
#' #Create a randomly generated raster layer r with a size of 250x250 grid cells.
#' #Then plot it as a hill shade layer using tools provided with the 'raster' package
#'
#' r <- loint(size=250)
#'
#' library(raster)
#' crs(r) <- CRS('+init=EPSG:6933')
#' slope <- terrain(r, opt='slope')
#' aspect <- terrain(r, opt='aspect')
#' hill <- hillShade(slope, aspect, 40, 270)
#' plot(hill, col=grey(0:100/100), legend=FALSE, main=NA)
#'
#' @export

loint <- function(size=250,
                  mid_z=5,
                  min_z=1,
                  max_z=100,
                  secondary_passes_num=1,
                  secondary_pass_type=2,
                  diagonal_pass_type=2,
                  primary_span=0.6,
                  secondary_span=0.4,
                  diagonal_span=0.2,
                  k=4,
                  noise_level=0.1) {


  mx_x_size <- size
  mx_y_size <- size

  mx <- matrix(NA , nrow = mx_x_size, ncol = mx_y_size)

  iterations_done <- NULL

  iterations_remaining <- c(1:mx_y_size)

  message("Run primary passes")

  while (k < mx_x_size) {

    intervals <- ceiling(seq(0,mx_x_size,mx_x_size/k))

    intervals[intervals < 1] <- 1

    if (is.null(iterations_done)) {

      intervals <- intervals

    } else {

      intervals <- intervals[!(intervals %in% iterations_done)]

    }


    iterations_done <- c(iterations_done, intervals)

    iterations_remaining <- iterations_remaining[!(iterations_remaining %in% intervals)]

    for (i_x in intervals) {

      starting_vals <- sample(c(min_z,mid_z,max_z),length(intervals), replace=TRUE)

      lo <- loess(starting_vals ~ intervals, span=primary_span,family = "gaussian")

      pred <- predict(lo, c(1:mx_x_size))

      mx[i_x, ] <- pred

    }


    for (i_y in intervals) {

      mx_df <-as.data.frame(mx)

      mx2 <- na.omit(mx_df[,i_y, drop = FALSE])

      lo <- loess(mx2[,1] ~ as.numeric(rownames(mx2)), span=primary_span,family = "gaussian")

      pred <- predict(lo, c(1:mx_x_size))

      mx[,i_y] <- pred


    }

    k <- k*2

  }

  #cicle through what is left

  for (k in iterations_remaining) {


    mx_df <-as.data.frame(mx)

    mx2 <- na.omit(mx_df[,k, drop = FALSE])

    lo <- loess(starting_vals ~ intervals, span=primary_span,family = "gaussian")

    pred <- predict(lo, c(1:mx_x_size))

    mx[,k] <- pred

  }

  #replace any NA resulting from the interpolation

  mx[is.na(mx)] <- 1

  for (z in c(1:secondary_passes_num)) {

    if (secondary_pass_type == 1 | secondary_pass_type == 2) {
      #rerun the interpolation through every row and column
      message("Run secondary pass rowwise")

      for (w in c(1:mx_x_size)) {


        lo <- loess(mx[w,] ~ c(1:mx_x_size), span=secondary_span,family = "gaussian")

        pred <- predict(lo, c(1:mx_x_size))

        mx[w,] <- pred

      }

      if (secondary_pass_type == 2) {

        message("Run secondary pass columnwise")

        for (w in c(1:mx_y_size)) {


          lo <- loess(mx[,w] ~ c(1:mx_y_size), span=secondary_span,family = "gaussian")

          pred <- predict(lo, c(1:mx_y_size))

          mx[,w] <- pred
        }

      }

    }

  }

  #diagonal pass. might be optional. apply only to vectors with length > number.

  #extract diagonal values

  if (diagonal_pass_type == 1 | diagonal_pass_type == 2) {

    message("Running diagonal pass")

    d <- row(mx) - col(mx)

    s <- split(mx, d)

    for (stack in c(min(d):max(d))) {

      stack_vals <- eval(parse(text=noquote(paste("s$`",(get("stack")),"`", sep=""))))

      if (length(stack_vals) >=10 ) {

        lo <- loess(stack_vals ~ c(1:length(stack_vals)), span=diagonal_span,family = "gaussian")

        pred <- predict(lo, c(1:length(stack_vals)))


        eval(parse(text=noquote(paste("s$`",(get("stack")),"` <- pred", sep=""))))  #do operations on selected diagonal

        mx <- matrix(unsplit(s,d), nrow = mx_x_size, ncol = mx_y_size) #rebuild matrix

      }
    }

    if (diagonal_pass_type == 2) {

      message("Running inverse diagonal pass")

      # reverse diagonal

      r.d <- row(mx) + col(mx)

      r.s <- split(mx, r.d)

      for (stack in c(min(r.d):max(r.d))) {

        stack_vals <- eval(parse(text=noquote(paste("r.s$`",(get("stack")),"`", sep=""))))

        if (length(stack_vals) >=10 ) {

          lo <- loess(stack_vals ~ c(1:length(stack_vals)), span=diagonal_span,family = "gaussian")

          pred <- predict(lo, c(1:length(stack_vals)))


          eval(parse(text=noquote(paste("r.s$`",(get("stack")),"` <- pred", sep=""))))  #do operations on selected diagonal

          mx <- matrix(unsplit(r.s,r.d), nrow = mx_x_size, ncol = mx_y_size) #rebuild matrix
        }
      }
    }

  }



  r <- raster(nrows = mx_x_size, ncols = mx_y_size, xmn = 0, xmx = mx_x_size, ymn = 0, ymx = mx_y_size, vals = mx)

  min_new <- 1
  max_new <- max_z

  min_old <- cellStats(r, "min")
  max_old <- cellStats(r, "max")

  r1 <- ((max_new-min_new)/(max_old-min_old))*(r-max_old)+max_new

  noise_level <- ifelse(noise_level==0, 0.00001, noise_level )

  noise_mx <- matrix(sample(seq((max_z*noise_level)/100,(max_z*noise_level)/10,(max_z*noise_level)/100), mx_x_size*mx_y_size, replace=TRUE) , nrow = mx_x_size, ncol = mx_y_size)
  noise_r <- raster(nrows = mx_x_size, ncols = mx_y_size, xmn = 0, xmx = mx_x_size, ymn = 0, ymx = mx_y_size, vals = noise_mx)

  final_r <- final_r <- sqrt(noise_r)+(r1)

  final_r <- (final_r-cellStats(final_r, "min"))/(cellStats(final_r, "max")-cellStats(final_r, "min"))*100

}
