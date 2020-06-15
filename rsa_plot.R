#' @title Plots a response surface of a polynomial equation of second degree
#'
#' @description
#' Plots an RSA object, or a response surface with specified parameters
#'
#' @details
#' Each plot type has its distinctive advantages. The two-dimensional contour plot gives a clear view of the position of the principal axes and the stationary point. The 3d plot gives a three dimensional impression of the surface, allows overplotting of the original data points (in case an RSA object is provided), and allows the interactive adjustment of regression weights in the \code{\link{RSA}} function. The interactive plot allows rotating and exploring a three-dimensional surface with the mouse (nice for demonstration purposes).
#' If you want to export publication-ready plots, it is recommended to export it with following commands:
#' \code{p1 <- plot(r1, bw=TRUE)
#' trellis.device(device="cairo_pdf", filename="RSA_plot.pdf")
#' print(p1)
#' dev.off()}
#'
#' @aliases plotRSA
#'
#' @importFrom aplpack compute.bagplot
#'
#' @export
#' @param x Either an RSA object (returned by the \code{RSA} function), or the coefficient for the X predictor
#' @param y Y coefficient
#' @param x2 X^2 coefficient
#' @param y2 Y^2 coefficient
#' @param xy XY interaction coefficient
#' @param w W coefficient (for (un)constrained absolute difference model)
#' @param wx WX coefficient (for (un)constrained absolute difference model)
#' @param wy WY coefficient (for (un)constrained absolute difference model)
#' @param y3 Y^3 coefficient
#' @param x2y X^2Y coefficient
#' @param xy2 XY^2 coefficient
#' @param x3 X^3 coefficient
#' @param b0 Intercept
#' @param xlim Limits of the x axis
#' @param ylim Limits of the y axis
#' @param zlim Limits of the z axis
#' @param xlab Label for x axis
#' @param ylab Label for y axis
#' @param zlab Label for z axis
#' @param main the main title of the plot
#' @param cex.main Factor for main title size
#' @param surface Method for the calculation of the surface z values. "predict" takes the predicted values from the model, "smooth" uses a thin plate smoother (function \code{Tps} from the \code{fields} package) of the raw data
#' @param lambda lambda parameter for the smoother. Default (NULL) means that it is estimated by the smoother function. Small lambdas around 1 lead to rugged surfaces, big lambdas to very smooth surfaces.
#' @param rotation Rotation of the 3d surface plot (when type == "3d")
#' @param label.rotation Rotation of the axis labls (when type == "3d")
#' @param gridsize Number of grid nodes in each dimension
#' @param bw Print surface in black and white instead of colors?
#' @param legend Print color legend for z values?
#' @param cex.tickLabel Font size factor for tick labels
#' @param cex.axesLabel Font size factor for axes labels
#' @param type \code{3d} for 3d surface plot, \code{contour} for 2d contour plot, "interactive" for interactive rotatable plot. Shortcuts (i.e., first letter of string) are sufficient
#' @param points A list of parameters which define the appearance of the raw scatter points: 
#'	\itemize{
#'		\item data: Data frame which contains the coordinates of the raw data points. First column = x, second = y, third = z. This data frame is automatically generated when the plot is based on a fitted RSA-object
#'		\item show = TRUE: Should the original data points be overplotted?
#'		\item color = "black": Color of the points
#' 		\item value="raw": Plot the original z value, "predicted": plot the predicted z value
#'		\item jitter = 0: Amount of jitter for the raw data points. For z values, a value of 0.005 is reasonable
#'		\item cex = .5: multiplication factor for point size
#' 		\item out.mark = FALSE: If set to TRUE, outliers according to Bollen & Jackman (1980) are printed as red X symbols, but only when they have been removed in the RSA function: \code{RSA(..., out.rm=TRUE)}.
#'			\itemize{
#'				\item If out.rm == TRUE (in RSA()) and out.mark == FALSE (in plotRSA()), the outlier is removed from the model and *not plotted* in plotRSA.
#'				\item If out.rm == TRUE (in RSA()) and out.mark == TRUE (in plotRSA()), the outlier is removed from the model but plotted and marked in plotRSA.
#'				\item If out.rm == FALSE (in RSA()): Outliers are not removed and cannot be plotted.
#'				\item Example syntax: \code{plotRSA(r1, points=list(show=TRUE, out.mark=TRUE))}
#'		}
#'	}
#' As a shortcut, you can also set \code{points=TRUE} to set the defaults.

#' @param model If x is an RSA object: from which model should the response surface be computed?
#' @param demo Do not change that parameter (internal use only)
#' @param fit Do not change that parameter (internal use only)
#' @param param Should the surface parameters a1 to a5 be shown on the plot? In case of a 3d plot a1 to a5 are printed on top of the plot; in case of a contour plot the principal axes are plotted. Surface parameters are not printed for cubic surfaces.
#' @param coefs Should the regression coefficients b1 to b5 (b1 to b9 for cubic models) be shown on the plot? (Only for 3d plot)
#' @param axes A vector of strings specifying the axes that should be plotted. Can be any combination of c("LOC", "LOIC", "PA1", "PA2", "E2", "K1", "K2"). LOC = line of congruence, LOIC = line of incongruence, PA1 = first principal axis, PA2 = second principal axis, E2 = second extremum line in the CA or RRCA model, K1, K2 = boundary lines of the regions of significance in the CL or RRCL model.
#' @param axesStyles Define the visual styles of the axes LOC, LOIC, PA1, PA2, E2, K1, and K2. Provide a named list: \code{axesStyles=list(LOC = list(lty="solid",  lwd=2, col=ifelse(bw==TRUE, "black", "blue"))}. It recognizes three parameters: \code{lty}, \code{lwd}, and \code{col}. If you define a style for an axis, you have to provide all three parameters, otherwise a warning will be shown.
#' @param project A vector of graphic elements that should be projected on the floor of the cube. Can include any combination of c("LOC", "LOIC", "PA1", "PA2", "contour", "points", "E2", "K1", "K2")
#' @param maxlines Should the maximum lines be plotted? (red: maximum X for a given Y, blue: maximum Y for a given X). Works only in type="3d"
#' @param link Link function to transform the z axes. Implemented are "identity" (no transformation; default), "probit", and "logit"
#' @param suppress.surface Should the surface be suppressed (only for \code{type="3d"})? Useful for only showing the data points, or for didactic purposes (e.g., first show the cube, then fade in the surface).
#' @param suppress.box Should the surrounding box be suppressed (only for \code{type="3d"})?
#' @param suppress.grid Should the grid lines be suppressed (only for \code{type="3d"})?
#' @param suppress.ticklabels Should the numbers on the axes be suppressed (only for \code{type="3d"})?
#' @param border Should a thicker border around the surface be plotted? Sometimes this border leaves the surrounding box, which does not look good. In this case the border can be suppressed by setting \code{border=FALSE}.
#' @param contour A list defining the appearance of contour lines (aka. height lines). show=TRUE: Should the contour lines be plotted on the 3d wireframe plot? (Parameter only relevant for \code{type="3d"}). color = "grey40": Color of the contour lines. highlight = c(): A vector of heights which should be highlighted (i.e., printed in bold). Be careful: the highlighted line is not necessarily exactly at the specified height; instead the nearest height line is selected.
#' @param hull Plot a bag plot on the surface (This is a bivariate extension of the boxplot. 50\% of points are in the inner bag, 50\% in the outer region). See Rousseeuw, Ruts, & Tukey (1999).
#' @param showSP Plot the stationary point? (only relevant for \code{type="contour"})
#' @param showSP.CI Plot the CI of the stationary point? (only relevant for \code{type="contour"})
#' @param distance A vector of three values defining the distance of labels to the axes
#' @param tck A vector of three values defining the position of labels to the axes (see ?wireframe)
#' @param pal A palette for shading. You can use \code{\link{colorRampPalette}} to construct a color ramp, e.g. \code{plot(r.m, pal=colorRampPalette(c("darkgreen", "yellow", "darkred"))(20))}. If \code{pal="flip"}, the default palette is used, but reversed (so that red is on top and green on the bottom).
#' @param pal.range Should the color range be scaled to the box (\code{pal.range = "box"}, default), or to the min and max of the surface (\code{pal.range = "surface"})? If set to "box", different surface plots can be compared along their color, as long as the zlim is the same for both.
#' @param pad Pad controls the margin around the figure (positive numbers: larger margin, negative numbers: smaller margin)
#' @param ... Additional parameters passed to the plotting function (e.g., sub="Title"). A useful title might be the R squared of the plotted model: \code{sub = as.expression(bquote(R^2==.(round(getPar(x, "r2", model="full"), 3))))}
#'
#' @references
#' Rousseeuw, P. J., Ruts, I., & Tukey, J. W. (1999). The Bagplot: A Bivariate Boxplot. The American Statistician, 53(4), 382-387. doi:10.1080/00031305.1999.10474494
#' @seealso \code{\link{demoRSA}}, \code{\link{RSA}}
#'
#' @examples
#' # Plot response surfaces from known parameters
#' # example of Edwards (2002), Figure 3
#' \dontrun{
#' # Default: 3d plot:
#' plotRSA(x=.314, y=-.118, x2=-.145, y2=-.102, xy=.299, b0=5.628)
#' # Contour plot:
#' plotRSA(x=.314, y=-.118, x2=-.145, y2=-.102, xy=.299, b0=5.628, type="c")
#' # Interactive plot (try the mouse!):
#' plotRSA(x=.314, y=-.118, x2=-.145, y2=-.102, xy=.299, b0=5.628, type="i")
#'
#' # Plot response surface from an RSA object
#' set.seed(0xBEEF)
#' n <- 300
#' err <- 2
#' x <- rnorm(n, 0, 5)
#' y <- rnorm(n, 0, 5)
#' df <- data.frame(x, y)
#' df <- within(df, {
#' 	diff <- x-y
#' 	absdiff <- abs(x-y)
#' 	SD <- (x-y)^2
#' 	z.diff <- diff + rnorm(n, 0, err)
#' 	z.abs <- absdiff + rnorm(n, 0, err)
#' 	z.sq <- SD + rnorm(n, 0, err)
#' 	z.add <- diff + 0.4*x + rnorm(n, 0, err)
#' 	z.complex <- 0.4*x + - 0.2*x*y + + 0.1*x^2 - 0.03*y^2 + rnorm(n, 0, err)
#' })
#' 
#' r1 <- RSA(z.sq~x*y, df, models=c("SQD", "full", "IA"))
#' plot(r1)	# default: model = "full"
#' plot(r1, model="SQD", points=list(show=TRUE, value="predicted"))
#' }



#b0=0; x=0; y=0; x2=0; y2=0; xy=0; w=0; wx=0; wy=0;  zlim=NULL; xlim=c(-2, 2); ylim=c(-2, 2); rotation=list(x=-45, y=45, z=35); legend=TRUE; cex=1.2; type="3d"; points=TRUE; demo=FALSE; model="full"; 

#b0=-9; x=0; y=0; x2=0; y2=0; xy=0; w=0; wx=1; wy=-1;  zlim=NULL; xlim=c(-2, 2); ylim=c(-2, 2); rotation=list(x=-45, y=45, z=35); legend=TRUE; cex=1.2; type="3d"; points=TRUE; demo=FALSE; model="full"; fit=NULL; link="identity"; param=TRUE; gridsize=21;bw=FALSE; pal=NULL; axes=c("LOC", "LOIC", "PA1", "PA2"); distance=c(1, 1, 1); tck=c(1, 1, 1); xlab="X"; ylab="Y"; zlab="Z"; border=TRUE;

## old rotation
# rotation=list(x=-45, y=45, z=35), label.rotation=list(x=45, y=-25, z=94)
# distance=c(1, 1, 1), tck=c(1, 1, 1)

plotRSA <- function(x=0, y=0, x2=0, y2=0, xy=0, w=0, wx=0, wy=0, x3=0, xy2=0, x2y=0, y3=0, b0=0, 
                    type="3d", model="full", 
                    xlim=NULL, ylim=NULL, zlim=NULL, 
                    xlab=NULL, ylab=NULL, zlab=NULL, main="",
                    surface="predict", lambda=NULL, 
                    suppress.surface=FALSE, suppress.box = FALSE, suppress.grid = FALSE,
                    suppress.ticklabels=FALSE,
                    rotation=list(x=-63, y=32, z=15), label.rotation=list(x=19, y=-40, z=92), 
                    gridsize=21, bw=FALSE, legend=TRUE, param=TRUE, coefs=FALSE,
                    axes=c("LOC", "LOIC", "PA1", "PA2"),
                    axesStyles=list(
                      LOC = list(lty="solid",  lwd=2, col=ifelse(bw==TRUE, "black", "blue")),
                      LOIC= list(lty="solid",  lwd=2, col=ifelse(bw==TRUE, "black", "blue")),
                      PA1 = list(lty="dotted", lwd=2, col=ifelse(bw==TRUE, "black", "gray30")),
                      PA2 = list(lty="dotted", lwd=2, col=ifelse(bw==TRUE, "black", "gray30"))
                    ),
                    project=c("contour"), maxlines=FALSE,
                    cex.tickLabel=1, cex.axesLabel=1, cex.main=1, 
                    points = list(data=NULL, show=NA, value="raw", jitter=0, color="black", cex=.5, out.mark=FALSE),
                    fit=NULL, link="identity", 
                    tck=c(1.5, 1.5, 1.5), distance=c(1.3, 1.3, 1.4), border=FALSE, 
                    contour = list(show=FALSE, color="grey20", highlight = c()),
                    hull=NA, showSP=FALSE, showSP.CI=FALSE, 
                    pal=NULL, pal.range="box", 
                    pad=0, demo=FALSE, ...) {
  
  
  # ---------------------------------------------------------------------
  # Warnings and error handling ...
  if (!identical(xlim, ylim)) {print("Note: Axes dimensions are not equal. The visual diagonal is *not* the line of numerical congruence! Consider choosing identical values for xlim and ylim.")}
  
  if (class(x) == "RSA") {
    stop("If you want to plot an RSA object, please use plot(...); plotRSA should be only used when you directly provide the regression coefficients.")
  }
  
  if (length(x) > 1 | length(y) > 1 | length(x2) > 1 | length(xy) > 1 | length(y2) > 1 | length(b0) > 1) {
    stop("Provide a single number to each regression coefficient")
  }
  
  
  # is the model a cubic model?
  is.cubicmodel <- model %in% c("cubic","CA","RRCA","CL","RRCL")
  
  # remove LOC, LOIC etc. when they do not make sense.
  if (any(c(w, wx, wy) != 0)) {
    axes <- ""
    project <- project[!project %in% c("PA1", "PA2", "LOC", "LOIC")]
  }
  if (is.cubicmodel) {
    axes <- axes[!axes %in% c("PA1", "PA2")]
    project <- project[!project %in% c("PA1", "PA2")]
  }
  if ((!model %in% c("CA","RRCA")) | x2 == 0 | x3 == 0) {
    axes <- axes[!axes %in% c("E2")]
    project <- project[!project %in% c("E2")]
  }
  if ((!model %in% c("CL","RRCL")) | x2 == 0 | x3 == 0 | is.null(fit)) {
    axes <- axes[!axes %in% c("K1", "K2")]
    project <- project[!project %in% c("K1", "K2")]
  }
  
  # define the defaults
  
  if (!is.null(points$data)) {
    points$data <- data.frame(points$data)	# a tibble causes an error ...
  }
  
  if (is.null(points) || (typeof(points) == "logical" && points == TRUE)) {
    points <- list(show=TRUE, value="raw", jitter=0, color="black", cex=.5, out.mark=FALSE)
  }
  if (is.null(points) || (typeof(points) == "logical" && points == FALSE)) {
    points <- list(show=FALSE)
  }
  
  if (is.null(points$show)) points$show <- TRUE
  if (is.na(points$show)) {
    if (is.null(points$data)) {
      points$show <- FALSE	
    } else {
      points$show <- TRUE
    }
  }
  if (is.null(points$value)) points$value <- "raw"
  if (is.null(points$color)) points$color <- "black"
  if (is.null(points$jitter)) points$jitter <- 0
  if (is.null(points$cex)) points$cex <- 0.5
  if (is.null(points$out.mark)) points$out.mark <- FALSE
  if (points$show==TRUE & is.null(points$data)) {
    warning("You must provide a data frame with the coordinates of the raw data points (points = list(show = TRUE, data = ???)). Points are not plotted.")
    points$show <- FALSE
  }
  if (!is.null(fit)) {
    if (points$out.mark==TRUE & fit$out.rm==FALSE) {
      warning("Outliers can only be marked in the plot when they were removed in the RSA function: RSA(..., out.rm=TRUE).")
      points$out.mark <- FALSE
    }
  }
  
  if (is.null(contour$show)) contour$show <- TRUE
  if (is.null(contour$color)) contour$color <- "grey30"
  if (is.null(contour$highlight)) contour$highlight <- c()
  
  # define default behavior for the "hull" parameter
  if (is.na(hull)) {
    if (!is.null(fit) | !is.null(points$data)) {
      hull <- TRUE
    } else {
      hull <- FALSE
    }
  }
  
  type <- match.arg(type, c("interactive", "3d", "contour"))
  surface <- match.arg(surface, c("predict", "smooth"))
  points[["value"]] <- match.arg(points[["value"]], c("raw", "predicted"))
  
  # define defaults of axes styles
  if (is.null(axesStyles[["LOC"]])) axesStyles[["LOC"]] <- list(lty="solid",  lwd=2, col=ifelse(bw==TRUE, "black", "blue"))
  if (is.null(axesStyles[["LOIC"]])) axesStyles[["LOIC"]] <- list(lty="solid",  lwd=2, col=ifelse(bw==TRUE, "black", "blue"))
  if (is.null(axesStyles[["PA1"]])) axesStyles[["PA1"]] <- list(lty="dotted", lwd=2, col=ifelse(bw==TRUE, "black", "gray30"))
  if (is.null(axesStyles[["PA2"]])) axesStyles[["PA2"]] <- list(lty="dotted", lwd=2, col=ifelse(bw==TRUE, "black", "gray30"))	
  if (is.null(axesStyles[["E2"]])) axesStyles[["E2"]] <- list(lty="solid", lwd=2, col=ifelse(bw==TRUE, "black", "deeppink"))	
  if (is.null(axesStyles[["K1"]])) axesStyles[["K1"]] <- list(lty="solid", lwd=2, col=ifelse(bw==TRUE, "black", "deeppink"))	
  if (is.null(axesStyles[["K2"]])) axesStyles[["K2"]] <- list(lty="solid", lwd=2, col=ifelse(bw==TRUE, "black", "deeppink"))	
  
  if (demo == FALSE) {
    if (is.null(xlab)) {
      if (!is.null(points$data)) {
        xlab <- colnames(points$data)[1]
      } else {
        xlab <- "X"
      }
    }
    if (is.null(ylab)) {
      if (!is.null(points$data)) {
        ylab <- colnames(points$data)[2]
      } else {
        ylab <- "Y"
      }
    }
    if (is.null(zlab)) {
      if (!is.null(points$data)) {
        zlab <- colnames(points$data)[3]
      } else {
        zlab <- "Z"
      }
    }
    
    if (is.null(xlim)) {xlim <- c(-2.1, 2.1)}
    if (is.null(ylim)) {ylim <- c(-2.1, 2.1)}
  }
  
  
  if (is.null(points$data) & surface == "smooth") {
    warning("Smoothing only works if data points are provided (points=list(data=???))! Reverting to surface = 'predict'")
    surface <- "predict"
  }
  
  C <- c(x, y, x2, y2, xy, w, wx, wy, x3, x2y, xy2, y3)
  
  if (!model %in% c("absunc", "absdiff")  & !is.cubicmodel) {
    if (!is.null(fit) & model != "null") {
      SP <- RSA.ST(fit, model=model)
      PAR <- getPar(fit, "coef", model=model)
      
      SP.text <- paste0("a", 1:5, ": ", f2(SP$SP$estimate, 2), p2star(SP$SP$p.value), collapse="   ")
      
      a4rs_par <- PAR[PAR$label == "a4.rescaled", ]
      if (nrow(a4rs_par) == 1) {
        SP.text <- paste0(SP.text, "/ a4(rescaled): ", f2(a4rs_par$est, 2), p2star(a4rs_par$pvalue))
      }
      SP.text <- paste0(SP.text, "\n")
      
      
      meanlevel <- PAR[PAR$label == "meaneffect", ]
      if (nrow(meanlevel) == 1) {
        SP.text <- paste0(SP.text, "mean-level effect = ", f2(meanlevel$est, 2), p2star(meanlevel$pvalue), "    ")
      }
      
      C_par <- PAR[PAR$label == "C", ]
      if (nrow(C_par) == 1) {
        SP.text <- paste0(SP.text, "C = ", f2(C_par$est, 2), p2star(C_par$pvalue), "    ")
      }
      
      S_par <- PAR[PAR$label == "S", ]
      if (nrow(S_par) == 1) {
        SP.text <- paste0(SP.text, "S = ", f2(S_par$est, 2), p2star(S_par$pvalue), "    ")
      }
      
    } else {
      SP <- RSA.ST(x=x, y=y, xy=xy, x2=x2, y2=y2)
      SP.text <- paste0("a", 1:5, ": ", f2(SP$SP$estimate, 2), p2star(SP$SP$p.value), collapse="    ")			
    }		
  } else {
    SP <- NULL
    param <- FALSE
    SP.text <- ""
  }
  
  # Print coefs in 3d plot?
  COEFS <- ""
  if (coefs == TRUE) {
    COEFS <- paste0("b1 = ", f2(x, 3), "\n", "b2 = ", f2(y, 3), "\n", "b3 = ", f2(x2, 3), "\n", "b4 = ", f2(xy, 3), "\n", "b5 = ", f2(y2, 3), "\n")
    if (is.cubicmodel){
      COEFS <- paste0(COEFS, "b6 = ", f2(x3, 3), "\n", "b7 = ", f2(x2y, 3), "\n", "b8 = ", f2(xy2, 3), "\n", "b9 = ", f2(y3, 3), "\n")
    }
  }
  
  
  # ---------------------------------------------------------------------
  # Calculate positions of raw points
  
  xpoints <- ypoints <- zpoints <- NA
  if (points$out.mark == TRUE & is.null(fit)) {
    warning("Outliers can only be marked if an RSA-object is provided. Points are not plotted.")
    points$show <- FALSE
  }
  if (points$show == TRUE | hull==TRUE) {
    if (is.null(points$data)) stop("You must provide raw data if you want to plot raw points or the bagplot.")
    if (points$out.mark == FALSE) {
      data.used <- points$data
    }
    if (points$out.mark == TRUE & !is.null(fit)) {
      data.used <- fit$data[, c(fit$IV1, fit$IV2, fit$DV)]	# this includes all data points
    }
    
    if (points$jitter > 0) {
      data.used[, 1] <- data.used[, 1] + rnorm(length(data.used[, 1]), 0, points$jitter)
      data.used[, 2] <- data.used[, 2] + rnorm(length(data.used[, 2]), 0, points$jitter)
    }
    
    if (points$value == "raw") {
      zpoints <- as.vector(data.used[, 3])
    } else if (points$value == "predicted") {
      N <- colnames(data.used)
      data.used2 <- add.variables(formula(paste0(N[3], " ~ ", N[1], "*", N[2])), data.used)
      
      # calculate predicted values
      zpoints <- b0 + colSums(C*t(data.used2[, c(
        N[1],
        N[2],
        paste0(N[1], "2"),
        paste0(N[2], "2"),
        paste0(N[1], "_", N[2]),
        "W",
        paste0("W_", N[1]),
        paste0("W_", N[2]),
        paste0(N[1], "3"),
        paste0(N[1], "2", "_", N[2]),
        paste0(N[1], "_", N[2], "2"),
        paste0(N[2], "3")
      )]))
      
      zpoints <- as.vector(zpoints)
    }
    
    xpoints <- as.vector(data.used[, 1])
    ypoints <- as.vector(data.used[, 2])	
  }
  
  
  
  # build data set
  grid <- gridsize
  new <- data.frame(x = rep(seq(xlim[1], xlim[2], length.out=grid), grid), y = rep(seq(ylim[1], ylim[2], length.out=grid), each=grid))
  new2 <- add.variables(z~x+y, new)
  
  # calculate z values of the surface
  if (surface == "predict") {
    new2$z <- b0 + colSums(C*t(new2[, c(1:5, 9:11, 15:18)]))
  }
  if (surface == "smooth") {
    
    if (!requireNamespace("fields", quietly = TRUE)) {
      stop('`fields` package needed for smooth surfaces to work. Please install it with install.packages("fields")', call. = FALSE)
    }
    
    tpsfit <- fields::Tps(points$data[, 1:2], points$data[, 3], scale.type="unscaled", lambda=lambda)
    new2$z <- fields::predict.Krig(tpsfit, new[, c("x", "y")])
    param <- FALSE
    axes <- ""
  }
  
  # impose link functions, both to the surface and the raw values
  logit <- function (x) {log(x/(1-x))}
  invlogit <- function (x) {1/(1+exp(-x))}
  link <- match.arg(link, c("identity", "logit", "probit"))
  if (link == "probit") {	
    # surface
    z.trans <- 1.7 * new2$z
    new2$z <- invlogit(z.trans)
    
    # raw data points
    if (points$value == "predicted") {
      zpoints.trans <- 1.7 * zpoints
      zpoints <- invlogit(zpoints.trans)
    }
  }
  if (link == "logit") {
    # surface
    new2$z <- invlogit(new2$z)
    
    # raw data points
    if (points$value == "predicted") {
      zpoints <- invlogit(zpoints)
    }
  }
  
  
  # determine zlim
  if (!is.null(points$data) & demo==FALSE & is.null(zlim)) {
    # old: set zlim according to fitted surface
    #zlim <- c(min(min(new2$z, na.rm=TRUE), min(points$data[, 3], na.rm=TRUE)), max(max(new2$z, na.rm=TRUE), max(points$data[, 3], na.rm=TRUE)))
    
    # new: set zlim according to actual data range
    zlim <- c(min(points$data[, 3], na.rm=TRUE), max(points$data[, 3], na.rm=TRUE))
  } else {
    if (is.null(zlim) & link != "probit") zlim <- c(min(new2$z), max(new2$z))
    if (is.null(zlim) & link == "probit") zlim <- c(0, 1)
  }
  zlim.final <- zlim
  
  
  # Catch border case: completely flat surface. Remove contour, redefine zlim
  if (var(new2$z) == 0) {
    contour$show <- FALSE
    project <- project[-which(project == "contour")]
    if (zlim[1] == 0 && zlim[2] == 0) zlim <- xlim
  }
  
  
  ## Define colors
  
  # flip palette?
  flip <- FALSE
  if (!is.null(pal) && pal=="flip") {
    flip <- TRUE
    pal <- NULL
  }
  if (bw == FALSE) {
    
    # RdYlGn palette
    if (is.null(pal)) {
      pal <- c("#A50026","#D73027","#F46D43","#FDAE61","#FEE08B","#FFFFBF","#D9EF8B","#A6D96A","#66BD63","#1A9850","#006837")
      if (flip==TRUE) {pal <- rev(pal)}
    }
    
    gridCol <- ifelse(contour$show == TRUE, "grey20", "grey20")
  } else {
    
    # B/W palette
    if (is.null(pal)) {
      pal <- colorRampPalette(c("#FFFFFF", "#AAAAAA", "#030303"), bias=2)(11)
      if (flip==TRUE) {pal <- rev(pal)}
    }
    gridCol <- ifelse(contour$show == TRUE, "grey20", "grey20")
  }
  if (length(pal) < 2) {legend <- FALSE}
  
  if (suppress.grid == TRUE) {
    gridCol <- "transparent"
  }
  
  # ---------------------------------------------------------------------
  #  calculate bag plot: bag = outer, loop = inner
  
  if (hull==TRUE) {
    BAG <- aplpack::compute.bagplot(xpoints, ypoints)
    
    # close the polygon
    h1 <- rbind(BAG$hull.bag, BAG$hull.bag[1, ])
    h2 <- rbind(BAG$hull.loop, BAG$hull.loop[1, ])
    
    # approx: interpolate the points of the bag (in order to get a more smooth fitting line on the z-axis)
    minDist <- min(diff(xlim)/gridsize, diff(ylim)/gridsize)/2
    
    h1 <- interpolatePolygon(h1[, 1], h1[, 2], minDist=minDist)	
    h2 <- interpolatePolygon(h2[, 1], h2[, 2], minDist=minDist)	
    
    # calculate predicted values
    bagpoints <- add.variables(z~x+y, data.frame(x=h1$x, y=h1$y))
    bagpoints$z <- b0 + colSums(C*t(bagpoints[, c(1:5, 9:11, 15:18)]))
    bag <- data.frame(X  = bagpoints$x, Y  = bagpoints$y, Z = bagpoints$z)
    
    looppoints <- add.variables(z~x+y, data.frame(x=h2$x, y=h2$y))
    looppoints$z <- b0 + colSums(C*t(looppoints[, c(1:5, 9:11, 15:18)]))
    loop <- data.frame(X  = looppoints$x, Y  = looppoints$y, Z = looppoints$z)			
  }
  
  
  ## ======================================================================
  ## Interactive plot
  ## ======================================================================
  
  if (type == "interactive") {
    if (!requireNamespace("rgl", quietly = TRUE)) {
      stop("`rgl` package needed for interactive plots. Please install it with install.packages('rgl').", call. = FALSE)
    }
    
    P <- list(x=seq(xlim[1], xlim[2], length.out=grid), y=seq(ylim[1], ylim[2], length.out=grid))
    DV2 <- matrix(new2$z, nrow=grid, ncol=grid, byrow=FALSE)
    R <- range(DV2)
    col2 <- as.character(cut(1:(R[2] - R[1] + 1), breaks=length(pal), labels=pal))
    
    rgl::open3d(cex=cex.main)
    rgl::rgl.viewpoint(-30, -90, fov=0)
    rgl::rgl.light(theta = 0, phi = 90, viewpoint.rel = TRUE, ambient = "#FF0000", diffuse = "#FFFFFF", specular = "#FFFFFF")
    rgl::persp3d(P$x, P$y, DV2, xlab = xlab, ylab = ylab, zlab = zlab, color=col2[DV2 - R[1] + 1], main=main, ...)
    
    if (contour$show == TRUE) {
      contours <- contourLines(P, z=DV2)
      for (i in 1:length(contours)) {
        with(contours[[i]], rgl::lines3d(x, y, level, col=contour$color))
      }
    }
    
    if (points$show == TRUE) {
      if (points$out.mark == FALSE) {
        rgl::points3d(data.frame(xpoints, ypoints, zpoints), col=points$color)
      }
      if (points$out.mark == TRUE) {
        if (!is.null(fit)) {
          colvec <- rep(points$color, nrow(fit$data))
          colvec[fit$outliers] <- "red"
          rgl::points3d(fit$data[, c(fit$IV1, fit$IV2, fit$DV)], col=colvec)
          rgl::text3d(fit$data[fit$outliers, c(fit$IV1, fit$IV2, fit$DV)], col="red", texts="X")
        } else {
          warning("Please provide an RSA-object to mark outliers.")
        }
      }
    }
    
    p1 <- NULL	# no plot object is returned	
  }	
  
  
  
  ## ======================================================================
  ## Wireframe plot
  ## ======================================================================
  
  
  if (type == "3d") {
    
    mypanel2 <- function(x, y, z, xlim, ylim, zlim, xlim.scaled, ylim.scaled, zlim.scaled, axes, axesList, x.points=NULL, y.points=NULL, z.points=NULL, SPs="", ...) {
      
      
      # rescale absolute x, y, and z values so that they fit into the box
      RESCALE.Z <- function(z1) {
        Z2 <- zlim.scaled[1] + diff(zlim.scaled) * (z1 - zlim[1]) / diff(zlim)
        return(Z2)
      }
      RESCALE <- function(n) {
        X2 <- xlim.scaled[1] + diff(xlim.scaled) * (n$X - xlim[1]) / diff(xlim)
        Y2 <- ylim.scaled[1] + diff(ylim.scaled) * (n$Y - ylim[1]) / diff(ylim)
        Z2 <- zlim.scaled[1] + diff(zlim.scaled) * (n$Z - zlim[1]) / diff(zlim)
        df <- data.frame(X=X2, Y=Y2, Z=Z2)
        df <- df[df$X >= min(xlim.scaled) & df$X <= max(xlim.scaled) & df$Y >= min(ylim.scaled) & df$Y <= max(ylim.scaled) &  df$Z >= min(zlim.scaled) & df$Z <= max(zlim.scaled), ]
        return(df)
      }
      
      
      # ---------------------------------------------------------------------
      # 1. Projection on bottom of cube
      if (length(project) > 0) {
        for (p in project) {
          if (p %in% c("LOC", "LOIC", "PA1", "PA2", "E2", "K1", "K2")) {
            if (is.null(axesList[[p]])) break;
            
            a0 <- RESCALE(getIntersect2(p0=axesList[[p]]$p0, p1=axesList[[p]]$p1))
            if (nrow(a0) <= 1) break;
            panel.3dscatter(x = a0$X, y = a0$Y, z = rep(RESCALE.Z(min(zlim.final) + .01), nrow(a0)), 
                            xlim = xlim, ylim = ylim, zlim = zlim,
                            xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, 
                            type="l", col.line=axesList[[p]]$style[["col"]], lty=axesList[[p]]$style[["lty"]], lwd=axesList[[p]]$style[["lwd"]], ...)
          }
          
          if (p == "hull") {
            bag.rescale <- RESCALE(bag)
            panel.3dscatter(x = bag.rescale$X, y = bag.rescale$Y, z = rep(RESCALE.Z(min(zlim.final) + .01), nrow(bag.rescale)), 
                            xlim = xlim, ylim = ylim, zlim = zlim,
                            xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, 
                            type="l", col.line="white", lty="dashed", lwd=2, ...)
            
            loop.rescale <- RESCALE(loop)
            panel.3dscatter(x = loop.rescale$X, y = loop.rescale$Y, z = rep(RESCALE.Z(min(zlim.final) + .01), nrow(loop.rescale)), 
                            xlim = xlim, ylim = ylim, zlim = zlim,
                            xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, 
                            type="l", col.line="black", lty="dashed", lwd=2, ...)			
            
          }
          
          if (p == "points") {
            x2 <- xlim.scaled[1] + diff(xlim.scaled) * (x.points - xlim[1]) / diff(xlim)
            y2 <- ylim.scaled[1] + diff(ylim.scaled) * (y.points - ylim[1]) / diff(ylim)
            z2 <- rep(RESCALE.Z(min(zlim.final) + .01), length(x2))
            
            panel.3dscatter(x = x2, y = y2, z = z2, 
                            xlim = xlim, ylim = ylim, zlim = zlim,
                            xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled,
                            pch=20, col=points$color, cex=points$cex, ...)
          }
        }
      }
      
      
      # project contour on bottom
      if (contour$show == TRUE | "contour" %in% project) {
        # "abuse" ggplot to compute the contour lines
        cs <- ggplot(new2, aes_string(x="x", y="y", fill="z", z="z")) + stat_contour(bins=ifelse(length(pal)>1, length(pal)+1, 8))
        cLines <- ggplot_build(cs)
        C0 <- cLines$data[[1]][, c("x", "y", "level", "group")]
        colnames(C0) <- c("X", "Y", "Z", "group")	# C0 keeps the contour lines
        
        
        if ("contour" %in% project) {	
          for (cL in C0$group) {
            C1 <- RESCALE(C0[C0$group==cL, c("X", "Y", "Z")])
            panel.3dscatter(x = C1$X, y = C1$Y, z = rep(RESCALE.Z(min(zlim.final) + .01), nrow(C1)), 
                            xlim = xlim, ylim = ylim, zlim = zlim,
                            xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled,
                            type="l", col.line=contour$color, lty="solid", lwd=1, ...)
          }	
        }
      }
      
      
      
      # ---------------------------------------------------------------------
      # 2. Borders, back part
      
      
      if (border==TRUE & suppress.surface==FALSE) {
        # Make boundary of grid a bit thicker
        box1 <- new2[new2$y == max(new2$y), ]
        box2 <- new2[new2$y == min(new2$y), ]
        box3 <- new2[new2$x == max(new2$x), ]
        box4 <- new2[new2$x == min(new2$x), ]
        box <- rbind(data.frame(box1, side=1), data.frame(box2, side=2), data.frame(box3, side=3), data.frame(box4, side=4))
        
        x.box <- xlim.scaled[1] + diff(xlim.scaled) * (box$x - xlim[1]) / diff(xlim)
        y.box <- ylim.scaled[1] + diff(ylim.scaled) * (box$y - ylim[1]) / diff(ylim)
        z.box <- zlim.scaled[1] + diff(zlim.scaled) * (box$z - zlim[1]) / diff(zlim)
        
        # plot the back lines of the border
        panel.3dscatter(x = x.box[box$side==1], y = y.box[box$side==1], z = z.box[box$side==1], xlim = xlim, ylim = ylim, zlim = zlim, xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, type="l", col.line=gridCol, lwd=4, ...)
        panel.3dscatter(x = x.box[box$side==3], y = y.box[box$side==3], z = z.box[box$side==3], xlim = xlim, ylim = ylim, zlim = zlim, xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, type="l", col.line=gridCol, lwd=4, ...)
      }
      
      # ---------------------------------------------------------------------
      # 3. the surface
      if (suppress.surface==FALSE) {
        panel.3dwire(x = x, y = y, z = z, xlim = xlim, ylim = ylim, zlim = zlim,
                     xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled,
                     col=gridCol, lwd=0.5, ...)
        
      }
      
      
      # ---------------------------------------------------------------------
      # 4. plot of LOC and LOIC, and other axes
      if (suppress.surface==FALSE) {
        for (a in axes) {
          if (!is.null(axesList[[a]])) {
            a0 <- RESCALE(getIntersect2(p0=axesList[[a]]$p0, p1=axesList[[a]]$p1))
            if (nrow(a0) <= 1) break;
            panel.3dscatter(x = a0$X, y = a0$Y, z = a0$Z, xlim = xlim, ylim = ylim, zlim = zlim,
                            xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, 
                            type="l", col.line=axesList[[a]]$style[["col"]], lty=axesList[[a]]$style[["lty"]], lwd=axesList[[a]]$style[["lwd"]], ...)
          }
        }   
      }
      
      # ---------------------------------------------------------------------
      # 4b. plot of maximum lines
      if (maxlines == TRUE & suppress.surface==FALSE) {
        # maximum X for a given Y
        a0 <- RESCALE(getIntersect2(p0=-(C[1]/C[5]), p1=-((2*C[3])/C[5])))
        panel.3dscatter(x = a0$X, y = a0$Y, z = a0$Z, xlim = xlim, ylim = ylim, zlim = zlim,
                        xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled,
                        type="l", col.line="red", lty="dashed", lwd=2, ...)
        
        
        
        a0 <- RESCALE(getIntersect2(p0=-(C[2]/(2*C[4])), p1=-((C[5])/(2*C[4]))))
        #a0 <- RESCALE(getIntersect2(p0=-(C[2]/C[5]), p1=-((2*C[4])/C[5])))
        panel.3dscatter(x = a0$X, y = a0$Y, z = a0$Z, xlim = xlim, ylim = ylim, zlim = zlim,
                        xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled,
                        type="l", col.line="blue", lty="dashed", lwd=2, ...)
        
      }
      
      
      
      # ---------------------------------------------------------------------
      # 5. Borders, front part	  
      
      if (border==TRUE & suppress.surface==FALSE) {
        # plot the front boundary lines
        panel.3dscatter(x = x.box[box$side==2], y = y.box[box$side==2], z = z.box[box$side==2], xlim = xlim, ylim = ylim, zlim = zlim, xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, type="l", col.line=gridCol, lwd=3, ...)
        panel.3dscatter(x = x.box[box$side==4], y = y.box[box$side==4], z = z.box[box$side==4], xlim = xlim, ylim = ylim, zlim = zlim, xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, type="l", col.line=gridCol, lwd=3, ...)
      }
      
      if (param == TRUE) {
        grid::grid.text(SPs, .02, .95, just="left", gp=grid::gpar(cex=cex.axesLabel*0.8))
      }  
      
      if (coefs == TRUE) {
        grid::grid.text(COEFS, .80, .87, just="left", gp=grid::gpar(cex=cex.axesLabel*0.8))
      }  
      
      
      
      # ---------------------------------------------------------------------
      # 6a: The bag plot, if requested
      
      if (hull==TRUE & suppress.surface==FALSE) {	
        
        # bag (= inner bag)
        if (any(bag$X < xlim[1] | bag$X > xlim[2] | bag$Y < ylim[1] | bag$Y > ylim[2])) {
          warning("The bag is partly outside the plotting region. Bag is not displayed, please adjust xlim and ylim to include the full range of raw data.")
        } else {
          bag.rescale <- RESCALE(bag)
          panel.3dscatter(x = bag.rescale$X, y = bag.rescale$Y, z = bag.rescale$Z, xlim = xlim, ylim = ylim, zlim = zlim, xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, type="l", col.line="grey30", lty="dashed", lwd=2, ...)
          
        }
        
        # loop (= outer bag)
        if (any(loop$X < xlim[1] | loop$X > xlim[2] | loop$Y < ylim[1] | loop$Y > ylim[2])) {
          warning("The loop is partly outside the plotting region. Loop is not displayed, please adjust xlim and ylim to include the full range of raw data.")
        } else {
          loop.rescale <- RESCALE(loop)
          panel.3dscatter(x = loop.rescale$X, y = loop.rescale$Y, z = loop.rescale$Z, xlim = xlim, ylim = ylim, zlim = zlim, xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled, type="l", col.line="black", lty="dashed", lwd=2, ...)
          
        }
      }	  	
      
      
      # ---------------------------------------------------------------------
      # 6b. Raw data points scatter plot	  
      if (points$show == TRUE) {
        
        x2 <- xlim.scaled[1] + diff(xlim.scaled) * (x.points - xlim[1]) / diff(xlim)
        y2 <- ylim.scaled[1] + diff(ylim.scaled) * (y.points - ylim[1]) / diff(ylim)
        z2 <- zlim.scaled[1] + diff(zlim.scaled) * (z.points - zlim[1]) / diff(zlim)
        
        panel.3dscatter(x = x2, y = y2, z = z2, xlim = xlim, ylim = ylim, zlim = zlim,
                        xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled,
                        pch=20, col=points$color, cex=points$cex, ...)
        # plot outliers
        if (points$out.mark==TRUE) {
          panel.3dscatter(x = x2[fit$outliers], y = y2[fit$outliers], z = z2[fit$outliers], xlim = xlim, ylim = ylim, zlim = zlim,
                          xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled,
                          pch=4, col="red", cex=points$cex, ...)
        }
        
      }	
      
      # ---------------------------------------------------------------------
      # 7. plot contour lines on surface:
      
      if (contour$show == TRUE & suppress.surface==FALSE) {
        # C0 keeps the contour lines and has been computed before
        for (cL in C0$group) {
          C1 <- RESCALE(C0[C0$group==cL, c("X", "Y", "Z")])
          
          if (contour$show == TRUE) {
            panel.3dscatter(x = C1$X, y = C1$Y, z = C1$Z, 
                            xlim = xlim, ylim = ylim, zlim = zlim,
                            xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled,
                            type="l", col.line=contour$color, lty="solid", lwd=1, ...)
          }								
        }
        
        # highlight specific contour lines?
        if (length(contour$highlight) > 0) {
          C2 <- C0[C0$Z %in% f0(unique(C0$Z), contour$highlight), ]
          for (cL in C2$group) {
            C3 <- RESCALE(C2[C2$group==cL, c("X", "Y", "Z")])
            panel.3dscatter(x = C3$X, y = C3$Y, z = C3$Z, xlim = xlim, ylim = ylim, zlim = zlim,
                            xlim.scaled = xlim.scaled, ylim.scaled = ylim.scaled, zlim.scaled = zlim.scaled,
                            type="l", col.line=contour$color, lty="solid", lwd=2, ...)
          }
          
        }	
      }
      
    }  # of mypanel2
    
    
    
    # local function: compute the surface line, defined by a line on the X-Y plane (p0 = intercept, p1=slope)
    getIntersect2 <- function(p0, p1, Z=NULL) {
      X <- seq(min(xlim), max(xlim), length.out=grid*2)
      Y <- p0 + p1*X
      n <- data.frame(X, Y)
      n2 <- add.variables(z~X+Y, n)
      n2$Z <- b0 + colSums(c(x, y, x2, y2, xy, x3, x2y, xy2, y3)*t(n2[, c("X","Y","X2","Y2","X_Y","X3","X2_Y","X_Y2","Y3")]))
      if (!is.null(Z)) n2$Z <- Z
      return(n2[, c("X", "Y", "Z")])
    }
    
    axesList <- list()
    axesList[["LOC"]]  <- list(p0=0, p1=1, style=axesStyles[["LOC"]])
    axesList[["LOIC"]] <- list(p0=0, p1=-1, style=axesStyles[["LOIC"]])
    
    if (x2 != y2 & !is.cubicmodel) {
      axesList[["PA1"]] <- list(p0=SP$p10, p1=SP$p11, style=axesStyles[["PA1"]])
      axesList[["PA2"]] <- list(p0=SP$p20, p1=SP$p21, style=axesStyles[["PA2"]])	
    }	
    
    axesList[["E2"]] <- list(p0=(2*x2/(3*x3)), p1=1, style=axesStyles[["E2"]])
    
    if ((model=="CL" | model=="RRCL") & !is.null(fit)){
      clrange <- clRange(fit, model=model)
      if (!is.na(clrange$k1)){
        axesList[["K1"]] <- list(p0=2*clrange$k1, p1=-1, style=axesStyles[["K1"]])
      } 
      if (!is.na(clrange$k2)){
        axesList[["K2"]] <- list(p0=2*clrange$k2, p1=-1, style=axesStyles[["K2"]])
      }
    }
    
    
    # Define color range: Relative to surface min/max, or relative to box (zlim)?
    if (pal.range == "box") {
      at <- seq(zlim[1], zlim[2], length.out=length(pal)-1)
    } else if (pal.range == "surface") {
      at <- seq(min(new2$z), max(new2$z), length.out=length(pal)-1)
    }
    
    # define the appearance of the color legend
    CK <- FALSE
    if (legend == TRUE) {
      CK <- list(labels=list(cex=cex.axesLabel))
    }
    
    # Define appearance of the surrounding box
    axesCol <- "black"
    boxCol <- "black"
    
    if (suppress.box == TRUE) {
      axesCol <- "transparent"
      boxCol <- NA
    }
    
    p1 <- wireframe(z ~ x*y, new2,  drape=TRUE, 
                    scales 	= list(arrows = FALSE, cex=cex.tickLabel, col = axesCol, font = 1, tck=tck, distance=distance), 
                    xlab	= list(cex=cex.axesLabel, label=xlab, rot=label.rotation[["x"]]), 
                    ylab	= list(cex=cex.axesLabel, label=ylab, rot=label.rotation[["y"]]), 
                    zlab	= list(cex=cex.axesLabel, label=zlab, rot=label.rotation[["z"]]), zlim=zlim, 
                    main	= list(cex=cex.main, label=main),
                    screen	= rotation,
                    at		= at, col.regions=pal, colorkey=CK, 
                    par.settings = list(
                      axis.line = list(col = "transparent"), 
                      layout.heights = list(top.padding=pad, bottom.padding=pad), 
                      layout.widths=list(left.padding=pad, right.padding=pad),
                      box.3d = list(col=boxCol)), 
                    axes	= axes,
                    axesList= axesList, 
                    SPs		= SP.text,
                    COEFS	= COEFS, 
                    panel.3d.wireframe = mypanel2,
                    x.points=xpoints, y.points=ypoints, z.points=zpoints)
    
  }  # of type == "3d"
  
  
  
  ## ======================================================================
  ## Contour plot
  ## ======================================================================
  if (type == "contour") {
    if (!all(C == 0)) {
      
      # Define color range: Relative to surface min/max, or relative to box (zlim)?
      if (pal.range == "box") {
        limits <- c(zlim[1], zlim[2])
      } else if (pal.range == "surface") {
        limits <- c(min(new2$z), max(new2$z))
      }
      
      p1 <- ggplot(new2, aes_string(x="x", y="y", z="z")) + geom_tile(aes_string(fill="z")) + 
        scale_fill_gradientn(zlab, colours=pal, limits=limits) + theme_bw() + 
        theme(aspect.ratio=1) + xlab(xlab) + ylab(ylab)
      
      if (legend==FALSE) {
        p1 <- p1 + guides(fill=FALSE)
      }
      
     # p1 <- p1 + geom_contour(bins = 40, alpha=.6, colour = "grey50")
      
      # highlight specific contour lines?
      if (length(contour$highlight) > 0) {
        cLines <- ggplot_build(p1)
        C0 <- cLines$data[[2]][, c("x", "y", "level", "group")]
        
        # Find closest values in contours
        C1 <- C0[C0$level %in% f0(unique(C0$level), contour$highlight), ]
        p1 <- p1 + geom_path(data=C1, aes_string(x="x", y="y", group="group", z="level"), size=1.1)
      }
      
      # (in)congruence lines
      if ("LOC" %in% axes) {
        p1 <- p1 + geom_abline(aes(intercept=0, slope=1), color="grey20")
      }
      if ("LOIC" %in% axes) {
        p1 <- p1 + geom_abline(aes(intercept=0, slope=-1), linetype="dotted", size=1, color="grey20")
      }
      
      if (!model %in% c("absunc", "absdiff")  & !is.cubicmodel){
        if (("PA1" %in% axes) & !any(is.na(SP[c("p10", "p11")]))) {
          p1 <- p1 + geom_abline(data=data.frame(SP[c("p10", "p11")]), aes_string(intercept="p10", slope="p11"), color="grey20")
        }
        if (("PA2" %in% axes) & !any(is.na(SP[c("p20", "p21")]))) {
          p1 <- p1 + geom_abline(data=data.frame(SP[c("p20", "p21")]), aes_string(intercept="p20", slope="p21"), linetype="dotted", color="grey20")
        }
      }
      
      if ("E2" %in% axes) {
        E20 <- 2*x2/(3*x3)
        p1 <- p1 + geom_abline(aes(intercept=E20, slope=1), color="deeppink")
      }
      
      if ("K1" %in% axes) {
        k1 <- clRange(fit, model=model)$k1
        if (!is.na(k1)){
          p1 <- p1 + geom_abline(aes(intercept=2*k1, slope=-1), color="deeppink")
        }
      }
      
      if ("K2" %in% axes) {
        k2 <- clRange(fit, model=model)$k2
        if (!is.na(k2)){
          p1 <- p1 + geom_abline(aes(intercept=2*k2, slope=-1), color="deeppink")
        }
      }
      
      if (!model %in% c("absunc", "absdiff")  & !is.cubicmodel){
        if (showSP==TRUE & !any(is.na(SP[c("X0", "Y0")])) & !model %in% c("RR", "SQD", "SSQD", "SRSQD", "SRR", "SRRR") & !is.cubicmodel) {
          p1 <- p1 + annotate("point", x=SP$X0, y=SP$Y0, z=max(new2$z))
        }
      }
      
      
      if (points$show == TRUE) {
        if (points$out.mark==FALSE) {
          p1 <- p1 + annotate("point", x=xpoints, y=ypoints, color=points$color, size=3*points$cex)
        }
        if (points$out.mark==TRUE) {
          colvec <- rep(points$color, nrow(fit$data))
          colvec[fit$outliers] <- "red"
          shapevec <- rep(19, nrow(fit$data))
          shapevec[fit$outliers] <- 4
          p1 <- p1 + annotate("point", x=fit$data[, fit$IV1], y=fit$data[, fit$IV2], color=colvec, size=3*points$cex, shape=shapevec)
        }
      }
      
      if (hull==TRUE & !is.null(points$data)) {
        p1 <- p1 + annotate("path", x=bag$X, y=bag$Y, linetype="solid", size=1, color="grey10")
        p1 <- p1 + annotate("path", x=loop$X, y=loop$Y, linetype="dashed", size=1, color="grey10")
      }
      
      # plot CI of SP
      if (showSP==TRUE & showSP.CI==TRUE & !is.null(fit) & !is.cubicmodel) {
        PAR <- getPar(fit, "coef", model=model)
        p1 <- p1 + annotate("errorbar", x=SP$X0, y=SP$Y0, ymin=PAR[PAR$label=="Y0", "ci.lower"], ymax=PAR[PAR$label=="Y0", "ci.upper"], z=max(new2$z), width=.3)
        p1 <- p1 + annotate("errorbarh", x=SP$X0, y=SP$Y0, xmin=PAR[PAR$label=="X0", "ci.lower"], xmax=PAR[PAR$label=="X0", "ci.upper"], z=max(new2$z), height=.3)
      }
      
      # Title
      if (main != "") p1 <- p1 + ggtitle(main)
      
      p1 <- p1 + coord_cartesian(xlim=xlim, ylim=ylim)
      
    }
  }
  
  
  return(p1)
}


#' @method plot RSA
#' @export

# Purpose: Extract the model parameters, xlim, xlab, etc. from the fitted object and give it to the plotRSA function
plot.RSA <- function(x, ...) {
  fit <- x
  
  extras <- match.call(expand.dots = FALSE)$...
  
  if (is.null(extras)) {extras <- list()}
  if (is.null(extras[["model"]])) {extras[["model"]] <- "full"}
  
  C <- coef(fit$models[[extras$model]])
  if (fit$models[[extras$model]]@Options$estimator != "DWLS") {
    extras[["b0"]] <- as.numeric(ifelse(is.na(C[paste0(fit$DV, "~1")]), 0, C[paste0(fit$DV, "~1")]))
  } else {			
    # the threshold is the negative of the intercept ...
    extras[["b0"]] <- -as.numeric(ifelse(is.na(C[paste0(fit$DV, "|t1")]), 0, C[paste0(fit$DV, "|t1")]))
  }
  extras$x <- as.numeric(ifelse(is.na(C["b1"]), 0, C["b1"]))
  extras$y <- as.numeric(ifelse(is.na(C["b2"]), 0, C["b2"]))
  extras$x2 <- as.numeric(ifelse(is.na(C["b3"]), 0, C["b3"]))
  extras$y2 <- as.numeric(ifelse(is.na(C["b5"]), 0, C["b5"]))
  extras$xy <- as.numeric(ifelse(is.na(C["b4"]), 0, C["b4"]))
  extras$w <- as.numeric(ifelse(is.na(C["w1"]), 0, C["w1"]))
  extras$wx <- as.numeric(ifelse(is.na(C["w2"]), 0, C["w2"]))
  extras$wy <- as.numeric(ifelse(is.na(C["w3"]), 0, C["w3"]))
  
  # cubic parameters
  extras$x3 <- as.numeric(ifelse(is.na(C["b6"]), 0, C["b6"]))
  extras$x2y <- as.numeric(ifelse(is.na(C["b7"]), 0, C["b7"]))
  extras$xy2 <- as.numeric(ifelse(is.na(C["b8"]), 0, C["b8"]))
  extras$y3 <- as.numeric(ifelse(is.na(C["b9"]), 0, C["b9"]))
  
  if (is.null(extras[["xlab"]])) {extras[["xlab"]] <- fit$IV1}
  if (is.null(extras[["ylab"]])) {extras[["ylab"]] <- fit$IV2}
  if (is.null(extras[["zlab"]])) {extras[["zlab"]] <- fit$DV}
  
  extras$fit <- fit
  
  # define the defaults
  if (is.null(extras$points) || (typeof(extras$points) == "logical" && extras$points == TRUE)) {
    extras$points <- list(show=TRUE, value="raw", jitter=0, color="black", cex=.5, out.mark=FALSE)
  }
  if (is.null(extras$points) || (typeof(extras$points) == "logical" && extras$points == FALSE)) {
    extras$points <- list(show=FALSE, value="raw", jitter=0, color="black", cex=.5, out.mark=FALSE)
  }
  if (is.null(extras$points$out.mark)) extras$points$out.mark <- FALSE
  
  if (extras$points$out.mark == FALSE) {
    data.used <- fit$data[fit$data$out==FALSE, ]
  }
  if (extras$points$out.mark == TRUE) {
    data.used <- fit$data
  }
  
  extras$points$data <- data.used[, c(fit$IV1, fit$IV2, fit$DV, colnames(fit$data)[which(!colnames(fit$data) %in% c(fit$IV1, fit$IV2, fit$DV))])]
  
  adjust <- FALSE
  if (is.null(extras$xlim)) {
    extras$xlim <- c(min(data.used[, fit$IV1], na.rm=TRUE), max(data.used[, fit$IV1], na.rm=TRUE))
    # expand range by 20% at each end
    extras$xlim[1] <- extras$xlim[1]*ifelse(extras$xlim[1]<0, 1.1, 0.9)
    extras$xlim[2] <- extras$xlim[2]*ifelse(extras$xlim[2]<0, 0.9, 1.1)
    adjust <- TRUE
  }
  
  if (is.null(extras$ylim)) {
    extras$ylim <- c(min(data.used[, fit$IV2], na.rm=TRUE), max(data.used[, fit$IV2], na.rm=TRUE))
    extras$ylim[1] <- extras$ylim[1]*ifelse(extras$ylim[1]<0, 1.1, 0.9)
    extras$ylim[2] <- extras$ylim[2]*ifelse(extras$ylim[2]<0, 0.9, 1.1)
    adjust <- TRUE
  }
  
  if (adjust == TRUE) {
    extras$xlim[1] <- extras$ylim[1] <- min(extras$xlim[1], extras$ylim[1])
    extras$xlim[2] <- extras$ylim[2] <- max(extras$xlim[2], extras$ylim[2])
  }
  
  do.call(plotRSA, as.list(extras), envir = parent.frame())
}


# helpers.R

# compute interaction, squared, cubic, dummy variables, etc. for RSA
add.variables <- function(formula, df) {
  IV1 <- all.vars(formula)[2]
  IV2 <- all.vars(formula)[3]
  
  IV12 <- paste0(IV1, "2")
  IV22 <- paste0(IV2, "2")
  IV13 <- paste0(IV1, "3")
  IV23 <- paste0(IV2, "3")
  IV_IA <- paste0(IV1, "_", IV2)
  IV_IA2 <- paste0(IV1, "2", "_", IV2)
  IV_IA3 <- paste0(IV1, "_", IV2, "2")
  
  
  
  df[, IV12] <- df[, IV1]^2
  df[, IV22] <- df[, IV2]^2
  df[, IV_IA] <- df[, IV1]*df[, IV2]
  
  # three new variables for piecewise regression (test absolute difference score) - Edwards (2002) model
  df$W.JRE <- ifelse(df[, IV1] >= df[, IV2], 0, 1)
  df[, paste0("W.JRE_", IV1)] <- df$W.JRE*df[, IV1]
  df[, paste0("W.JRE_", IV2)] <- df$W.JRE*df[, IV2]
  
  # three new variables for piecewise regression (test absolute difference score) - new model Schoenbrodt 2012
  df$W <- ifelse(df[, IV1] >= df[, IV2], 1, -1)
  df$W[df[, IV1] == df[, IV2]] <- 0
  df[, paste0("W_", IV1)] <- df$W*df[, IV1]
  df[, paste0("W_", IV2)] <- df$W*df[, IV2]
  
  df$diff <- df[, IV2] - df[, IV1]
  df$SD <- df$diff^2
  df$absdiff <- abs(df$diff)
  
  # cubic terms
  df[, IV13] <- df[, IV1]^3
  df[, IV_IA2] <- df[, IV1]^2*df[, IV2]
  df[, IV_IA3] <- df[, IV1]*df[, IV2]^2
  df[, IV23] <- df[, IV2]^3
  
  return(df)
}


# helper function: takes a list of lavaan models (can include NULLs), and returns the usual anova object
anovaList <- function(modellist) {
  mods <- modellist[!sapply(modellist, function(x) is.null(x))]
  mods <- mods[!sapply(mods, function(x) !inspect(x, "converged"))]
  
  if (length(mods) == 0) {
    return(list(n.mods=0))
  }
  
  # put them in order (using df)
  DF <- sapply(mods, fitmeasures, "df")
  mods <- mods[order(DF, decreasing = FALSE)]
  
  
  # prevent lavaan error in case that some models (but not all) have scaled test statistics
  
  # detect such cases (condition copy-pasted from lav_test_LRT.R in lavaan)
  mods.scaled <- unlist( lapply(mods, function(x) {
    any(c("satorra.bentler", "yuan.bentler", "yuan.bentler.mplus", "mean.var.adjusted", "scaled.shifted") 
        %in% unlist(sapply(slot(x, "test"), "[", "test")) ) }))
  
  # change internal label
  if( ! ( all(mods.scaled) | !any(mods.scaled) ) ) {
    mods[[ which(sapply(mods, fitmeasures, "df") == 0) ]]@test[[2]]$test <- mods[[ which(mods.scaled)[1] ]]@test[[2]]$test
  } 
  
  
  pStr <- sapply(1:length(mods), function(x){ 
    if(x==1) {
      paste("mods[[",x,"]]",sep = "")
    } else {
      paste("force(mods[[",x,"]])",sep = "")
    }
  })
  #pStr2 <- paste0("lavTestLRT(", paste(pStr, collapse=", "), ", method='satorra.bentler.2010')")
  pStr2 <- paste0("lavTestLRT(", paste(pStr, collapse=", "), ", method='default')")
  
  a1 <- eval(parse(text = pStr2))
  
  if (length(mods) > 1) {
    rownames(a1) <- names(mods)
  }
  
  attr(a1, "n.mods") <- length(mods)
  return(list(ANOVA=a1, models=mods, n.mods=length(mods)))
}


## internal helper function: compare models
# mL = model list
# set = label that is attached to the results
cModels <- function(mL, set, free.max) {
  aL1 <- anovaList(mL)
  if (aL1$n.mods > 1) {
    N <- lavaan::nobs(aL1$models[[1]])
    a1 <- cbind(aL1$ANOVA[, c(1, 4:7)], plyr::ldply(aL1$models, function(X) {
      F <- fitmeasures(X)
      R <- inspect(X, "r2")
      names(R) <- "R2"
      n <- lavaan::nobs(X)
      k <- free.max - F["df"] # number of parameters, including coefficients of control variables (if there are any)
      
      suppressWarnings({		
        R2.p <- ifelse(k==0,
                       NA,
                       pf(((n-k-1)*R)/(k*(1-R)), k, n-k-1, lower.tail=FALSE))
      })
      
      names(R2.p) <- "R2.p"
      
      # compute AICc
      K <- k + 2  # number of parameters, including coefficients of control variables (if there are any), intercept and residual variance (thus the +2)
      AICc <- -2*F["logl"] + 2*K + 2*(K*(K+1))/(n-K-1)
      names(AICc) <- NULL
      
      return(c(AICc=AICc, F[c("cfi", "srmr")], R, R2.p))
    }))
    a1 <- a1[, !grepl(".id", colnames(a1))]
    a1$k <- free.max - a1$Df
    a1$R2.adj <- 1 - ((1-a1$R2))*((N-1)/(N-a1$k-1))
    a1$delta.R2 <- c(NA, a1$R2[1:(nrow(a1)-1)] - a1$R2[2:(nrow(a1))])			
    a1$model <- rownames(a1)
    a1$set <- set
    return(a1)
  }
}


# internal helper function: F-test to compare R^2 difference between two nested models
# x = output object of RSA()
# unrestricted = name of unrestricted model (model with more estimated parameters)
# restricted = name of restricted model (less estimated parameters); with default setting "interceptonly", the function will return R^2 and respective pvalue of the unrestricted model
R2difftest <- function(x, unrestricted="", restricted="interceptonly"){
  
  n <- lavaan::nobs(x$models[[unrestricted]])
  free.max <- getFreeParameters(x$models[[unrestricted]])
  
  Fu <- fitmeasures(x$models[[unrestricted]])
  Ru <- inspect(x$models[[unrestricted]], "r2")
  ku <- free.max - Fu["df"]
  
  if ( restricted=="interceptonly" ){
    Rr <- 0
    kr <- 0
  } else {
    Fr <- fitmeasures(x$models[[restricted]])
    Rr <- inspect(x$models[[restricted]], "r2")
    kr <- free.max - Fr["df"]
  }
  
  suppressWarnings({		
    R2.p <- ifelse(ku==kr,
                   NA,
                   pf( ( ( n - ku - 1 ) * ( Ru - Rr ) ) / ( (ku-kr) * (1 - Ru ) ), ku-kr, n-ku-1, lower.tail=FALSE)
    )
  })
  
  delta.R2 <- Ru-Rr
  names(delta.R2) <- "R2"
  names(R2.p) <- "R2.p"
  
  return(list(
    delta.R2=delta.R2,
    R2.p=R2.p, 
    ku=ku
  ))
  
}


# simple wrapper: formats a number in f.2 format
f2 <- function(x, digits=2, prepoint=0, skipZero=FALSE) {
  
  if (skipZero == TRUE) {zero <- "."} else {zero <- "0."}
  
  if (length(dim(x)) == 2) {
    apply(x, 2, function(x2) {gsub("0.", zero, sprintf(paste("%",prepoint,".",digits,"f",sep=""), x2) , fixed=TRUE)})
  } else {
    gsub("0.", zero, sprintf(paste("%",prepoint,".",digits,"f",sep=""), x) , fixed=TRUE)
  }
}

# converts p values in stars
p2star <- function(val) {
  
  res <- val
  
  for (i in 1:length(val)) {
    res[i] <- ""
    if (is.na(val[i])) next();
    if (val[i] <= 0.1) res[i] <- "\U2020"
    if (val[i] <= 0.05) res[i] <- "*"
    if (val[i] <= 0.01) res[i] <- "**"
    if (val[i] <= 0.001) res[i] <- "***"
  }
  
  return(res)
}

# nicely formats a p-value
p0 <- function(x) {
  if (is.na(x)) return("NA")
  if (x >= .001) return(paste0("p = ", f2(x, 3, skipZero=TRUE)))
  if (x <  .001) return("p <.001")	
}
p <- Vectorize(p0)

# returns number of maximum free parameters of a regression model
getFreeParameters <- function(model) {
  VARS <- nrow(inspect(model, "free")$beta)	# number of variables
  df.max <- (VARS*(VARS+1))/2		# maximum df
  df.pred <- ((VARS-1)*(VARS))/2 + 1 # df bound in the predictors (i.e., (co)variances of the predictors & variance of DV)
  free.max <- df.max - df.pred	# maximum of free parameters
  return(free.max)
}




# computes the coordinates of an arbitrary intersection of the surface,
# defined by a line on the X-Y plane (p0 = intercept, p1=slope)
getIntersect <- function(b0=0, x=0, y=0, x2=0, xy=0, y2=0, p0, p1, xlim=c(-2, 2), grid=21) {
  X <- seq(min(xlim), max(xlim), length.out=grid)
  Y <- p0 + p1*X
  n <- data.frame(X, Y)
  n2 <- add.variables(z~X+Y, n)
  n2$Z <- b0 + colSums(c(x, y, x2, y2, xy)*t(n2[, c(1:5)]))
  return(n2[, c("X", "Y", "Z")])
}


model <- function(x, model="full") x$models[[model]]

syntax <- function(x, model="full") cat(x$models[[model]]@Options$model)




# transforms p-values to colors
pRamp <- function(p, sig=.05, borderline=.10, bias=.8) {
  # calculate bias that the color transition is at the borderline value
  bias2 <- .33/(borderline/(1 - sig))
  cR1 <- colorRamp(c("red", "red", "orange"), bias=bias, space="Lab")
  cR2 <- colorRamp(c("orange", "green", "green"), bias=bias2, space="Lab")
  
  p2 <- rep("#FFFFFF", length(p))
  if (length(p[p < sig])>0) {
    p2[p < sig] <- rgb(cR1(p[p < sig]/sig), maxColorValue=255)
  }
  if (length(p[p >= sig])>0) {
    p2[p >= sig] <- rgb(cR2((p[p >= sig] - sig) / (1 - sig)), maxColorValue=255)
  }
  return(p2)
}


# helper function: find closest value in vector
f0 <- function (vec, target, unique = TRUE) {
  ret <- vec[sapply(target, function(x) which.min(abs(x - vec)))]
  if (unique) { ret <- unique(ret) }
  ret
}

# compute the predicted value from a single pair of predictors
predictRSA <- function(object, X, Y, model="full") {
  C <- coef(object$models[[model]])
  if (object$models[[model]]@Options$estimator != "DWLS") {
    b0 <- as.numeric(ifelse(is.na(C[paste0(object$DV, "~1")]), b0, C[paste0(object$DV, "~1")]))
  } else {
    # the threshold is the negative of the intercept ...
    b0 <- -as.numeric(ifelse(is.na(C[paste0(object$DV, "|t1")]), b0, C[paste0(object$DV, "|t1")]))
  }
  x <- as.numeric(ifelse(is.na(C["b1"]), 0, C["b1"]))
  y <- as.numeric(ifelse(is.na(C["b2"]), 0, C["b2"]))
  x2 <- as.numeric(ifelse(is.na(C["b3"]), 0, C["b3"]))
  y2 <- as.numeric(ifelse(is.na(C["b5"]), 0, C["b5"]))
  xy <- as.numeric(ifelse(is.na(C["b4"]), 0, C["b4"]))
  w <- as.numeric(ifelse(is.na(C["w1"]), 0, C["w1"]))
  wx <- as.numeric(ifelse(is.na(C["w2"]), 0, C["w2"]))
  wy <- as.numeric(ifelse(is.na(C["w3"]), 0, C["w3"]))
  
  # cubic parameters
  x3 <- as.numeric(ifelse(is.na(C["b6"]), 0, C["b6"]))
  x2y <- as.numeric(ifelse(is.na(C["b7"]), 0, C["b7"]))
  xy2 <- as.numeric(ifelse(is.na(C["b8"]), 0, C["b8"]))
  y3 <- as.numeric(ifelse(is.na(C["b9"]), 0, C["b9"]))
  
  
  C <- c(x, y, x2, y2, xy, w, wx, wy,x3, x2y, xy2, y3)
  
  # compute predicted value
  Z <- b0 + colSums(C*t(cbind(X, Y, X^2, Y^2, X*Y, 0, 0, 0, X^3, X^2*Y, X*Y^2, Y^3)))
  return(Z)
}


# fills up the long edges of a polygon with intermediate points
# If an edge is longer than minDist, new points re inserted.
# @param x Vector of x values
# @param y Vector of y values
interpolatePolygon <- function(x, y, minDist, plot=FALSE) {
  minDist <- minDist^2	# compare with squared x^2 + y^2 (faster)
  interp <- data.frame()
  pol <- data.frame(x, y)
  colnames(pol) <- c("x", "y")
  for (i in 1:(nrow(pol)-1)) {
    # get distance
    D <- (pol[i, 1] - pol[i+1, 1])^2 + (pol[i, 2] - pol[i+1, 2])^2
    if (D > minDist) {
      N <- ceiling(sqrt(D)/sqrt(minDist)) # number of interpolations
      APPROX <- data.frame(
        x = seq(pol[i, 1], pol[i+1, 1], length.out=N),
        y = seq(pol[i, 2], pol[i+1, 2], length.out=N)
      )
      interp <- rbind(interp, APPROX)
    } else if (D>0 & D <= minDist){
      interp <- rbind(interp, pol[i, ])
      if (i==1) colnames(interp) <- c("x", "y")
    }
  }
  interp <- rbind(interp, pol[nrow(pol), ])
  if (plot==TRUE) {
    plot(pol, col="red")
    points(interp[, 1], interp[, 2], col="green", pch=20)
    lines(interp[, 1], interp[, 2], col="darkgreen")
    text(pol, label=1:nrow(pol))
  }
  return(interp)
}
