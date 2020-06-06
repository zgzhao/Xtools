bnw_barplot <- function(height, shade=1:7, ...) {
    args <- list(...)
    angle1 <- c(0, 45, 135, 88, 0, 45, 180)[shade]
    angle2 <- c(0, 45, 135, 88, 0, 135, 90)[shade]
    dens <- c(0, rep(10, 6))[shade]
    barplot(height, angle=angle1, density=dens, col='gray40', ...)
    barplot(height, angle=angle2, density=dens, col='gray40', add=TRUE)
}

bnw_legend <- function(x, y=NULL, legend, shade=1:7, ...) {
    angle1 <- c(0, 45, 135, 88, 0, 45, 180)[shade]
    angle2 <- c(0, 45, 135, 88, 0, 135, 90)[shade]
    dens <- c(0, rep(30, 6))[shade]
    legend(x, y, legend, angle=angle1, density=dens, ...)
    legend(x, y, legend, angle=angle2, density=dens, ...)
}

## shade <- c(2, 4, 5, 1, 3)
## bnw_barplot(1:5, shade=shade)
## bnw_legend("topleft", legend=LETTERS[1:6], shade=shade)
