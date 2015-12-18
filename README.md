# INSTALL
The `Xtools` package contains some R functions often use in our lab. You can use it freely.

The package can be install with `devtools` from `github` by execute the following commands in R:

```
## install.packages('devtools')
require(devtools)
install_github('zgzhao/Xtools')
```
# LT50s shiny app

After `Xtools` installation, you can use the shiny app as:

```
require(shiny)
xpp <- system.file("examples", "LT50s", package="Xtools")
runApp(xpp)
```
