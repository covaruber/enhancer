##################################################################################################
#Startup function
#this function is executed once the library is loaded
.onAttach = function(library, pkg)
{
  Rv = R.Version()
  if(!exists("getRversion", baseenv()) || (getRversion() < "3.5.0"))
    stop("This package requires R 3.5.0 or later")
  if(interactive()) {
    packageStartupMessage(black(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(black(paste("[] Package enhancer 1.1.0 (2025-12)                                 []",sep="")),appendLF=TRUE)
    packageStartupMessage(paste0(black("[] Author: Giovanny Covarrubias-Pazaran",paste0(bgGreen(white(" ")), bgWhite(magenta("M")), bgRed(white(" ")),"  ", bgRed(bold(yellow(" (") )),bgRed(bold(white("W"))), bgRed(bold(yellow(") "))) ) ,"                 []")),appendLF=TRUE)
    packageStartupMessage(black(paste("[]==================================================================[]")),appendLF=TRUE)
  }
  invisible()
}
