.First.lib <- function(libname, pkgname) {

    library.dynam("Harshlight", pkgname, libname)

}

.Last.lib <- function(libpath) {

	library.dynam.unload("Harshlight", libpath)

}
