pkgLoad <- function( packages = "favourites" ) {

    if( length( packages ) == 1L && packages == "favourites" ) {
        packages <- c( "keras", "stringr", "scclusteval", "Seurat",
                       "ggplot2", "Signac", "scater", "gridExtra", 
                       "biomaRt", "scran", "cowplot", "Matrix",
                       "data.table", "GenomeInfoDb", "EnsDb.Hsapiens.v75", "patchwork", "rhdf5",
                       "Rcpp","rdist","dplyr","ChIPpeakAnno","hypeR"
        )
    }

    packagecheck <- match( packages, utils::installed.packages()[,1] )

    packagestoinstall <- packages[ !is.na( packagecheck ) ]
    packagestoinstall_2 <- packages[is.na( packagecheck ) ]

    if( length( packagestoinstall ) > 0L ) {
        utils::install.packages( packagestoinstall,
                             repos = "http://cran.csiro.au"
        )
    } else {
        print( "All requested packages from CRAN already installed" )
    }
    if( length( packagestoinstall_2 ) > 0L ) {
	if (!requireNamespace("BiocManager", quietly = TRUE))
    	install.packages("BiocManager")
	BiocManager::install(packagestoinstall_2)	
    } else {
        print( "All requested packages from BIOCONDUCTOR already installed" )
    }
    for( package in packages ) {
        suppressPackageStartupMessages(
            library( package, character.only = TRUE, quietly = TRUE )
        )
    }
}
pkgLoad("favourites")