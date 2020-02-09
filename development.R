require(devtools)
load_all()
check(document = FALSE)

build_vignettes()

# build with: 
# bash R-devel CMD build LaplacesDemon --compact-vignettes="gs+qpdf"
devtools::build(args = '--compact-vignettes=gs+qpdf')

#
devtools::revdep_maintainers()
