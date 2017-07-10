# MIxT App R Package (mixtApp)
R package that provides the statistical analyses used in the compute service for
the [MIxT](http://mixt-blood-tumor.bci.mcgill.ca) web application. 

# Install 
Using R 
```
> devtoools::install_github("vdumeaux/mixtApp")
```

or with the shell 

```
$ git clone https://github.com/vdumeaux/mixtApp.git
$ R CMD INSTALL mixtApp
```

if you plan on using it with your data later. 


# Data? 
If you want to use the package to analyze your own data, you'll have to replace 
the data in the `data/` folder with your own data and rebuild the package. We 
use the `data-raw/datasets.R` import script to retrieve data and place it in the
`data/` folder in the R package. Modify this script to load your data and
install the package. Below is a short code snippet that shows how to do it
yourself.

```
$ git clone git@github.com:vdumeaux/mixtApp.git 
$ cd mixtApp
# modify the data-raw/datasets.R file to load your data. 
$ R -f data-raw/datasets.R
$ R CMD INSTALL .
```
