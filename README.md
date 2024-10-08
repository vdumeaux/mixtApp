# MIxT App R Package (mixtApp) 
The mixtApp R package provides the data and statistical analyses used in the compute service for
the MIxT web application (Matched Interactions Across Tissues - [blood:tumor example](http://mixt-blood-tumor.bci.mcgill.ca)). 

The MIxT web application is designed for exploring the results from the MIxT
analysis comparing transcriptional profiles from two matched tissues across
individuals. 

For a more detailed description of the ideas and design of the MIxT web
application please refer to
["Building Applications For Interactive Data Exploration In Systems Biology" by Fjukstad et al.](biorxiv.org/content/early/2017/05/24/141630) 
and the source code [here](https://github.com/fjukstad/mixt). 

If you want to modify the data being used in the web application you'll have to
rebuild the mixtApp package yourself, and start the different services
manually.

First clone down the [mixtApp](https://github.com/vdumeaux/mixtApp)
repository: 

```
$ git clone https://github.com/vdumeaux/mixtApp.git
```
Install dependencies in R
```
install.packages("devtools")
install.packages("ggplot2")

install.packages("BiocManager")
BiocManager::install("Biobase")
install.packages("DBI")
install.packages("RSQLite")
BiocManager::install("AnnotationDbi")
BiocManager::install("GO.db")
install.packages(c("latticeExtra","Hmisc", "reshape"))
BiocManager::install(c("illuminaHumanv3.db", "illuminaHumanv4.db"))
BiocManager::install(c("limma"))
install.packages(c("WGCNA","animation", "plyr", "dplyr", "igraph", "GGally", "network", "sna"))
BiocManager::install(c("breastCancerVDX", "hgu133a.db"))
BiocManager::install(c("genefu", "iC10"))
BiocManager::install("preprocessCore")
```

```
$ cd mixtApp
# modify the data-raw/datasets.R file to load your data. 
$ R -f data-raw/datasets.R
$ R CMD INSTALL .
```
Data must be formatted as described in the [mixtR package](https://github.com/vdumeaux/mixtR). 

Next up is building the compute service
container with your new mixtApp package.

Still in the `mixtApp/` directory, run: 

```
docker build -t compute-service .
```

which build the container for you. Now you can start it up to accept requests on
port `8787` by running 

```
docker run --name=compute-service -t compute-service
```

and it should appear with the `docker ps` command: 

```
 docker ps
CONTAINER ID        IMAGE               COMMAND                  CREATED             STATUS              PORTS               NAMES
08c9889b705e        compute-service     "/bin/sh -c 'go ru..."   4 seconds ago       Up 4 seconds        80/tcp              compute-service
```

The compute service is now running, so next up is starting the web application
container. This is lucily one liner: 

```
docker run -p 8000:80 --link compute-service -e COMPUTE_SERVICE=compute-service:80 --name=mixt -t fjukstad/mixt
```

That's it!  You can now visit the application running on
[localhost:8000](http://localhost:8000). 

If you need more details on the docker commands you can have a look
at[docs.docker.com](https://docs.docker.com).
