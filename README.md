# WNV_temp_dep_trait_RS_validation
Code and data for "Remote sensing of temperature-dependent mosquito and viral traits predicts field surveillance-based disease risk" (DOI: 10.5281/zenodo.13344415).

This dataset contains a bird species list for California that was used to calculate Shannon's diversity within the study region in this manuscript, to use as a covariate in models predicting West Nile virus infection rates in mosquitos in the northern Central Valley. The dataset also contains shapefiles that were used to 1) request ECOSTRESS data for the region of interest ("Northern_Central_Valley.shp") to be used in modeling temperature-dependent mosquito and viral traits for West Nile virus risk, and 2) summarize environmental and vertebrate host data to include as covariates in models predicting mosquito abundance and West Nile virus infection rates as described in the manuscript ("trap_clusters.shp"). This dataset additionally includes the final summarized datasets, associated data dictionary describing each variable in the datasets, and code files required to replicate the analyses and produce the figures in the associated manuscript.

Description of the data and file structure
The species list is a simple list of the species present in California, minus rare species and species not modeled by ebird status and trends. It includes common and scientific name.
The polygon files are simple shapefiles that define regions of interest in the manuscript for the purpose of data processing and analysis as described in the manuscript.
The csv data files "Abundance_data.csv" and "WNV_data.csv" contain columns of variables, with rows associated with unique trap station clusters. The data dictionary csv files "abundance_cols.csv" and "WNV_cols.csv" describe what each variable is measuring.

Sharing/Access information
Bird species data came from the CBRC: https://www.californiabirds.org/default.asp Central Valley boundary file came from the USGS: https://water.usgs.gov/GIS/metadata/usgswrd/XML/pp1766\_Alluvial\_Bnd.xml

Code/Software
Code files "WNV_TempR0_abundance_clean.R" and "WNV_TempR0_transmission_clean.R" contain the code used to replicate the main results and produce figures and tables that appear in the associated manuscript and supplementary information.

Methods
The list of California bird species was downloaded from the California Bird Records Committee, CBRC (https://www.californiabirds.org/default.asp), and filtered to remove rare species (as listed by the CBRC) as well as remove species that were not modeled in the ebird status and trends data for 2019.

The region of interest polygon for the northern Central Valley was obtained by intersecting the alluvial boundary of the Central Valley polygon from USGS (https://water.usgs.gov/GIS/metadata/usgswrd/XML/pp1766_Alluvial_Bnd.xml), with a polygon defining our region of interest for requesting ECOSTRESS data from the NASA APPEEARS portal, as described in the associated manuscript.

The mosquito trap station cluster polygon was created by buffering trap station clusters calculated from vector surveillance data as described in the associated manuscript with a 1500 meter radius buffer.

The csv data files contain estimates of Culex tarsalis abundance and infection rates with West Nile virus, by trap station cluster polygon. The dataset also includes environmental conditions (inclduing remotely sensed temperature-dependent mosquito and viral traits) and bird host abundance and diversity used to model their effects on mosquito and West Nile virus surveillance, as described in detail in the associated manuscript. The data sources and processing are also described in detail in the associated manuscript.

Funding
National Science Foundation, Award: 2011147, DEB
United States Department of Agriculture, Award: 2023-68016-40683, National Institute of Food and Agriculture
National Science Foundation, Award: 2339209, DEB
