SITS - Satellite Image Time Series Analysis for Earth Observation Data
Cubes
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="inst/extdata/sticker/sits_sticker.png" alt="SITS icon" align="right" height="150" width="150"/>

<!-- badges: start -->
<!-- [![Build Status](https://drone.dpi.inpe.br/api/badges/e-sensing/sits/status.svg)](https://drone.dpi.inpe.br/e-sensing/sits) -->

[![CRAN
status](https://www.r-pkg.org/badges/version/sits)](https://cran.r-project.org/package=sits)
[![Build
status](https://cloud.drone.io/api/badges/e-sensing/sits/status.svg)](https://cloud.drone.io/e-sensing/sits)
[![codecov](https://codecov.io/gh/e-sensing/sits/branch/dev/graph/badge.svg?token=hZxdJgKGcE)](https://codecov.io/gh/e-sensing/sits)
[![Documentation](https://img.shields.io/badge/docs-online-blueviolet)](https://e-sensing.github.io/sitsbook/)
[![Life
cycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Software
License](https://img.shields.io/badge/license-GPL--2-green)](https://github.com/e-sensing/sits/blob/master/LICENSE)

<!-- badges: end -->

## Overview

`sits` is an open source R package for satellite image time series
analysis. It enables users to apply machine learning techniques for
classifying image time series obtained from earth observation data
cubes. The basic workflow in `sits` is:

1.  Select an image collection available on cloud providers AWS,
    Microsoft Planetary Computer, Digital Earth Africa and Brazil Data
    Cube.
2.  Build a regular data cube using analysis-ready image collections.
3.  Extract labelled time series from data cubes to be used as training
    samples.
4.  Perform quality control using self-organised maps.
5.  Filtering time series samples for noise reduction.
6.  Use the samples to train machine learning models.
7.  Tune machine learning models for improved accuracy.
8.  Classify data cubes using machine learning models.
9.  Post-process classified images with Bayesian smoothing to remove
    outliers.
10. Estimate uncertainty values of classified images.
11. Evaluate classification accuracy using best practices.
12. Improve results with active learning and self-supervised learning
    methods.

<img src="inst/extdata/markdown/figures/sits_general_view.jpg" title="Conceptual view of data cubes (source: authors)" alt="Conceptual view of data cubes (source: authors)" width="60%" height="60%" style="display: block; margin: auto;" />

## Documentation

Detailed documentation on how to use `sits` is available in the e-book
[“Satellite Image Time Series Analysis on Earth Observation Data
Cubes”](https://e-sensing.github.io/sitsbook/).

## `sits` on Kaggle

Those that want to evaluate the `sits` package before installing are
invited to run the examples available on
[Kaggle](https://www.kaggle.com/esensing/code). These examples provide a
fast-track introduction to the package. We recommend running them in the
following order:

1.  [Introduction to
    SITS](https://www.kaggle.com/esensing/introduction-to-sits)
2.  [Working with time series in
    SITS](https://www.kaggle.com/esensing/working-with-time-series-in-sits)
3.  [Creating data cubes in
    SITS](https://www.kaggle.com/esensing/creating-data-cubes-in-sits)
4.  [Raster classification in
    SITS](https://www.kaggle.com/esensing/raster-classification-in-sits)
5.  [Using SOM for sample quality control in
    SITS](https://www.kaggle.com/esensing/using-som-for-sample-quality-control-in-sits)

## Installation

### Pre-Requisites

The `sits` package relies on the geospatial packages `sf`, `stars`,
`gdalcubes` and `terra`, which depend on the external libraries GDAL and
PROJ. Please follow the instructions for installing `sf` together with
GDAL available at the [RSpatial sf github
repository](https://github.com/r-spatial/sf).

### Obtaining `sits`

`sits` can be installed from CRAN:

``` r
install.packages("sits")
```

The development version is available on github.

``` r
devtools::install_github("e-sensing/sits", dependencies = TRUE)
```

``` r
# load the sits library
library(sits)
#> Using configuration file: /home/sits/R/x86_64-pc-linux-gnu-library/4.2/sits/extdata/config.yml
#> Color configurations found in /home/sits/R/x86_64-pc-linux-gnu-library/4.2/sits/extdata/config_colors.yml
#> To provide additional configurations, create an YAML file and inform its path to environment variable 'SITS_CONFIG_USER_FILE'.
#> Using raster package: terra
#> SITS - satellite image time series analysis.
#> Loaded sits v1.1.0-5.
#>         See ?sits for help, citation("sits") for use in publication.
#>         See demo(package = "sits") for examples.
```

## Building Earth Observation Data Cubes

### Image Collections Accessible by `sits`

The `sits` package allows users to created data cubes from
analysis-ready data (ARD) image collections available in cloud services.
The collections accessible in `sits` 1.1.0.5 are:

1.  Brazil Data Cube
    ([BDC](http://brazildatacube.org/en/home-page-2/#dataproducts)):
    Open data collections of Sentinel-2, Landsat-8 and CBERS-4 images.
2.  Microsoft Planetary Computer
    ([MPC](https://planetarycomputer.microsoft.com/catalog)): Open data
    collection of Sentinel-2/2A and Landsat-8
3.  Earth on AWS ([AWS](https://aws.amazon.com/earth/)): Sentinel-2/2A
    level 2A collections.
4.  Digital Earth Africa
    ([DEAFRICA](https://www.digitalearthafrica.org/)): Open data
    collection of Sentinel-2/2A and Landsat-8 for Africa.
5.  [USGS](https://landsatlook.usgs.gov/stac-browser): Landsat-4/5/7/8
    collections, which are not open data.
6.  Swiss Data Cube ([SDC](https://www.swissdatacube.org/)): Open data
    collection of Sentinel-2/2A and Landsat-8.

Open data collections do not require payment of access fees. Except for
those in the Brazil Data Cube, these collections are not regular.
Irregular collections require further processing before they can be used
for classification using machine learning models.

### Building a Data Cube from an ARD Image Collection

The following code defines an irregular data cube of Sentinel-2/2A
images available in the Microsoft Planetary Computer, using the open
data collection `"SENTINEL-2-L2A"`. The geographical area of the data
cube is defined by the tiles `"20LKP"` and `"20LLKP"`, and the temporal
extent by a start and end date. Access to other cloud services works in
similar ways.

``` r
s2_cube <- sits_cube(
    source = "MPC",
    collection = "SENTINEL-2-L2A",
    tiles = c("20LKP", "20LLP"),
    bands = c("B03", "B08", "B11", "SCL"),
    start_date = as.Date("2018-07-01"),
    end_date = as.Date("2019-06-30"),
    progress = FALSE
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%
```

This cube is irregular. The timelines of tiles `"20LKP"` and `"20LLKP"`
and the resolutions of the bands are different. Sentinel-2 bands `"B03"`
and `"B08"` have 10-meters resolution, while band `"B11"` and the cloud
band `"SCL"` have 20-meters resolution. Irregular collections need an
additional processing step to be converted to regular data cubes, as
described below.

<img src="inst/extdata/markdown/figures/datacube_conception.jpg" title="Conceptual view of data cubes (source: authors)" alt="Conceptual view of data cubes (source: authors)" width="90%" height="90%" style="display: block; margin: auto;" />

After defining an irregular ARD image collection from a cloud service
using `sits_cube()`, users should run `sits_regularize()` to build a
regular data cube. This function uses the [gdalcubes R
package](https://github.com/appelmar/gdalcubes), described in [Appel and
Pebesma, 2019](https://www.mdpi.com/2306-5729/4/3/92).

``` r
gc_cube <- sits_regularize(
    cube          = s2_cube,
    output_dir    = tempdir(),
    period        = "P15D",
    res           = 60, 
    multicores    = 4
)
```

The above command builds a regular data cube with all bands interpolated
to 60 m spatial resolution and 15-days temporal resolution. Regular data
cubes are the input to the `sits` functions for time series retrieval,
building machine learning models, and classification of raster images
and time series.

The cube can be shown in a leaflet using `sits_view()`.

``` r
# View a color composite on a leaflet
sits_view(s2_cube[1,], green = "B08", blue = "B03", red = "B04")
```

## Working with Time Series in `sits`

### Accessing Time Series in Data Cubes

`sits` has been designed to use satellite image time series to derive
machine learning models. After the data cube has been created, time
series can be retrieved individually or by using CSV or SHP files, as in
the following example. The example below uses a data cube in a local
directory, whose images have been obtained from the `"MOD13Q1-6"`
collection of the Brazil Data Cube.

``` r
library(sits)
# this data cube uses images from the Brazil Data Cube that have 
# downloaded to a local directory
data_dir <- system.file("extdata/raster/mod13q1", package = "sits")
# create a cube from downloaded files
raster_cube <- sits_cube(
    source = "BDC",
    collection = "MOD13Q1-6",
    data_dir = data_dir,
    delim = "_",
    parse_info = c("X1", "X2", "tile", "band", "date"),
    progress = FALSE
)
# obtain a set of samples defined by a CSV file
csv_file <- system.file("extdata/samples/samples_sinop_crop.csv",
                        package = "sits")
# retrieve the time series associated with the samples from the data cube
points <- sits_get_data(raster_cube, samples = csv_file)
#> All points have been retrieved
# show the time series
points[1:3,]
#> # A tibble: 3 × 7
#>   longitude latitude start_date end_date   label    cube      time_series      
#>       <dbl>    <dbl> <date>     <date>     <chr>    <chr>     <list>           
#> 1     -55.8    -11.7 2013-09-14 2014-08-29 Cerrado  MOD13Q1-6 <tibble [23 × 2]>
#> 2     -55.8    -11.7 2013-09-14 2014-08-29 Cerrado  MOD13Q1-6 <tibble [23 × 2]>
#> 3     -55.7    -11.7 2013-09-14 2014-08-29 Soy_Corn MOD13Q1-6 <tibble [23 × 2]>
```

After a time series has been obtained, it is loaded in a tibble. The
first six columns contain the metadata: spatial and temporal location,
label assigned to the sample, and coverage from where the data has been
extracted. The spatial location is given in longitude and latitude
coordinates. The first sample has been labelled “Pasture”, at location
(-55.65931, -11.76267), and is considered valid for the period
(2013-09-14, 2014-08-29). To display the time series, use the `plot()`
function.

``` r
plot(points[1,])
```

<img src="man/figures/README-unnamed-chunk-9-1.png" title="Plot of point at location (-55.65931, -11.76267) labelled as Pasture" alt="Plot of point at location (-55.65931, -11.76267) labelled as Pasture" style="display: block; margin: auto;" />

For a large number of samples, where the amount of individual plots
would be substantial, the default visualization combines all samples
together in a single temporal interval.

``` r
# select the "ndvi" band
samples_ndvi <- sits_select(samples_modis_4bands, "NDVI")
# select only the samples with the cerrado label
samples_cerrado <- dplyr::filter(samples_ndvi, 
                  label == "Cerrado")
plot(samples_cerrado)
```

<img src="man/figures/README-unnamed-chunk-10-1.png" title="Samples for NDVI band for Cerrado class" alt="Samples for NDVI band for Cerrado class" style="display: block; margin: auto;" />

## Time Series Clustering and Filtering

### Clustering for sample quality control

Clustering methods in `sits` improve the quality of the samples and to
remove those that might have been wrongly labeled or that have low
discriminatory power. Good samples lead to good classification maps.
`sits` provides support for sample quality control using Self-organizing
Maps (SOM). The process of clustering with SOM is done by
`sits_som_map()`, which creates a self-organizing map and assesses the
quality of the samples.

``` r
# load the kohonen library
library(kohonen)
# create a SOM map from the samples
som_map <- sits_som_map(samples_modis_4bands,
                        grid_xdim = 6,
                        grid_ydim = 6)
# plot the map
plot(som_map)
#> Warning in par(opar): argument 1 does not name a graphical parameter
```

<img src="man/figures/README-unnamed-chunk-11-1.png" title="Samples analysis using SOM (grid 6x6)" alt="Samples analysis using SOM (grid 6x6)" style="display: block; margin: auto;" />

This function uses the [“kohonen” R
package](https://www.jstatsoft.org/article/view/v087i07) to compute a
SOM grid \[7\]. Each sample is assigned to a neuron, and neurons are
placed in the grid based on similarity. Each neuron will be associated
with a discrete probability distribution. Homogeneous neurons (those
with a single class) are assumed to be composed of good quality samples.
Heterogeneous neurons (those with two or more classes with significant
probability) are likely to contain noisy samples. Noisy samples can then
be identified and removed from the sample set using
`sits_som_clean_samples()`.

``` r
# create a new sample set removing noisy points
new_samples <- sits_som_clean_samples(som_map)
```

### Filtering

Satellite image time series are contaminated by atmospheric influence
and directional effects. To make the best use of available satellite
data archives, methods for satellite image time series analysis need to
deal with data sets that are *noisy* and *non-homogeneous*. For data
filtering, `sits` supports Savitzky–Golay (`sits_sgolay()`) and
Whittaker (`sits_whittaker()`) filters. As an example, we show how to
apply the Whittaker smoother to a 16-year NDVI time series.

``` r
# apply Whitaker filter to a time series sample for the NDVI band from 2000 to 2016
# merge with the original data
# plot the original and the modified series
point_ndvi <- sits_select(point_mt_6bands, bands = "NDVI")
point_ndvi %>% 
    sits_filter(sits_whittaker(lambda = 10)) %>% 
    sits_merge(point_ndvi) %>% 
    plot()
```

<img src="man/figures/README-unnamed-chunk-13-1.png" title="Whittaker filter of NDVI time series" alt="Whittaker filter of NDVI time series" style="display: block; margin: auto;" />

## Time Series Classification

### Training Machine Learning Models

`sits` provides support for the classification of both individual time
series as well as data cubes. The following machine learning methods are
available in `sits`:

-   Support vector machines (`sits_svm()`)
-   Random forests (`sits_rfor()`)
-   Extreme gradient boosting (`sits_xgboost()`)
-   Multi-layer perceptrons (`sits_mlp()`)
-   Deep Residual Networks (`sits_resnet()`) (see ref. \[8\])
-   1D convolution neural networks (`sits_tempcnn()`) (see ref. \[9\])
-   Temporal self-attention encoder (`sits_tae()`) (see ref. \[10\])
-   Lightweight temporal attention encoder (`sits_lighttae()`) (see ref.
    \[11\] and \[12\])

The following example illustrate how to train a dataset and classify an
individual time series. First we use the `sits_train()` function with
two parameters: the training dataset (described above) and the chosen
machine learning model (in this case, TempCNN). The trained model is
then used to classify a time series from Mato Grosso Brazilian state,
using `sits_classify()`. The results can be shown in text format using
the function `sits_show_prediction()` or graphically using `plot`.

``` r
# training data set
data("samples_modis_4bands")
# point to be classified
data("point_mt_6bands")
# Select the NDVI and EVI bands 
# Filter the band to reduce noise
# Train a deep learning model
tempcnn_model <- samples_modis_4bands %>% 
    sits_select(bands = "NDVI") %>% 
    sits_train(ml_method = sits_tempcnn()) 
# Select NDVI and EVI bands of the  point to be classified
# Filter the point 
# Classify using TempCNN model
# Plot the result
point_mt_6bands %>% 
  sits_select(bands = "NDVI") %>% 
  sits_classify(tempcnn_model) %>% 
  plot()
```

<img src="man/figures/README-unnamed-chunk-14-1.png" title="Classification of NDVI time series using TempCNN" alt="Classification of NDVI time series using TempCNN" style="display: block; margin: auto;" />

The following example shows how to classify a data cube organized as a
set of raster images. The result can also be visualized interactively
using `sits_view()`.

``` r
# Create a data cube to be classified
# Cube is composed of MOD13Q1 images from the Sinop region in Mato Grosso (Brazil)
data_dir <- system.file("extdata/raster/mod13q1", package = "sits")
sinop <- sits_cube(
    source = "BDC",
    collection = "MOD13Q1-6",
    data_dir = data_dir,
    delim = "_",
    parse_info = c("X1", "X2", "tile", "band", "date"),
    progress = FALSE
)
# Classify the raster cube, generating a probability file
# Filter the pixels in the cube to remove noise
probs_cube <- sits_classify(sinop, ml_model = tempcnn_model)
# apply a bayesian smoothing to remove outliers
bayes_cube <- sits_smooth(probs_cube)
# generate a thematic map
label_cube <- sits_label_classification(bayes_cube)
# plot the the labelled cube
plot(label_cube, title = "Land use and Land cover in Sinop, MT, Brazil in 2018")
```

<img src="man/figures/README-unnamed-chunk-15-1.png" title="Land use and Land cover in Sinop, MT, Brazil in 2018" alt="Land use and Land cover in Sinop, MT, Brazil in 2018" style="display: block; margin: auto;" />

## Additional information

For more information, please see the on-line book [“SITS: Data analysis
and machine learning for data cubes using satellite image time
series”](https://e-sensing.github.io/sitsbook/).

### References

#### Citable papers for sits

If you use `sits`, please cite the following paper:

-   \[1\] Rolf Simoes, Gilberto Camara, Gilberto Queiroz, Felipe Souza,
    Pedro R. Andrade, Lorena Santos, Alexandre Carvalho, and Karine
    Ferreira. “Satellite Image Time Series Analysis for Big Earth
    Observation Data”. Remote Sensing, 13: 2428, 2021.
    <doi:10.3390/rs13132428>.

Additionally, the sample quality control methods that use self-organized
maps are described in the following reference:

-   \[2\] Lorena Santos, Karine Ferreira, Gilberto Camara, Michelle
    Picoli, Rolf Simoes, “Quality control and class noise reduction of
    satellite image time series”. ISPRS Journal of Photogrammetry and
    Remote Sensing, 177:75-88, 2021.
    <doi:10.1016/j.isprsjprs.2021.04.014>.

#### Papers that use sits to produce LUCC maps

-   \[3\] Rolf Simoes, Michelle Picoli, et al., “Land use and cover maps
    for Mato Grosso State in Brazil from 2001 to 2017”. Sci Data
    7(34), 2020. <doi:10.1038/s41597-020-0371-4>.

-   \[4\] Michelle Picoli, Gilberto Camara, et al., “Big Earth
    Observation Time Series Analysis for Monitoring Brazilian
    Agriculture”. ISPRS Journal of Photogrammetry and Remote
    Sensing, 2018. <doi:10.1016/j.isprsjprs.2018.08.007>.

-   \[5\] Karine Ferreira, Gilberto Queiroz et al., “Earth Observation
    Data Cubes for Brazil: Requirements, Methodology and Products”.
    Remote Sens. 12:4033, 2020. <doi:10.3390/rs12244033>.

#### Papers that describe software used in sits

We thank the authors of these papers for making their code available to
be used in connection with sits.

-   \[6\] Marius Appel and Edzer Pebesma, “On-Demand Processing of Data
    Cubes from Satellite Image Collections with the Gdalcubes Library.”
    Data 4 (3): 1–16, 2020. <doi:10.3390/data4030092>.

-   \[7\] Ron Wehrens and Johannes Kruisselbrink, “Flexible
    Self-Organising Maps in kohonen 3.0”. Journal of Statistical
    Software, 87(7), 2018. <doi:10.18637/jss.v087.i07>.

-   \[8\] Hassan Fawaz, Germain Forestier, Jonathan Weber, Lhassane
    Idoumghar, and Pierre-Alain Muller, “Deep learning for time series
    classification: a review”. Data Mining and Knowledge Discovery,
    33(4): 917–963, 2019. \<arxiv:1809.04356\>.

-   \[9\] Charlotte Pelletier, Geoffrey I. Webb, and Francois Petitjean.
    “Temporal Convolutional Neural Network for the Classification of
    Satellite Image Time Series.” Remote Sensing 11 (5), 2019.
    <doi:10.3390/rs11050523>.

-   \[10\] Vivien Garnot, Loic Landrieu, Sebastien Giordano, and Nesrine
    Chehata, “Satellite Image Time Series Classification with Pixel-Set
    Encoders and Temporal Self-Attention”, Conference on Computer Vision
    and Pattern Recognition, 2020. \<doi:
    10.1109/CVPR42600.2020.01234\>.

-   \[11\] Vivien Garnot, Loic Landrieu, “Lightweight Temporal
    Self-Attention for Classifying Satellite Images Time Series”, 2020.
    \<arXiv:2007.00586\>.

-   \[12\] Maja Schneider, Marco Körner, “\[Re\] Satellite Image Time
    Series Classification with Pixel-Set Encoders and Temporal
    Self-Attention.” ReScience C 7 (2), 2021.
    <doi:10.5281/zenodo.4835356>.

#### R packages used in sits

The authors are thankful for the contributions of Marius Appel, Tim
Appelhans, Henrik Bengtsson, Robert Hijmans, Edzer Pebesma, and Ron
Wehrens, respectively chief developers of the packages `gdalcubes`,
`leafem`, `data.table`, `terra/raster`, `sf`/`stars`, and `kohonen`. The
`sits` package is also much indebted to the work of the RStudio team,
including the `tidyverse`. We are indepted to Daniel Falbel for his and
the `torch` packages. We thank Charlotte Pelletier and Hassan Fawaz for
sharing the python code that has been reused for the TempCNN and ResNet
machine learning models. We would like to thank Maja Schneider for
sharing the python code that helped the implementation of the
`sits_lighttae()` and `sits_tae()` model. We recognise the importance of
the work by Chris Holmes and Mattias Mohr on the STAC specification and
API.

## Acknowledgements for Financial and Material Support

We acknowledge and thank the project funders that provided financial and
material support:

1.  Amazon Fund, established by the Brazilian government with financial
    contribution from Norway, through the project contract between the
    Brazilian Development Bank (BNDES) and the Foundation for Science,
    Technology and Space Applications (FUNCATE), for the establishment
    of the Brazil Data Cube, process 17.2.0536.1.

2.  Coordenação de Aperfeiçoamento de Pessoal de Nível Superior-Brasil
    (CAPES) and from the Conselho Nacional de Desenvolvimento Científico
    e Tecnológico (CNPq), for providing MSc and PhD scholarships.

3.  Sao Paulo Research Foundation (FAPESP) under eScience Program grant
    2014/08398-6, for for providing MSc, PhD and post-doc scholarships,
    equipment, and travel support.

4.  International Climate Initiative of the Germany Federal Ministry for
    the Environment, Nature Conservation, Building and Nuclear Safety
    (IKI) under grant 17-III-084- Global-A-RESTORE+ (“RESTORE+:
    Addressing Landscape Restoration on Degraded Land in Indonesia and
    Brazil”).

5.  Microsoft Planetary Computer under the GEO-Microsoft Cloud Computer
    Grants Programme.

## How to contribute

The `sits` project is released with a [Contributor Code of
Conduct](https://github.com/e-sensing/sits/blob/master/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.
