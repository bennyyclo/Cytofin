
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cytofin

CytofIn (CyTOF integration) is an R package for homogenizing and
normalizing heterogeneous [mass cytometry
(CyTOF)](https://pubmed.ncbi.nlm.nih.gov/21551058/) data from diverse
data sources. Specifically, CytofIn provides functions that perform the
following tasks:

-   **Dataset homogenization** - CyTOF datasets that were collected
    separately may differ in which markers were included in their
    antibody panels; in addition, they may use different naming
    conventions for their panels’ shared markers. Thus, data mining
    across multiple CyTOF datasets requires **homogenization,** the
    process of aligning each dataset’s antibody panels so that they can
    be analyzed together. In CytofIN, data homogenization (i.e. panel
    alignment) is performed with the `cytofin_homogenize()` function
    that leverages user-provided panel information to combine datasets.
-   **Dataset normalization** - Combined analysis of multiple CyTOF
    datasets is likely to be confounded by dataset-to-dataset batch
    effects due to differences in instrumentation and experimental
    protocols between groups. To normalize multiple CyTOF datasets with
    respect to these batch effects, CytofIN provides 3 functions:
    `cytofin_prep_anchors()`, `cytofin_normalize()`, and
    `cytofin_normalize_nrs()`.

The general CytofIn workflow unfolds in 3 steps. First, users align the
panels of the CyTOF datasets being integrated using
`cytofin_homogenize()`. Second, users generate reference statistics from
“generalized anchors” identified on each CyTOF plate (see below) using
`cytofin_prep_anchors()`. Finally, users can then normalize/batch
correct the datasets relative to one another using their choice of
`cytofin_normalize()` or `cytofin_normalize_nrs()`, each of which
performs the normalization procedure differently (see below).

## Installation

To install CytofIn, run the following code:

``` r
library(devtools)
install_github("bennyyclo/Cytofin")
```

To attach the CytofIn package to your current R session, run the
following line:

``` r
library(cytofin)
```

## Usage

### CyTOF data homogenization

Here, the term “homogenization” refers to the process of aligning the
antigen panels of multiple CyTOF experiments by (1) removing all
channels that are not shared across all cohorts and (2) standardizing
the antigen names used to refer to each channel so that existing
analysis tools (like the `flowCore` and `tidyverse` packages) can be
applied in later analytical steps. In CytofIn, dataset homogenization is
performed using the `cytofin_homogenize()` function.

The `cytofin_homogenize()` function takes several arguments. The first
of these is `metadata_path`, a string that specifies the file path to a
.csv or .xlsx metadata file containing information about each of the
.fcs files being analyzed. Specifically, the metadata file will have one
row for each .fcs file being analyzed and must contain the following
columns (all of which will be converted to character vectors):

-   **filename -** Required. The name of the .fcs file within its local
    directory.
-   **cohort -** Required. The name of the cohort (i.e. experimental
    source) of each .fcs file.
-   **plate\_number -** Required. The name of the CyTOF plate (e.g.
    “plate1”, “plate2”, etc.) on which the sample corresponding to each
    .fcs file was analyzed during data acquisition.
-   **patient\_id -** Optional. The name of the patient to whom each
    .fcs file corresponds.
-   **condition -** Optional. The stimulation condition corresponding to
    each .fcs file (i.e. “basal”, “IL-3”, etc.).
-   **is\_anchor -** Required. A numeric column indicating whether or
    not each sample should be used as an “anchor” for the batch
    correction procedure (1 if yes; 0 if no). Exactly one anchor should
    be identified for each CyTOF plate being analyzed.
-   **validation -** Optional. The name of the
    [bead-normalized](https://pubmed.ncbi.nlm.nih.gov/23512433/) .fcs
    file corresponding to each input file listed in the `filename`
    column (per the gold standard in CyTOF batch correction). Most users
    will ignore this column (because bead-normalized data will not be
    available), but it can be used to validate the results of the
    CytofIn batch normalization algorithms if bead-normalized data are
    available.

Importantly, only the fields marked as “required” are needed for
`cytofin_homogenize()` to work; “NA” can be recorded for any/all
optional columns that don’t apply to the experimental design of the
files being analyzed (for example, if no stimulation conditions were
used in the studies being integrated, enter “NA” for each element of the
`condition` column). Alternatively, these columns can be omitted from
the metadata table entirely. The `cytofin_generate_metadata_template`
function is provided to generate an example metadata .csv file filled
with dummy example data in a location specified by the user:

``` r
# specify the path where you'd like to store the template file
my_path <- file.path("~", "Desktop", "template_folder")

# generate the template file, which then can be edited manually 
cytofin_generate_metadata_template(template_path = my_path)
```

The second argument for `cytofin_homogenize` is `panel_path`, a string
that specifies the file path to a .csv or .xlsx file containing
information about the panel(s) of each of the .fcs files being analyzed.
Each row represents a channel (i.e. a protein measurement) to be
included in the final, homogenized panel. This file must contain the
following columns:

-   **metal\_name -** A character vector representing the name of the
    metal isotope measured by each channel.
-   **antigen\_name -** A character vector representing the name of the
    antigen associated with a given metal isotope in the consensus panel
    (the final antigen name to assign to a given channel during
    homogenization).
-   **antigen\_pattern -** A regular expression used to match antigen
    names that may differ slightly across different .fcs files. For
    example, the regular expression “(C\|c)(D\|d)45” will detect all of
    the following channel names: “cd45”, “CD45”, “Cd45”, “cD45”.
-   **lineage -** A numeric vector representing whether or not a marker
    is a lineage marker (1 if yes; 0 otherwise).
-   **functional -** A numeric vector representing whether or not a
    marker is a functional marker (1 if yes; 0 otherwise).
-   **general -** A numeric vector representing whether or not a marker
    is a “general” (i.e. neither a lineage nor a functional) marker (1
    if yes; 0 otherwise).

The layout of this antigen panel is displayed graphically below.

![](./images/Slide1.png?raw=true "Title")

As above, the `cytofin_generate_panel_template` function is provided to
generate an example metadata .csv file filled with dummy example data:

``` r
# specify the path where you'd like to store the template file
my_path <- 
  file.path("~", "Desktop", "template_folder")

# generate the template file, which then can be edited manually 
cytofin_generate_panel_template(template_path = my_path)
```

For many users, the most difficult part of filling out the consensus
panel information table will be designing the regular expressions for
the `antigen_pattern` column. However, in most cases the required
regular expressions will be quite simple: For a primer on regular
expressions (and their use in the
[`stringr`](https://stringr.tidyverse.org/) package) written by
[RStudio](https://www.rstudio.com/about/), install the `stringr` package
and read the following vignette:

``` r
vignette(topic = "regular-expressions", package = "stringr")
```

The next two arguments for `cytofin_homogenize` are `input_data_path`
and `output_data_path`, two strings that indicate which directory input
.fcs files should be read from and which directory homogenized .fcs
files should be written to, respectively. Lastly, the final two
arguments are optional: `prefix` allows the user to specify the prefix
appended to each input .fcs file name to get the name of the
corresponding output (i.e. homogenized) .fcs file name, and `verbose` is
a boolean value (default = FALSE) specifying if chatty print statements
should be made while the homogenization is performed.

Using these arguments, `cytofin_homogenize` can homogenize a set of
CyTOF files with distinct antigen naming conventions. Specifically, the
program performs a regular expression search to match the synonymous
term in the panel and correct the antigen name with standardized names
in the panel.

Example function call:

``` r
# define input paths 
metadata_path <- 
  "/Users/tkeyes/GitHub/cytofin/inst/extdata/test_metadata_raw.csv"
panel_path <- 
  "/Users/tkeyes/GitHub/cytofin/inst/extdata/test_panel.csv"
input_data_path <- 
  "/Users/tkeyes/GitHub/cytofin/inst/extdata/test_raw_fcs_files"

# define output path
# --Change this line to wherever you want the output files saved!--
output_data_path <- file.path("~", "temp", "out_test")

# call homogenization function
cytofin_homogenize(
  metadata_path = metadata_path, 
  panel_path = panel_path, 
  input_data_path = input_data_path, 
  output_data_path = output_data_path
)
```

This function call will save homogenized .fcs files to the directory
located at `output_data_path`. These files will be different from the
input .fcs files in the `input_data_path` directory in that they will
only contain channels whose antigen names match the `antigen_pattern`
column of the reference panel located at `panel_path`. All other
channels will be removed, and the names of the channels with matches in
`antigen_pattern` will be standardized to the names given in the
`antigen_name` column of the reference panel.

The input files for this homogenization run were as follows:

``` r
list.files(input_data_path, pattern = ".fcs$")
#>  [1] "ALL05v2_Plate2_healthy basal1.fcs" "ALL05v2_Plate2_UPN94 das.fcs"     
#>  [3] "ALL08_Plate8_Healthy03 basal.fcs"  "ALL08_Plate8_UPN26 basal.fcs"     
#>  [5] "CRLF2_Plate1_Healthy 04 BCR.fcs"   "CRLF2_Plate1_UPN53 das + TSLP.fcs"
#>  [7] "MS_Plate5_Healthy BM.fcs"          "MS_Plate5_SU978 Basal.fcs"        
#>  [9] "SJ_Plate2_Healthy_BM.fcs"          "SJ_Plate2_TB010950_Basal.fcs"
```

…and the corresponding output file saved in the `output_data_path`
directory are now as follows:

``` r
list.files(output_data_path, pattern = ".fcs$")
#>  [1] "final_outputnormalized_ALL05v2_Plate2_UPN94 das.fcs"
#>  [2] "homogenized_ALL05v2_Plate2_healthy basal1.fcs"      
#>  [3] "homogenized_ALL05v2_Plate2_UPN94 das.fcs"           
#>  [4] "homogenized_ALL08_Plate8_Healthy03 basal.fcs"       
#>  [5] "homogenized_ALL08_Plate8_UPN26 basal.fcs"           
#>  [6] "homogenized_CRLF2_Plate1_Healthy 04 BCR.fcs"        
#>  [7] "homogenized_CRLF2_Plate1_UPN53 das + TSLP.fcs"      
#>  [8] "homogenized_MS_Plate5_Healthy BM.fcs"               
#>  [9] "homogenized_MS_Plate5_SU978 Basal.fcs"              
#> [10] "homogenized_SJ_Plate2_Healthy_BM.fcs"               
#> [11] "homogenized_SJ_Plate2_TB010950_Basal.fcs"
```

### CyTOF batch normalization using external anchors

After dataset homogenization, **batch correction** (or **batch
normalization**) can be performed across datasets. The following
schematic diagram illustrates how normalization is performed in
`CytofIn`:

![Alt text](./images/Slide2.PNG?raw=true "Title")

TO DO: Fix a few things about this schematic diagram

In words, CytofIn performs batch normalization though the use of
**generalized anchors -** or \[definition of generalized anchors\] **-**
that users identify on each CyTOF plate. To identify a sample as a
generalized anchor in CytofIn, users must mark its row in the metadata
table with a “1” in the `is_anchor` column. These anchors are used to
define a **universal mean** and/or **universal variance** that represent
the central tendency and dispersion, respectively, of the target
distribution to which all samples will be batch corrected using the
user’s choice from one of five batch correction functions.

In other words, CytofIn’s batch normalization procedure has two steps:

1.  Preparation of external anchors  
2.  Application of a transformation function that performs the batch
    correction (of which `CytofIn` provides 5 options)

We detail each of these steps below.

#### 1. Anchor preparation:

The `cytofin_prep_anchors` function concatenates the identified anchor
files (one file per plate) and then calculates summary statistics that
are used for batch correction in later steps of the pipeline. First,
CytofIn calculates the mean and standard deviation of each channel in
the homogenized dataset across all cells from samples identified as
generalized anchors. These values represent the overall central tendency
and dispersion, respectively, of each channel among the healthy control
samples on each CyTOF plate and are thus termed the **universal means**
and **universal variances** of the analysis. Accordingly, the universal
mean and universal variance vectors will each have *g* elements, where
*g* is the number of channels in the consensus antigen panel in the
panel information table. The universal mean and universal variance
vectors are used in the `meanshift`, `variance`, and `z-score` methods
of batch correction (see below).

In addition, the mean of all of the elements of the universal mean
vector (i.e. the mean of all channel means) and the mean of all of the
elements of the universal variance vector (i.e. the mean of all channel
variances) are calculated. These values represent the central tendency
and dispersion of antigen measurements in general among the healthy
control samples on each CyTOF plate and are no longer channel-specific.
Thus, we call them the *bulk mean* and *bulk variance*, as they are used
in the `meanshift_bulk` batch correction method implemented in
`cytofin_homogenize` (see below), which is not channel-specific.

`cytofin_prep_anchors` returns the universal mean vector, universal
standard deviation vector, bulk mean, and bulk standard deviation as a
`list()`. In addition, users are given an option to save these
statistics as an .rds file in a specified directory in order to avoid
performing redundant calculations in future analyses.

Specifically, `cytofin_prep_anchors` takes 4 required arguments:

-   `metadata_path`: A connection leading to an .xlsx or .csv file
    containing a metadata table with information about each file to be
    analyzed. This file should be identical to that used for
    `cytofin_homogenize`.
-   `panel_path`: A connection leading to an .xlsx or .csv file
    containing information about the standardized antigen panel in the
    homogenized dataset. This file should be identical to that used for
    `cytofin_homogenize`.
-   `input_data_path`: A connection to a directory containing the input
    .FCS files from which to draw summary statistics
-   `output_path`: A connection to a directory where the output .rds and
    .FCS files will be saved. The default is “none”, in which case no
    output files will be stored (and the only effect of the function
    will be to return the calculated statistics as a `list()`).

In addition, `cytofin_prep_anchors` also takes 2 optional arguments
relating to the conventional arcsinh transformation performed on the raw
ion counts of the input data. These optional arguments are as follows:

-   `shift_factor`: The scalar value `a` in the following equation used
    to transform CyTOF raw data ion counts using the hyperbolic arcsinh
    function: `new_x <- asinh(a + b * x)`. Defaults to 0.

-   `scale_factor`: The scalar value `b` in the following equation used
    to transform CyTOF raw data ion counts using the hyperbolic arcsinh
    function: `new_x <- asinh(a + b * x)`. Defaults to 0.2.

Finally, here is an example functional call of `cytofin_prep_anchors`:

``` r
input_data_path <- file.path("~", "temp", "out_test")

anchor_statistics <- 
  cytofin_prep_anchors(
    metadata_path = metadata_path, 
    panel_path = panel_path, 
    input_data_path = file.path("~", "temp", "out_test"), 
    output_path = file.path("~", "temp", "out_test", "anchor_prep")
  )

print(anchor_statistics)
#> $universal_var
#>         Time Event_length    (Pd102)Di    (Pd104)Di    (Pd105)Di    (Pd106)Di 
#>   1.28235792   0.16399756   6.78770451   0.89290897   5.74351522   4.00916670 
#>    (Pd108)Di    (Pd110)Di    (In113)Di    (In115)Di    (La139)Di    (Pr141)Di 
#>   6.47944462   6.14839951   3.14291787   3.69776978   0.31651260   0.20067263 
#>    (Nd142)Di    (Nd143)Di    (Nd144)Di    (Nd145)Di    (Nd146)Di    (Sm147)Di 
#>   0.88280840   0.50837979   0.18512779   0.27893442   0.79089548   1.30174061 
#>    (Nd148)Di    (Sm149)Di    (Nd150)Di    (Sm152)Di    (Eu153)Di    (Sm154)Di 
#>   1.53148051   0.24234410   0.19237185   0.78984151   3.36668746   0.64687396 
#>    (Gd156)Di    (Gd158)Di    (Gd160)Di    (Dy161)Di    (Dy162)Di    (Dy163)Di 
#>   0.62963342   0.21865740   2.88801028   0.07940630   0.12194444   0.07128214 
#>    (Dy164)Di    (Ho165)Di    (Er166)Di    (Er167)Di    (Er168)Di    (Er170)Di 
#>   0.44285804   1.04235848   0.28206380   4.31831331   3.59089444   3.35406088 
#>    (Yb171)Di    (Yb172)Di    (Yb173)Di    (Yb174)Di    (Lu175)Di    (Yb176)Di 
#>   1.95310084   0.67905696   0.13911985   6.12832312   1.77734024   0.53625671 
#>    (Ir191)Di    (Ir193)Di 
#>   3.21574811   3.27089639 
#> 
#> $universal_mean
#>         Time Event_length    (Pd102)Di    (Pd104)Di    (Pd105)Di    (Pd106)Di 
#>  14.50995327   2.30820954   3.48055714   1.06062913   4.08199057   4.77092034 
#>    (Pd108)Di    (Pd110)Di    (In113)Di    (In115)Di    (La139)Di    (Pr141)Di 
#>   2.69248853   3.31279576   1.34656332   2.31588156   0.35046633   0.19319399 
#>    (Nd142)Di    (Nd143)Di    (Nd144)Di    (Nd145)Di    (Nd146)Di    (Sm147)Di 
#>   0.57791130   0.34730008   0.20086489   0.34646560   0.61382685   0.56851774 
#>    (Nd148)Di    (Sm149)Di    (Nd150)Di    (Sm152)Di    (Eu153)Di    (Sm154)Di 
#>   1.13302732   0.15299272   0.19208744   0.43406391   2.13362865   0.45270859 
#>    (Gd156)Di    (Gd158)Di    (Gd160)Di    (Dy161)Di    (Dy162)Di    (Dy163)Di 
#>   0.34711746   0.17472376   1.30261426   0.11212254   0.13257570   0.07266354 
#>    (Dy164)Di    (Ho165)Di    (Er166)Di    (Er167)Di    (Er168)Di    (Er170)Di 
#>   0.22465161   0.48758658   0.28522175   2.63843957   2.43044297   0.80540655 
#>    (Yb171)Di    (Yb172)Di    (Yb173)Di    (Yb174)Di    (Lu175)Di    (Yb176)Di 
#>   1.30095098   0.65077576   0.15830507   2.43474419   1.14821570   0.54578885 
#>    (Ir191)Di    (Ir193)Di 
#>   3.80031272   4.48577210 
#> 
#> $bulk_var
#> [1] 1.467075
#> 
#> $bulk_mean
#> [1] 0.969387
```

As shown above, the returned value is a list with 4 items in it: the
universal variance vector (`universal_var`), the universal mean vector
(`universal_mean`), the bulk variance (`bulk_var`) and the bulk mean
(`bulk_mean`). Note that the elements of `universal_var` and
`universal_mean` are named with their corresponding metal names (not
antigen names), as this interfaces a bit more conveniently with the
`flowCore` functions that CytofIn uses under-the-hood.

Importantly, you only need to use `cytofin_prep_anchors` if you plan to
batch normalize your .fcs files using external anchors identified on
each plate (using `cytofin_normalize`). If you plan to batch normalize
your .fcs files using non-redundancy scores from each sample’s most
stable channels (using `cytofin_normalize_nrs`),

#### 2. Batch correction

After the generalized anchors’ summary statistics are computed, batch
correction can be performed using either `cytofin_normalize` or
`cytofin_normalize_nrs`. These two functions perform batch correction in
different ways, and which of them is most applicable to a given analysis
will differ from user to user. We recommended that users try using both
and then manually inspect/visualize the batch-corrected data in order to
determine which method they prefer.

##### cytofin\_normalize: Batch correction using external anchors on each plate

To perform batch correction using external anchors identified on each
plate, use `cytofin_normalize`. This batch normalization strategy
assumes that the anchors on each plate are relatively similar to one
another, and it uses this similarity to adjust the marker expression
measurements on each plate based on how much each plate’s anchor differs
from the other anchors. The `cytofin_normalize` function takes several
required arguments:

-   Start

-   Start

-   Start

In addition to these required arguments, `cytofin_normalize` takes
several optional arguments:

-   Start

-   Start

-   Start

Using these arguments, a call to `cytofin_normalize` will perform the
batch correction and save the output (i.e. normalized) .fcs files to the
directory specified by `output_data_path`. An example function call is
given here:

``` r
output_data_path <- 
  file.path("~", "temp", "out_test", "anchor_prep")

norm_result <- 
  cytofin_normalize(
    metadata_path = metadata_path, 
    panel_path = panel_path, 
    anchor_statistics = anchor_statistics, 
    input_data_path = input_data_path, 
    output_data_path = output_data_path, 
    mode = "meanshift"
  )
#> Warning in dir.create(output_data_path): '/Users/tkeyes/temp/out_test/
#> anchor_prep' already exists
```

\[Some description of the results that are returned by this function
call\]

\[some description of how these results can be used to make plots of the
normalization results using `cytofin_make_plots`:

``` r
# we make just the first plot for illustrative purposes
cytofin_make_plots(
  normalization_result = dplyr::slice(norm_result, 1), 
  val_path = "none"
)
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

##### cytofin\_normalize\_nrs: Batch correction using stable channels (internal anchors) within each .fcs file

![Alt text](./images/Slide3.PNG?raw=true "Title")

Description: In the event that the external references are not
available, internal anchors can be used. Here, we identified the most
stable channels as internal anchors using a PCA-based non-redundancy
score (NRS). A minimal of three channels should be selected to establish
an internal reference from which signals can be calibrated between CyTOF
files.

Function definition:

``` r
# local paths
metadata_path <- 
  file.path("~", "GitHub", "cytofin", "inst", "extdata", "test_metadata_raw.csv")
panel_path <- panel_path <- 
  here::here("inst", "extdata", "test_panel.csv")
input_data_path <- 
  file.path("~", "temp", "out_test")
output_data_path <- 
  file.path("~", "temp", "out_test", "final_output")

# call function
norm_result_nrs <- 
  cytofin_normalize_nrs(
    metadata_path = metadata_path, 
    panel_path = panel_path, 
    input_data_path = input_data_path, 
    output_data_path = output_data_path, 
    nchannels = 5, 
    make_plot = FALSE
  )
#> Warning in dir.create(output_data_path): '/Users/tkeyes/temp/out_test/
#> final_output' already exists
```
