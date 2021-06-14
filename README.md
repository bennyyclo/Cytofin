
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

-   **filename -** The name of the .fcs file within its local directory.
-   **cohort -** The name of the cohort (i.e. experimental source) of
    each .fcs file.
-   **plate\_number -** The name of the CyTOF plate (e.g. “plate1”,
    “plate2”, etc.) on which the sample corresponding to each .fcs file
    was analyzed during data acquisition.
-   **patient\_id -** The name of the patient to whom each .fcs file
    corresponds.
-   **condition -** The stimulation condition corresponding to each .fcs
    file (i.e. “basal”, “IL-3”, etc.).
-   **population -** The cell population contained in each .fcs file
    (i.e. “T-cells”, “B-cells”, etc.).
-   **validation -** TO DO - I’m not sure what this column is?

All of these fields are required for `cytofin_homogenize()` to work;
“NA” can be recorded for any/all columns that don’t apply to the
experimental design of the files being analyzed (for example, if no
stimulation conditions were used in the studies being integrated, enter
“NA” for each element of the `condition` column). The
`cytofin_generate_metadata_template` function is provided to generate an
example metadata .csv file filled with dummy example data:

``` r
# specify the path where you'd like to store the template file
my_path <- file.path("~", "Desktop", "template_folder", "metadata_template.csv")

# generate the template file, which then can be edited manually 
cytofin_generate_metadata_template(template_path = my_path)
```

The second argument for `cytofin_homogenize` is `panel_path`, a string
that specifies the file path to a .csv or .xlsx file containing
information about the panel(s) of each of the .fcs files being analyzed.
Each row represents a channel (i.e. a protein measurement) to be
included in the final, homogenized panel. This file must contain the
following columns:

-   **desc -** TO DO ???
-   **range -** TO DO ???
-   **metal\_pattern -** A regular expression used to \_\_\_\_\_. I
    DON’T THINK THIS IS EVER USED IN THE ACTUAL HOMOGENIZATION
    FUNCTION???
-   **antigen\_pattern -** A regular expression used to \_\_\_\_\_.
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
my_path <- file.path("~", "Desktop", "template_folder", "panel_template.csv")

# generate the template file, which then can be edited manually 
cytofin_generate_panel_template(template_path = my_path)
```

The final two arguments for `cytofin_homogenize` are `input_data_path`
and `output_data_path`, two strings that indicate which directory input
.fcs files should be read from and which directory homogenized .fcs
files should be written to, respectively.

Using these arguments, `cytofin_homogenize` can homogenize a set of
CyTOF files with distinct antigen naming conventions. Specifically, the
program performs a regular expression search to match the synonymous
term in the panel and correct the antigen name with standardized names
in the panel.

Example function call:

``` r
# define input and output paths 
# TO DO: We need to find a better way of getting these data to the user 
#        referencing this vignette. Or we need to make the data smaller. 
metadata_path <- 
  here::here("inst", "extdata", "test_metadata_raw.csv")
panel_path <- 
  here::here("inst", "extdata", "test_panel.csv")
input_data_path <- 
  here::here("inst", "extdata", "test_raw_fcs_files")
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
`antigen_pattern` will be standardized to the names given in the `range`
column of the reference panel.

### CyTOF data normalization using external anchors

After dataset homogenization, batch correction (or “normalization”) can
be performed across datasets. The following schematic diagram
illustrates how normalization is performed:

![Alt text](./images/Slide2.PNG?raw=true "Title")

TO DO: Fix a few things about this schematic diagram

The external anchor normalization proceeds with two steps:

1.  Preparation of external anchors and
2.  Application of transformation function.

We detail each of these steps below.

#### 1. Anchor preparation:

The `cytofin_prep_anchors` function concatenates the identified anchor
files (one file per plate/batch) and subsequently generates summary
statistics used for batch correction by later steps of the pipeline.
These summary statistics include the univeral (i.e. overall) means and
variances of all channels in the homogenized dataset as well as the
means of both of these values (to be used in bulk, non-channel-specific
batch correction). It can both return these statistics as a `list()` and
save them as an .rds file in a specified directory.

`cytofin_prep_anchors` takes 4 arguments: \* `metadata_path`: A
connection leading to an .xlsx or .csv file containing a metadata table
about each file to be analyzed. This file should be identical to that
used for `cytofin_homogenize`. \* `panel_path`: A connection leading to
an .xlsx or .csv file containing information about the standardized
antigen panel in the homogenized dataset. This file should be identical
to that used for `cytofin_homogenize`. \* `input_data_path`: A
connection to a directory containing the input .FCS files from which to
draw summary statistics \* `output_path`: A connection to a directory
where the output .rds and .FCS files will be saved. The default is
“none”, in which case no output files will be stored (and the only
effect of the function will be to return the calculated statistics as a
`list()`).

An example:

``` r
anchor_statistics <- 
  cytofin_prep_anchors(
    metadata_path = metadata_path, 
    panel_path = panel_path, 
    input_data_path = file.path("~", "temp", "out_test"), 
    output_path = file.path("~", "temp", "out_test", "anchor_prep")
  )

print(anchor_statistics)
#> $var_uni
#>         Time Event_length    (Pd102)Di    (Pd104)Di    (Pd105)Di    (Pd106)Di 
#>    1.2162778    0.1558900    4.3739226    5.2718740    2.6143893    3.6480251 
#>    (Pd108)Di    (Pd110)Di    (In113)Di    (In115)Di    (La139)Di    (Pr141)Di 
#>    4.4172818    4.0906834    1.5580572    3.0209456    0.3526106    0.4121863 
#>    (Nd142)Di    (Nd143)Di    (Nd144)Di    (Nd145)Di    (Nd146)Di    (Sm147)Di 
#>    1.4156903    0.5137832    1.0072002    0.2515747    1.8266198    1.2210086 
#>    (Nd148)Di    (Sm149)Di    (Nd150)Di    (Sm152)Di    (Eu153)Di    (Sm154)Di 
#>    2.4468164    2.7853649    0.4828961    1.5855550    2.6799599    0.8680412 
#>    (Gd156)Di    (Gd158)Di    (Gd160)Di    (Dy161)Di    (Dy162)Di    (Dy163)Di 
#>    3.0887633    1.8616815    3.0766629    0.9811344    0.1051780    0.0623276 
#>    (Dy164)Di    (Ho165)Di    (Er166)Di    (Er167)Di    (Er168)Di    (Er170)Di 
#>    1.0875415    4.1854272    0.7091670    5.0975416    4.2993913    2.0480149 
#>    (Yb171)Di    (Yb172)Di    (Yb173)Di    (Yb174)Di    (Lu175)Di    (Yb176)Di 
#>    0.9510328    3.6103713    0.3994399    4.4140730    1.0460688    1.8253507 
#>    (Ir191)Di    (Ir193)Di 
#>    2.8895530    2.7402837 
#> 
#> $mean_uni
#>         Time Event_length    (Pd102)Di    (Pd104)Di    (Pd105)Di    (Pd106)Di 
#>   14.7167663    2.2581173    4.2323231    4.2997551    5.4608738    3.1754406 
#>    (Pd108)Di    (Pd110)Di    (In113)Di    (In115)Di    (La139)Di    (Pr141)Di 
#>    2.4091200    1.8342731    0.9310619    2.3723671    0.4888809    0.3658603 
#>    (Nd142)Di    (Nd143)Di    (Nd144)Di    (Nd145)Di    (Nd146)Di    (Sm147)Di 
#>    1.4420595    0.5468615    0.8490954    0.4168375    1.5083718    0.5970459 
#>    (Nd148)Di    (Sm149)Di    (Nd150)Di    (Sm152)Di    (Eu153)Di    (Sm154)Di 
#>    1.7862250    1.5815316    0.4839106    0.9445937    2.1161911    0.5480379 
#>    (Gd156)Di    (Gd158)Di    (Gd160)Di    (Dy161)Di    (Dy162)Di    (Dy163)Di 
#>    1.5095232    1.2012469    2.6788510    0.4609075    0.1496267    0.1014609 
#>    (Dy164)Di    (Ho165)Di    (Er166)Di    (Er167)Di    (Er168)Di    (Er170)Di 
#>    0.8211192    2.2591851    0.7192399    3.8276557    3.7954262    0.5737639 
#>    (Yb171)Di    (Yb172)Di    (Yb173)Di    (Yb174)Di    (Lu175)Di    (Yb176)Di 
#>    0.6453832    2.0596321    0.3955951    3.8735154    0.8437249    1.5804540 
#>    (Ir191)Di    (Ir193)Di 
#>    4.5485968    5.2404215 
#> 
#> $var_uni_mean
#> [1] 1.858537
#> 
#> $mean_uni_mean
#> [1] 1.507341
```

\[TO DO: Some description of these values\].

#### 2. Batch correction

Description: the annorm function applied different transformation
functions (modes) to normalize each anchor to the reference statistics
generated by the anprep function.

Function definition:

`annorm (control_metadata_filename, control_data_filename, sample_metadata_filename, panel_filename, input_file_dir, val_file_dir="none" ,output_file_dir, mode)`

Input:

`control_metadata_file`: metadata table of anchor CyTOF files
(.fcs)(must be in the current directory).

`control_data_filename`: RData object containing anchor referene
statistics (must be in the current directory).

`sample_metadata_filename`: metadata table of homogenized CyTOF files
(.fcs)(must be in the current directory).

`panel_filename`: standardized antigen panel table file
(.xlsx/.csv)(must be in the current directory).

`input_file_dir`: folder directory containing input homogenized CyTOF
data file.

`val_file_dir`: folder directory containing validation homogenized CyTOF
data file (optional).

`output_file_dir`: folder directory containing output normalized CyTOF
data file.

`mode`: transformation function used for normaliztion (*meanshift*,
*meanshift\_bulk*, *variance*, *z\_score*, *beadlike*).

Output: normalized CyTOF files.

**CyTOF data normalization using internal anchors**

![Alt text](./images/Slide3.PNG?raw=true "Title")

Description: In the event that the external references are not
available, internal anchors can be used. Here, we identified the most
stable channels as internal anchors using a PCA-based non-redundnacy
score (NRS). A minimal of three channels should be selected to establish
an internal refernece from which signal can be calibrated between CyTOF
files.

Function definition:

`annorm_nrs(sample_metadata_filename, panel_filename, input_file_dir, val_file_dir="none", output_file_dir, nchannels)`

Input:

`sample_meta_filename`: metadata table of homogenized CyTOF files
(.fcs)(must be in the current directory).

`panel_filename`: standardized antigen panel table file
(.xlsx/.csv)(must be in the current directory).

`val_file_dir`: folder directory containing validation homogenized CyTOF
data file (optional).

`output_file_dir`: folder directory containing output normalized CyTOF
data file.

`nchannels`: number of stabilized channels used for normalization.

Output: normalized CyTOF files.

**Computational pipeline for CyTOF data integration**

Below is a demo Rscript using CytofIn package for CyTOF data
integration.

    #import cytofin R package
    library(cytofin)

    #homogenization antigen panel, use the demo data supplied with the package
    metadata_filename <- paste0(path.package("cytofin"),"/extdata/test_metadata_raw.csv")
    panel_filename <- paste0(path.package("cytofin"),"/extdata/test_panel.csv")
    input_file_dir <- paste0(path.package("cytofin"),"/extdata/test_raw_fcs_files/")
    output_file_dir <- "out_test/"
    homogenize(metadata_filename, panel_filename, input_file_dir, output_file_dir)

    #prep external anchor 
    anchor_metadata_filename <- paste0(path.package("cytofin"),"/extdata/test_anchor_metadata_raw.csv")
    input_file_dir <- output_file_dir #use the homogenized files
    anprep(anchor_metadata_filename, panel_filename, input_file_dir)

    #data normalization using external anchors and meanshift transofmration function
    val_file_dir <- paste0(path.package("cytofin"),"/extdata/test_batch_fcs_files/")
    anchor_data_filename <- "./Prep_control.RData"
    output_file_dir <- "norm_test/"
    mode <- "meanshift"
    annorm(anchor_metadata_filename, anchor_data_filename, metadata_filename, panel_filename, 
    input_file_dir, val_file_dir, output_file_dir, mode)
    annorm(anchor_metadata_filename, anchor_data_filename, metadata_filename, panel_filename, 
    input_file_dir, "none", output_file_dir, mode)

    #data normalization using 4 internal channels and meanshift_bulk transformation function
    nchannels <- 4
    output_file_dir <- "norm_test2/"
    annorm_nrs(metadata_filename, panel_filename, input_file_dir, val_file_dir, 
    output_file_dir, nchannels)
    annorm_nrs(metadata_filename, panel_filename, input_file_dir, "none", 
    output_file_dir, nchannels)
