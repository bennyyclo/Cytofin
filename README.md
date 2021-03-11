# CytofIN

CytofIN (CyTOF integration) is an R package for homogenizing and integrating heterogenous CyTOF data from diverse data sources.

Before CyTOF data integration, all CyTOF files need to be homogenized to have consistent channels. CytofIN requires that all input CyTOF files be homogenized based on a user-provided standardized panel with user defined search pattern. To normalize the CyTOF data, CytofIN uses a novel generalized anchor strategy that defines the based line of the signal between batch to correct for batch effects. One anchor needs to be identified by the user from each plate (batch). A reference anchor is generated based on the mean expression of all identified anchors from each plate (batch). Next, a user-specified transformation function is applied to fit each plate-specific anchor to the reference data distribution and the same transformation is then applied to correct the sample data signal on each plate.  

CytofIN provided three functions for CyTOF data integration:

1. homogenize-this function performs batch homogenization of cytof data based on a user-defined panel and search pattern. 
2. anprep-this function generates reference statistics from anchors identified from each plate (batch).
3. annorm-this function performs signal normalization using transformation function based on the anchor from anprep.
4. annorm_nrs-this function performs signal normalization using stabilized channels as internal anchors. 

**Installation**

1. Go to https://github.com/bennyyclo/Cytofin and click on clone or download button to download the package to the working directory. Alternatvely, from the command line,

```git clone https://github.com/bennyyclo/Cytofin```

2. Once unpacked, go inside the cytofin folder and initiate the R environment.

3. Install the CytofIN package using the following R commands:

```
library(devtools)

devtools::install()
```

**CyTOF data homogenization**

![Alt text](./images/Slide1.png?raw=true "Title")

Description:
the homognize function takes a user input antigen panel table, which includes standardized antigen name and associated antigen search pattern. Given two CyTOF files with distinct antigen naming, the program performs a regular expression search to match the synonymous term in the panel and correct the antigen name with standardized names in the panel.

Function Definition: 

```homogenize(metadata_filename, panel_filename, input_file_dir, output_file_dir)```

Input: 

```metadata_file```: metadata table of raw CyTOF files (.fcs)(must be in the current directory).

```panel_filename```: standardized antigen panel table file (.xlsx/.csv)(must be in the current directory).

```input_file_dir```: folder directory containing input raw CyTOF files.

```output_file_dir```: folder directory containing output homogenized files.

Output: homogenized CyTOF file with user-defined channels presented in the standardized antigen table.  



**CyTOF data normalization using external anchors**


![Alt text](./images/Slide2.PNG?raw=true "Title")

The external anchor normalization steps include: 1. preparation of external anchors and 2. application of transformation function.

1. Anchors preparation:

Description: 
the anprep function concatenates the identified anchor file, one file per plate/batch, and subsequently generates summary statistics including mean and variance which will be used for batch correction. 

Function definition: 

```anprep(metadata_filename, panel_filename, input_file_dir)```

Input: 

```metadata_filename```: metadata table of anchor CyTOF files (.fcs)(must be in the current directory).

```panel_filename```: standardized antigen panel table file (.xlsx/.csv)(must be in the current directory).

```input_file_dir```: folder directory containing output data.


Output: an RData object containig reference statistics and concatenated anchor FCS files.

2. Data transformation:

Description: the annorm function applied different transformation functions (modes) to normalize each anchor to the referenece statistcs generated by the anprep function.

Function definition:

```annorm (control_metadata_filename, control_data_filename, sample_metadata_filename, panel_filename, input_file_dir, val_file_dir="none" ,output_file_dir, mode)```

Input: 

```control_metadata_file```: metadata table of anchor CyTOF files (.fcs)(must be in the current directory).

```control_data_filename```: RData object containing anchor referene statistics (must be in the current directory).

```sample_metadata_filename```: metadata table of homogenized CyTOF files (.fcs)(must be in the current directory).

```panel_filename```: standardized antigen panel table file (.xlsx/.csv)(must be in the current directory).

```input_file_dir```: folder directory containing input homogenized CyTOF data file.

```val_file_dir```: folder directory containing validation homogenized CyTOF data file (optional).

```output_file_dir```: folder directory containing output normalized CyTOF data file.

```mode```: transformation function used for normaliztion (_meanshift_, _meanshift_bulk_, _variacne_, _z_score_, _beadlike_).

 
Output: normalized CyTOF files.

**CyTOF data normalization using internal anchors**


![Alt text](./images/Slide3.PNG?raw=true "Title")

Description:
In the event that the external references are not available, internal anchors can be used. Here, we identifed the most stable channels as internal anchors using a PCA-based non-redundnacy score. A minimal of three channels should be selected to establish an internal refernece from which signal can be calibrated between CyTOF files.

Function definition:

```annorm_nrs(sample_metadata_filename, panel_filename, input_file_dir, val_file_dir="none", output_file_dir, nchannels)```

Input: 

```sample_meta_filename```: metadata table of homogenized CyTOF files (.fcs)(must be in the current directory).

```panel_filename```: standardized antigen panel table file (.xlsx/.csv)(must be in the current directory).

```val_file_dir```: folder directory containing validation homogenized CyTOF data file (optional).

```output_file_dir```: folder directory containing output normalized CyTOF data file.

```nchannels```: number of stabilized channels used for normalization.


Output: normalized CyTOF files.

**Computational pipeline for CyTOF data integration**

Below is an demo Rscript using Cytofin package for CyTOF data integration.

```
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

```
