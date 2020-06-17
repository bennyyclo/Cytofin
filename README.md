# Cytofin

Cytofin (Cytof integration) is an R package for homogenizing and integrating hetereogenous Cytof data from diverse data source.

Before CyTOF data inegration, all cytof files need to be homogenized to have consistent channels. Cytofin requires that all input cytof be homogenized based on a user provided standarized panel with user defined search pattern. To normalize the cytof data, Cytofin use a novel generalized anhor strategy that define the based line of the signal between batch to correct for batch effects. One anchor need to be identified by the user from each plate (batch). A reference anchor is generated based on the mean expression of all identified anchors from each plate (batch). Next, different transformation functions is applied to fit each plate-specific anchor to the reference data distribution and the same transformation is then appleid to correct the sample data signal on each plate.  

Cytofin provided three functions for Cytof data integration:

1. homogenize-this function performs batch homogenization of cytof data based on a user-defined panel and search pattern. 
2. anprep-this generate reference anchors from anchors identified from each plate (batch).
3. annorm-this function performed signal normalization using transformation function based on the anchor from anprep.

1. CyTOF data homogenization

![Alt text](./images/Slide1.png?raw=true "Title")

Requred Input: raw CyTOF files (.fcs), standardized antigen panel table file (.xlsx/.csv), metadata table (.xlsx/.csv) 
Output: The homognize function take a user input antigen panel table, which include standardized antigen name and associated antigen search pattern. Given two Cytof files with distinct antigen naming, the program perform a regular expression search to match the synomymous term in the panel and correct the antigen name with standardized name in the panel. The output of this function is homogenized CyTOF file with user defined channel defined by the standardized antigen table.  

2. CyTOF data normalization

![Alt text](./images/Slide2.png?raw=true "Title")

The normalization include two steps: 1. preparation of external anchors and 2. application of transformation function.

1. Anchors preparation:

Requred Input: homogenized CyTOF files (.fcs), and anchor table file (.xlsx/.csv).
Output: the anprep function concatenated the identified anchor file, one file per plate/batch, and subsequently generated summary statistics including mean and variacne which will be used for batch correction. The values will be store in the RData object. The function also outputs the concatenated FCS files.









