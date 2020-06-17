# Cytofin

Cytofin (Cytof integration) is an R package for homogenizing and integrating hetereogenous Cytof data from diverse data source.

Before CyTOF data inegration, all cytof files need to be homogenized to have consistent channels. Cytofin requires that all input cytof be homogenized based on a user provided standarized panel with user defined search pattern. To normalize the cytof data, Cytofin use a novel generalized anhor strategy that define the based line of the signal between batch to correct for batch effects. One anchor need to be identified by the user from each plate (batch). A reference anchor is generated based on the mean expression of all identified anchors from each plate (batch). Next, different transformation functions is applied to fit each plate-specific anchor to the reference data distribution and the same transformation is then appleid to correct the sample data signal on each plate.  

Cytofin provided three functions for Cytof data integration:

1. homogenize-this function performs batch homogenization of cytof data based on a user-defined panel and search pattern. 
2. anprep-this generate reference anchors from anchors identified from each plate (batch).
3. annorm-this function performed signal normalization using transformation function based on the anchor from anprep.

