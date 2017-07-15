# Intensity

This script is to quantify cell-cell borders fluorescence intensities. 
## Input data

**INPUT FILES** should be localised to two subfolders:
* *tifs_original*: unmodified average projections of original files 
numbered sequencially starting with 1 (i.e. 1,2,3,4...)

* *borders*: the corresponding 8-bit images segmented with 
Packing Analyser. The subfolders (1,2,3,4...) should contain
tracked_bd.png and vertices.png files from Packing Analyser.

## Running the script
1. Start by running Intensity.m

1. It gives choice to analyse embryos or wings.

1. After selecting current sample type, the script asks you whether 
to include borders with intensities below background. I would
recommend use "No" as default as it helps to exclude artefacts, but in
case of very weak signal or high cytoplasmic signal, choose "Yes".

1. Then you can choose location of the files to be analysed


1. Then the script automatically runs either Intensity_wing or
Intensity_embryo depending on the choice.

## Output data
1. Folder **Cells** contains .csv files with information about cell parameters
in each individual image.

1. Folder **Distribution no BG** contains graphs of border intensities vs. angle
after background subtraction and fitted with a line

1. Folder **Distribution with BG** contains graphs of border intensities vs. angle
without background subtraction and fitted with a line

1. Folder **Intensity** contains information about all borders intensity in each image

1. Folder **Summary** contains pulled information about all images:

* "cell_number.csv" - number of cells that are completely within image in each image

* "Distribution_all_no_BG.tif" - graphs of pulled border intensities vs. angle
after background subtraction and fitted with a line

* "Distribution_all_with_BG.tif" - graphs of pulled border intensities vs. angle
without background subtraction and fitted with a line

<<<<<<< HEAD
* "Intensity_angles.scv" - pulled data about all cell-cell borders

* "Intensity_wing.scv" or "Intensity_embryo.scv" - averaged data on fluorescence
=======
* "Intensity_angles.csv" - pulled data about all cell-cell borders

* "Intensity_wing.csv" or "Intensity_embryo.csv" - averaged data on fluorescence
>>>>>>> origin/master
intensity from each image. In case of "Intensity_embryo.scv" data is subdivided 
into 0-10, 10-40 and 40-90 angles relative to image orientation. Image orientation
is calculated as mean orientation of all cells.

* "vertices_all" contains pulled information about all cells in all images.

## References
Bulgakova N.A., Brown N.H. (2016) Drosophila p120-catenin is crucial for endocytosis of the dynamic E-cadherin-Bazooka complex. Journal of Cell Science, 129(3), 477-482, [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/26698216). 
