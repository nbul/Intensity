# Intensity

This script is to quantify cell-cell borders fluorescence intensities. 

1. Start by running Intensity.m

2. It gives choice to analyse embryos or wings.

3. After selecting current sample type, the script asks you whether 
to include borders with intensities below background. I would
recommend use "No" as default as it helps to exclude artefacts, but in
case of very weak signal or high cytoplasmic signal, choose "Yes".

3. Then you can choose location of the files to be analysed
INPUT FILES should be localised to two subfolders:
tifs_original: unmodified average projections of original files 
numbered sequencially starting with 1 (i.e. 1,2,3,4...)

tifs_8bit: the corresponding 8-bit images segmented with 
Packing Analyser. The subfolders (1,2,3,4...) should contain
tracked_bd.png and vertices.png files from Packing Analyser.

4. Then the script automatically runs either Intensity_wing or
Intensity_embryo depending on the choice.

OUTPUT FILES:
Folder "Cells" contains .csv files with information about cell parameters
in each individual image.

Folder "Distribution no BG" contains graphs of border intensities vs. angle
after background subtraction and fitted with a line

Folder "Distribution with BG" contains graphs of border intensities vs. angle
without background subtraction and fitted with a line

Folder "Intensity" contains information about all borders intensity in each image

Folder "Summary" contains pulled information about all images:

"cell_number.csv" - number of cells that are completely within image in each image

"Distribution_all_no_BG.tif" - graphs of pulled border intensities vs. angle
after background subtraction and fitted with a line

"Distribution_all_with_BG.tif" - graphs of pulled border intensities vs. angle
without background subtraction and fitted with a line

"Intensity_angles.scv" - pulled data about all cell-cell borders

"Intensity_wing.scv" or "Intensity_embryo.scv" - averaged data on fluorescence
intensity from each image. In case of "Intensity_embryo.scv" data is subdivided 
into 0-10, 10-40 and 40-90 angles relative to image orientation. Image orientation
is calculated as mean orientation of all cells.

"vertices_all" contains pulled information about all cells in all images.
