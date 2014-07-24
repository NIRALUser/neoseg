
Quick Info
=================

-*- Compiling the program:
Run cmake on ./mu/brainseg for the command line version
Run cmake on ./mu/brainseg_gui for the GUI version, needs FLTK 1.1

-*- Example build sequence for UNIX:
mkdir bin-brainseg
cd bin-brainseg
ccmake /some/path/mu/brainseg
make

-*- Directory contents:

aniso = gradient anisotropic filter for a set of 3D MR images
bias = bias correction using polynomial least squares fit
brainseg = Expectation Maximization segmentation for brain MRI
brainseg_gui = GUI wrapper for the segmentation tool
common = basic classes and functions
conn = connected components
mireg = registration using Mutual Information and affine transformation
robust = robust estimation / clustering

-*- Brief description of segmentation parameters:

Atlas Directory
Directory where you can find the brain atlas, it must contain the following
files: white.gipl, gray.gipl, csf.gipl, rest.gipl, and template.gipl

Output Directory
Directory where the program will write all the output files (transforms and
output images)

Suffix
Tag that is used for all output files, every file will have the following
form: xxx_suffix.yyy. The string should only contain alphanumeric characters.
Can leave this blank, intended use is for tagging separate sessions
(ex: "test", date/time).

Input Images
List of images to be used for segmentation (ex: T1, T2). The choice of the
first image in the list is important as all output image will be in the space
of the first image.

Anisotropic Diffusion
Filtering to remove some textures, set iterations to 0 if you want to skip this
step. The time step should be low enough so it remains stable (<= 0.1).

Prior Weights
Scale the priors for wm, gm, and csf globally. Useful adjustment when
segmentations consistently generate too much of a particular label.

Maximum Bias Degree
Maximum polynomial degree for bias correction. Degree 4 should be enough for
most cases. Set to zero to disable bias correction.

-*- GUI help:

Add
Adds an image to the list of input images (ex. T1, T2, PD).

Remove
Delete the last entry from the list of input images.

File menu > Load params
Read a .seg parameter file and updates current input values.

File menu > Save params
Write the parameters to a .seg file. Useful if you want to run batch jobs
with the command line version. Ex: use the GUI to write a bunch of .seg files
and write a script that runs the segmentations one at a time.

Run
Perform segmentation with current parameters, has limited error checking.
