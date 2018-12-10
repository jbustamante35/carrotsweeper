# CarrotSweeper
MATLAB-based functions to decompose the shape of objects in an images. This is
a pipeline that extracts contours from a set of images, run the contours
through Principal Components Analysis, then iteratively sweep through those
PCs. More features to come. As of now, this is optimized for images of carrot
roots on a black background, where filenames contain data on that root's
Location, Curvature, Unique Identifier (UID), etc.

## Background
This is simply a subset of files from HypoQuantyl needed to run a full analysis
for Scott Brainard's corn project. Scott images a large collection of carrot
roots of different genotypes and is trying to figure out what genes define
particular traits/phenotypes. He is particularly interested in being able to
predict the offspring phenotype when crossing different genotypes.

He came to Nathan first to talk about the best way to analyze his large
datasets, but at the time he was juggling 20+ projects and so handed him off
to me. The ways he wants to analyze his data is going to give me pretty good
practice for what I'm doing with my hypocotyls, so I should get a lot out of
this.

## Summary
Generally what we want to do is decompose the 2D shape of carrot roots in order
to generate a model in which Scott can define what controls root shape/size/etc.

Our main approach has been to perform a Principal Components Analysis on the
contours of these roots to reduce the dimensionality of the contour shape to
just a small set of features.

As of now we've just been doing PCA on the x-coordinates and y-coordinates of
each carrot's contour. Common sense told us that the x-coordinate components
would control the length (distance from shoulder to tip), whereas the major
y-coordinate components would control width (mean distance from each shoulder
to both sides of the tip). Additional components would likely reflect more
subtle root features.

### CarrotSweeper Pipeline
The general pipeline for our analyses (as of September 9, 2018) is as follows:

Scott:
- Collect color and black-and-white (bw) images of carrot roots
  laid on black backgrounds <br />
- Cut off the head and leave a small gap from the edge of the 
  window <br />
    * This was necessary because the bw image combines the 
      edge of the photo with the contour. <br />
    * This considerably weakened the accuracy of the data <br />
- Generate filenames that contain plenty of identifiable 
  information <br />
    * Location, Genotype, UID, Curvature, etc. <br />

Julian:
- Read the bw images from the directory/subdirectories <br />
- Generate contours and normalize their coordinates (more detail belowbr />
- Rasterize x-/y-coordinates and run them through Principal Components Analysis <br />
- Iteratively sweep up/down by 1 standard deviation through the resulting PCs <br />
- Concatenate the PC scores with each carrot's corresponding UID in a csv table <br />

[to be continued]

## Getting Started
### Dependencies
- MATLAB R2017B (recommended, but most likely works on most older versions)
- Image Processing Toolbox
- Statistics and Learning Toolbox

This repository uses many functions from my other GitHub repo,
[HypoQuantyl](https://github.com/jbustamante35/hypoquantyl) (HQ). The full list
of functions taken from HQ are found in the [dependencies](./dependencies) text
file in the main directory of this repo.

### Installation
Clone this repository to your desired location, then simply add
carrotsweeper/scripts to your MATLAB search path.

```
git clone https://github.com/jbustamante35/carrotsweeper
```

## Author
**Julian Bustamante**, Cellular and Molecular Biology (<jbustamante@wisc.edu>) <br />
	University of Wisconsin - Madison <br />
	Department of Botany <br />


