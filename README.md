# CarrotSweeper
Digital image-based phenotyping platform to analyze carrot root shape
attributes. <br />

**[TODO]** <br />
Redo all these blank sections

## Background

## Summary

### CarrotSweeper Pipeline
There are two main sub-directories to perform the major image analysis tasks:

The [carrot-straightener](./carrot-straightener) directory contains the
functions required to generate a midline from the image masks, and then run a
straightening algorithm of that mask.

The [pca-tools](./pca-tools) directory contains the analysis functions to
measure and quantify the straightened masks.


---
## Getting Started
### Dependencies
- MATLAB R2020A (recommended, but likely works on most older versions)
- Image Processing Toolbox
- Statistics and Learning Toolbox

This repository uses many functions from my other GitHub repo,
[HypoQuantyl](https://github.com/jbustamante35/hypoquantyl) (HQ). The full list
of functions taken from HQ are found in the [dependencies](./dependencies) text
file in the main directory of this repo.

### Installation
Clone this repository to your file system, then simply add carrotsweeper to
your MATLAB search path.

```
git clone https://github.com/jbustamante35/carrotsweeper
```

---
## Author
**Julian Bustamante**, Cellular and Molecular Biology (<jbustamante@wisc.edu>) <br />
	University of Wisconsin - Madison <br />
	Department of Botany <br />


