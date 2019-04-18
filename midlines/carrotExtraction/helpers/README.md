# CarrotSweeper Helper Functions
This is a set of functions from the _PhytoMorph_ codebase that are called in
the main CarrotSweeper code. To filter out the unnecessary functions, I brought
in the entire 1.2 GB _PhytoMorph_ repository and went down the rabbit hole of
figuring out each dependency's dependencies. There was probably a better way to
do this, but this method was the one I started with and just mindlessly
continued until I isolated the minimum number of functions needed to get the
code working. <br />

## Mex-compiled C code
This directory also has a C code that was compiled with the `mex` command. As
of 04.18.2019, this includes `ba_interp2`, `MtimesX`, and `DGradient`. I
compiled it for Linux, Windows, and MacOS and stored them in separate
directories. <br />


