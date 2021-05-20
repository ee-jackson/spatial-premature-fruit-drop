# spatial-premature-fruit-drop

This repository contains code for the project: Spatial patterns of premature fruit drop in a tropical forest plant community

## Data sets used in this project

The seed rain data sets: `BCI_TRAP200_20190215_spcorrected.txt`, `trapLocations.csv`, `pairedTraps.csv`, and seed trait data set: `20120227_seedsMassForTraits.csv`, were provided by Joe Wright (wrightj@si.edu) and are not shared publicly.

The Barro Colorado Island 50-ha plot census data, `bci.tree*`, are in the [Dryad Digital Repository](https://doi.org/10.15146/5xcp-0d46). 

## Contents:

### code/
The `code/` directory contains two subdirectories, `exploration/` which contains R Markdown exploratory analyses, and `scripts/` which contains all the code for cleaning, combining, and analysing the data. All paths in the scripts are relative to the root directory (where the .Rproj file lives).

Each .R script has a summary at the top of what it does. The scripts are numbered in the order in which they would typically be run. All R Markdown files in `exploration/` are knitted to `github_documents` to make the GitHub repo browsable.
