## tree_connect.csv
Calculated neighbourhood density indices for each trap x year x species observation.
Output of script `07_calculate-neighbourhood-densities.R`

| variable | description |
|:---------|:------------|
| `year` | The phenological year. The start of the phenological year was estimated per species as described in the manuscript. If the date of the weekly seed trap census was after the start of the phenological year then `year` = calendar year. If the date of census was before the start of the phenological year then `year` = calendar year - 1.|
| `trap` | Seed trap ID. Labelled from "trap_001" to "trap_450", n = 450 |
| `sp4` |  4-letter ID code of plant species. This code does not change and can be used to track earlier botanical names to compare across publications. |
| `conn_RH` | The neighbourhood density of reproductive-sized heterospecifics |
| `conn_RC` | The neighbourhood density of reproductive-sized conspecifics | 
| `conn_NRC` | The neighbourhood density of non-reproductive-sized conspecifics |
| `alpha` | The chosen alpha value for the calculation of neighbourhood density (see manuscript) |

## species_abundance.csv 
Relative abundance of species in the 50-ha plot, 
summarised from 50-ha plot censuses conducted in 1990, 1995, 2000, 2005, 2010, 2010, 2015 and 2022.
Output of script `08_calculate-abundance.R`

| variable | description |
|:---------|:------------|
| `sp4` | 4-letter ID code of plant species. This code does not change and can be used to track earlier botanical names to compare across publications. |
| `genus` |Genus name. |
| `species` | Species name. Species names were matched to a static copy of The World Flora Online (http://www.worldfloraonline.org) v.2024.12 valid at 2025-04-23 |
| `median_abundance` | The median number of reproductive-sized individuals in the 50-ha plot across years | 
| `mean_abundance` | The mean number of reproductive-sized individuals in the 50-ha plot across years |
| `max_abundance` | The max number of reproductive-sized individuals in the 50-ha plot across years |

## trap_locations.csv
Co-ordinates for seed traps in the 50-ha plot.

| variable | description |
|:---------|:------------|
| `trap` | Seed trap ID. Labelled from "trap_001" to "trap_450", n = 450 |
| `quadrat` | Quadrat ID; created by overlaying a 50 ✕ 50 m grid over the 50-ha plot. Referred to in the manuscript as "subplot". Labelled from "quadrat_001" to "quadrat_100", n = 100 |
| `x` | The x coordinate within the plot, meters from the west border of the plot, always in [0,1000] |
| `y` | The y coordinate, meters from the south border, always in [0,500] |

## BCI_Plot_50ha.shp
A shapefile showing the location of the 50-ha plot on Barro Colorado Island.
The geometry is a polygon and uses the EPSG:4326 Coordinate Reference System (CRS).
