# Visium Codebase

The included files are for reproducing the analysis and figures generated in the manuscript, starting from the output of SpaceRanger. Acompaning data from the manuscript might be necessary, such as snRNAseq atlas to map. Downloading the final Visium object from KPMP would allow the user to skip the data preparation scripts and move directly to the figures. The data preparation scripts should be run in the following order:
1. start_spatial_objects.R
2. preliminary_plots.R
3. merge_and_remove_razor_edge.R
4. tal_neighborhoods.R

The remaining scripts could be run in any order.

