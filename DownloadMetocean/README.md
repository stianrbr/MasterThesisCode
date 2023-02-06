# DownloadMetocean

This project defines the code needed for downloading hindcast data from several dataset avaiable from MET Norway.
The project is written as a part of a Master Thesis in Marine Hydrodynamics.
## Use:

* inputs.py: Where you should include all the inputs to the code. Only file that should be edited if one only is to download data from a single coordinate.

* misc.py: Definition of helper functions, placed in separate file for readability.

* download_functions.py: Definition of the actual functions. Separate functions are written for each data set due to different folder structures and nomenclature in the different projects.

* main.py: Main file of project, which is the one that must be run to download the data. The option for parallel downloads are included, where the number of parallels are specified in inputs.py.

* combine_yearly_data.py: Separate scripts for gathering the yearly data downloaded by main.py. Takes all the files in the savefolder specified by inputs.py, and gathers them in a separate file