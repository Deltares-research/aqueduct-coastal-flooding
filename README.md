Scripts to simulate coastal flood hazard at the global scale based on near-shore still water levels using a simple GIS routing.  
These scripts were initially developed for [Aqueduct Floods](http://floods.wri.org/).
See publication list for more background.

Quick start guide
-----------------
1. Install python using the miniconda distribution from [conda](https://docs.conda.io/en/latest/miniconda.html). 
2. Install all required dependencies using `conda env create -f environment.yml` with the provide [environment.yml file](environment.yml). Note that the scripts requires python 2.7.
3. Install PCRaster 4.2, see [quick-start-guide](http://pcraster.geo.uu.nl/quick-start-guide/)
4. Start the inun conda environment with `conda activate inun`
5. Prepare the ini file with references to the a DEM and forcing datasets, see example [configuration](coastal_inun.ini).
6. Run the model from CLI with e.g. `python scripts\coastal_inun.py -i coastal_inun.ini -d <output.nc> -b <waterlevels.nc> -v <waterlevel_variable_name> -s <sealevel_map.tif>`


Publications
------------
- Haer, T., et. al.: Coastal and river flood risk analyses for guiding economically optimal flood adaptation policies: A country-scale study for Mexico, Philos. Trans. R. Soc. A Math. Phys. Eng. Sci., 376(2121), 20170329, doi:10.1098/rsta.2017.0329, 2018.
- Tiggeloven, T., et. al.: Global-scale benefit–cost analysis of coastal flood adaptation to different flood risk drivers using structural measures, Nat. Hazards Earth Syst. Sci., 20(4), 1025–1044, doi:10.5194/nhess-20-1025-2020, 2020.
- Ward, P. J., et. al.: Aqueduct Floods Methodology, Washington D.C. Available from: https://www.wri.org/publication/aqueduct-floods-methodology, 2020.


Authors
-------
- Dirk Eilande (@DirkEilander)
- Hessel C. Winsemius (@hcwinsemius)
- Gundula Winter


License
-------
Copyright (c) Deltares

Distributed under GPL V3 Licencse, see [LICENSE](LICENSE.txt)