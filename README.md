# Well-doublet-optimization
Supporting code for the publication: Spatial analysis of thermal groundwater use based on optimal sizing and placement of well doublets

Installation: 
Download all the files from the repository and place them into one folder.

The files are the following:
- 'optimization_thermal_GW_use.py': main Python script for optimizing the well doublets.
- 'support_functions.py': Python script with the implementation of supporting functions. 
- 'input_data' : folder that contains the required input data
- 'results': folder that containts an exemplary result (this folder can be also empty, but should exist in the main directory before running the script)
    
Usage:
1) Install the required Python packages: "numpy", "csv", "pandas" and "geopandas".
2) Install Python-MIP package (https://www.python-mip.com/) and Gurobi solver (https://www.gurobi.com/documentation/quickstart.html) for the best performance.
3) Open the script 'optimization_thermal_GW_use.py' and choose optimization scenario by setting the values for: dist_ratio (distance ratio $r_\Delta$) and min_V_l (minimum pump rate of a doublet $q_{\mathrm{min}}$).
   Default is: dist_ratio=2 and min_V_l=5 [L/s].
4) Run the script 'optimization_thermal_GW_use.py'. 
5) The result files will be saved in the folder 'results\scenario_folder'. The 'scenario_folder' has the format: 'scenario_rA_VB', where A and B are the values of $r_\Delta$ and $q_{\mathrm{min}}$, respectively.  
6) Two results will be saved: 
    - optimization_result.csv - csv file containing summary of optimization results
    - optimal_wells.shp - shape file containing optimal well locations. 
