# -*- coding: utf-8 -*-
"""
Optimization of thermal groundwater use

@author: Smajil Halilovic
"""

import pandas as pd
import geopandas as gpd
import os
import csv
from mip import Model, xsum, maximize, BINARY
import time
from support_functions import *

# define paths for input files:
file_well_positions = 'input_data/potential_wells.gpkg'   # file containing potential well locations for all city blocks
file_gw_direction = 'input_data/flow_dir_blocks.gpkg'     # file containing groundwater flow direction for each city block

# input data as dataframes:
data_df = gpd.read_file(file_well_positions)
data_df_dir = gpd.read_file(file_gw_direction)

# list of city block IDs - optimization is done for each city block separately:
block_IDs = data_df['baublock'].unique()

## For breaking the loop after certain number of blocks are optimized (and not optimizing all the blocks):
# block counter:
block_count = 0
# starting block:
start_block_number = 0
# ending block:
stop_block_number = 11

## Define optimization scenario:
# distance ratio: r = external_distance/internal_distance:
dist_ratio = 2
# minimum pumping rate of a doublet if installed (L/s):
min_V_l = 5
# regulatory minimum (internal) distance between extraction and injection wells of one doublet (10 [m]):
dist_min = 10
# minimum (external) distance between neighboring doublets:
dist_ext_min = dist_ratio*dist_min

# name of the scenario (for saving the results):
scenario_name = 'scenario_r' + str(dist_ratio) + '_V' + str(min_V_l)
scenario_folder = 'results/' + scenario_name
# generate directory for the scenario results:
if not os.path.exists(scenario_folder):
    os.mkdir(scenario_folder)

# define paths for result files:
file_well_result = scenario_folder + '/optimal_wells.shp'       # shape file
csv_file_result = scenario_folder + '/optimization_result.csv'  # CSV file
results_csv_row_list = [] # list of results to be saved in CSV file

# get the start time
start_time = time.time()

# initialize result data frame - for saving the results later:
result_df = None

## Iterate over city blocks - i.e. optimize each block seperately:
for b_ID in block_IDs:
    # optimization for this block:
    print("Optimizing the block:")
    print(b_ID)
    # increase block counter:
    block_count += 1

    # For skipping the certain number of blocks (and not starting from the 1st block):
    if block_count < start_block_number:
        continue

    # filter data for this block:
    in_block = data_df['baublock'] == b_ID
    data_block = data_df[in_block]
    # copy of data frame (for this block) - for saving the results later:
    result_data_block = data_block.copy()

    print("Number of potential well locations (extraction or injection) in the block:")
    print(len(data_block))

    # Groundwater flow direction for this block:
    # filter data for this block:
    in_block_dir = data_df_dir['baublock'] == b_ID
    data_block_dir = data_df_dir[in_block_dir]
    dir_flow = float(data_block_dir.iloc[0].at["dir"])
    direction = dir_flow - 90  # deg
    # rotation angle [rad] - needed in "upstream-downstream" constraint
    rotation_rad = direction * np.pi / 180  # rad

    # Ranges for extraction and injection wells - needed for optimization variables and constraints:
    Ext = range(len(data_block))
    Inj = range(len(data_block))

    # Dictionary that relates lines and wells - needed for constraints
    line_ext_wells = {}
    for well_ID in Ext:
        line_ID = data_block.iloc[well_ID]['line']
        if line_ID in line_ext_wells:
            line_ext_wells[line_ID].append(well_ID)
        else:
            line_ext_wells[line_ID] = []
            line_ext_wells[line_ID].append(well_ID)

    # List of line IDs:
    Lines_list = list(line_ext_wells.keys())
    Lines = range(len(Lines_list))

    # Generate optimization model
    optim_model = Model(solver_name="GRB")

    # Add optimization variables
    di = [optim_model.add_var(var_type=BINARY) for i in Inj]        # decision/optimization vars. for inj. wells
    de = [optim_model.add_var(var_type=BINARY) for e in Ext]        # decision/optimization vars. for ext. wells
    V_l = [optim_model.add_var() for l in Lines]                    # optimization vars. for pump rates of lines
    d_V = [optim_model.add_var(var_type=BINARY) for l in Lines]     # decision/optimization vars. for lines

    ## Objective function
    # maximize energy extracted from groundwater, i.e. the amount of pumped water:
    optim_model.objective = maximize(xsum(V_l[l] for l in Lines))

    ## Constraints

    ## Preprocessing some data needed for constraints:
    # Calculate for each line:
    median_alpha_lines = {}     # median 'alpha'
    max_V_lines = {}            # maximal pump rate - based on drawdown and upconing (flooding)
    max_Dx_lines = {}           # line length - max. distance between two points on a line
    # iterate over lines:
    for line_ID in line_ext_wells:
        # filter data for this line:
        on_line = data_block['line'] == line_ID
        data_line = data_block[on_line]
        median_alpha_lines[line_ID] = data_line['pump_rate_breakthrough_1m_l_s'].median()
        max_V_drawdown_line = data_line['pump_rate_drawdown_l_s'].max()
        max_V_flood_line = data_line['pump_rate_flooding_l_s'].max()
        max_V_lines[line_ID] = min(max_V_drawdown_line, max_V_flood_line)
        # find line length:
        line_len = 0
        for i in line_ext_wells[line_ID]:       # injection wells
            for j in line_ext_wells[line_ID]:   # extraction wells
                dist_ij = data_block.iloc[i]['geometry'].distance(data_block.iloc[j]['geometry'])
                line_len = max(line_len, dist_ij)
        max_Dx_lines[line_ID] = line_len

    # Constraint: Number of ext/inj wells on each line is 0 or 1, depending on if the doublet (on that line) is installed:
    for key in line_ext_wells:
        optim_model += d_V[Lines_list.index(key)] == xsum(de[i] for i in line_ext_wells[key])
    for key in line_ext_wells:
        optim_model += d_V[Lines_list.index(key)] == xsum(di[i] for i in line_ext_wells[key])

    # Constraint: Doublet can pump only if it is installed:
    for line_ID in line_ext_wells:  # iterate over lines
        optim_model += V_l[Lines_list.index(line_ID)] <= max_V_lines[line_ID] * d_V[Lines_list.index(line_ID)]
        # Constraint: If installed, the doublet must pump at least with min_V_l (e.g. 1 L/s):
        optim_model += min_V_l * d_V[Lines_list.index(line_ID)] <= V_l[Lines_list.index(line_ID)]

    # Constraint about regulatory 10m minimum distance between installed extraction and injection wells of one doublet
    # AND
    # Injection well must be downstream compared to extraction well:
    for key in line_ext_wells:
        for i in line_ext_wells[key]:       # injection wells
            for j in line_ext_wells[key]:   # extraction wells
                dist_ij = data_block.iloc[i]['geometry'].distance(data_block.iloc[j]['geometry'])
                x_inj = data_block.iloc[i]['geometry'].x
                y_inj = data_block.iloc[i]['geometry'].y
                x_ext = data_block.iloc[j]['geometry'].x
                y_ext = data_block.iloc[j]['geometry'].y
                x_r_inj, y_r_inj = rotation(x_inj, y_inj, rotation_rad)
                x_r_ext, y_r_ext = rotation(x_ext, y_ext, rotation_rad)
                if dist_ij < dist_min or x_r_inj <= x_r_ext:
                    optim_model += de[j] + di[i] <= 1
                #if dist_ij < dist_min or data_block.iloc[i]['geometry'].y <= data_block.iloc[j]['geometry'].y:
                #    optim_model += de[j] + di[i] <= 1

    # Drawdown constraint: V_l <= V_drawdown
    for line in Lines:  # iterate over lines
        line_ID = Lines_list[line]
        for j in line_ext_wells[line_ID]:  # extraction wells on that line
            optim_model += V_l[line] <= (de[j] * data_block.iloc[j]['pump_rate_drawdown_l_s'] + max_V_lines[line_ID] * (1 - de[j]))

    # Upconing (Flooding) constraint: V_l <= V_flooding
    for line in Lines:  # iterate over lines
        line_ID = Lines_list[line]
        for i in line_ext_wells[line_ID]:  # injection wells on that line
            optim_model += V_l[line] <= (di[i] * data_block.iloc[i]['pump_rate_flooding_l_s'] + max_V_lines[line_ID] * (1 - di[i]))

    # Internal breakthrough constraint: V_l <= V_breakthrough = alpha * Dist(ext,inj)
    # between installed extraction and injection wells of one doublet:
    for key in line_ext_wells:
        for i in line_ext_wells[key]:       # injection wells
            for j in line_ext_wells[key]:   # extraction wells
                if i!=j:
                    dist_ij = data_block.iloc[i]['geometry'].distance(data_block.iloc[j]['geometry'])
                    alpha_ij = (data_block.iloc[i]['pump_rate_breakthrough_1m_l_s'] + data_block.iloc[j]['pump_rate_breakthrough_1m_l_s'])/2
                    optim_model += V_l[Lines_list.index(key)] <= (alpha_ij * dist_ij + max_V_lines[key] * (2 - de[j] - di[i]))

    # Constraint: Distance between neighboring doublets
    # iterate over lines:
    for line_l in Lines:
        l = Lines_list[line_l]
        if len(line_ext_wells[l])>1:
            for line_p in Lines[line_l+1:]: # iterate over neighboring lines
                p = Lines_list[line_p]
                if len(line_ext_wells[p]) > 1:
                    # Calculate distance between two lines (can be also done in preprocessing):
                    # use "first" two points on each line to define the line functions:
                    well_i = line_ext_wells[l][0]
                    well_j = line_ext_wells[l][1]
                    m_l, c_l = lin_equ([data_block.iloc[well_i]['geometry'].x, data_block.iloc[well_i]['geometry'].y],
                                       [data_block.iloc[well_j]['geometry'].x, data_block.iloc[well_j]['geometry'].y])
                    well_p_i = line_ext_wells[p][0]
                    well_p_j = line_ext_wells[p][1]
                    m_p, c_p = lin_equ([data_block.iloc[well_p_i]['geometry'].x, data_block.iloc[well_p_i]['geometry'].y],
                                       [data_block.iloc[well_p_j]['geometry'].x, data_block.iloc[well_p_j]['geometry'].y])
                    # distance between two parallel lines l and p:
                    dist_l_p = distance_lines(m_l, c_l, c_p)
                    # check if these lines are 'relatively close' and if it makes sense to check the constraint:
                    if (max_Dx_lines[l] + max_Dx_lines[p]) > (2/dist_ratio)*dist_l_p:
                        param_M4 = (max_V_lines[l] * median_alpha_lines[p]) + (max_V_lines[p] * median_alpha_lines[l])
                        optim_model += (V_l[line_l] * median_alpha_lines[p]) + (V_l[line_p] * median_alpha_lines[l]) \
                                       <= ((2 / dist_ratio) * dist_l_p * median_alpha_lines[l] * median_alpha_lines[p] + param_M4 * (2 - d_V[line_l] - d_V[line_p]))
                    # Additional Constraint about min distance between doublets based on regulatory internal distance and distance ratio:
                    # e.g. dist_ext_min = 20[m]
                    if dist_l_p < dist_ext_min:
                        optim_model += d_V[line_l] + d_V[line_p] <= 1

    # Solve optimization problem:
    optim_model.optimize()
    print('optimal solution cost {} found'.format(optim_model.objective_value))

    # All optimization variables (solution):
    opt_vars = [v.x for v in optim_model.vars]
    print("Total number of optimization variables:")
    print(len(opt_vars))
    #print(opt_vars)

    # Optimization variables di for injection wells (solution):
    opt_i = opt_vars[:len(data_block)]
    opt_di = [int(round(k)) for k in opt_i]
    print("optimal solution - di:")
    print(opt_di)
    print("nr. of installed injection wells:")
    print(sum(opt_di))

    # Add these variables as a new column to result dataframe:
    result_data_block['optimal_inj_wells'] = opt_di

    # Optimization variables di for injection wells (solution):
    opt_e = opt_vars[len(data_block):2*len(data_block)]
    opt_de = [int(round(k)) for k in opt_e]
    print("optimal solution - de:")
    print(opt_de)
    print("nr. of installed extraction wells:")
    print(sum(opt_de))

    # Add these variables as a new column to result dataframe:
    result_data_block['optimal_ext_wells'] = opt_de

    # Optimization variables Vl for lines/doublets (solution):
    opt_Vl = opt_vars[(2*len(data_block)):(2*len(data_block) + len(Lines_list))]
    print("optimal solution - Vl - pump rates of doublets [L/s]:")
    print(opt_Vl)
    print("optimal solution - total Vl:")
    opt_total_Vl = sum(opt_Vl)
    print(opt_total_Vl)

    # Optimization variables dV for lines/doublets (solution):
    opt_d_V = opt_vars[-len(Lines_list):]
    opt_dV = [int(round(k)) for k in opt_d_V]
    print("optimal solution - dV - decision for doublets (lines):")
    print(opt_dV)
    print("nr. of installed doublets (lines) in the current city block:")
    opt_n_HP = sum(opt_dV)
    print(opt_n_HP)

    # Summary list of results for the current city block (one row in CSV results file):
    results_csv_row = [b_ID, opt_n_HP, opt_total_Vl]
    results_csv_row.extend(opt_Vl)
    results_csv_row_list.append(results_csv_row)

    # Updating (appending) results dataframe with results from this city block :
    if result_df is None:
        result_df = result_data_block.copy()
    else:
        result_df = pd.concat([result_df, result_data_block])

    # For breaking the loop after certain number of city blocks are optimized:
    if block_count > stop_block_number:
        break

# get the end time
end_time = time.time()

# get and print the execution time
elapsed_time = end_time - start_time
print('Execution time:', elapsed_time, 'seconds')

# Saving results in CSV file:
with open(csv_file_result, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(results_csv_row_list)

# Saving the results to shape file:
result_df.to_file(file_well_result)
