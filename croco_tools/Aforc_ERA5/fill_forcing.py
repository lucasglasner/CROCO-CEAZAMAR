'''
 # @ Author: lucasg
 # @ Create Time: 2023-02-27 11:16:33
 # @ Modified by: lucasg
 # @ Modified time: 2023-02-27 11:16:39
 # @ Description: 
 # ERA5 dataset doesnt have daily sea surface temperatures nor sea surface salinities.
 # This script creates a fake ERA5 SST and SSS from a different dataset, all in the context 
 # of croco_tools. With this ERA5-named files the interpolation routine from croco tools should
 # work without problems.
  
 *** The SST and SSS in the forcing files is only for nudging ***
 lucas.glasner@ceaza.cl
 '''
 
import xarray as xr
import numpy as np


