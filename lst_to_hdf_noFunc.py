# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 16:56:40 2023

@author: Mehdi Heris, Andrew Kittredge
"""
# main raster address
# This is the NLCD LC reference raster
mainRstPath=r"G:\\NLCD\\nlcd_2021_land_cover_l48_20230630.img"

# landsat folder address
lstDir =r"G:/Landsat"

# hdf files addresses
hdfLstPath=r"H:\\lst-national-mosaic-july2024.hdf5"
hdfDatePath=r"E:\\date-national-mosaic-july2024.hdf5"

# Landsat 8-9 clear sky value
# Update: Adds clear sky over water value
clearQaVal=[21824, 21952]   

# output directory
resultDir = r"H:\\lst-national-mosaic-july2024-result.csv" 

from datetime import date
import os
import numpy as np
import rasterio
from glob import glob
from rasterio.warp import calculate_default_transform, reproject, Resampling
import h5py
import json
import tarfile
import calendar
import random
import time
import pandas as pd

# set the timer
t1=time.time()

# this function finds the dataset (level) that is empty
## gets the hdf file and the window to get the data
def findLevelToWrite(f,startRow,endRow,startCol,endCol):
    level=-1
    for key in f.items():
        #print (key[0])
        d=f[key[0]]
        t=d[startRow:endRow,startCol:endCol]
        if np.nanmax(t)==0:
            level=key[0]
            break
    return (level)

# This function converts a yyyymmdd date to days since 2013-01-01
def dateToInteger (string):
    originDate = date(2013, 1, 1)
    y = int(string[0:4])
    m = int(string[4:6])
    d = int(string[6:])
    stringDate = date(y, m, d)
    intDate = abs(stringDate - originDate)
    return(intDate.days)

# this function gets a mother raster and a child raster
   # we want to write the child raster in the mother raster
   # to make sure that the bounaries fit and we find the overlap
   # we use this function to find the overlapping boundaries
def calculatableBounds2(mainRstBounds,fillingRstBounds):
    # match the top side
    if fillingRstBounds[3]>mainRstBounds[3]:
        endY = mainRstBounds[3]
    else:
        endY = fillingRstBounds[3]
    # match the bottom side  
    if mainRstBounds[1]>fillingRstBounds[1]:
        startY = mainRstBounds[1]
    else:
        startY = fillingRstBounds[1]
    # match the left side
    if mainRstBounds[0]>fillingRstBounds[0]:
        startX = mainRstBounds[0]
    else:
        startX = fillingRstBounds[0]
    # match right side
    if fillingRstBounds[2]>mainRstBounds[2]:
        endX = mainRstBounds[2]
    else:
        endX = fillingRstBounds[2] 
    return (startX,startY,endX,endY)

# this function recieves the calculated boundary and rasters
# then calculates the start and end rows and cols that will be used to fill arrays
def sliceFillingRst (calBounds,rstBounds,rstShape, cellSizeAr):
    # calBounds is (startX,startY,endX,endY)
    startRow = int((rstBounds[3] - calBounds[3])/cellSizeAr)
    endRow   = int(rstShape[0] - ((calBounds[1] - rstBounds[1])/cellSizeAr))
    startCol = int((calBounds[0] - rstBounds[0])/cellSizeAr)
    endCol   = int(rstShape[1] - ((rstBounds[2]-calBounds[2])/cellSizeAr))
    return (startRow,endRow,startCol,endCol)

def sliceFillingRst_2 (calBounds,rstBounds,rstShape, cellSizeAr):
    # calBounds is (startX,startY,endX,endY)
    startRow = int(np.floor((rstBounds[3] - calBounds[3])/cellSizeAr))
    endRow   = int(np.ceil(rstShape[0] - ((calBounds[1] - rstBounds[1])/cellSizeAr)))
    startCol = int(np.floor((calBounds[0] - rstBounds[0])/cellSizeAr))
    endCol   = int(np.ceil(rstShape[1] - ((rstBounds[2]-calBounds[2])/cellSizeAr)))
    return (startRow,endRow,startCol,endCol)

def sliceFillingRst_3 (calBounds,rstBounds,rstShape, cellSizeAr):
    startRow = int(np.floor((rstBounds[3] - calBounds[3])/cellSizeAr))
    endRow   = int(np.floor(rstShape[0] - ((calBounds[1] - rstBounds[1])/cellSizeAr)))
    startCol = int(np.floor((calBounds[0] - rstBounds[0])/cellSizeAr))
    endCol   = int(np.ceil(rstShape[1] - ((rstBounds[2]-calBounds[2])/cellSizeAr)))
    return (startRow,endRow,startCol,endCol)

def sliceFillingRst_4 (calBounds,rstBounds,rstShape, cellSizeAr):
    startRow = int(np.floor((rstBounds[3] - calBounds[3])/cellSizeAr))
    endRow   = int(np.ceil(rstShape[0] - ((calBounds[1] - rstBounds[1])/cellSizeAr)))
    startCol = int(np.floor((calBounds[0] - rstBounds[0])/cellSizeAr))
    endCol   = int(np.floor(rstShape[1] - ((rstBounds[2]-calBounds[2])/cellSizeAr)))
    return (startRow,endRow,startCol,endCol)

def lstWrite (path, mS, kS):
    with h5py.File(path,'r+') as fst:
        try:
            # let's find which level is empty to write
            level=(findLevelToWrite(fst,mS[0],mS[1],mS[2],mS[3]))
            # let's defind the slice from the kid that should be written
            kidToWrite=reprojected_array[kS[0]:kS[1],kS[2]:kS[3]]
            # dataset to write
            fst[level][mS[0]:mS[1],mS[2]:mS[3]]=kidToWrite
            lstSuccessFlag = True
        except(TypeError):
            print("Broadcast error: ", stAdr)
            #lstErrorList.append(stAdr)
            lstSuccessFlag = False
        return(level, lstSuccessFlag)


def dateWrite (path, mS, kS):
    with h5py.File(path,'r+') as fst_date:
        try:
            # let's find which level is empty to write
            level_date=(findLevelToWrite(fst_date,mS[0],mS[1],mS[2],mS[3]))
            # let's defind the slice from the kid that should be written
            kidToWrite_date=reprojected_array_datetime[kS[0]:kS[1],kS[2]:kS[3]]
            # dataset to write
            fst_date[level_date][mS[0]:mS[1],mS[2]:mS[3]]=kidToWrite_date
            dateSuccessFlag = True
        except(TypeError):
            print("Broadcast datetime error: ", stAdr)
            #dateErrorList.append(stAdr)
            dateSuccessFlag = False
        return(level_date, dateSuccessFlag)


np.seterr(divide='ignore', invalid='ignore')

# get the projection
with rasterio.open(mainRstPath) as dst:
    dst_crs = dst.crs

# get raster metadata, shape, bounds, cell size
with rasterio.open(mainRstPath) as dst:
    kwdsMainRaster = dst.meta.copy()
    mainRstShape=dst.shape
    mainRstBounds = dst.bounds
    mainRasterCellSize=kwdsMainRaster['transform'][0]

# this block creates 2 hdf5 files and creates 52 (13*4) datasets in them. 
#lst hdf

### Commented out to append existing files ###
# with h5py.File(hdfLstPath,'w') as fst:
#     for i in range(1,53):
#         dsetLevels = fst.create_dataset(name=f'level{i:02}',shape=mainRstShape,
#                                       dtype='uint16')
# # date hdf, uint32
# with h5py.File(hdfDatePath,'w') as date_fst:
#     for i in range(1,53):
#         dsetLevels = date_fst.create_dataset(name=f'level{i:02}',shape=mainRstShape,
#                                       dtype='uint16')

# EarthExplorer API returns C2L2 scenes as tar files
# Get list of tar files
tarPaths =[]
for t in glob(os.path.join(lstDir, '*.tar')):
    tarPaths.append(t)
    
tarFiles = [os.path.basename(tarPath) for tarPath in tarPaths]

# #Sample of 500 scenes
# random.seed(902)
# index500 = random.sample(range(len(tarFiles)), 500)   
# tarFiles500 = [tarFiles[i] for i in index500]

# Generate paths of B10 and QA pixel rasters in each tar file
lstSTs = []
lstQAs = []
lstDates = []
#for f in tarFiles:
for f in tarFiles:
    display_id = f.split(".")[0]
    st_path = "tar://" + lstDir + "\\" + f + "!" + display_id + "_ST_B10.TIF"
    qa_path = "tar://" + lstDir + "\\" + f + "!" + display_id + "_QA_PIXEL.TIF"
    lstDate = display_id.split("_")[3]
    lstSTs.append(st_path)
    lstQAs.append(qa_path)
    lstDates.append(lstDate)

# Mask B10 raster, reproject, find H5 level to write, write to H5
# for i in range(len(tarFiles)):
hdf_write_results = []
 
#for i in range(len(tarFiles)):
for i in range(4819 ,len(tarFiles)):
    current_file = tarFiles[i]
    stAdr = lstSTs[i]
    qaAdr = lstQAs[i]
    strDate = lstDates[i]
    print('started working on tarFiles[', i, ']')    
    with rasterio.open(stAdr,'r') as rstST, rasterio.open(qaAdr, 'r') as rstQA:
            arST=rstST.read(1)
            arQA=rstQA.read(1)
            
            intDate = dateToInteger(strDate)
            
            kwargs = rstST.meta.copy()
            #print(arST.shape)
            #print(rstST.meta.copy())
            
            ######## mask the raster by multiplying non-clear pixel values by 0
            arMask = np.where(np.logical_or(arQA == 21824, arQA == 21952),
                              1, 0)
            arST_masked = arST * arMask
            ######## reproject the raster
            dst_transform, dst_width, dst_height = calculate_default_transform(
                                                                        rstST.crs,
                                                                        dst_crs,
                                                                        rstST.width,
                                                                        rstST.height,
                                                                        resolution=(mainRasterCellSize,
                                                                                    mainRasterCellSize),
                                                                        *rstST.bounds)
            # we update the profile that we got from the input raster
            # with the calculated value. This profile will be used for the destination raster
            kwargs.update({
                'crs': dst_crs,
                'transform': dst_transform,
                'width': dst_width,
                'height': dst_height,
                'dtype':arST.dtype
            })
            #print(kwargs)
            # here we create a zeros raster with the shape that is calculated
            destination =( np.zeros((dst_height,dst_width), arST.dtype))
            # to know what reproject does try reproject?. this will return a corrected array
            projected_raster = reproject(
                # input masked array
                arST_masked,
                # destination array that is zeros now
                destination,
                src_transform=rstST.transform,
                src_crs=rstST.crs,
                dst_transform=dst_transform,
                dst_crs=dst_crs,
                dst_nodata=rstST.nodata,
                resampling=Resampling.nearest)
            ##### raster reprojected
            kidShape=projected_raster[0].shape
    
            kidCellSize=kwargs['transform'][0]
            left = kwargs['transform'][2]
            top =  kwargs['transform'][5]
            bottom = top - (kidCellSize*kidShape[0])
            right = left + (kidCellSize*kidShape[1])
            kidBounds=[left,bottom,right,top]
            reprojected_array = projected_raster[0]
            
            # Date array
            reprojected_array_datetime = np.zeros_like(reprojected_array, dtype=np.uint16)
            reprojected_array_datetime[:] = intDate
            # Create date array
            
            # this returns (startX,startY,endX,endY)
            calBoundsKidMother=calculatableBounds2(mainRstBounds,kidBounds)
            # lets find the start and end or rows and cols for each array: both the mother and the kid
            # this funciton returns (startRow,endRow,startCol,endCol)
            currentSliceFn = 1
            motherSlice=sliceFillingRst(calBoundsKidMother,mainRstBounds,mainRstShape,mainRasterCellSize)
            kidSlice= sliceFillingRst(calBoundsKidMother,kidBounds,kidShape,kidCellSize)
            motherSliceShape = (motherSlice[1] - motherSlice[0], motherSlice[3] - motherSlice[2])
            kidSliceShape = (kidSlice[1] - kidSlice[0], kidSlice[3] - kidSlice[2])
            
            # If mother slice and kid slice shapes, run 2nd sliceFillingRst function with rint
            # Clean solution would check for out of bounds errors
            if kidSliceShape != motherSliceShape:
                # Reset kidShape
                print("Trying with sliceFillingRst_2")
                currentSliceFn = 2
                kidShape=projected_raster[0].shape
                motherSlice=sliceFillingRst_2(calBoundsKidMother,mainRstBounds,mainRstShape,mainRasterCellSize)
                kidSlice= sliceFillingRst_2(calBoundsKidMother,kidBounds,kidShape,kidCellSize)
                motherShape_2 = (motherSlice[1] - motherSlice[0], motherSlice[3] - motherSlice[2])
                kidShape_2 = (kidSlice[1] - kidSlice[0], kidSlice[3] - kidSlice[2])
                
                if kidShape_2 != motherShape_2:
                    print("Trying with sliceFillingRst_3")
                    currentSliceFn = 3
                    kidShape=projected_raster[0].shape
                    motherSlice=sliceFillingRst_3(calBoundsKidMother,mainRstBounds,mainRstShape,mainRasterCellSize)
                    kidSlice= sliceFillingRst_3(calBoundsKidMother,kidBounds,kidShape,kidCellSize)
                    motherShape_3 = (motherSlice[1] - motherSlice[0], motherSlice[3] - motherSlice[2])
                    kidShape_3 = (kidSlice[1] - kidSlice[0], kidSlice[3] - kidSlice[2])
                
                    if kidShape_3 != motherShape_3:
                        print("Trying with sliceFillingRst_4")
                        currentSliceFn = 4
                        kidShape=projected_raster[0].shape
                        motherSlice=sliceFillingRst_4(calBoundsKidMother,mainRstBounds,mainRstShape,mainRasterCellSize)
                        kidSlice= sliceFillingRst_4(calBoundsKidMother,kidBounds,kidShape,kidCellSize)
                        motherShape_4 = (motherSlice[1] - motherSlice[0], motherSlice[3] - motherSlice[2])
                        kidShape_4 = (kidSlice[1] - kidSlice[0], kidSlice[3] - kidSlice[2])
            
                        
            l, lf = lstWrite(hdfLstPath, motherSlice, kidSlice)
           
            dl, df = dateWrite(hdfDatePath, motherSlice, kidSlice)
            
            sceneWriteDict = {'i' : i,
                              'scene' : current_file,
                              'slice_function' : currentSliceFn,
                              'level_lst' : l,
                              'level_date' : dl,
                              'lst_write_success' : lf,
                              'date_write_success' : df,
                              'motherSlice' : motherSlice,
                              'kidSlice' : kidSlice}
           
            hdf_write_results.append(sceneWriteDict)
            
            result_df = pd.DataFrame(hdf_write_results)
            
            result_df.to_csv(resultDir)
            
          
t2=time.time()
print ('time passed: ',t2-t1)
#df_sceneSuccess = pd.DataFrame(hdf_write_results)
#df_sceneSuccess.to_csv(resultDir)




