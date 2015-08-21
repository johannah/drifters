
# This post will take a closer look at magnetic data collected by ocean drifters. The drifters were deployed off of the northwestern United States in the summer of 2012 by the Ocean Bottom Magenetology Laboratory. Five drifters collected data and reported back via an iridium link. The data is described in the following cell.
#
# More information about this project and data can be found at: http://deeptow.whoi.edu/research/magdrifter.html

# The APS magnetic sensor made a measurement for 128 seconds every 20 minutes. Those 128 samples were averaged into one value of X,Y and Z magnetic field value. Every 4 hours the data were sent via Iridium to the Argo drifter database website
#
# Once a day (every 24 hrs) a calibration cycle was run where 256 sequential raw data values were recorded. The data fields are ascii with an identifying first character to describe the type of data.
#
# ----------------------------------------------------
#
# S: status line
#
# S 000000000066760 123 270 2011.08.13 23:03:06 47.651360 -129.042221 8
#
# S drifter_id rec_# day date time latitude longitude ?
#
# ----------------------------------------------------
#
# M: measurement
#
# M 1 -62475.9 -32540.4 -10721.9 19.39 47.9019 -128.9508 1.6 2011 224.80556 17.49
#
# M meas_# xf yf zf temperature latitude longitude ? yr decimal_julian _day
#
# ----------------------------------------------------
#
# E:
#
# E 12.7 0
#
# E battery_voltage ?
#
# ----------------------------------------------------
#
# C: Calibration header
#
# C 8 2011 225.12518
#
# C id yr decimal_julian_day
#
# c: calibration measurement

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import sys
import math
from datetime import datetime
from glob import glob
from datetime import timedelta
plt.style.use('ggplot')
from mpl_toolkits.basemap import Basemap

#Convert julian day described in the data to datetime format
def convert_julian_frac(julian_frac, year):
    """
    julian_frac is string in the form of a float
    """
    frac, julian_day = math.modf(float(julian_frac)+1)
    #The drifters reported both 0 and 356 for julian days in the calendar year
    #When I get access to source code, I will try to determine which is correct
    if int(julian_day) > 365:
        julian_day = julian_day-365
        year = int(year) + 1
    mins, hrs = math.modf(frac*24.)
    secs, mins = math.modf(mins*60)
    usecs, secs = math.modf(secs*60)
    dval= '%s %s %s %s %s' %(year, int(julian_day), int(hrs), int(mins), int(secs))
    dtval = datetime.strptime(dval, '%Y %j %H %M %S')
    return dtval

def load_data(fname):
    """Input the name of the drifter file downloaded from the website. This function parses the two types of data,
    averaged measurements, M, and calibration measurements, C
    """
    dval = open(fname, 'r')
    mheaders = ["S_drifter_id", "S_record_num", "S_day", "S_date", "S_time", "S_lat", "S_lon", "S_?",
               "M_meas_num", "M_xf", "M_yf", "M_zf", "M_temp", "M_lat", "M_lon", "M_?",
               "M_datetime", "battery_voltage"]
    cheaders = [
                "C_id", "C_start_datetime","C_sample_datetime",
                "c_num", "c_1", "c_2", "c_3", "c_4",
               ]
    mdata = []
    cdata = []
    #initialize batter voltage
    bvalue = -1

    for line in dval:
        line_vals = line.split(' ')
        line_vals = [x for x in line_vals if x!='']
        line_vals[-1] = line_vals[-1].strip()
        if line_vals[0] == 'S':
            #S drifter_id rec_num day date time latitude longitude ?
            mstatus = line_vals[1:]
        if line_vals[0] == 'M':
            #M meas_num xf yf zf temperature latitude longitude ? yr decimal_julian _day
            mdt = convert_julian_frac(line_vals[10], line_vals[9])
            mval = mstatus + line_vals[1:9] + [mdt, bvalue]
            mdata.append(mval)
        if line_vals[0] == 'C':
            #C id yr decimal_julian_day
            jdt = convert_julian_frac(line_vals[3], line_vals[2])
            Cid = line_vals[1]
            Cdf = jdt
        if line_vals[0] == 'c':
            #calibration measurement
            cdata.append([Cid, Cdf, Cdf+timedelta(0,int(line_vals[1]))] + line_vals[1:])
        if line_vals[0] == 'E':
            #battery_voltage ?
            bvalue = line_vals[1]

    mddict = {}
    mdatanp = np.asarray(mdata)
    for nn, name in enumerate(mheaders):
        mddict[name] = mdatanp[:,nn]
    mddpd = pd.DataFrame(mddict)
    mddpd['M_datetime'] = mddpd['M_datetime']
    mddpd['S_date'] = pd.to_datetime(mddpd['S_date'])
    mddpd.index = mddpd['M_datetime']
    mddpd = mddpd.convert_objects(convert_numeric=True)
    #remove rows with invalid battery voltage
    mddpd = mddpd[mddpd['battery_voltage'] > 3]
    #remove bad latitude data
    mddpd = mddpd[mddpd['M_lat'] != -90.0]
    mddpd = mddpd[mddpd['M_lat'] != 0.0]
    mddpd = mddpd[mddpd.index > pd.to_datetime('2012-06-01 00:00:00')]
    cddict = {}
    cdatanp = np.asarray(cdata)
    for nn, name in enumerate(cheaders):
        cddict[name] = cdatanp[:,nn]
    cddpd = pd.DataFrame(cddict)
    cddpd['C_sample_datetime'] = pd.to_datetime(cddpd['C_sample_datetime'])
    cddpd.index = cddpd['C_sample_datetime']
    cddpd = cddpd.convert_objects(convert_numeric=True)
    return mddpd, cddpd




def plot_drifter_dict(drifter_dict):
    """
    drifter dict is a dictionary containing eachthe drifter name with a subdictionary
    that contains a measurement pd and a calibration pd.

    This function plots measpds
    The function returns the basemap map class and the bounds of the dataset
    as a list [ minlat, maxlat, minlon, maxlon ]
    """
    #plot the tracks of the drifters
    #first we need to determine the extents of the lat/lon
    minlon, minlat = 999, 999
    maxlon, maxlat = -999, -999
    for drifter in drifter_dict.keys():
        minlond = drifter_dict[drifter]['meas']['M_lon'].min()
        if minlond < minlon:
            minlon = minlond
        maxlond = drifter_dict[drifter]['meas']['M_lon'].max()
        if maxlond > maxlon:
            maxlon = maxlond
        minlatd = drifter_dict[drifter]['meas']['M_lat'].min()
        if minlatd < minlat:
            minlat = minlatd
        maxlatd = drifter_dict[drifter]['meas']['M_lat'].max()
        if maxlatd > maxlat:
            maxlat = maxlatd
    print("Extents:", minlat, maxlat, minlon, maxlon)
    #build a map
    m = Basemap(projection='merc', ellps='WGS84',
               llcrnrlat=minlat,
               urcrnrlat=maxlat,
               llcrnrlon=minlon,
               urcrnrlon=maxlon,
               resolution='h')



    #plot map
    plt.figure()
    colors = ['r', 'b', 'g', 'k', 'c']
    for nn, drifter in enumerate(drifter_dict.keys()):
        dpd = drifter_dict[drifter]['meas']
        x, y = m(np.asarray(dpd['M_lon']),
                 np.asarray(dpd['M_lat']))
        m.scatter(x, y, c=colors[nn], edgecolor="None", label=drifter, s=3)
    m.drawcoastlines()
    plt.legend()
    return m, [minlat, maxlat, minlon, maxlon]

if __name__ == '__main__':
    #My raw files are located at data and are names as: drifter_<name>.txt
    drifter_data_dir = 'data'
    drifter_files = glob(os.path.join(drifter_data_dir, 'drifter*.txt'))
    drifter_dict = {}
    for df in drifter_files:
        dname = os.path.split(df)[1].replace('.txt', '')
        print("Loading %s and writing measured and calibration data files" %dname)
        measpd, calpd = load_data(df)
        mdname = df.replace('drifter_', 'measure_')
        measpd.to_csv(mdname, header=True, sep=' ', index=True)
        cdname = df.replace('drifter_', 'cal_')
        calpd.to_csv(cdname, header=True, sep=' ', index=True)
        drifter_dict[dname] = {'meas':measpd, 'cal':calpd}
    plot_drifter_dict(drifter_dict)
    plt.show()
