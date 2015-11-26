import math
import requests
import os
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
from igrf12py.igrf12fun import runigrf12, plotigrf 


base_url = 'http://www.ndbc.noaa.gov/view_text_file.php?filename=%s&dir=data/historical/stdmet/'
# these are buoys within the drifter region that were active in 2012/2013
buoy_list = {
      46002:(42.614, -130.490),
      46005:(45.958, -131.000),
      46011:(34.956, -121.019),
      46012:(37.363, -122.881),
      46013:(38.242, -123.301),
      46015:(42.764, -124.832),
      46029:(46.159, -124.514),
      46047:(32.403, -119.536),
      46061:(47.353, -124.731),
      46087:(48.494, -124.728),
      46089:(45.893, -125.819),
      46211:(46.858, -124.244),
      46229:(43.767, -124.549),
      46239:(36.342, -122.102), 
      46246:(49.904, -145.243),
      46089:(45.893, -125.819),
     'cdea2':(56.388, -134.637)}

def compute_distance(lat1, long1, lat2, long2):
 
    # Convert latitude and longitude to 
    # spherical coordinates in radians.
    degrees_to_radians = math.pi/180.0
         
    # phi = 90 - latitude
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians
         
    # theta = longitude
    theta1 = long1*degrees_to_radians
    theta2 = long2*degrees_to_radians
         
    # Compute spherical distance from spherical coordinates.
         
    # For two locations in spherical coordinates 
    # (1, theta, phi) and (1, theta', phi')
    # cosine( arc length ) = 
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length
     
    cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) + 
           math.cos(phi1)*math.cos(phi2))
    arc = math.acos( cos )
    # multiply by radius of earth in km

 
    # Remember to multiply arc by the radius of the earth 
    # in your favorite set of units to get length.
    return arc*6371

def get_request(stationid, year, dirp=''):
    fname = "%sh%s.txt.gz" %(stationid, year)
    path = os.path.join(dirp, fname[:-3])
    if not os.path.exists(path):
        print("downloading from %s" %path)
        rr = requests.get(base_url %(fname))
        response = rr.text
        if 'Unable to access' in response:
            print("Unable to access data at %s" %rr.url)
            return ''
        fp = open(path, 'w')
        fp.write(rr.text)
        fp.close()
    return path

def parse_request(path):
    if not os.path.exists(path):
        return ''
    else:
        ss = pd.read_csv(path, 
                   delim_whitespace=True, 
                   header=0,
                   skiprows=[1],
                   parse_dates={'buoy_datetime':['#YY', 'MM', 
                                                 'DD', 'hh', 'mm']}, 
                   )
        ss['buoy_datetime'] = pd.to_datetime(ss['buoy_datetime'], 
                                  format="%Y %m %d %H %M")
        return ss

def get_all_buoys():
    buoys = pd.DataFrame()
    for buoyid in buoy_list.keys():
        for yr in [2012, 2013]:
            path = get_request(buoyid, yr, 'buoys')
            bpd = parse_request(path)
            if type(bpd) != str:
                bpd['buoy_lat'] = buoy_list[buoyid][0]
                bpd['buoy_lon'] = buoy_list[buoyid][1]
                bpd['buoy_id'] = buoyid
                buoys = buoys.append(bpd)
    return buoys

def find_nearest_buoy(lat, lon):
    buoy_dists = []
    for buoy in buoy_list.keys():
        blat, blon = buoy_list[buoy]
        dd = compute_distance(blat, blon, lat, lon)
        buoy_dists.append((buoy, dd))
    buoys_sorted = sorted(buoy_dists, key=lambda x: x[1])
    return buoys_sorted

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

def load_data(fname, drifter_type, launch_datetime='2012-06-01 00:00:00', end_datetime='2014-06-01 00:00:00'):
    """Input the name of the drifter file downloaded from the website. This function parses the two types of data,
    averaged measurements, M, and calibration measurements, C
    """
    min_bat = 11.5
    dval = open(fname, 'r')
    #initialize battery voltage
    bvalue = -1
    grf12 = -1
    cd = {"id":[], "cal_start_datetime":[], "sample_datetime":[], "num":[], 
           "x":[], "y":[], "z":[], "f":[], "temp":[], "lat":[], "lon":[], "bat":[],
          }
          
    md = {"id":[], "sample_datetime":[], "num":[], 
           "x":[], "y":[], "z":[], "f":[], "temp":[], "lat":[], "lon":[], "bat":[], 
          }
    nsams = 250
    calsx = [0]*nsams
    calsy = [0]*nsams
    calsz = [0]*nsams
    do_return = False
    cl = {"id":[],  "cal_start_datetime":[], "calsx":[], "calsy":[], "calsz":[],
          "temp":[], "lat":[], "lon":[], "bat":[]}
    
    for line in dval:
        line_vals = line.split(' ')
        
        line_vals = [x for x in line_vals if x!='']
        line_vals[-1] = line_vals[-1].strip()
      
        if line_vals[0] == 'S':
            # S: status line
            # S 000000000066760 123 270 2011.08.13 23:03:06 47.651360 -129.042221 8
            # S drifter_id rec_# day date time latitude longitude ?
            # remove the S character
            # this message is sent when the data is uploaded to the server
            #mstatus = line_vals[1:]
            d = pd.to_datetime("%s %s" %(line_vals[4], line_vals[5]))
            S_drifter_id = line_vals[1]
            S_record_num = line_vals[2]
            S_datetime = d
            S_lat = float(line_vals[6])
            S_lon = float(line_vals[7])
            
            
        if line_vals[0] == 'M':
            # M: measurement
            # M 1 -62475.9 -32540.4 -10721.9 19.39 47.9019 -128.9508 1.6 2011 224.80556 17.49
            # M meas_# xf yf zf temperature latitude longitude ? yr decimal_julian _day
            # convert julian day to something that is more easily readable
            mdt = convert_julian_frac(line_vals[10], line_vals[9])
    
            # now we have the measurement value, datetime, and battery values. This is the averaged vals
        
            # Always use the lat/lon, temp included here for the averaged data
            M_lat = float(line_vals[6])
            M_lon = float(line_vals[7])
            x, y, z, f = get_xyz(line_vals[2:5], drifter_type)

            md['lat'].append(M_lat)
            md['lon'].append(M_lon)
            md['id'].append(S_drifter_id)
            md["sample_datetime"].append(mdt)
            
            md['x'].append(z)
            md['y'].append(y)
            md['z'].append(z)
            md['f'].append(f)
            
            md['temp'].append(float(line_vals[5]))
            md['num'].append(int(line_vals[1]))
            md['bat'].append(bvalue)
            
            
        if line_vals[0] == 'C':
            # The date reported here is always the start time of the sample period
            # C: Calibration header
            # C 8 2011 225.12518
            # C id yr decimal_julian_day
            jdt = convert_julian_frac(line_vals[3], line_vals[2])
            # store the calibration value
            Cid = line_vals[1]
            # store the datetime
            Cdf = jdt
            
        if line_vals[0] == 'c':
            # calibration measurement, add this to the header value for each value
            # offset the time with the measurement frequency
            # the first few and last few of this cal value seem to be bad
            C_count = int(line_vals[1])
      
            cdt = Cdf + timedelta(0, C_count)
            x, y, z, f = get_xyz(line_vals[2:5], drifter_type)
            ctemp = float(line_vals[5])
            cd["sample_datetime"].append(cdt)
            cd['x'].append(x)
            cd['y'].append(y)
            cd['z'].append(z)

            cd['temp'].append(ctemp)
            cd['num'].append(C_count)
            cd['cal_start_datetime'].append(Cdf)
            cd['bat'].append(bvalue)
            cd['lat'].append(S_lat)
            cd['lon'].append(S_lon)
            cd['id'].append(Cid)
            cd['f'].append(f)

        if line_vals[0] == 'E':
            # E:
            # E 12.7 0
            # E battery_voltage ?
            bvalue = float(line_vals[1])
    
    cpd = pd.DataFrame.from_dict(cd, orient='columns')
    cpd = cpd.set_index('sample_datetime', drop=False)
    # ignore data before june 1, this will need to be tailored for each drifter
    cpd = cpd[cpd.index > pd.to_datetime(launch_datetime)]
    cpd = cpd[cpd.index < pd.to_datetime(end_datetime)]
    cpd = cpd[cpd['bat'] > min_bat]
    #remove bad latitude data
    cpd = cpd[cpd['lat'] != -90.0]
    cpd = cpd[cpd['lat'] != 0.0]

    mpd = pd.DataFrame.from_dict(md, orient='columns')
    mpd = mpd.set_index('sample_datetime', drop=False)
    # ignore data before june 1, this will need to be tailored for each drifter
    mpd = mpd[mpd.index > pd.to_datetime(launch_datetime)]
    mpd = mpd[mpd.index < pd.to_datetime(end_datetime)]
    
    # theoretically, APS should be good down to 4.95 V (+4.95V to +12V. Input current is 40mA.)
    mpd = mpd[mpd['bat'] > min_bat]
    #remove bad latitude data
    mpd = mpd[mpd['lat'] != -90.0]
    mpd = mpd[mpd['lat'] != 0.0]
    
    return mpd, cpd

def get_model(dt, lats, lons):
    alt = 0
    # set isv to 0 (main field) because that is what it is in example ??
    isv = 0
    # set itype to 1 (geodectic)
    itype = 1
    mx,my,mz,mf,yeardec = runigrf12(dt, isv, itype, [alt], lats, lons)
    return mx[0], my[0], mz[0], mf[0], yeardec

def get_model_df(df):
    xs = []
    ys = []
    zs = []
    fs = []    
    alt = 0
    # set isv to 0 (main field) because that is what it is in example ??
    isv = 0
    # set itype to 1 (geodectic)
    itype = 1
    for i in df.index:
        x,y,z,f,yr = mx,my,mz,mf,yeardec = runigrf12(i, isv, itype, alt, 
                                                     df.loc[i,'lat'], df.loc[i,'lon'])
        xs.append(x[0])
        ys.append(y[0])
        zs.append(z[0])
        fs.append(f[0])
    df.loc[:,'igrfx'] = xs
    df.loc[:,'igrfy'] = ys
    df.loc[:,'igrfz'] = zs
    df.loc[:,'igrff'] = fs
    return df

def to_total_field(x, y, z):
    """convert to total magnetic field"""
    return np.sqrt(x**2 + y**2 + z**2)

def fix_hmr(d):
    k = np.sign(d)
    dd = abs(int(d))*10
    # convert to hex and discard last byte
    hexd = "{0:#0{1}x}".format(dd,8)[:-2]
    dd = int(hexd, 16)
    # convert to nT
    dd = (dd*6.6667)
    return dd

def get_xyz(xyz, drifter_type="APS"):
    x,y,z = [float(a) for a in xyz]
    if drifter_type == 'HMR':
        # the APS sensors were written in the wrong format, convert to correct
        x = fix_hmr(x)
        y = fix_hmr(y)
        z = fix_hmr(z)
    
    f = to_total_field(x, y, z)
    return x, y, z, f

def get_buoy_data(buoypd, datetimes, lats, lons):
    """
    :param datetimes: np array of panda datetimes
    :param lats: np array of latitudes
    :param lons: np array of longitudes
    """
    
    datetimes = pd.to_datetime(datetimes)
    data = pd.DataFrame()
    
    for dt, lat, lon in zip(datetimes, lats, lons):
        buoys_sorted = find_nearest_buoy(lat, lon)
    
        for (buoy, dis) in buoys_sorted:
            # if within 300 km, use this buoy
            if dis < 500: 
                md = buoypd.loc[buoypd.loc[:,'buoy_id'] == buoy].copy(deep=True)
                md['sample_datetime'] = dt
                md['mins_diff'] = abs(md['sample_datetime']-md['buoy_datetime']).astype('timedelta64[m]')
                min_time = md['mins_diff'].min()
                
                ## if nearest val within x minutes, use it, otherwise, look further away
                if min_time < 180:
                    closest = md.loc[md.loc[:,'mins_diff'] == min_time]
                    data = data.append(closest)
                    # don't bother searching other distances
                    break
            else:
                # don't search farther away in sorted list
                #print('no buoys close enough to', lat, lon)
                break
    data.index = data['sample_datetime']
    return data
        
def parse_raw_files(drifter_data_dir, drifter_dict):
    """parse raw drifter files into provided dictionary and write to meas/cal/list files. This function is really slow, but only needs to be called once to write the parsed files.
    :param drifter_data_dir: relative path where meas/cal/list_name.txt files are stored
    :param drifter_dict: dictionary with names of drifter (ie sleepy) as key
    """
    buoypd = get_all_buoys()
    for dname in drifter_dict.keys():
        dpath = os.path.join(drifter_data_dir, 'drifter_' + dname + '.txt')
        print("Loading %s and writing measured and calibration data files" %dname)
        drifter_type = drifter_dict[dname]['type']
        mpd, cpd, lpd = load_data(dpath, drifter_type, drifter_dict[dname]['launch'], drifter_dict[dname]['end'])
        mpd = get_model_df(mpd)
        bpd = get_buoy_data(buoypd, mpd.index, mpd['lat'], mpd['lon'])
        # join buoy data with drifter data
        mbpd = pd.merge(mpd, bpd, left_index=True, right_index=True)
        mpath = os.path.join(drifter_data_dir, 'meas_' + dname + '.txt')
        mbpd.to_csv(mpath, header=True, sep=' ', index=True)
        drifter_dict[dname]['meas'] = mbpd
        
        cpd = get_model_df(cpd)
        bpd = get_buoy_data(buoypd, cpd.index, cpd['lat'], cpd['lon'])
        # join buoy data with drifter data
        cbpd = pd.merge(cpd, bpd, left_index=True, right_index=True)
        cpath = os.path.join(drifter_data_dir, 'cal_' + dname + '.txt')
        cbpd.to_csv(cpath, header=True, sep=' ', index=True)
        drifter_dict[dname]['cal'] = cbpd

    return drifter_dict

def parse_txt_files(drifter_data_dir, drifter_dict):
    """load files previously produced by parse_raw_files
    :param drifter_data_dir: relative path where meas/cal/list_name.txt files are stored
    :param drifter_dict: dictionary with names of drifter (ie sleepy) as key
    """
    for dname in drifter_dict.keys():
        print("Loading %s meas, cal, and list data files" %dname)
        mpath = os.path.join(drifter_data_dir, 'meas_' + dname + '.txt')
        mpd = pd.read_csv(mpath, header=0, sep=' ', parse_dates=0, index_col=0, low_memory=False)
        drifter_dict[dname]['meas'] = mpd
        
        cpath = os.path.join(drifter_data_dir, 'cal_' + dname + '.txt')
        cpd = pd.read_csv(cpath, header=0, sep=' ', parse_dates=0, index_col=0, low_memory=False)
        drifter_dict[dname]['cal'] = cpd
        
        lpath = os.path.join(drifter_data_dir, 'list_' + dname + '.txt')
        lpd = pd.read_csv(lpath, header=0, sep=' ', parse_dates=0, index_col=0, low_memory=False)
        drifter_dict[dname]['list'] = lpd
    return drifter_dict
