import pandas as pd
import xarray as xr
import os
import time

def download_year_NORA10(year, variables, lat, lon, savefolder):
    """
    Function for downloading NORA10 data.
    Download specified variables and saves it as a H5-file.
    """
    print(f"Downloading: {year} from NORA10EI")
    df = pd.DataFrame()
    start_date = '{}-01-01'.format(year)
    end_date = '{}-12-31'.format(year)
    date_list = pd.date_range(start=start_date, end=end_date, freq='m')
    for i in range(len(date_list)):
        for ii, d in enumerate(["01", "11", "21"]):
            print("\tDownloading: {}-{}-{} with process {}".format(date_list[i].strftime('%Y'), date_list[i].strftime('%m'), d, os.getpid()))
            infile = "https://thredds.met.no/thredds/dodsC/nora10ei/wam10ei_{}.nc".format(date_list[i].strftime('%Y%m'+d))
            for count in range(6):
                try:
                    ds = xr.open_dataset(infile)
                except:
                    print(f"Failed reading data set: {infile}, trying again")
                    time.sleep(10)
                else:
                    break
            if i == 0 and ii == 0:
                x0, y0, lat0, lon0 = find_nearest(ds.longitude, ds.latitude, lat, lon, product="NORA10")
            ds_selected = ds.sel(Xc=x0, Yc=y0)
            temp = ds_selected[variables].to_dataframe().reset_index()
            df = pd.concat([df, temp], ignore_index=True)
            del ds, ds_selected, temp
        print("\tDownloaded: {}-{}-{} with process {}".format(date_list[i].strftime('%Y'), date_list[i].strftime('%m'),
                                                              date_list[i].strftime('%d'), os.getpid()))
    df.to_hdf(savefolder + "\\NORA10_{}.h5".format(year), key='df', mode='w')


def download_year_NORA3(year, variables, lat, lon, savefolder):
    """
    Function for downloading NORA3 wind data.
    Download specified variables and saves it as a H5-file.
    """
    print(f"Downloading: {year} from NORA3")
    df = pd.DataFrame()
    start_date = '{}-01-01'.format(year)
    end_date = '{}-12-31'.format(year)
    date_list = pd.date_range(start=start_date, end=end_date, freq='d')
    for i in range(len(date_list)):
        print("\tDownloading: {}-{}-{} with process {}".format(date_list[i].strftime('%Y'), date_list[i].strftime('%m'),
                                             date_list[i].strftime('%d'), os.getpid()))
        for ii, h in enumerate(["00", "06", "12", "18"]):
            for iii, hh in enumerate(["003", "004", "005", "006", "007", "008"]):
                infile = f"https://thredds.met.no/thredds/dodsC/nora3/{date_list[i].strftime('%Y')}/{date_list[i].strftime('%m')}/{date_list[i].strftime('%d')}/{h}/fc{date_list[i].strftime('%Y')}{date_list[i].strftime('%m')}{date_list[i].strftime('%d')}{h}_{hh}_fp.nc"
                for count in range(6):
                    try:
                        ds = xr.open_dataset(infile)
                    except:
                        print(f"Failed reading data set: {infile}, trying again")
                        time.sleep(10)
                    else:
                        break
                if i == 0 and ii == 0 and iii == 0:
                    x0, y0, lat0, lon0 = find_nearest(ds.longitude, ds.latitude, lat, lon, product="NORA3")
                ds_selected = ds.sel(x=x0.values, y=y0.values, drop=True)
                temp = ds_selected[variables].squeeze(dim=["height4", "x", "y"], drop=True).to_dataframe().reset_index()
                df = pd.concat([df, temp], ignore_index=True)
                del ds, ds_selected, temp
        print("\tDownloaded: {}-{}-{} with process {}".format(date_list[i].strftime('%Y'), date_list[i].strftime('%m'),
                                                             date_list[i].strftime('%d'), os.getpid()))
    df.to_hdf(savefolder + "\\NORA3_{}_lat_{}_lon_{}.h5".format(year,str(lat).replace(".","_"),str(lon).replace(".", "_")), key='df', mode='w')

def download_year_NORA3_wave(year, variables, lat, lon, savefolder):
    """
    Function for downloading NORA10 data.
    Download specified variables and saves it as a H5-file.
    """
    print(f"Downloading: {year} from NORA3 - Windsurfer")
    df = pd.DataFrame()
    start_date = '{}-01-01'.format(year)
    end_date = '{}-12-31'.format(year)
    date_list = pd.date_range(start=start_date, end=end_date, freq='d')
    for i in range(len(date_list)):
        print("\tDownloading: {}-{}-{} with process {}".format(date_list[i].strftime('%Y'), date_list[i].strftime('%m'),
                                             date_list[i].strftime('%d'), os.getpid()))
        infile = f"https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_files/{date_list[i].strftime('%Y')}/{date_list[i].strftime('%m')}/{date_list[i].strftime('%Y')}{date_list[i].strftime('%m')}{date_list[i].strftime('%d')}_MyWam3km_hindcast.nc"
        for count in range(6):
            try:
                ds = xr.open_dataset(infile)
            except:
                print(f"Failed reading data set: {infile}, trying again")
                time.sleep(10)
            else:
                break
        if i == 0:
            x0, y0, lat0, lon0 = find_nearest(ds.longitude, ds.latitude, lat, lon, product="NORA3_wave")
        ds_selected = ds.sel(rlat=x0.values, rlon=y0.values, drop=True)
        temp = ds_selected[variables].squeeze(dim=["rlat", "rlon"], drop=True).to_dataframe(dim_order=["time"]).reset_index()
        df = pd.concat([df, temp], ignore_index=True)
        del ds, ds_selected, temp
        print("\tDownloaded: {}-{}-{} with process {}".format(date_list[i].strftime('%Y'), date_list[i].strftime('%m'),
                                                             date_list[i].strftime('%d'), os.getpid()))
    df.to_hdf(savefolder + "\\NORA3_{}_lat_{}_lon_{}.h5".format(year,str(lat).replace(".","_"),str(lon).replace(".", "_")), key='df', mode='w')