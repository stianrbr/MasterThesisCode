import numpy as np

def distance_2points(lat1, lon1, lat2, lon2):
    """
    Taken from: https://github.com/KonstantinChri/MET_waves

    Calculates the distance between two points specified by coordinates
    """
    R = 6371.0
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    distance = R * c  # in km
    return distance

def find_nearest(lon_model, lat_model, lat0, lon0, product):
    """
    Adopted from: https://github.com/KonstantinChri/MET_waves

    Finds the nearest sets of coordinates in the dataset
    """
    print('\t\tFinding nearest point...')
    dx = distance_2points(lat0, lon0, lat_model, lon_model)
    temp = dx.where(dx == dx.min(), drop=True)

    if product == "NORA3":
        lat0 = temp.longitude
        lon0 = temp.latitude
        x0 = temp.x
        y0 = temp.y
    elif product == "NORA10":
        lat0 = temp.longitude
        lon0 = temp.latitude
        x0 = temp.Xc
        y0 = temp.Yc
    elif product == "NORA3_wave":
        x0 = temp.rlat
        y0 = temp.rlon
        lat0 = 0
        lon0 = 0
    print("\t\tFound nearest")
    return x0, y0, lat0, lon0