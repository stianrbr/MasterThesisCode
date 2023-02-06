import pandas as pd
import os
import seaborn as sns
from inputs import savefolder
sns.set(rc={'text.usetex' : True})
sns.set_theme(style="whitegrid", context="paper")

datafolder = savefolder
files = os.listdir(datafolder)

data_files = [f for f in files if f.endswith(".h5")]
df = pd.DataFrame()

for d in data_files:
    temp = pd.read_hdf(datafolder+"\\"+d, key="df")
    df = pd.concat([df, temp], ignore_index=True)

try:
    os.mkdir(datafolder+"\\Combined_data\\")
except FileExistsError:
    pass

df.to_hdf(datafolder+"\\Combined_data\\combined_dataset.h5", key="df")
