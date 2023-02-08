#Not yet implemented

import pandas as pd

import Process as Pro



def save_data(k,t,xy):
    dataframe1 = pd.DataFrame(t,columns = ("Time",))
    dataframe1["t"] = pd.Series(t)  
    for j in range(n_particles):
        dataframe1["X"+str(j)] = pd.Series(xy[:,0,j])
        dataframe1["Y"+str(j)] = pd.Series(xy[:,1,j])
    kfilled = str(k).zfill(5)
    dataframeCombined.to_csv(f"{folder}/{filename}Data{kfilled}.csv")
    dataframeGrid = pd.DataFrame(grid)
    dataframeGrid.to_csv(f"{folder}/{filename}Grid{kfilled}.csv")
    return
print("Done")