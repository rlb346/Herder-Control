#Not yet implemented

import pandas as pd

import Parameters as P

t = P.t

def save_data(k):
    dataframe1 = pd.DataFrame(t,columns = ("Time",))
    for i in range(len(ustoreelec[0])):
        dataframe1["U"+str(i)] = pd.Series(ustoreelec[:,i])  
    for j in range(n_particles):
        dataframe1["X"+str(j)] = pd.Series(xy[:,0,j])
        dataframe1["Y"+str(j)] = pd.Series(xy[:,1,j])
    dataframe2 = pd.DataFrame({"length": length, "rd": rd, "v0":v0, "W_1": W_1, "mu":mu,"D":D, "R": radius}, index = [0])
    dataframeCombined = pd.concat([dataframe1,dataframe2], axis = 1)
    kfilled = str(k).zfill(5)
    dataframeCombined.to_csv(f"{folder}/{filename}Data{kfilled}.csv")
    dataframeGrid = pd.DataFrame(grid)
    dataframeGrid.to_csv(f"{folder}/{filename}Grid{kfilled}.csv")
    return
print("Done")