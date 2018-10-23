import os
import glob
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# Obtain the path to the images. Print the number of images in total.
path = str(os.getcwd())
ext = "fit"  # Extension of the files to read
folders = ["20171030", "20171113"]  # Folders with the data to analyze [20171030, 20171113]
pf = "AST0111_Fits/"  # Parent folder which contains the other folders
nf = len(folders)

# Add the paths from every folder into one list and the number of files in each folder in another list
paths = []
npaths = []
for f in folders:
    paths.append(path + "/" + pf + f + "/*." + ext)
    npaths.append(len(path) + len(f))

# Add the files from every folder into one list and the number of files in each folder in another list
files = []
nfiles = []
for p in paths:
    pt = glob.glob(p)
    n = len(pt)
    files.append(pt)
    nfiles.append(n)
    print(p)
    print("Files:", n)

# Print the files obtained
for f in range(nf):
    print("")
    for n in range(nfiles[f]):
        print(n, " ", files[f][n])
    print("")

# Import the data -> data[folder][file]
datas = []
exposures = []
gains = []
filters = []
for f in range(nf):
    file_datas = []
    file_exposures = []
    file_gains = []
    file_filters = []

    for n in range(nfiles[f]):
        filename = files[f][n]

        with fits.open(filename) as hdul:
            hdr = hdul[0].header
            info = hdul.info()
            data = hdul[0].data
            data_norm = data/np.sum(data)

            file_datas.append(data)
            file_gains.append(hdr['EGAIN'])
            file_exposures.append(hdr['EXPOSURE'])
            file_filters.append(hdr['FILTER'])

    datas.append(file_datas)
    gains.append(file_gains)
    exposures.append(file_exposures)
    filters.append(file_filters)

# Identify the index for each filter in each folder
nRl, nGl, nBl = [], [], []
for f in range(nf):
    list_R = []
    list_G = []
    list_B = []

    for n in range(nfiles[f]):
        color = filters[f][n]

        if color == 'Red':
            list_R.append(n)

        if color == 'Green':
            list_G.append(n)

        if color == 'Blue':
            list_B.append(n)

    nRl.append(list_R)
    nGl.append(list_G)
    nBl.append(list_B)


# TESTING
indx = 0
print("\nFolder selected:", folders[indx])
nR, nG, nB = nRl[indx][0], nGl[indx][0], nBl[indx][0]
dataR = datas[indx][nR]
dataG = datas[indx][nG]
dataB = datas[indx][nB]

plt.imshow(-dataR, cmap='Reds')
plt.show()
plt.imshow(-dataG, cmap='Greens')
plt.show()
plt.imshow(-dataB, cmap='Blues')
plt.show()
plt.imshow(dataR + dataG + dataB, cmap='nipy_spectral')
plt.show()

dataR_log = -2.5*np.log(dataR/exposures[indx][nR])
dataG_log = -2.5*np.log(dataG/exposures[indx][nG])
dataB_log = -2.5*np.log(dataB/exposures[indx][nB])

plt.imshow(dataR_log, cmap='nipy_spectral')
plt.show()
plt.imshow(dataG_log, cmap='nipy_spectral')
plt.show()
plt.imshow(dataB_log, cmap='nipy_spectral')
plt.show()

print(np.mean(dataR_log))
print(np.max(dataR_log))
print(np.min(dataR_log))

plt.plot(dataR_log)
plt.show()

dataBG_log = dataB_log - dataG_log
plt.imshow(dataBG_log, cmap='inferno')
plt.show()
