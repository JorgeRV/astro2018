import os
import glob
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


# Returns the velocity (in the units of c) given the observed and base frequency
def doppler(nu_obs, nu_0):
    return (nu_obs - nu_0)*c/nu_0


# Calculate the distance between the LSR and the measured point. It gives 2 values
# l is the galactic latitude, vel is the velocity of the observed point
# ang_vel is the angular velocity of the LSR and r0 its distance to the galaxy center
def rdist(l, vel, gvel, r0):
    rad = l*np.pi/180
    res1 = r0*(np.cos(rad) + np.sqrt(np.cos(rad)**2 + 1/(vel/(gvel*np.sin(rad)) + 1)**2 - 1))
    res2 = r0*(np.cos(rad) - np.sqrt(np.cos(rad)**2 + 1/(vel/(gvel*np.sin(rad)) + 1)**2 - 1))
    return res1, res2


# Obtain the path to the images. Print the number of images in total.
ext = "rad"  # Extension of the files to read
path = str(os.getcwd())
data_folder = Path(path + "/datos/")
path = path + "/datos/*." + ext
archivos = glob.glob(path)
arch = len(archivos)
print("Files:", arch)

# Print the path to the images
for i in range(arch):
    print(i, " ", archivos[i])

# Parameters used in the calculation
freq_max = []
vlsr_list = []
local_vel = 220  # km/s
local_dist = 8*3.086*10**16  # km
hydrogen_freq = 1420.406  # MHz
c = 3*10**5  # km/s

# Cycle to use all of the images
for nfile in range(len(archivos)):
    # print("---------------------------------------------------------------")
    # Filename: the path to a specified image which is going to be analyzed.
    # Arname: the name of the image to be analyzed
    filename = archivos[nfile]  # The value can change to select other pictures. To use all see bottom of the code
    arname = filename[106:]  # Mod the number to have cut to the name of the file only
    angle = arname[1:-4]
    print("\n\nFile selected:", arname, "\n")
    print("Galactic angle:", angle)

    # Import the file to be analyzed
    file_to_open = data_folder / str(arname)
    raw_data = open(file_to_open)
    rdata = raw_data.read().splitlines()

    # Extracting the useful data (intensities and the frequency used)
    data = []
    freq, diff, vlsr = 0, 0, 0
    for i in range(2, len(rdata) - 3):
        line = rdata[i].split("\t")
        fline = line[9:len(line) - 2]
        fline = list(map(float, fline))
        data.append(fline)

        # Take the frequency used and the interval between frequency
        freq, diff = float(line[5]), float(line[6])

        # Take the VLSR given in the data in m/s
        vlsr = float(line[-1])
        vlsr_list.append(vlsr)

    data = np.array(data)

    # Take the vertical mean of the data
    m_data = np.mean(data, axis=0)

    # Calculating the median and the standard deviation
    median = np.median(m_data)
    sigma = np.std(m_data)

    # We need to eliminate the noise of the spectrum
    del_indx = []
    for i in range(len(m_data)):
        if m_data[i] < median - 1*sigma:
            del_indx.append(i)

    m_data = np.delete(m_data, del_indx)

    # The data is reversed, so we ned to invert the list
    m_data = m_data[::-1]

    # Range of frequencies measured (x-axis of the spectrum)
    freq_list = []
    for i in range(len(m_data)):
        freq_list.append(freq + (i + int(len(del_indx)/2))*diff)

    plt.plot(freq_list, m_data)
    plt.title("Spectrum at " + angle)
    plt.xlabel("Frequency (MHz)")
    plt.ticklabel_format(style='plain', axis='x', scilimits=(0, 0))
    plt.savefig("spectrum_g" + angle + ".png")
    plt.show()

    # Take the index of the maximum value of each spectrum
    imax = int(np.argmax(m_data))
    # Take the frequency corresponding with the index
    fmax = freq_list[imax]
    freq_max.append(fmax)
    print("Frequency of max intensity: ", freq, "MHz")

    obs_vel = doppler(fmax, hydrogen_freq)
    corr_vel = doppler(fmax, hydrogen_freq) - vlsr
    print("Observed relative velocity:", obs_vel, "km/s")
    print("Corrected velocity with the VLSR:", corr_vel, "km/s")

    print("VLSR:", vlsr, "km/s")

    s1, s2 = rdist(int(angle), corr_vel, local_vel, local_dist)

    print("Distances obtained:", s1/(3.086*10**16), s2/(3.086*10**16), "kpc")
