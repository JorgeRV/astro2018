import os
import glob
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from numpy.lib import scimath as sm

def km_to_pc(dist):
    return dist/(3.086*10**13)


def smooth(x, y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')

    plt.plot(x, y, 'o')
    plt.plot(x, y_smooth, 'r-', lw=2)

    return y_smooth


# Returns the velocity (in the units of c) given the observed and base frequency
def doppler(nu_obs, nu_0):
    return (nu_obs - nu_0)*c/nu_0


# Calculate the distance between the LSR and the measured point. It gives 2 values
# l is the galactic latitude, vel is the velocity of the observed point
# ang_vel is the angular velocity of the LSR and r0 its distance to the galaxy center
def rdist(l, vel, gvel, r0):
    rad = l*np.pi/180
    # res1 = r0*(np.cos(rad) + np.sqrt(np.cos(rad)**2 + 1/(vel/(gvel*np.sin(rad)) + 1)**2 - 1))
    # res2 = r0*(np.cos(rad) - np.sqrt(np.cos(rad)**2 + 1/(vel/(gvel*np.sin(rad)) + 1)**2 - 1))

    r = (r0*gvel*np.sin(rad))/(gvel*np.sin(rad) + vel)
    print("r:", km_to_pc(r))
    print("Delta r:", km_to_pc(r - r0), "pc")

    res1 = r0 * np.cos(rad) + sm.sqrt(r ** 2 - (r0 * np.sin(rad)) ** 2)
    res2 = r0 * np.cos(rad) - sm.sqrt(r ** 2 - (r0 * np.sin(rad)) ** 2)

    print(km_to_pc(res1), km_to_pc(res2))

    if np.iscomplex(res1) or np.iscomplex(res2):
        res1, res2 = np.abs(res1), np.abs(res2)

    return res1, res2


# Obtain the path to the images. Print the number of images in total.
ext = "rad"  # Extension of the files to read
folder = "/datos/"
path = str(os.getcwd())
data_folder = Path(path + "/datos/")
npath = len(path) + len(folder)
path = path + folder + "*." + ext
archivos = glob.glob(path)
arch = len(archivos)
print("Files:", arch)

# Print the path to the images
for i in range(arch):
    print(i, " ", archivos[i])

# Parameters used in the calculation
freq_max = []  # List for max frequency found per angle
dist_list = []  # List of distances between the measured peak and the Solar System
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
    arname = filename[npath:]  # Mod the number to have cut to the name of the file only
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

    # Now we move the plot to the vertical 0 axis
    m_data = m_data - np.min(m_data)

    # Recalculating the median and the standard deviation
    median = np.median(m_data)
    sigma = np.std(m_data)

    # Range of frequencies measured (x-axis of the spectrum)
    freq_list = []
    for i in range(len(m_data)):
        freq_list.append(freq + (i + int(len(del_indx)/2))*diff)

    # Smooth the data
    m_data = smooth(freq_list, m_data, 7)

    # Plot the data
    plt.plot(freq_list, m_data)
    plt.title("Spectrum at " + angle)
    plt.xlabel("Frequency (MHz)")
    plt.ticklabel_format(style='plain', axis='x', scilimits=(0, 0))
    plt.savefig("spectrum_g" + angle + ".png")
    plt.show()

    # Take the peak's frequencies with an intensity greater than the threshold
    threshold = sigma
    peak_freqs = []

    # For the next angles, a threshold of 2 sigmas is used
    if angle in ["20", "160", "310", "320", "330"]:
        threshold = 2*sigma

    for hindx in range(1, len(m_data) - 1):
        ang_freqs = []
        height = m_data[hindx]
        if height > m_data[hindx - 1]:
            if height > m_data[hindx + 1]:
                if height > threshold:
                    peak_freqs.append(freq_list[hindx])

    print("Peaks found:", peak_freqs)

    # Calculate the distance s using positive and negative solutions
    dlist = [int(angle)]
    for f in peak_freqs:
        obs_vel = doppler(f, hydrogen_freq)
        corr_vel = obs_vel - vlsr
        print("\nFrequency:", f)
        print("Observed relative velocity:", obs_vel, "km/s")
        print("Corrected velocity with the VLSR:", corr_vel, "km/s")
        print("VLSR:", vlsr, "km/s")

        s1, s2 = rdist(int(angle), corr_vel, local_vel, local_dist)
        print("Distances obtained:", km_to_pc(s1), km_to_pc(s2), "pc")

        # Skip the values with complex numbers
        if np.isnan(s1) and np.isnan(s2):
            continue

        dlist.append([f, km_to_pc(s1), km_to_pc(s2)])

    dist_list.append(dlist)

    # Take the index of the maximum value of each spectrum
    imax = int(np.argmax(m_data))
    # Take the frequency corresponding with the index
    fmax = freq_list[imax]
    freq_max.append(fmax)
    print("Frequency of max intensity: ", freq, "MHz")

np_dist_list = np.array(dist_list)
np_dist_list.sort()
for i in np_dist_list:
    print(i)
