import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt


# Flux in erg cm–2 s–1 Hz–1 given an AB magnitude
def flux_ab(mag_ab):
    return 10**(-0.4*(mag_ab + 48.6))


# Gives the frequency from a given wavelength (Angstroms) in Hz
def frecuency(wlen):
    return c*(10**10)/wlen


# Calculates the exposure time needed given the parameters
def calculate_exp(mv, sn, wlen):
    v_ab = mv - 0.044  # Conversion to AB magnitude from Johnson/Cousins +/- 0.004
    flux_nu = flux_ab(v_ab)  # This flux is in (erg cm–2 s–1 Hz–1)
    flux_lambda = flux_nu*c/(wlen**2)  # This flux is in (erg cm–2 s–1 Å–1)

    print("\nThe flux is: ", flux_lambda, "(erg cm–2 s–1 Å–1)")

    nr_ph = flux_lambda*(wave/(h*c))  # Number of photons in (Nr_ph cm−2 s−1 Å−1)
    print("Number of photons: ", nr_ph, "(Nr_ph cm−2 s−1 Å−1)")

    prim_area = np.pi*(r_telescope**2)  # Area of the primary telescope
    area = prim_area*(1 - 1/10)  # Effective area of the telescope. The secondary mirror covers 10% of the area

    nr_ph = nr_ph*area  # Number of photons collected by the telescope
    print("Number of photons captured in the mirror: ", nr_ph, "(Nr_ph s−1 Å−1)")

    slit_loss = sp.erf(seeing/(slit_width*np.sqrt(2)))  # Fraction of the photons remaining after passing the slit
    nr_ph = nr_ph*slit_loss  # Number of photons that pass the slit of the spectrograph
    print("Number of photons after passing the slit of the spectrograph: ", nr_ph, "(Nr_ph s−1 Å−1) (Slit loss:")

    nr_ph = nr_ph*eff_f4_600_5500  # Number of photons that make into the CCD
    print("Number of photons that makes into the CCD: ", nr_ph, "(Nr_ph s−1 Å−1)")

    nr_ph = nr_ph * disp_f4_600_5500  # Number of photons per pixel
    print("Number of photons captured per pixel: ", nr_ph, "(Nr_ph s−1 pix−1)")

    flux_sky = flux_ab(sky_mag)*c/(wlen**2)  # This flux is in (erg cm–2 s–1 Å−1 ☐”–1)
    nr_pix = int(3*pix_seeing*pix_width)  # Number of pixels in 3" for seeing (3 sigmas) and 1" width of slit (area)
    flux_sky = flux_sky*nr_pix*ccd_scale  # This flux is in (erg cm–2 s–1 Å−1)
    n_sky = flux_sky*(wlen/(h*c))*area*eff_f4_600_5500*disp_f4_600_5500  # Number of photons of the sky per pixel in the CCD
    print("Number of photons captured from the sky per pixel: ", n_sky, "(Nr_ph s−1 pix−1)")

    time = sn**2*(nr_ph + n_sky + n_dark**2) + np.sqrt(sn**4*(nr_ph + n_sky + n_dark**2)**2 + 4*nr_ph**2*n_ccd**2)
    time = time/(2*nr_ph**2)

    print("\nThe observation will require approximately", time, "seconds to complete for a S/N of", sn)
    print("Noise contributions: Star:", np.sqrt(nr_ph*time), " Sky:", np.sqrt(n_sky*time), " Dark:", np.sqrt(n_dark**2*time), " CCD:", n_ccd**2)

    times = np.array([i for i in range(3000)])
    sn_as_time = sn_time(nr_ph, n_sky, n_dark, n_ccd, times)

    plt.plot(times, sn_as_time)
    plt.xlabel("Time (s)")
    plt.ylabel("S/N")
    plt.savefig("Grafico_sn_vs_tiempo.png")
    plt.show()

    slits = np.linspace(0.01, 5, 100)
    nr_ph_slit = number_photons(flux_lambda, area, seeing, slits, eff_f4_600_5500, disp_f4_600_5500)
    sn_as_slit = sn_time(nr_ph_slit, sky_photons(flux_sky, seeing, slits, area), n_dark, n_ccd, time)

    plt.plot(slits, sn_as_slit)
    plt.xlabel("Slit width (\")")
    plt.ylabel("S/N")
    plt.savefig("Grafico_sn_vs_slit.png")
    plt.show()

    seeings = np.linspace(0.01, 5, 100)
    nr_ph_seeing = number_photons(flux_lambda, area, seeings, slit_width, eff_f4_600_5500, disp_f4_600_5500)
    sn_as_seeing = sn_time(nr_ph_seeing, sky_photons(flux_sky, seeings, slit_width, area), n_dark, n_ccd, time)

    plt.plot(seeings, sn_as_seeing)
    plt.xlabel("Seeing (\")")
    plt.ylabel("S/N")
    plt.savefig("Grafico_sn_vs_seeing.png")
    plt.show()

    plt.clf()

    return exit()


def number_photons(flux, area_tel, seeing, slit_w, efficiency, dispersion):
    return flux*(wave/(h*c))*area_tel*sp.erf(seeing/(slit_w*np.sqrt(2)))*efficiency*dispersion


def sky_photons(flux_sky, seeing, slit_w, area_tel):
    nr_pix = 3 * seeing * slit_w/ccd_scale**2
    return flux_sky * nr_pix * ccd_scale * (wave / (h * c)) * area_tel * eff_f4_600_5500 * disp_f4_600_5500


def sn_time(nr_ph, noise_sky, noise_dark, ron, t):
    return nr_ph*t/np.sqrt(nr_ph*t + noise_sky*t + noise_dark**2*t + ron**2)


def exposure_arg(boolarg):
    global seeing
    global slit_width
    if boolarg:
        print("Please select the new parameters:")
        mv = float(input("Magnitude V: "))
        sn = float(input("Signal To Noise Ratio: "))
        seeing = float(input("Seeing range (\"): "))
        slit_width = float(input("Slit width (\"): "))
        calculate_exp(mv, sn, wave)
    if not boolarg:
        print("The parameters are set as the given conditions.")
        mv = mag_v
        sn = expc_sn
        calculate_exp(mv, sn, wave)
    return


def menu(op):
    if op == "n":
        return exposure_arg(False)
    if op == "y":
        return exposure_arg(True)
    else:
        return exit()


def choose():
    op = str.lower(input("Do you want to change parameters? Press t to terminate the program. (y/n/t) "))
    if op != "y" and op != "n" and op != "t":
        print("Please select Y, N or T. (Yes/No/Terminate)")
        op = choose()
    return op


# Parameters and constants values
expc_sn = 70  # SN ratio = 70/pix
mag_v = 14  # Johnson/Cousins magnitude green V = 14
mag_vb = 0.8  # Johnson/Cousins magnitude V-B = 0.8
mag_b = mag_v + mag_vb  # Johnson/Cousins magnitude blue B = 14.8
wave = 5500  # Wavelength is in Angstroms Å
h = 6.626*10**(-27)  # Planck constant in erg*s
c = 3*10**18  # Speed of light in Angstroms Å
r_telescope = 6500  # Radius of the primary telescope in cm
eff_f4_600_5500 = 0.24  # Throughput of f/4 channel at 5500 Å with 600 lines/mm and 8.6º blaze angle (no units)
disp_f4_600_5500 = 0.378  # Dispersion of f/4 channel at 5500 Å with 600 lines/mm and 8.6º blaze angle in Å/pixel
seeing = 1  # Seeing of the telescope in ". For LCO is 1". FWHM = seeing = sigma = 1" in this case (standard deviation)
slit_width = 1  # Width of the slit in " considered to see the star
ccd_scale = 0.111  # Scale of the f/4 detector in "/pix
pix_width = slit_width/ccd_scale  # Number of pixels corresponding to the width of the slit
pix_seeing = seeing/ccd_scale  # Number of pixels corresponding to the seeing in the f/4 channel
n_dark = 27.10/3600  # Dark noise contribution for f/4 in (e- s-1 pix-1)
n_ccd = 4.7375  # RON contribution (CCD read noise). Taken from the mean for the 8 CCD of f/4 (2x2, fast) in e-/pix
n_ccd_slow = 2.65  # RON contribution (CCD read noise). Taken from the mean for the 8 CCD of f/4 (2x2, slow) in e-/pix
sky_mag = 21.8  # AB magnitude of the sky for V filter (24.7 for results like LCO)
airmass = 1  # Mass of air that the signal travels from space to telescope (1 to 3) in "

# Initial configuration and program initialization
author = "Jorge Riveros"
name_target = "Globular Cluster 47 Tucanae"
print("Exposure Time Calculator made by", author)
print("This calculator is going to determine the exposure time needed to have a Signal To Noise Ratio (SN) of", expc_sn)
print("The target is the", name_target, "with magnitude V of", mag_v, "and magnitude V-B of", mag_vb, "at", wave, "Å")

option = choose()
menu(option)
