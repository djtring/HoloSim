# All mathimatical calculations are pulled from the file "HoloSimBackend.ipynb"

import numpy as np 
from scipy.integrate import dblquad

class Calculator:

    def calculate_step1(self, D, nu):
        c = 3 * 10**8  # Speed of light in m/s
        lam = c / (nu * 1e9)  # Wavelength in meters
        R_F = (2 * D**2) / lam  # Far-field cutoff in meters
        R_min = 5 * D  # Minimum transmitter distance in meters
        R_react = 0.62 * np.sqrt(D**3 / lam)  # Reactive near-field cutoff in meters
        return R_F, R_min, R_react

    def calculate_step2(self, a):
        delta_d = a * np.sqrt(2) / 2 * 1e2  # Spatial resolution in cm
        return delta_d

    def calculate_step3(self, delta_d, D, nu, f_1, f_apo, f_osr, dtheta):
        t_int = 6.2e4 * f_osr * f_1 * f_apo**2 / (dtheta * nu * D)  # Integration time in seconds
        t_map = 171768 * f_osr * f_1 * f_apo**2 * D / (dtheta * nu * delta_d**2)  # Total map time in hours
        return t_int, t_map

    def calculate_step4(self, delta_d, nu, D, f_oss):
        f_1 = 1.13  # Primary beam taper factor
        f_apo = 1.3  # Apodization smoothing factor
        f_osr = 2.2  # Oversampling factor between rows

        theta_ext = 1717.7 * f_1 * f_apo / (delta_d * nu)  # Angular extent of map in degrees
        theta_b = 61836.6 * f_1 / (nu * D)  # Beamwidth in arcsec
        theta_sr = theta_b / f_osr  # Sampling interval between rows in arcsec
        theta_ss = theta_b / f_oss  # Sampling interval along scan in arcsec
        return theta_ext, theta_b, theta_sr, theta_ss

    def calculate_step5(self, delta_d, delta_z, D, nu, f_apo, f_osr, f_oss):
        c = 3e8
        lam = c / (nu * 1e9)
        theta_point = 6e-4 * delta_z * nu  # Pointing accuracy
        N_row = D / (delta_d * 1e-2) * f_apo * f_osr  # Total number of rows
        SNR = 10 * np.log10(
            0.044 * lam * np.sqrt(N_row**2 / (f_oss * f_osr)) * 1 / delta_z * 1 / np.sqrt(f_apo) * 1e6
        )  # SNR calculation
        return theta_point, N_row, SNR

    def calculate_step6(self, SNR, D_t, D, z, T_sys, B, nu):
        k = 1.38e-23        # Boltzmann constant (J/K)
        G = 33              # Gain (dB)
        c = 3e8             # Speed of light (m/s)

        B_hz = B * 1e6      # Bandwidth in Hz
        nu_hz = nu * 1e9    # Frequency in Hz
        lam = c / nu_hz     # Wavelength (m)

        noise_floor = 10 * np.log10(k * B_hz * T_sys / 1e-3)

        w_0 = D_t / 2
        z_R = np.pi * w_0**2 / lam
        w_z = w_0 * np.sqrt(1 + (z / z_R)**2)

        def integrand(r, theta):
            return r * np.exp(-r**2 / w_z**2)

        overlap, _ = dblquad(
            integrand,
            0,
            D / 2,
            lambda r: 0,
            lambda r: 2 * np.pi
        )

        total_power = np.pi * w_z**2 / 2
        eta = (overlap / total_power)**2
        eta_dB = 10 * np.log10(eta)

        P_dB = noise_floor + SNR + eta_dB + G
        P = 1e-3 * 10**(P_dB / 10)
        return P