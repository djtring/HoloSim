# All mathimatical calculations are pulled from the file "ALMA 3.0" which refrences the ALMA paper. 

import numpy as np 
from scipy.integrate import dblquad
from tkinter import messagebox

class Calculator:

    def calculate_step1(D, nu):
        c = 3 * 10**8  # Speed of light in m/s
        lam = c / (nu * 10**9)  # Wavelength in meters
        R_F = (2 * D**2) / lam  # Far-field cutoff in meters
        R_min = 5 * D  # Minimum transmitter distance in meters
        R_react = 0.62 * np.sqrt(D**3 / lam)  # Reactive near-field cutoff in meters
        return R_F, R_min, R_react
    
    def calculate_step2(a):
        delta_d = a * np.sqrt(2) / 2 * 10**2  # Calculate spatial resolution in cm
        return delta_d
    
    def calculate_step3(delta_d): 
        try: 
            # Retrieve inputs from the correct entry boxes
            D = float(entry_2.get().strip())  # Diameter input (D) in meters
            nu = float(entry_5.get().strip())  # Frequency input (nu) in GHz
            f_1 = float(entry_3.get().strip())  # Primary beam taper factor (entry_3)
            f_apo = float(entry_6.get().strip())  # Apodization smoothing factor (entry_6, default 1.3)
            f_osr = float(entry_7.get().strip())  # Oversampling factor between rows (entry_7)
            dtheta = float(entry_19.get().strip())  # Rotation rate of the dish antenna (arcsec/sec) (entry_19)

            # Integration time (t_int) calculation from the paper
            t_int = 6.2 * 10**4 * f_osr * f_1 * f_apo**2 / (dtheta * nu * D)  # seconds
            
            # Total map time (t_map) calculation from the paper
            t_map = 171768 * f_osr * f_1 * f_apo**2 * D / (dtheta * nu * delta_d**2)  # hours

            return t_int, t_map

        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numbers or check previous step.")
            return None, None

    def calculate_step4(delta_d):
        try: 
            # Retrieve inputs from the relevant entry boxes
            nu = float(entry_5.get().strip())  # Frequency input from entry_5 (in GHz)
            D = float(entry_2.get().strip())  # Diameter input from entry_2
            f_oss = float(entry_12.get().strip())  # User input for f_oss from entry_12
        
            # Constants from the paper
            f_1 = 1.13  # Primary beam taper factor
            f_apo = 1.3  # Apodization smoothing factor
            f_osr = 2.2  # Oversampling factor between rows

            # Ensure that delta_d is correctly pulled from Step 2
            # Here, delta_d comes from Step 2, already calculated, and is passed into this function
            print(f"delta_d (cm) from Step 2: {delta_d}")

            # Perform calculations according to the paper's (Alma) backend logic
            theta_ext = 1717.7 * f_1 * f_apo / (delta_d * nu)  # Angular extent of map in degrees
            print(f"Calculated theta_ext: {theta_ext}")

            theta_b = 61836.6 * f_1 / (nu * D)  # Beamwidth in arcsec
            

            theta_sr = theta_b / f_osr  # Sampling interval between rows in arcsec
            

            theta_ss = theta_b / f_oss  # Sampling interval along scan in arcsec
            

            # Return the results
            return theta_ext, theta_b, theta_sr, theta_ss

        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numbers or check previous step.")
            return None, None, None, None
        
    def claculate_step5(delta_d):
        try:
            # Retrieve the input for delta_z from entry_16
            delta_z = float(entry_16.get().strip())  # Surface Deformation (μm)

            # Retrieve D and nu from Step 1
            D = float(entry_2.get().strip())  # Diameter from Step 1
            nu = float(entry_5.get().strip())  # Frequency from Step 1

            # Constants
            f_apo = 1.3  # Example constant from previous steps
            f_osr = 2.2  # Example constant
            f_oss = float(entry_12.get().strip())  # User input for f_oss

            # Calculate wavelength (lam) using nu
            c = 3 * 10**8  # Speed of light in m/s
            lam = c / (nu * 10**9)  # Wavelength in meters (from Step 1)

            # delta_d is now passed from Step 2
            print(f"delta_d from Step 2: {delta_d} cm")

            # Perform the Step 5 calculations
            theta_point = 6 * 10**-4 * delta_z * nu  # Pointing accuracy
            
            N_row = D / (delta_d * 10**-2) * f_apo * f_osr  # Total number of rows
            
            SNR = 10 * np.log10(0.044 * lam * np.sqrt(N_row**2 / (f_oss * f_osr)) * 1 / delta_z * 1 / np.sqrt(f_apo) * 10**6)  # SNR calculation
            
            return theta_point, N_row, SNR

        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numbers or check previous step.")
            return None, None, None
        
    def calculate_step6(SNR):
         try:
            # Retrieve inputs from entry boxes
            D_t = float(entry_17.get().strip())     # Transmitter diameter (m)
            D = float(entry_2.get().strip())        # Receiver diameter (m)
            z = float(entry_8.get().strip())        # Distance between transmitter and receiver (m)
            T_sys = float(entry_9.get().strip())    # System temperature (K)
            B = float(entry_10.get().strip())       # Bandwidth in MHz
            nu = float(entry_5.get().strip())       # Frequency in GHz

            # Constants
            k = 1.38e-23        # Boltzmann constant (J/K)
            G = 33              # Gain (dB)
            c = 3e8             # Speed of light (m/s)

            # Convert units
            B_hz = B * 1e6      # Bandwidth in Hz
            nu_hz = nu * 1e9    # Frequency in Hz
            lam = c / nu_hz     # Wavelength (m)

            # Noise floor in dBm
            noise_floor = 10 * np.log10(k * B_hz * T_sys / 1e-3)

            # Beam waist at the transmitter
            w_0 = D_t / 2
            z_R = np.pi * w_0**2 / lam
            w_z = w_0 * np.sqrt(1 + (z / z_R)**2)

            # Define the integrand for the overlap integral
            def integrand(r, theta):
                return r * np.exp(-r**2 / w_z**2)

            # Calculate the overlap integral over the aperture
            overlap, _ = dblquad(
                integrand,
                0,
                D / 2,
                lambda r: 0,
                lambda r: 2 * np.pi
            )

            # Total power over the aperture
            total_power = np.pi * w_z**2 / 2

            # Beam-coupling efficiency
            eta = (overlap / total_power)**2
            eta_dB = 10 * np.log10(eta)

            # Calculate transmitter output power in dBm
            P_dB = noise_floor + SNR + eta_dB + G

            # Convert to watts
            P = 1e-3 * 10**(P_dB / 10)

            return P
         
         except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numbers or check previous step.")
            return None
        
        


