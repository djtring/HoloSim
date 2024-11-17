#Holosim w/o graphs

import os
import tkinter as tk
from pathlib import Path
from tkinter import Tk, Canvas, Entry, Button, PhotoImage, messagebox, Toplevel, Label
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from scipy.integrate import dblquad
from tkinter import Toplevel, messagebox, StringVar, ttk, LabelFrame, Label, Entry, Button, filedialog
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import fitz 
from PIL import Image, ImageTk

OUTPUT_PATH = Path(__file__).resolve().parent
ASSETS_PATH = OUTPUT_PATH

def relative_to_assets(path: str) -> Path:
    return ASSETS_PATH / Path(path)

def calculate_step1(D, nu):
    c = 3 * 10**8  # Speed of light in m/s
    lam = c / (nu * 10**9)  # Wavelength in meters
    R_F = (2 * D**2) / lam  # Far-field cutoff in meters
    R_min = 5 * D  # Minimum transmitter distance in meters
    R_react = 0.62 * np.sqrt(D**3 / lam)  # Reactive near-field cutoff in meters
    return R_F, R_min, R_react

def calculate_step2(a):
    delta_d = a * np.sqrt(2) / 2 * 10**2  # Calculate spatial resolution in cm
    print(f"Calculated delta_d: {delta_d}")  # Debugging print
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

        # Debugging print statements to check intermediate values
        print(f"D = {D}, nu = {nu}, f_1 = {f_1}, f_apo = {f_apo}, f_osr = {f_osr}, dtheta = {dtheta}, delta_d = {delta_d}")

        # Integration time (t_int) calculation from the paper
        t_int = 6.2 * 10**4 * f_osr * f_1 * f_apo**2 / (dtheta * nu * D)  # seconds
        print(f"t_int = {t_int}")  # Debugging output

        # Total map time (t_map) calculation from the paper
        t_map = 171768 * f_osr * f_1 * f_apo**2 * D / (dtheta * nu * delta_d**2)  # hours
        print(f"t_map = {t_map}")  # Debugging output

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

        # Perform calculations according to the paper's backend logic
        theta_ext = 1717.7 * f_1 * f_apo / (delta_d * nu)  # Angular extent of map in degrees
        print(f"Calculated theta_ext: {theta_ext}")

        theta_b = 61836.6 * f_1 / (nu * D)  # Beamwidth in arcsec
        print(f"Calculated theta_b: {theta_b}")

        theta_sr = theta_b / f_osr  # Sampling interval between rows in arcsec
        print(f"Calculated theta_sr: {theta_sr}")

        theta_ss = theta_b / f_oss  # Sampling interval along scan in arcsec
        print(f"Calculated theta_ss: {theta_ss}")

        # Return the results
        return theta_ext, theta_b, theta_sr, theta_ss

    except ValueError:
        messagebox.showerror("Input Error", "Please enter valid numbers or check previous step.")
        return None, None, None, None
    
def calculate_step5(delta_d):
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
        print("theta_point < " + str(theta_point))  # deg

        N_row = D / (delta_d * 10**-2) * f_apo * f_osr  # Total number of rows
        print('N_row = ' + str(N_row))

        SNR = 10 * np.log10(0.044 * lam * np.sqrt(N_row**2 / (f_oss * f_osr)) * 1 / delta_z * 1 / np.sqrt(f_apo) * 10**6)  # SNR calculation
        print('SNR > ' + str(SNR))

        return theta_point, N_row, SNR

    except ValueError:
        messagebox.showerror("Input Error", "Please enter valid numbers or check previous step.")
        return None, None, None

def calculate_step6():
    try:
        # Retrieve inputs from entry boxes
        D_t = float(entry_17.get().strip())  # Transmitter diameter
        D = float(entry_2.get().strip())     # Diameter from Step 1 (used as Receiver diameter)
        z = float(entry_8.get().strip())     # Distance between transmitter and receiver
        T_sys = float(entry_9.get().strip()) # System temperature
        B = float(entry_10.get().strip())    # Bandwidth in MHz

        # SNR from previous step (or set a default value)
        SNR = 28.23  # This should match your value from Step 5

        # Constants from the paper
        k = 1.38 * 10**-23  # Boltzmann constant in J/K
        G = 33              # Gain in dB

        # Noise floor
        noise_floor = 10 * np.log10(k * B * 1e6 * T_sys / 1e-3)
        print(f"Noise floor: {noise_floor}")

        # Wavelength (lam)
        lam = 3e8 / (float(entry_5.get()) * 1e9)  # Frequency from Step 1 (GHz)
        print(f"Wavelength (lambda): {lam}")

        # Beam waist at the transmitter
        w_0 = D_t / 2                  # Beam waist
        z_R = np.pi * w_0**2 / lam     # Rayleigh range
        w_z = w_0 * np.sqrt(1 + (z / z_R)**2)  # Beam radius at the receiver location
        print(f"Beam waist (w_0): {w_0}, Rayleigh range (z_R): {z_R}, Beam radius (w_z): {w_z}")

        # Define integrand for the overlap integral
        def integrand(r, theta):
            return r * np.exp(-2 * r**2 / w_z**2)

        # Calculate the overlap integral over the aperture
        overlap, _ = dblquad(
            integrand,
            0,
            D / 2,                       # Use D from Step 1
            lambda r: 0,
            lambda r: 2 * np.pi
        )
        print(f"Overlap integral: {overlap}")

        # Calculate beam-coupling efficiency
        total_power = np.pi * w_z**2 / 2
        eta = (overlap / total_power)**2
        eta_dB = 10 * np.log10(eta)
        print(f"Beam-coupling efficiency (η): {eta}, Beam-coupling efficiency in dB (η_dB): {eta_dB}")

        # Calculate transmitter output power in dBm
        P_dB = noise_floor + SNR - eta_dB + G
        print(f"Transmitter output power (P_dB): {P_dB}")

        # Convert to watts
        P = 1e-3 * 10**(P_dB / 10)
        print(f"Transmitter output power: {P} W")

        return P

    except ValueError:
        messagebox.showerror("Input Error", "Please enter valid numbers for Step 6 inputs.")
        return None
    


def calculate_all_steps():
    global result_text_items

    try:
        # Step 1: Retrieve and process inputs for Diameter (D) and Frequency (ν)
        D_input = entry_2.get().strip()    # Diameter input
        nu_input = entry_5.get().strip()   # Frequency input

        print(f"Retrieved D_input: {D_input}")
        print(f"Retrieved nu_input: {nu_input}")

        if not D_input or not nu_input:
            raise ValueError("Empty input for Step 1")

        D = float(D_input)
        nu = float(nu_input)

        # Step 1 calculations
        R_F, R_min, R_react = calculate_step1(D, nu)

        # Step 2: Retrieve and process input for parameter 'a'
        a_input = entry_1.get().strip()
        print(f"Retrieved a_input: {a_input}")

        if not a_input:
            raise ValueError("Empty input for Step 2")

        a = float(a_input)

        # Step 2 calculations
        delta_d = calculate_step2(a)

        # Step 3: Retrieve delta_d from Step 2, and use it in Step 3 calculations
        t_int, t_map = calculate_step3(delta_d)

        # Step 4: Retrieve inputs for f_oss (Frequency Offset)
        f_oss_input = entry_12.get().strip()
        print(f"Retrieved f_oss_input: {f_oss_input}")

        if not f_oss_input:
            raise ValueError("Empty input for Step 4")

        f_oss = float(f_oss_input)

        # Step 4 calculations using delta_d from Step 2
        theta_ext, theta_b, theta_sr, theta_ss = calculate_step4(delta_d)

        # Step 5: Retrieve input for surface deformation (delta_z)
        delta_z_input = entry_16.get().strip()
        print(f"Retrieved delta_z_input: {delta_z_input}")

        if not delta_z_input:
            raise ValueError("Empty input for Step 5")

        delta_z = float(delta_z_input)

        # Step 5 calculations using delta_d from Step 2
        theta_point, N_row, SNR = calculate_step5(delta_d)

        # Step 6: Retrieve inputs for transmitter calculations
        D_t_input = entry_17.get().strip()
        z_input = entry_8.get().strip()
        T_sys_input = entry_9.get().strip()
        B_input = entry_10.get().strip()

        print(f"Retrieved D_t_input: {D_t_input}")
        print(f"Retrieved z_input: {z_input}")
        print(f"Retrieved T_sys_input: {T_sys_input}")
        print(f"Retrieved B_input: {B_input}")

        if not D_t_input or not z_input or not T_sys_input or not B_input:
            raise ValueError("Empty input for Step 6")

        D_t = float(D_t_input)
        z = float(z_input)
        T_sys = float(T_sys_input)
        B = float(B_input)

        # Step 6 calculations (transmitter output power)
        transmitter_power = calculate_step6()

        # Clear previous results if they exist
        canvas.delete("results")
        # Display Step 1 Results
        canvas.create_text(645.0, 810.0, anchor="nw", text=f"{R_min:.2f} m", fill="#000000", font=("Inter Medium", 15), tags="results")

        # Display Step 2 Result
        canvas.create_text(645, 840, anchor="nw", text=f"{delta_d:.2f} cm", fill="#000000", font=("Inter Medium", 15), tags="results")

        # Display Step 3 Results
        canvas.create_text(645, 870, anchor="nw", text=f"{t_int:.2f} s", fill="#000000", font=("Inter Medium", 15), tags="results")
        canvas.create_text(645, 900, anchor="nw", text=f"{t_map:.2f} hr", fill="#000000", font=("Inter Medium", 15), tags="results")

        # Display Step 4 Results
        canvas.create_text(860, 810, anchor="nw", text=f"{theta_ext:.2f} deg", fill="#000000", font=("Inter Medium", 15), tags="results")
        canvas.create_text(860, 840, anchor="nw", text=f"{theta_b:.2f} arcsec", fill="#000000", font=("Inter Medium", 15), tags="results")
        canvas.create_text(860, 870, anchor="nw", text=f"{theta_sr:.2f} arcsec", fill="#000000", font=("Inter Medium", 15), tags="results")
        canvas.create_text(860, 900, anchor="nw", text=f"{theta_ss:.2f} arcsec", fill="#000000", font=("Inter Medium", 15), tags="results")

        # Display Step 5 Results
        canvas.create_text(1080, 810, anchor="nw", text=f"{theta_point:.2f} deg", fill="#000000", font=("Inter Medium", 15), tags="results")
        canvas.create_text(1080, 840, anchor="nw", text=f"{N_row:.2f}", fill="#000000", font=("Inter Medium", 15), tags="results")
        canvas.create_text(1080, 870, anchor="nw", text=f"{SNR:.2f} dB", fill="#000000", font=("Inter Medium", 15), tags="results")

        # Display Step 6 Result (Transmitter Power)
        canvas.create_text(1080, 900, anchor="nw", text=f"{transmitter_power:.4e} W", fill="#000000", font=("Inter Medium", 15), tags="results")
        # Now, display the inputs at the bottom
        print(f"Starting to display inputs")

        display_inputs()
    
    except ValueError as e:
        messagebox.showerror("Calculation Error", str(e))


# Store references to the canvas text items so they can be deleted later
result_text_items = []
def calculate_and_display_step1():
    global result_text_items

    try:
        # Retrieve and process inputs from the correct entry boxes
        D_input = entry_2.get().strip()  # For Diameter (D)
        nu_input = entry_5.get().strip()  # For Frequency (ν)

        # Debug print to check if inputs are correctly retrieved
        print(f"D_input: {D_input}, nu_input: {nu_input}")

        if not D_input or not nu_input:
            raise ValueError("Empty input")

        D = float(D_input)
        nu = float(nu_input)

        # Perform the calculations
        R_F, R_min, R_react = calculate_step1(D, nu)

        # Clear previous results if they exist
        for item in result_text_items:
            canvas.delete(item)
        result_text_items.clear()

        # Display the new results to the right of the specified text boxes
        result_text_items.append(canvas.create_text(145.0, 210.0, anchor="nw", text=f"{R_F:.2f} m", fill="#000000", font=("Inter Medium", 13)))
        result_text_items.append(canvas.create_text(145.0, 230.0, anchor="nw", text=f"{R_min:.2f} m", fill="#000000", font=("Inter Medium", 13)))
        result_text_items.append(canvas.create_text(145.0, 250.0, anchor="nw", text=f"{R_react:.2f} m", fill="#000000", font=("Inter Medium", 13)))

    except ValueError:
        # Show a warning message if inputs are invalid
        messagebox.showwarning(
            "Input Error",
            "Please enter valid numerical values for D and ν (Frequency)."
        )

def open_Reference_Window():
    
    base_path = os.path.dirname(__file__)  

    # Construct the relative path to the PDF
    pdf_path = os.path.join(base_path, 'HoloSim_Quick_Reference.pdf')

    # Open the PDF in the default viewer, depending on the system
    if os.name == 'posix':  # For macOS and Linux
        os.system(f'open "{pdf_path}"')  # For macOS
        os.system(f'xdg-open "{pdf_path}"')  # For Linux

step2_text_items = []
def calculate_and_display_step2():
    global step2_text_items

    try:
        # Retrieve and process the input from the correct entry box
        a_input = entry_1.get().strip()  # For parameter a

        # Debug print to check if input is correctly retrieved
        print(f"a_input: {a_input}")

        if not a_input:
            raise ValueError("Empty input")

        a = float(a_input)

        # Perform the calculation for spatial resolution
        delta_d = calculate_step2(a)

        # Clear previous results if they exist
        for item in step2_text_items:
            canvas.delete(item)
        step2_text_items.clear()

        # Display the result on the canvas to the right of the "<" symbol
        # Adjusted output to print to the right of the text
        step2_text_items.append(canvas.create_text(590.0, 230.0, anchor="nw", text=f"{delta_d:.2f} cm", fill="#000000", font=("Inter Medium", 15)))


    except ValueError:
        # Show a warning message if the input is invalid
        messagebox.showwarning(
            "Input Error",
            "Please enter a valid numerical value for a."
        )



step3_text_items = []
def calculate_and_display_step3():
    global step3_text_items

    try:
        # Retrieve and process inputs from the correct entry boxes
        D_input = entry_2.get().strip()
        nu_input = entry_5.get().strip()

        if not D_input or not nu_input:
            raise ValueError("Empty input from previous step")

        # Get the correct delta_d from Step 2 (pass it directly)
        delta_d = calculate_step2(float(entry_1.get().strip()))  # Assuming entry_1 is for Step 2 input
        
        # Perform the calculations and pass delta_d into calculate_step3
        t_int, t_map = calculate_step3(delta_d)

        # Clear previous results if they exist
        for item in step3_text_items:
            canvas.delete(item)
        step3_text_items.clear()

        # Display the results on the canvas to the right of the "=" symbol
        step3_text_items.append(canvas.create_text(950.0, 240.0, anchor="nw", text=f"{t_int:.2f} s", fill="#000000", font=("Inter SemiBold", 13)))
        step3_text_items.append(canvas.create_text(950.0, 265.0, anchor="nw", text=f"{t_map:.2f} hr", fill="#000000", font=("Inter SemiBold", 13)))

    except ValueError:
        # Show a warning message if the input is invalid
        messagebox.showwarning(
            "Input Error",
            "Please enter valid numerical values for Step 3."
        )


step4_text_items = []
def calculate_and_display_step4():
    global step4_text_items

    try:
        # Retrieve and process inputs from the correct entry boxes
        nu_input = entry_5.get().strip()  # For Frequency (ν)
        D_input = entry_2.get().strip()  # For Diameter (D)
        f_oss_input = entry_12.get().strip()  # For f_oss
        
        # Check if any input is missing
        if not nu_input or not D_input or not f_oss_input:
            raise ValueError("Empty input")

        # Convert inputs to floats
        nu = float(nu_input)
        D = float(D_input)
        f_oss = float(f_oss_input)

        # Retrieve the delta_d value from Step 2 (which has already been calculated)
        delta_d = calculate_step2(float(entry_1.get().strip()))  # Assuming 'a' is in entry_1
        print(f"delta_d from Step 2: {delta_d}")

        # Perform the calculations using the updated calculate_step4 function
        theta_ext, theta_b, theta_sr, theta_ss = calculate_step4(delta_d)

        # Clear previous results if they exist
        for item in step4_text_items:
            canvas.delete(item)
        step4_text_items.clear()

        
        step4_text_items.append(canvas.create_text(130.0, 520.0, anchor="nw", text=f"{theta_ext:.2f} deg", fill="#000000", font=("Inter Medium", 14)))
        step4_text_items.append(canvas.create_text(130.0, 543.0, anchor="nw", text=f"{theta_b:.2f} arcsec", fill="#000000", font=("Inter Medium", 14)))
        step4_text_items.append(canvas.create_text(130.0, 568.0, anchor="nw", text=f"{theta_sr:.2f} arcsec", fill="#000000", font=("Inter Medium", 14)))
        step4_text_items.append(canvas.create_text(130.0, 590.0, anchor="nw", text=f"{theta_ss:.2f} arcsec", fill="#000000", font=("Inter Medium", 14)))

    except ValueError:
        # Show a warning message if inputs are invalid
        messagebox.showwarning(
            "Input Error",
            "Please enter valid numerical values for frequency (ν), D, and f_oss in Step 4."
        )


step5_text_items = []
def calculate_and_display_step5():
    global step5_text_items

    try:
        # Retrieve and process inputs from the correct entry box
        delta_z_input = entry_16.get().strip()  # For surface deformation (delta_z)

        if not delta_z_input:
            raise ValueError("Empty input")

        # Retrieve delta_d from Step 2
        delta_d = calculate_step2(float(entry_1.get().strip()))  # Pull delta_d from Step 2

        # Perform the calculations using the updated calculate_step5
        theta_point, N_row, SNR = calculate_step5(delta_d)

        # Clear previous results if they exist
        for item in step5_text_items:
            canvas.delete(item)
        step5_text_items.clear()

        # Display the new results and store the item IDs in the list
        step5_text_items.append(canvas.create_text(555.0, 528.0, anchor="nw", text=f"{theta_point:.2f} deg", fill="#000000", font=("Inter Medium", 15)))
        step5_text_items.append(canvas.create_text(555.0, 552.0, anchor="nw", text=f"{N_row:.2f}", fill="#000000", font=("Inter Medium", 15)))
        step5_text_items.append(canvas.create_text(555.0, 575.0, anchor="nw", text=f"{SNR:.2f} dB", fill="#000000", font=("Inter Medium", 15)))

    except ValueError:
        # Show a warning message if inputs are invalid
        messagebox.showwarning(
            "Input Error",
            "Please enter a valid numerical value for delta_z in Step 5."
        )


step6_text_items = []
def calculate_and_display_step6():
    global step6_text_items

    try:
        # Perform the Step 6 calculation using the updated calculate_step6 function
        result = calculate_step6()

        # Ensure result is valid
        if result is None:
            return

        # Clear previous results if they exist
        for item in step6_text_items:
            canvas.delete(item)
        step6_text_items.clear()

        # Display the result to the right of the "Transmitter output power = " text
        step6_text_items.append(
            canvas.create_text(
                960.0, 600.0, anchor="nw", text=f"{result:.4e} W", fill="#000000", font=("Inter Medium", 14)
            )
        )

    except ValueError:
        messagebox.showwarning(
            "Input Error",
            "Please enter valid numerical values for all required inputs in Step 6, and ensure 'D' from Step 1 is correctly entered."
        )

#This file will be HoloSim final v1 w/o any graphs- HoloSimModern has the graphs.  

window = Tk()

icon_image = PhotoImage(file=relative_to_assets("HS_logo.png"))
window.iconphoto(False, icon_image)

window.geometry("1200x940")
window.configure(bg="#808080")
window.title("HoloSim")

canvas = Canvas(
    window,
    bg="#808080",
    height=940,
    width=1200,
    bd=0,
    highlightthickness=0,
    relief="ridge"
)

canvas.place(x=0, y=0)
canvas.create_text(
    13.0,
    15.0,
    anchor="nw",
    text="HoloSim.",
    fill="#FFFFFF",
    font=("Inter", 42, "bold")
)

def load_and_display_image(canvas, image_path, x, y):
    image = PhotoImage(file=relative_to_assets(image_path))
    canvas.create_image(x, y, image=image)
    return image

def create_entry(canvas, image_path, x, y, width, height, default_text=""):
    entry_image = PhotoImage(file=relative_to_assets(image_path))
    entry_bg = canvas.create_image(x, y, image=entry_image)
    entry = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
    entry.place(x=x - width / 2, y=y - height / 2, width=width, height=height)
    entry.insert(0, default_text)
    return entry

def create_button(image_path, x, y, width, height, command):
    button_image = PhotoImage(file=relative_to_assets(image_path))
    button = Button(image=button_image, borderwidth=0, highlightthickness=0, command=command, relief="flat")
    button.place(x=x, y=y, width=width, height=height)
    return button

def open_help_window():
   
    base_path = os.path.dirname(__file__)  

    # Construct the relative path to the PDF
    pdf_path = os.path.join(base_path, 'Submillimeter_Wave_Holography_for_Large_Dish_Antennas.pdf')

    # Open the PDF in the default viewer, depending on the system
    if os.name == 'posix':  # For macOS and Linux
        os.system(f'open "{pdf_path}"')  # For macOS
        os.system(f'xdg-open "{pdf_path}"')  # For Linux
    else:
        messagebox.showerror("Error", "Unsupported operating system.")


def display_inputs():
    
    try:
        # Ensure the canvas is cleared before drawing new text
        canvas.delete("inputs")

        # Retrieve and display each entry's value
        inputs = [
            {"value": entry_2.get().strip(), "unit": "m", "x": 60.0, "y": 810.0},
            {"value": entry_5.get().strip(), "unit": "GHz", "x": 60.0, "y": 840.0},
            {"value": entry_1.get().strip(), "unit": "m", "x": 60.0, "y": 870.0},
            {"value": entry_3.get().strip(), "unit": "", "x": 60.0, "y": 905.0},
            {"value": entry_6.get().strip(), "unit": "", "x": 210.0, "y": 810.0},
            {"value": entry_7.get().strip(), "unit": "", "x": 210.0, "y": 840.0},
            {"value": entry_19.get().strip(), "unit": "sec", "x": 210.0, "y": 875.0},
            {"value": entry_12.get().strip(), "unit": "", "x": 210.0, "y": 903.0},
            {"value": entry_16.get().strip(), "unit": "μm", "x": 360.0, "y": 810.0},
            {"value": entry_17.get().strip(), "unit": "m", "x": 360.0, "y": 840.0},
            {"value": entry_8.get().strip(), "unit": "m", "x": 360.0, "y": 870.0},
            {"value": entry_9.get().strip(), "unit": "K", "x": 360.0, "y": 900.0},
            {"value": entry_10.get().strip(), "unit": "MHz", "x": 475.0, "y": 810.0},
        ]

        # Print statements for debugging
        print("Displaying inputs:")
        for item in inputs:
            value_str = item["value"]
            unit = item["unit"]
            x = item["x"]
            y = item["y"]

            # Check if the entry has a value
            if value_str == "":
                print(f"Warning: Empty input detected at coordinates ({x}, {y})")

            # Prepare the text to display (just the value and unit)
            text = f"{value_str} {unit}".strip()

            # Display the text on the canvas
            canvas.create_text(x, y, anchor="nw", text=text, fill="#000000", font=("Inter Medium", 15), tags="inputs")
            print(f"Displayed: '{text}' at ({x}, {y})")

    except Exception as e:
        print(f"Error displaying inputs: {e}")



# images used to make each box for each step and seperating lines
image_1 = load_and_display_image(canvas, "Image_1.png", 200.0, 230.0)
image_2 = load_and_display_image(canvas, "Image_2.png", 600.0, 230.0)
image_3 = load_and_display_image(canvas, "Image_3.png", 1000.0, 230.0)
image_4 = load_and_display_image(canvas, "Image_4.png", 200.0, 560.0)
image_5 = load_and_display_image(canvas, "Image_5.png", 600.0, 560.0)
image_6 = load_and_display_image(canvas, "Image_6.png", 1000.0, 560.0)
image_Line_1 = load_and_display_image(canvas, "Line_1.png", 550.0, 920.0)
image_Line_2 = load_and_display_image(canvas, "Line_2.png", 600.0, 790.0)

# Bottom rectangle used for "Calculate all" portion of the GUI
rect = canvas.create_rectangle(0.0, 725.0, 1200.0, 976.0, fill="#B6B9C9", outline="")
canvas.tag_lower(rect)

# Step 1 latex
image_D = load_and_display_image(canvas, "D.png", 185.0, 159.0)
image_nu = load_and_display_image(canvas, "nu.png", 205.0, 180.0)
image_R_F = load_and_display_image(canvas, "R_F.png", 100.0, 219.0)
image_R_min = load_and_display_image(canvas, "R_min.png", 105.0, 239.0)
image_R_react = load_and_display_image(canvas, "R_react.png", 107.0, 259.0)

# Step 2 latex
image_a = load_and_display_image(canvas, "a.png", 654.0, 159.0)
image_small_delt_d = load_and_display_image(canvas, "small_delt_d.png", 550.0, 239.0)

# Step 3 latex
image_f_1 = load_and_display_image(canvas, "f_1.png", 1008.0, 158.0)
image_f_apo = load_and_display_image(canvas, "f_apo.png", 1037.0, 178.0)
image_f_osr = load_and_display_image(canvas, "f_osr.png", 980.0, 203.0)
image_Theta_dot = load_and_display_image(canvas, "Theta_dot.png", 1018.0, 220.0)
image_t_int = load_and_display_image(canvas, "t_int.png", 918.0, 250.0)
image_t_map = load_and_display_image(canvas, "t_map.png", 920.0, 275.0)

# Step 4 latex
image_f_oss = load_and_display_image(canvas, "f_oss.png", 234.0, 490.0)
image_step4_outputs = load_and_display_image(canvas, "step4_outputs.png", 100.0, 560.0)

# Step 5 latex
image_small_delta_z = load_and_display_image(canvas, "small_delta_z.png", 621.0, 490.0)
image_step_5_output = load_and_display_image(canvas, "step_5_output.png", 520.0, 560.0)

# Step 6 latex
image_D_t = load_and_display_image(canvas, "D_t.png", 1000.0, 491.0)
image_T_sys = load_and_display_image(canvas, "T_sys.png", 997.0, 530.0)
image_T_R_z = load_and_display_image(canvas, "T_R_z.png", 1000.0, 509.0)
image_B = load_and_display_image(canvas, "B.png", 990.0, 551.0)

# Calculate all portion's latex

#inputs
image_Input_D = load_and_display_image(canvas, "Input_D.png", 30.0, 820.0)
image_Input_nu = load_and_display_image(canvas, "Input_nu.png", 30.0, 850.0)
image_Input_a = load_and_display_image(canvas, "Input_a.png", 30.0, 880.0)
image_Input_f_1 = load_and_display_image(canvas, "Input_f_1.png", 30.0, 910.0)
image_Input_f_apo = load_and_display_image(canvas, "Input_f_apo.png", 180.0, 820.0)
image_Input_f_osr = load_and_display_image(canvas, "Input_f_osr.png", 180.0, 850.0)
image_Input_theta_dot = load_and_display_image(canvas, "Input_theta_dot.png", 170.0, 880.0)
image_Input_f_oss = load_and_display_image(canvas, "Input_f_oss.png", 180.0, 910.0)
image_Input_delta_z = load_and_display_image(canvas, "Input_delta_z.png", 330.0, 820.0)
image_Input_D_t = load_and_display_image(canvas, "Input_D_t.png", 330.0, 850.0)
image_Input_z = load_and_display_image(canvas, "Input_z.png", 324.0, 880.0)
image_Input_T_sys = load_and_display_image(canvas, "Input_T_sys.png", 330.0, 910.0)
image_Input_B = load_and_display_image(canvas, "Input_B.png", 445.0, 820.0)

canvas.create_text(10.0, 730.0, anchor="nw", text="Final Parameters", fill="#000000", font=("Inter", 20, "bold"))
canvas.create_text(15.0, 765.0, anchor="nw", text="Inputs:", fill="#000000", font=("Inter Medium", 17))
canvas.create_text(570.0, 765.0, anchor="nw", text="Outputs:", fill="#000000", font=("Inter Medium", 17))


#outputs 
image_R_min_Equals = load_and_display_image(canvas, "R_min_Equals.png", 607.0, 820.0)
image_Delta_lower_d_less = load_and_display_image(canvas, "Delta_lower_d_less.png", 595.0, 850.0)

image_t_int_equals = load_and_display_image(canvas, "t_int_equals.png", 596.0, 880.0)
image_t_map_equals = load_and_display_image(canvas, "t_map_equals.png", 601.0, 910.0)
image_Theta_ext_Equals = load_and_display_image(canvas, "Theta_ext_Equals.png", 830.0, 820.0)
image_Theta_b_equals = load_and_display_image(canvas, "Theta_b_equals.png", 823.0, 850.0)

image_Theta_sr_equals = load_and_display_image(canvas, "Theta_sr_equals.png", 825.0, 880.0)
image_Theta_ss_equals = load_and_display_image(canvas, "Theta_ss_equals.png", 825.0, 910.0)
image_Theta_point_less = load_and_display_image(canvas, "Theta_point_less.png", 1040.0, 820.0)
image_N_row_equals = load_and_display_image(canvas, "N_row_equals.png", 1035.0, 850.0)

image_SNR = load_and_display_image(canvas, "SNR.png", 1031.0, 880.0)

canvas.create_text(1007.0, 895.0, anchor="nw", text="P =", fill="#000000", font=("Inter Medium", 20, "italic"))

def create_text_elements(canvas):
    # Steps 1-6 text
    canvas.create_text(55.0, 100.0, anchor="nw", text="Step 1: Transmitter Placement", fill="#000000", font=("Inter Medium", 20))
    canvas.create_text(50.0, 150.0, anchor="nw", text="Aperture Diameter  \n\n", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(315, 150, anchor="nw", text="m", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(50.0, 170.0, anchor="nw", text="Transmitter frequency  ", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(315, 170, anchor="nw", text="GHz", fill="#000000", font=("Inter SemiBold", 14))

    canvas.create_text(485.0, 100.0, anchor="nw", text="Step 2: Spatial Resolution   ", fill="#000000", font=("Inter Medium", 20))
    canvas.create_text(420.0, 150.0, anchor="nw", text="Distance between corner adjustors    :", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(770, 150, anchor="nw", text="m", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(570.0, 230.0, anchor="nw", text="< ", fill="#000000", font=("Inter Medium", 15))

    canvas.create_text(830.0, 85.0, anchor="nw", text="Step 3: Grid Point Integration Time and", fill="#000000", font=("Inter Medium", 20))
    canvas.create_text(930, 110, anchor="nw", text="Total Map Time", fill="#000000", font=("Inter Medium", 20))
    canvas.create_text(840.0, 150.0, anchor="nw", text="Primary beam taper factor      : ", fill="#000000", font=("Inter SemiBold", 13))
    canvas.create_text(840.0, 170.0, anchor="nw", text="Apodization smoothing factor          :", fill="#000000", font=("Inter SemiBold", 13))
    canvas.create_text(840.0, 193.0, anchor="nw", text="Oversampling factor         :", fill="#000000", font=("Inter SemiBold", 13))
    canvas.create_text(840.0, 215.0, anchor="nw", text="Rotation rate of the antenna     :", fill="#000000", font=("Inter SemiBold", 13))
    canvas.create_text(1160, 215.0, anchor="nw", text="sec", fill="#000000", font=("Inter SemiBold", 13))

    canvas.create_text(45.0, 410.0, anchor="nw", text="Step 4: Map Angular Extent and", fill="#000000", font=("Inter Medium", 20))
    canvas.create_text(105, 435, anchor="nw", text="Sampling Intervals", fill="#000000", font=("Inter Medium", 20))
    canvas.create_text(20.0, 480.0, anchor="nw", text="Oversampling factor along row        :", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(355, 480, anchor="nw", text="sec", fill="#000000", font=("Inter SemiBold", 14))

    canvas.create_text(460.0, 410.0, anchor="nw", text="Step 5: Pointing Accuracy and", fill="#000000", font=("Inter Medium", 20))
    canvas.create_text(520, 435, anchor = 'nw', text= "SNR Requirement", fill="#000000", font=("Inter Medium", 20))
    canvas.create_text(463.0, 480.0, anchor="nw", text="Surface deformation       :", fill="#000000", font=("Inter Medium", 15)) 
    canvas.create_text(740, 480, anchor="nw", text="μm", fill="#000000", font=("Inter Medium", 15))

    canvas.create_text(850.0, 430.0, anchor="nw", text="Step 6: Transmitter output power", fill="#000000", font=("Inter Medium", 20))
    canvas.create_text(850.0, 482.0, anchor="nw", text="Transmitter diameter      : ", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(1140, 482.0, anchor="nw", text="m", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(850.0, 500.0, anchor="nw", text="Distance between", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(1140, 500.0, anchor="nw", text="m", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(1030.0, 500.0, anchor="nw", text=":", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(850.0, 520.0, anchor="nw", text="System temperature        : ", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(1140, 520.0, anchor="nw", text="K", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(847.0, 542.0, anchor="nw", text=" Detector Bandwidth     :", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(1140, 542.0, anchor="nw", text="MHz", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(930.0, 600.0, anchor="nw", text="P", fill="#000000", font=("Inter SemiBold", 15, "italic"))
    canvas.create_text(940.0, 600.0, anchor="nw", text=" = ", fill="#000000", font=("Inter SemiBold", 15))


# Call the function to create text elements
create_text_elements(canvas)

entry_image_1 = load_and_display_image(canvas, "entry_1.png", 715.0, 160.0)
entry_1 = create_entry(canvas, "entry_1.png", 715.0, 160.0, 85.0, 15.0)

entry_image_2 = load_and_display_image(canvas, "entry_2.png", 265.0, 158.0)
entry_2 = create_entry(canvas, "entry_2.png", 265.0, 158.0, 78.0, 16.0)

entry_image_3 = load_and_display_image(canvas, "entry_3.png", 1110.0, 156.5)
entry_3 = create_entry(canvas, "entry_3.png", 1110.0, 156.5, 77.0, 17.0)

entry_image_5 = load_and_display_image(canvas, "entry_5.png", 265.0, 180.0)
entry_5 = create_entry(canvas, "entry_5.png", 265.0, 180.0, 77.0, 17.0)

entry_image_6 = load_and_display_image(canvas, "entry_6.png", 1110.0, 180.0)
entry_6 = create_entry(canvas, "entry_6.png", 1110.0, 180.0, 78.0, 16.0)
entry_6.insert(0, "1.3")  # Sets the default value to 1.3 on startup

entry_image_7 = load_and_display_image(canvas, "entry_7.png", 1110.0, 202.0)
entry_7 = create_entry(canvas, "entry_7.png", 1110.0, 202.0, 77.0, 17.0)

entry_image_8 = load_and_display_image(canvas, "entry_8.png", 1090.0, 510.5)
entry_8 = create_entry(canvas, "entry_8.png", 1090.0, 510.5, 77.0, 17.0)

entry_image_9 = load_and_display_image(canvas, "entry_9.png", 1090.0, 531.5)
entry_9 = create_entry(canvas, "entry_9.png", 1090.0, 531.5, 77.0, 17.0)

entry_image_10 = load_and_display_image(canvas, "entry_10.png", 1090.0, 552.5)
entry_10 = create_entry(canvas, "entry_10.png", 1090.0, 552.5, 77.0, 17.0)

entry_image_12 = load_and_display_image(canvas, "entry_12.png", 305.0, 489.0)
entry_12 = create_entry(canvas, "entry_12.png", 305.0, 489.0, 77.0, 17.0)

entry_image_16 = load_and_display_image(canvas, "entry_16.png", 690.0, 490.0)
entry_16 = create_entry(canvas, "entry_16.png", 690.0, 490.0, 77.0, 17.0)

entry_image_17 = load_and_display_image(canvas, "entry_17.png", 1090.0, 490.0)
entry_17 = create_entry(canvas, "entry_17.png", 1090.0, 490.0, 77.0, 17.0)

entry_image_19 = load_and_display_image(canvas, "entry_19.png", 1110.0, 224.0)
entry_19 = create_entry(canvas, "entry_19.png", 1110.0, 224.0, 78.0, 16.0)


button_image_1 = PhotoImage(file=relative_to_assets("button1.png"))
button_1 = Button(image=button_image_1, borderwidth=0, highlightthickness=0, command=calculate_all_steps, relief="flat")
button_1.place(x=1040, y=730, width=148, height=46)

button_image_2 = PhotoImage(file=relative_to_assets("button_2.png")) # Connect the updated handler
button_2 = Button(image=button_image_2, borderwidth=0, highlightthickness=0, command=calculate_and_display_step1, relief="flat")
button_2.place(x=110.0, y=300.0, width=135, height=50)

button_image_3 = PhotoImage(file=relative_to_assets("button_3.png")) # Connect Button 3 to the calculate_and_display_step4 function
button_3 = Button(image=button_image_3, borderwidth=0, highlightthickness=0, command=calculate_and_display_step4, relief="flat")
button_3.place(x=110, y=630, width=135, height=50)

button_image_4 = PhotoImage(file=relative_to_assets("button_4.png"))
button_4 = Button(image=button_image_4, borderwidth=0, highlightthickness=0, command=calculate_and_display_step5, relief="flat")
button_4.place(x=530, y=630, width=135, height=50)

button_image_5 = PhotoImage(file=relative_to_assets("button_5.png"))
button_5 = Button(image=button_image_5, borderwidth=0, highlightthickness=0, command=calculate_and_display_step2,  relief="flat")
button_5.place(x=530.0, y=300.0, width=135, height=50)

button_image_6 = PhotoImage(file=relative_to_assets("button_6.png"))
button_6 = Button(image=button_image_6,borderwidth=0,highlightthickness=0,command=calculate_and_display_step3, relief="flat") # Connect Button 6 to the calculate_and_display_step3 function
button_6.place(x=940.0, y=300, width=135, height=50)

button_image_7 = PhotoImage(file=relative_to_assets("button_7.png"))
button_7 = Button(image=button_image_7, borderwidth=0, highlightthickness=0, command=calculate_and_display_step6,  relief="flat") # Connect to Step 6 function
button_7.place(x=940.0, y=630.0, width=135, height=50)

Guide_Button_image = PhotoImage(file=relative_to_assets("Massive_Help_Button.png"))
Guide_Button = Button(image=Guide_Button_image, borderwidth=0, highlightthickness=0, command=open_help_window, relief="flat")
Guide_Button.place(x=1050.0, y=20.0, width=141, height=36)

button_image_14 = PhotoImage(file=relative_to_assets("button_14.png")) # help button 1 
button_14 = Button(image=button_image_14,borderwidth=0,highlightthickness=0,command=open_Reference_Window, relief="flat")
button_14.place(x=350, y=350, width=30, height=32)

button_image_15 = PhotoImage(file=relative_to_assets("button_15.png"))
button_15 = Button(image=button_image_15, borderwidth=0, highlightthickness=0, command=open_Reference_Window, relief="flat")
button_15.place(x=1150.3448486328125, y=680.2237548828125, width=30, height=32)

button_image_16 = PhotoImage(file=relative_to_assets("button_16.png"))
button_16 = Button(image=button_image_14, borderwidth = 0,highlightthickness=0, command=open_Reference_Window, relief= "flat")
button_16.place(x=1150.3448486328125, y=350, width=30, height=32)

button_image_17 = PhotoImage(file=relative_to_assets("button_17.png"))
button_17 = Button(image=button_image_17, borderwidth=0, highlightthickness=0, command=open_Reference_Window, relief="flat")
button_17.place(x=740, y=680, width=30, height=32)

button_image_18 = PhotoImage(file=relative_to_assets("button_18.png"))
button_18 = Button(image=button_image_18, borderwidth=0, highlightthickness=0, command=open_Reference_Window, relief="flat")
button_18.place(x=350, y=680, width=30, height=32)

button_image_19 = PhotoImage(file=relative_to_assets("button_19.png"))
button_19 = Button(image=button_image_19, borderwidth=0, highlightthickness=0, command=open_Reference_Window, relief="flat")
button_19.place(x=740, y=350, width=30, height=32)
window.resizable(False, False)

window.mainloop()