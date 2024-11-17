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

def calculate_step5_with_delta_z(delta_d, delta_z):
    # with the provided delta_z parameter in this specific function
    try:
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

        # Perform the Step 5 calculations with varying delta_z
        theta_point = 6 * 10**-4 * delta_z * nu  # Pointing accuracy
        N_row = D / (delta_d * 10**-2) * f_apo * f_osr  # Total number of rows
        SNR = 10 * np.log10(0.044 * lam * np.sqrt(N_row**2 / (f_oss * f_osr)) * 1 / delta_z * 1 / np.sqrt(f_apo) * 10**6)

        return theta_point, N_row, SNR

    except ValueError:
        return None, None, None

def open_graph_window(step, data_dict, x_label_text, y_label_text):
    graph_window = Toplevel()
    graph_window.title(f"Graph for Step {step}")
    graph_window.geometry("900x600")
    graph_window.configure(bg="#ecf0f1")

    # Graph frame
    graph_frame = tk.Frame(graph_window, width=450, height=600, bg="#ecf0f1")
    graph_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=10, pady=10)

    # Control frame
    control_frame = tk.Frame(graph_window, width=450, height=600, bg="#ecf0f1")
    control_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)

    # Set up Matplotlib figure and axes
    fig = Figure(figsize=(5, 5), dpi=100)
    ax = fig.add_subplot(111)

    # Initial plot with default X and Y axis selections
    ax.scatter(data_dict[x_label_text], data_dict[y_label_text], color='blue', label=f'Step {step} Data', s=10)
    ax.set_xlabel(x_label_text)
    ax.set_ylabel(y_label_text)
    ax.legend()

    # Canvas for displaying the graph
    canvas = FigureCanvasTkAgg(fig, master=graph_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # Dropdown for X-axis selection
    x_var = tk.StringVar(value=x_label_text)
    x_label = tk.Label(control_frame, text="Select X-axis:", font=("Arial", 12), bg="#ecf0f1", fg="#2c3e50")
    x_label.pack(pady=10)
    x_menu = ttk.Combobox(control_frame, textvariable=x_var, values=list(data_dict.keys()))
    x_menu.pack(pady=5)

    # Dropdown for Y-axis selection
    y_var = tk.StringVar(value=y_label_text)
    y_label = tk.Label(control_frame, text="Select Y-axis:", font=("Arial", 12), bg="#ecf0f1", fg="#2c3e50")
    y_label.pack(pady=10)
    y_menu = ttk.Combobox(control_frame, textvariable=y_var, values=list(data_dict.keys()))
    y_menu.pack(pady=5)

    # Input fields for domain (X-axis limits)
    domain_frame = tk.LabelFrame(control_frame, text="Domain (X-axis limits)", bg="#ecf0f1", fg="#2c3e50", padx=10, pady=10)
    domain_frame.pack(pady=10, fill=tk.X)

    x_min_label = tk.Label(domain_frame, text="X Min:", bg="#ecf0f1", fg="#2c3e50")
    x_min_label.pack(side=tk.LEFT, padx=5)
    x_min_entry = tk.Entry(domain_frame, width=10)  # Set width
    x_min_entry.pack(side=tk.LEFT, padx=5)

    x_max_label = tk.Label(domain_frame, text="X Max:", bg="#ecf0f1", fg="#2c3e50")
    x_max_label.pack(side=tk.LEFT, padx=5)
    x_max_entry = tk.Entry(domain_frame, width=10)  # Set width
    x_max_entry.pack(side=tk.LEFT, padx=5)

    # Input fields for range (Y-axis limits)
    range_frame = tk.LabelFrame(control_frame, text="Range (Y-axis limits)", bg="#ecf0f1", fg="#2c3e50", padx=10, pady=10)
    range_frame.pack(pady=10, fill=tk.X)

    y_min_label = tk.Label(range_frame, text="Y Min:", bg="#ecf0f1", fg="#2c3e50")
    y_min_label.pack(side=tk.LEFT, padx=5)
    y_min_entry = tk.Entry(range_frame, width=10)  # Set width
    y_min_entry.pack(side=tk.LEFT, padx=5)

    y_max_label = tk.Label(range_frame, text="Y Max:", bg="#ecf0f1", fg="#2c3e50")
    y_max_label.pack(side=tk.LEFT, padx=5)
    y_max_entry = tk.Entry(range_frame, width=10)  # Set width
    y_max_entry.pack(side=tk.LEFT, padx=5)

    # Save Graph Button
    save_button = tk.Button(control_frame, text="Save Graph", command=lambda: save_graph_as_image(fig))
    save_button.pack(pady=20)

    # Apply domain and range button
    apply_button = tk.Button(control_frame, text="Apply Domain & Range", command=lambda: update_domain_and_range())
    apply_button.pack(pady=10)

    # Update plot function to apply new selections
    def update_plot(*args):
        current_x_label = x_var.get()
        current_y_label = y_var.get()

        # Clear and replot with new data
        ax.clear()
        ax.scatter(data_dict[current_x_label], data_dict[current_y_label], color='blue', label=f'{current_x_label} vs {current_y_label}', s=10)
        ax.set_xlabel(current_x_label)
        ax.set_ylabel(current_y_label)
        ax.legend()
        canvas.draw()

    # Function to update domain and range
    def update_domain_and_range():
        try:
            # Get domain limits (X-axis)
            x_min = float(x_min_entry.get()) if x_min_entry.get() else None
            x_max = float(x_max_entry.get()) if x_max_entry.get() else None

            # Get range limits (Y-axis)
            y_min = float(y_min_entry.get()) if y_min_entry.get() else None
            y_max = float(y_max_entry.get()) if y_max_entry.get() else None

            # Ensure that x_min < x_max and y_min < y_max
            if x_min is not None and x_max is not None and x_min >= x_max:
                messagebox.showwarning("Invalid Domain", "X Min should be less than X Max.")
                return
            if y_min is not None and y_max is not None and y_min >= y_max:
                messagebox.showwarning("Invalid Range", "Y Min should be less than Y Max.")
                return

            # Apply the new limits
            if x_min is not None and x_max is not None:
                ax.set_xlim([x_min, x_max])
            if y_min is not None and y_max is not None:
                ax.set_ylim([y_min, y_max])

            canvas.draw()

        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter valid numerical values for domain and range.")

    # Bind the dropdowns to update the plot dynamically
    x_var.trace("w", update_plot)
    y_var.trace("w", update_plot)

    graph_window.mainloop()

def save_graph_as_image(fig):
    # Ask the user where they want to save the image
    file_path = filedialog.asksaveasfilename(defaultextension='.png', filetypes=[("PNG files", "*.png"), ("JPEG files", "*.jpg"), ("All files", "*.*")])
    if file_path:
        fig.savefig(file_path)
        messagebox.showinfo("Save Successful", f"Graph saved as {file_path}")

def plot_step1_graph():
    try:
        # Retrieve inputs from the GUI
        D = float(entry_2.get())
        nu = float(entry_5.get())

        # Generate a range of values for plotting
        D_values = np.linspace(D * 0.8, D * 1.2, 300)  # Vary D by ±20%
        nu_values = np.linspace(nu * 0.8, nu * 1.2, 300)  # Vary nu by ±20%
        R_F_values = []
        R_min_values = []
        R_react_values = []

        for D_val in D_values:
            R_F, R_min, R_react = calculate_step1(D_val, nu)
            R_F_values.append(R_F)
            R_min_values.append(R_min)
            R_react_values.append(R_react)

        # Create a new window for plotting
        graph_window = Toplevel()
        graph_window.title("Step 1 Graph")
        graph_window.geometry("900x600")

        # Create a figure for plotting
        fig = Figure(figsize=(6, 4), dpi=100)
        ax = fig.add_subplot(111)

        # Plot the default Y axis (Far-field Cutoff) with Diameter on the X-axis
        ax.scatter(D_values, R_F_values, label="R_F (m)", color='blue', s=10)
        ax.set_xlabel("D (m)")  # Label with variable and unit only
        ax.set_ylabel("R_F (m)")  # Label with variable and unit only
        ax.set_title("Step 1: Transmitter Placement")
        ax.legend()

        # Canvas for displaying the graph
        canvas = FigureCanvasTkAgg(fig, master=graph_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill='both', expand=True)

        # Control frame for selecting X and Y axes
        control_frame = LabelFrame(graph_window, text="Select X and Y Axis", padx=10, pady=10)
        control_frame.pack(fill='x', padx=10, pady=10)

        # Options for X and Y axes, including all calculated and retrieved values
        axis_options = ["Diameter (D)", "Frequency (nu)", "Far-field Cutoff (R_F)", "Minimum Distance (R_min)", "Reactive Near-field Cutoff (R_react)"]

        # Dropdown for X-axis selection
        x_var = StringVar(value="Diameter (D)")
        x_label = Label(control_frame, text="Select X-axis:")
        x_label.pack(side='left', padx=5)
        x_menu = ttk.Combobox(control_frame, textvariable=x_var, values=axis_options)
        x_menu.pack(side='left', padx=5)

        # Dropdown for Y-axis selection
        y_var = StringVar(value="Far-field Cutoff (R_F)")
        y_label = Label(control_frame, text="Select Y-axis:")
        y_label.pack(side='left', padx=5)
        y_menu = ttk.Combobox(control_frame, textvariable=y_var, values=axis_options)
        y_menu.pack(side='left', padx=5)

        # Entry fields for setting X-axis min and max
        x_min_label = Label(control_frame, text="X-axis Min:")
        x_min_label.pack(side='left', padx=5)
        x_min_entry = Entry(control_frame, width=5)
        x_min_entry.pack(side='left', padx=5)

        x_max_label = Label(control_frame, text="X-axis Max:")
        x_max_label.pack(side='left', padx=5)
        x_max_entry = Entry(control_frame, width=5)
        x_max_entry.pack(side='left', padx=5)

        # Save original limits for reset
        original_x_min = min(D_values)
        original_x_max = max(D_values)

        def update_plot(*args):
            current_x_label = x_var.get()
            current_y_label = y_var.get()

            # Retrieve new X-axis limits
            try:
                x_min = float(x_min_entry.get())
                x_max = float(x_max_entry.get())
            except ValueError:
                x_min = original_x_min
                x_max = original_x_max

            # Recalculate X and Y values for the new range
            x_values = np.linspace(x_min, x_max, 300)
            y_values = []

            if current_x_label == "Diameter (D)":
                for x_val in x_values:
                    if current_y_label == "Far-field Cutoff (R_F)":
                        y_val, _, _ = calculate_step1(x_val, nu)
                    elif current_y_label == "Minimum Distance (R_min)":
                        _, y_val, _ = calculate_step1(x_val, nu)
                    elif current_y_label == "Reactive Near-field Cutoff (R_react)":
                        _, _, y_val = calculate_step1(x_val, nu)
                    else:
                        y_val = x_val  # Default case for matching X to Y
                    y_values.append(y_val)
                ax.set_xlabel("D (m)")
            elif current_x_label == "Frequency (nu)":
                x_values = nu_values
                y_values = R_F_values if current_y_label == "Far-field Cutoff (R_F)" else R_min_values
                ax.set_xlabel("nu (GHz)")

            # Clear and replot the graph with updated values
            ax.clear()
            ax.scatter(x_values, y_values, label=f"{current_y_label} vs {current_x_label}", color='blue', s=10)
            ax.set_ylabel(f"{current_y_label} (unit)")  # Update Y-axis label dynamically
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(min(y_values), max(y_values))
            ax.set_title("Step 1: Transmitter Placement")
            ax.legend()
            canvas.draw()

        # Button to trigger plot update
        update_button = Button(control_frame, text="Update X-axis Limits", command=update_plot)
        update_button.pack(side='left', padx=10)

        # Reset button to restore original limits
        def reset_limits():
            x_min_entry.delete(0, 'end')
            x_max_entry.delete(0, 'end')
            ax.set_xlim(original_x_min, original_x_max)
            canvas.draw()

        reset_button = Button(control_frame, text="Reset X-axis Limits", command=reset_limits)
        reset_button.pack(side='left', padx=10)

        # Trigger plot update when X or Y axis selection changes
        x_var.trace("w", update_plot)
        y_var.trace("w", update_plot)
        x_min_entry.bind("<Return>", lambda event: update_plot())
        x_max_entry.bind("<Return>", lambda event: update_plot())

        # Save button to save the graph as an image
        save_button = Button(control_frame, text="Save Graph", command=lambda: save_graph_as_image(fig))
        save_button.pack(pady=10)

    except ValueError:
        messagebox.showwarning("Input Error", "Please enter valid numbers for Step 1 inputs.")

def plot_step2_graph():
    try:
        # Retrieve inputs from the GUI
        a = float(entry_1.get())  # Aperture
        nu = float(entry_5.get())  # Frequency
        D = float(entry_2.get())  # Diameter

        # Generate initial ranges of values for plotting
        a_values = np.linspace(a * 0.8, a * 1.2, 300)  # Vary a by ±20%
        nu_values = [nu] * 300
        D_values = np.linspace(D * 0.8, D * 1.2, 300)  # Vary D by ±20%

        # Initialize lists for spatial resolution and other calculations
        delta_d_values = []
        R_F_values = []
        R_min_values = []
        R_react_values = []

        # Precompute data for the initial range
        for a_val in a_values:
            R_F, R_min, R_react = calculate_step1(D, nu)  # Using constant D and nu
            delta_d = calculate_step2(a_val)  # Calculate spatial resolution
            delta_d_values.append(delta_d)
            R_F_values.append(R_F)
            R_min_values.append(R_min)
            R_react_values.append(R_react)

        # Create a new window for plotting
        graph_window = Toplevel()
        graph_window.title("Step 2 Graph")
        graph_window.geometry("900x600")

        # Create a figure for plotting
        fig = Figure(figsize=(6, 4), dpi=100)
        ax = fig.add_subplot(111)

        # Initial plot with default X and Y axis selections
        ax.scatter(nu_values, delta_d_values, color='green', s=10)
        ax.set_xlabel("nu (GHz)")  # Frequency variable and unit
        ax.set_ylabel("delta_d (cm)")  # Spatial Resolution variable and unit
        ax.set_title("Step 2: Frequency vs Spatial Resolution")

        # Canvas for displaying the graph
        canvas = FigureCanvasTkAgg(fig, master=graph_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill='both', expand=True)

        # Control frame for selecting X and Y axes
        control_frame = LabelFrame(graph_window, text="Select X and Y Axis", padx=10, pady=10)
        control_frame.pack(fill='x', padx=10, pady=10)

        # All available options for both X and Y axes
        axis_options = ["Aperture (a)", "Frequency (nu)", "Diameter (D)", "Spatial Resolution (delta_d)",
                        "Far-field Cutoff (R_F)", "Minimum Distance (R_min)", "Reactive Near-field Cutoff (R_react)"]

        # Dropdown for X-axis selection
        x_var = StringVar(value="Frequency (nu)")
        x_label = Label(control_frame, text="Select X-axis:")
        x_label.pack(side='left', padx=5)
        x_menu = ttk.Combobox(control_frame, textvariable=x_var, values=axis_options)
        x_menu.pack(side='left', padx=5)

        # Dropdown for Y-axis selection
        y_var = StringVar(value="Spatial Resolution (delta_d)")
        y_label = Label(control_frame, text="Select Y-axis:")
        y_label.pack(side='left', padx=5)
        y_menu = ttk.Combobox(control_frame, textvariable=y_var, values=axis_options)
        y_menu.pack(side='left', padx=5)

        # Save original limits for reset
        original_x_min = min(a_values)
        original_x_max = max(a_values)

        # Entry fields for setting X-axis min and max
        x_min_label = Label(control_frame, text="X-axis Min:")
        x_min_label.pack(side='left', padx=5)
        x_min_entry = Entry(control_frame, width=5)
        x_min_entry.pack(side='left', padx=5)

        x_max_label = Label(control_frame, text="X-axis Max:")
        x_max_label.pack(side='left', padx=5)
        x_max_entry = Entry(control_frame, width=5)
        x_max_entry.pack(side='left', padx=5)

        # Function to update the plot
        def update_plot(*args):
            current_x_label = x_var.get()
            current_y_label = y_var.get()

            # Retrieve new X-axis limits
            try:
                x_min = float(x_min_entry.get())
                x_max = float(x_max_entry.get())
            except ValueError:
                x_min = original_x_min
                x_max = original_x_max

            # Generate new data for the updated range
            x_values = np.linspace(x_min, x_max, 300)
            y_values = []

            if current_x_label == "Aperture (a)":
                y_values = [calculate_step2(x_val) if current_y_label == "Spatial Resolution (delta_d)"
                            else calculate_step1(D, nu)[0] for x_val in x_values]
                ax.set_xlabel("a (m)")
            elif current_x_label == "Frequency (nu)":
                x_values = [nu] * 300
                y_values = delta_d_values
                ax.set_xlabel("nu (GHz)")
            else:
                y_values = delta_d_values  # Default to precomputed Y values

            # Clear and replot the graph
            ax.clear()
            ax.scatter(x_values, y_values, color='blue', s=10)
            ax.set_ylabel(f"{current_y_label} (unit)")
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(min(y_values), max(y_values))
            ax.set_title(f"Step 2: {current_y_label} vs {current_x_label}")
            ax.legend()
            canvas.draw()

        # Button to trigger plot update
        update_button = Button(control_frame, text="Update X-axis Limits", command=update_plot)
        update_button.pack(side='left', padx=10)

        # Reset button to restore original limits
        def reset_limits():
            x_min_entry.delete(0, 'end')
            x_max_entry.delete(0, 'end')
            ax.set_xlim(original_x_min, original_x_max)
            canvas.draw()

        reset_button = Button(control_frame, text="Reset X-axis Limits", command=reset_limits)
        reset_button.pack(side='left', padx=10)

        # Trigger plot update when X or Y axis selection changes
        x_var.trace("w", update_plot)
        y_var.trace("w", update_plot)
        x_min_entry.bind("<Return>", lambda event: update_plot())
        x_max_entry.bind("<Return>", lambda event: update_plot())

        # Save button to save the graph as an image
        save_button = Button(control_frame, text="Save Graph", command=lambda: save_graph_as_image(fig))
        save_button.pack(pady=10)

    except ValueError:
        messagebox.showwarning("Input Error", "Please enter valid numbers for Step 2 inputs.")

def plot_step3_graph():
    try:
        # Retrieve inputs from the GUI
        D = float(entry_2.get())
        nu = float(entry_5.get())
        f_1 = float(entry_3.get())
        f_osr = float(entry_7.get())
        dtheta = float(entry_19.get())
        delta_d = calculate_step2(float(entry_1.get()))

        # Generate initial ranges for the parameters
        delta_d_values = np.linspace(delta_d * 0.8, delta_d * 1.2, 300)  # Vary delta_d by ±20%
        f_apo_values = np.linspace(float(entry_6.get()) * 0.5, float(entry_6.get()) * 1.5, 300)  # Vary f_apo by ±50%
        f_osr_values = np.linspace(f_osr * 0.8, f_osr * 1.2, 300)  # Vary f_osr by ±20%
        D_values = np.linspace(D * 0.8, D * 1.2, 300)  # Vary D by ±20%
        nu_values = np.linspace(nu * 0.8, nu * 1.2, 300)  # Vary nu by ±20%
        dtheta_values = np.linspace(dtheta * 0.8, dtheta * 1.2, 300)  # Vary dtheta by ±20%

        # Prepare lists for t_int and t_map calculations
        t_int_values = []
        t_map_values = []

        # Calculate t_int and t_map for each f_apo value
        for f_apo in f_apo_values:
            t_int = 6.2 * 10**4 * f_osr * f_1 * f_apo**2 / (dtheta * nu * D)
            t_map = 171768 * f_osr * f_1 * f_apo**2 * D / (dtheta * nu * delta_d**2)
            t_int_values.append(t_int)
            t_map_values.append(t_map)

        # Create a new window for plotting
        graph_window = Toplevel()
        graph_window.title("Step 3 Graph")
        graph_window.geometry("900x600")

        # Create a figure for plotting
        fig = Figure(figsize=(6, 4), dpi=100)
        ax = fig.add_subplot(111)

        # Initial plot with f_apo vs t_int
        ax.scatter(f_apo_values, t_int_values, label="Integration Time (t_int)", color="orange", s=10)
        ax.set_xlabel("f_apo ")  # Apodization Smoothing Factor variable and unit
        ax.set_ylabel("t_int (s)")  # Integration Time variable and unit
        ax.set_title("Integration Time vs Apodization Smoothing Factor")
        ax.legend()

        # Canvas for displaying the graph
        canvas = FigureCanvasTkAgg(fig, master=graph_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

        # Control frame for selecting X and Y axes
        control_frame = LabelFrame(graph_window, text="Select X and Y Axis", padx=10, pady=10)
        control_frame.pack(fill="x", padx=10, pady=10)

        # All available options for both X and Y axes
        axis_options = [
            "Spatial Resolution (delta_d)", "f_apo", "f_osr", "Integration Time (t_int)",
            "Total Map Time (t_map)", "Frequency (nu)", "Diameter (D)", "Rotation Rate (dtheta)"
        ]

        # Dropdown for X-axis selection
        x_var = StringVar(value="f_apo")
        x_label = Label(control_frame, text="Select X-axis:")
        x_label.pack(side="left", padx=5)
        x_menu = ttk.Combobox(control_frame, textvariable=x_var, values=axis_options)
        x_menu.pack(side="left", padx=5)

        # Dropdown for Y-axis selection
        y_var = StringVar(value="Integration Time (t_int)")
        y_label = Label(control_frame, text="Select Y-axis:")
        y_label.pack(side="left", padx=5)
        y_menu = ttk.Combobox(control_frame, textvariable=y_var, values=axis_options)
        y_menu.pack(side="left", padx=5)

        # Save original limits for reset
        original_x_min = min(f_apo_values)
        original_x_max = max(f_apo_values)

        # Entry fields for setting X-axis min and max
        x_min_label = Label(control_frame, text="X-axis Min:")
        x_min_label.pack(side="left", padx=5)
        x_min_entry = Entry(control_frame, width=5)
        x_min_entry.pack(side="left", padx=5)

        x_max_label = Label(control_frame, text="X-axis Max:")
        x_max_label.pack(side="left", padx=5)
        x_max_entry = Entry(control_frame, width=5)
        x_max_entry.pack(side="left", padx=5)

        # Update plot function to dynamically recalculate Y values for the new X range
        def update_plot(*args):
            current_x_label = x_var.get()
            current_y_label = y_var.get()

            try:
                x_min = float(x_min_entry.get())
                x_max = float(x_max_entry.get())
            except ValueError:
                x_min = original_x_min
                x_max = original_x_max

            # Generate a new range of X values
            x_values = np.linspace(x_min, x_max, 300)
            y_values = []

            # Recalculate Y values based on the selected axis
            for x_val in x_values:
                if current_x_label == "f_apo":
                    if current_y_label == "Integration Time (t_int)":
                        y_val = 6.2 * 10**4 * f_osr * f_1 * x_val**2 / (dtheta * nu * D)
                    elif current_y_label == "Total Map Time (t_map)":
                        y_val = 171768 * f_osr * f_1 * x_val**2 * D / (dtheta * nu * delta_d**2)
                    else:
                        y_val = x_val
                else:
                    y_val = x_val  # Default for non-dependent Y-axis
                y_values.append(y_val)

            # Clear and replot the graph
            ax.clear()
            ax.scatter(x_values, y_values, color="blue", s=10)
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(min(y_values), max(y_values))
            ax.set_xlabel(f"{current_x_label} (unit)")
            ax.set_ylabel(f"{current_y_label} (unit)")
            ax.set_title(f"Step 3: {current_y_label} vs {current_x_label}")
            ax.legend()
            canvas.draw()

        # Button to trigger plot update
        update_button = Button(control_frame, text="Update X-axis Limits", command=update_plot)
        update_button.pack(side="left", padx=10)

        # Reset button to restore original limits
        def reset_limits():
            x_min_entry.delete(0, "end")
            x_max_entry.delete(0, "end")
            ax.set_xlim(original_x_min, original_x_max)
            canvas.draw()

        reset_button = Button(control_frame, text="Reset X-axis Limits", command=reset_limits)
        reset_button.pack(side="left", padx=10)

        # Trigger plot update when X or Y axis selection changes
        x_var.trace("w", update_plot)
        y_var.trace("w", update_plot)
        x_min_entry.bind("<Return>", lambda event: update_plot())
        x_max_entry.bind("<Return>", lambda event: update_plot())

        # Save button to save the graph as an image
        save_button = Button(control_frame, text="Save Graph", command=lambda: save_graph_as_image(fig))
        save_button.pack(pady=10)

    except ValueError as e:
        messagebox.showwarning("Input Error", str(e))

def plot_step4_graph():
    try:
        # Retrieve inputs from the GUI
        nu = float(entry_5.get())
        D = float(entry_2.get())
        f_osr = float(entry_7.get())
        f_1 = float(entry_3.get())
        f_apo = float(entry_6.get())
        delta_d = calculate_step2(float(entry_1.get()))
        dtheta = float(entry_19.get())

        # Generate initial ranges for values
        f_osr_values = np.linspace(f_osr * 0.8, f_osr * 1.2, 300)
        f_1_values = np.linspace(f_1 * 0.8, f_1 * 1.2, 300)
        f_apo_values = np.linspace(f_apo * 0.8, f_apo * 1.2, 300)
        delta_d_values = np.linspace(delta_d * 0.8, delta_d * 1.2, 300)
        nu_values = np.linspace(nu * 0.8, nu * 1.2, 300)
        D_values = np.linspace(D * 0.8, D * 1.2, 300)
        dtheta_values = np.linspace(dtheta * 0.8, dtheta * 1.2, 300)

        theta_ext_values = []
        theta_b_values = []
        theta_sr_values = []
        theta_ss_values = []

        for f_osr_val in f_osr_values:
            theta_ext, theta_b, theta_sr, theta_ss = calculate_step4(f_osr_val)
            theta_ext_values.append(theta_ext)
            theta_b_values.append(theta_b)
            theta_sr_values.append(theta_sr)
            theta_ss_values.append(theta_ss)

        # Create a new window for plotting
        graph_window = Toplevel()
        graph_window.title("Step 4 Graph")
        graph_window.geometry("900x600")

        # Create a figure for plotting
        fig = Figure(figsize=(6, 4), dpi=100)
        ax = fig.add_subplot(111)

        # Initial plot with default selections
        ax.scatter(f_osr_values, theta_ext_values, label="Angular Extent (theta_ext)", color="purple", s=15)
        ax.set_xlabel("f_osr")
        ax.set_ylabel("theta_ext (degrees)")
        ax.set_title("Holography Map Sampling Intervals")
        ax.legend()

        # Canvas for displaying the graph
        canvas = FigureCanvasTkAgg(fig, master=graph_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

        # Control frame for selecting X and Y axes
        control_frame = LabelFrame(graph_window, text="Select X and Y Axis", padx=10, pady=10)
        control_frame.pack(fill="x", padx=10, pady=10)

        # All available options for both X and Y axes
        axis_options = [
            "Oversampling (f_osr)", "f_1", "f_apo", "Spatial Resolution (delta_d)",
            "Frequency (nu)", "Diameter (D)", "Rotation Rate (dtheta)",
            "Angular Extent (theta_ext)", "Beamwidth (theta_b)",
            "Sampling Interval (theta_sr)", "Sampling Interval (theta_ss)"
        ]

        # Dropdown for X-axis selection
        x_var = StringVar(value="Oversampling (f_osr)")
        x_label = Label(control_frame, text="Select X-axis:")
        x_label.pack(side="left", padx=5)
        x_menu = ttk.Combobox(control_frame, textvariable=x_var, values=axis_options)
        x_menu.pack(side="left", padx=5)

        # Dropdown for Y-axis selection
        y_var = StringVar(value="Angular Extent (theta_ext)")
        y_label = Label(control_frame, text="Select Y-axis:")
        y_label.pack(side="left", padx=5)
        y_menu = ttk.Combobox(control_frame, textvariable=y_var, values=axis_options)
        y_menu.pack(side="left", padx=5)

        # Save original limits for reset
        original_x_min = min(f_osr_values)
        original_x_max = max(f_osr_values)

        # Entry fields for setting X-axis min and max
        x_min_label = Label(control_frame, text="X-axis Min:")
        x_min_label.pack(side="left", padx=5)
        x_min_entry = Entry(control_frame, width=5)
        x_min_entry.pack(side="left", padx=5)

        x_max_label = Label(control_frame, text="X-axis Max:")
        x_max_label.pack(side="left", padx=5)
        x_max_entry = Entry(control_frame, width=5)
        x_max_entry.pack(side="left", padx=5)

        # Update plot function to dynamically recalculate Y values for the new X range
        def update_plot(*args):
            current_x_label = x_var.get()
            current_y_label = y_var.get()

            try:
                x_min = float(x_min_entry.get())
                x_max = float(x_max_entry.get())
            except ValueError:
                x_min = original_x_min
                x_max = original_x_max

            # Generate a new range of X values
            x_values = np.linspace(x_min, x_max, 300)
            y_values = []

            # Dynamically recalculate Y values based on the selected axis
            for x_val in x_values:
                if current_x_label == "Oversampling (f_osr)":
                    theta_ext, theta_b, theta_sr, theta_ss = calculate_step4(x_val)
                    if current_y_label == "Angular Extent (theta_ext)":
                        y_val = theta_ext
                    elif current_y_label == "Beamwidth (theta_b)":
                        y_val = theta_b
                    elif current_y_label == "Sampling Interval (theta_sr)":
                        y_val = theta_sr
                    elif current_y_label == "Sampling Interval (theta_ss)":
                        y_val = theta_ss
                    else:
                        y_val = x_val
                else:
                    y_val = x_val  # Default fallback for independent Y-axis

                y_values.append(y_val)

            # Clear and replot the graph
            ax.clear()
            ax.scatter(x_values, y_values, color="blue", s=10)
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(min(y_values), max(y_values))
            ax.set_xlabel(f"{current_x_label} (unit)")
            ax.set_ylabel(f"{current_y_label} (unit)")
            ax.set_title(f"{current_y_label} vs {current_x_label}")
            ax.legend()
            canvas.draw()

        # Button to trigger plot update
        update_button = Button(control_frame, text="Update X-axis Limits", command=update_plot)
        update_button.pack(side="left", padx=10)

        # Reset button to restore original limits
        def reset_limits():
            x_min_entry.delete(0, "end")
            x_max_entry.delete(0, "end")
            ax.set_xlim(original_x_min, original_x_max)
            canvas.draw()

        reset_button = Button(control_frame, text="Reset X-axis Limits", command=reset_limits)
        reset_button.pack(side="left", padx=10)

        # Trigger plot update when X or Y axis selection changes
        x_var.trace("w", update_plot)
        y_var.trace("w", update_plot)
        x_min_entry.bind("<Return>", lambda event: update_plot())
        x_max_entry.bind("<Return>", lambda event: update_plot())

        # Save button to save the graph as an image
        save_button = Button(control_frame, text="Save Graph", command=lambda: save_graph_as_image(fig))
        save_button.pack(pady=10)

    except ValueError as e:
        messagebox.showwarning("Input Error", str(e))

def plot_step5_graph():
    try:
        # Retrieve inputs from the GUI
        delta_z = float(entry_16.get())  # Surface deformation input
        delta_d = calculate_step2(float(entry_1.get()))  # Spatial resolution from Step 2
        f_apo = float(entry_6.get())  # Apodization smoothing factor from entry_6
        f_osr = float(entry_7.get())  # Oversampling factor from entry_7
        f_oss = float(entry_12.get())  # Oversampling factor from entry_12

        # Generate a range of values for plotting
        delta_z_values = np.linspace(delta_z * 0.8, delta_z * 1.2, 300)
        delta_d_values = np.linspace(delta_d * 0.8, delta_d * 1.2, 300)
        f_apo_values = np.linspace(f_apo * 0.8, f_apo * 1.2, 300)
        f_osr_values = np.linspace(f_osr * 0.8, f_osr * 1.2, 300)
        f_oss_values = np.linspace(f_oss * 0.8, f_oss * 1.2, 300)

        theta_point_values = []
        N_row_values = []
        SNR_values = []

        for delta_z_val in delta_z_values:
            theta_point, N_row, SNR = calculate_step5_with_delta_z(delta_d, delta_z_val)
            theta_point_values.append(theta_point)
            N_row_values.append(N_row)
            SNR_values.append(SNR)

        # Create a new window for plotting
        graph_window = Toplevel()
        graph_window.title("Step 5 Graph")
        graph_window.geometry("900x600")

        # Create a figure for plotting
        fig = Figure(figsize=(6, 4), dpi=100)
        ax = fig.add_subplot(111)

        # Initial plot with default selections
        ax.scatter(delta_z_values, theta_point_values, label="Pointing Accuracy (theta_point)", color="brown", s=10)
        ax.set_xlabel("delta_z (µm)")
        ax.set_ylabel("theta_point (degrees)")
        ax.set_title("Step 5: Pointing Accuracy and SNR")
        ax.legend()

        # Canvas for displaying the graph
        canvas = FigureCanvasTkAgg(fig, master=graph_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

        # Control frame for selecting X and Y axes
        control_frame = LabelFrame(graph_window, text="Select X and Y Axis", padx=10, pady=10)
        control_frame.pack(fill="x", padx=10, pady=10)

        # All available options for both X and Y axes
        axis_options = ["Surface Deformation (delta_z)", "Pointing Accuracy (theta_point)", "Number of Rows (N_row)",
                        "Spatial Resolution (delta_d)", "f_apo", "f_osr", "f_oss", "Signal-to-Noise Ratio (SNR)"]

        # Dropdown for X-axis selection
        x_var = StringVar(value="Surface Deformation (delta_z)")
        x_label = Label(control_frame, text="Select X-axis:")
        x_label.pack(side="left", padx=5)
        x_menu = ttk.Combobox(control_frame, textvariable=x_var, values=axis_options)
        x_menu.pack(side="left", padx=5)

        # Dropdown for Y-axis selection
        y_var = StringVar(value="Pointing Accuracy (theta_point)")
        y_label = Label(control_frame, text="Select Y-axis:")
        y_label.pack(side="left", padx=5)
        y_menu = ttk.Combobox(control_frame, textvariable=y_var, values=axis_options)
        y_menu.pack(side="left", padx=5)

        # Save original limits for reset
        original_x_min = min(delta_z_values)
        original_x_max = max(delta_z_values)

        # Entry fields for setting X-axis min and max
        x_min_label = Label(control_frame, text="X-axis Min:")
        x_min_label.pack(side="left", padx=5)
        x_min_entry = Entry(control_frame, width=5)
        x_min_entry.pack(side="left", padx=5)

        x_max_label = Label(control_frame, text="X-axis Max:")
        x_max_label.pack(side="left", padx=5)
        x_max_entry = Entry(control_frame, width=5)
        x_max_entry.pack(side="left", padx=5)

        # Update plot function to recalculate Y values dynamically
        def update_plot(*args):
            current_x_label = x_var.get()
            current_y_label = y_var.get()

            try:
                x_min = float(x_min_entry.get())
                x_max = float(x_max_entry.get())
            except ValueError:
                x_min = original_x_min
                x_max = original_x_max

            # Generate new X-axis values based on the user-defined range
            x_values = np.linspace(x_min, x_max, 300)
            y_values = []

            # Recalculate Y values dynamically for the new X range
            for x_val in x_values:
                if current_x_label == "Surface Deformation (delta_z)":
                    delta_z = x_val
                    theta_point, N_row, SNR = calculate_step5_with_delta_z(delta_d, delta_z)
                    if current_y_label == "Pointing Accuracy (theta_point)":
                        y_values.append(theta_point)
                    elif current_y_label == "Number of Rows (N_row)":
                        y_values.append(N_row)
                    elif current_y_label == "Signal-to-Noise Ratio (SNR)":
                        y_values.append(SNR)
                elif current_x_label == "Spatial Resolution (delta_d)":
                    y_values.append(delta_d)  # Placeholder; extend logic for other variables
                elif current_x_label == "f_apo":
                    y_values.append(f_apo)  # Placeholder
                elif current_x_label == "f_osr":
                    y_values.append(f_osr)  # Placeholder
                elif current_x_label == "f_oss":
                    y_values.append(f_oss)  # Placeholder

            # Clear and replot the graph
            ax.clear()
            ax.scatter(x_values, y_values, color="blue", s=10)
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(min(y_values), max(y_values))
            ax.set_xlabel(f"{current_x_label} (unit)")
            ax.set_ylabel(f"{current_y_label} (unit)")
            ax.set_title(f"{current_y_label} vs {current_x_label}")
            ax.legend()
            canvas.draw()

        # Button to trigger plot update
        update_button = Button(control_frame, text="Update X-axis Limits", command=update_plot)
        update_button.pack(side="left", padx=10)

        # Reset button to restore original limits
        def reset_limits():
            x_min_entry.delete(0, "end")
            x_max_entry.delete(0, "end")
            ax.set_xlim(original_x_min, original_x_max)
            canvas.draw()

        reset_button = Button(control_frame, text="Reset X-axis Limits", command=reset_limits)
        reset_button.pack(side="left", padx=10)

        # Trigger plot update when X or Y axis selection changes
        x_var.trace("w", update_plot)
        y_var.trace("w", update_plot)
        x_min_entry.bind("<Return>", lambda event: update_plot())
        x_max_entry.bind("<Return>", lambda event: update_plot())

        # Save button to save the graph as an image
        save_button = Button(control_frame, text="Save Graph", command=lambda: save_graph_as_image(fig))
        save_button.pack(pady=10)

    except ValueError:
        messagebox.showwarning("Input Error", "Please enter valid numbers for Step 5 inputs.")

def plot_step6_graph():
    try:
        # Retrieve inputs from the GUI
        D_t = float(entry_17.get())  # Transmitter Diameter
        D = float(entry_2.get())  # Antenna Diameter
        z = float(entry_8.get())  # Distance
        T_sys = float(entry_9.get())  # System Temperature
        B = float(entry_10.get())  # Bandwidth
        nu = float(entry_5.get())  # Frequency

        # Generate ranges of values for each selected variable
        D_t_values = np.linspace(D_t * 0.8, D_t * 1.2, 300)  # Vary D_t by ±20%
        D_values = np.linspace(D * 0.8, D * 1.2, 300)  # Vary D by ±20%
        z_values = np.linspace(z * 0.8, z * 1.2, 300)  # Vary z by ±20%
        T_sys_values = np.linspace(T_sys * 0.8, T_sys * 1.2, 300)  # Vary T_sys by ±20%
        B_values = np.linspace(B * 0.8, B * 1.2, 300)  # Vary B by ±20%
        nu_values = np.linspace(nu * 0.8, nu * 1.2, 300)  # Vary nu by ±20%

        # Use calculate_step6 to generate power values for D_t
        power_values_D_t = []
        for D_t_val in D_t_values:
            # Temporarily update Entry values to use calculate_step6 as is
            entry_17.delete(0, 'end')
            entry_17.insert(0, str(D_t_val))

            entry_2.delete(0, 'end')
            entry_2.insert(0, str(D))

            entry_8.delete(0, 'end')
            entry_8.insert(0, str(z))

            entry_9.delete(0, 'end')
            entry_9.insert(0, str(T_sys))

            entry_10.delete(0, 'end')
            entry_10.insert(0, str(B))

            entry_5.delete(0, 'end')
            entry_5.insert(0, str(nu))

            power_values_D_t.append(calculate_step6())

        # Create a new window for plotting
        graph_window = Toplevel()
        graph_window.title("Step 6 Graph")
        graph_window.geometry("900x600")

        # Create a figure for plotting
        fig = Figure(figsize=(6, 4), dpi=100)
        ax = fig.add_subplot(111)

        # Default plot with D_t on the x-axis and output power on the y-axis
        ax.scatter(D_t_values, power_values_D_t, label='Power vs D_t', color='blue', s=10)
        ax.set_xlabel("D_t (m)")  # Variable and unit for Transmitter Diameter
        ax.set_ylabel("P (W)")  # Variable and unit for Power
        ax.set_title("Step 6: Transmitter Output Power")
        ax.legend()

        # Canvas for displaying the graph
        canvas = FigureCanvasTkAgg(fig, master=graph_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill='both', expand=True)

        # Control frame for selecting X and Y axes
        control_frame = LabelFrame(graph_window, text="Select X and Y Axis", padx=10, pady=10)
        control_frame.pack(fill='x', padx=10, pady=10)

        # Limited axis options for both X and Y
        axis_options = ["Bandwidth (B)", "System Temperature (T_sys)", "Distance (z)",
                        "Transmitter Diameter (D_t)", "Power (P)", "Diameter (D)", "Frequency (nu)"]

        # Dropdown for X-axis selection
        x_var = StringVar(value="Transmitter Diameter (D_t)")
        x_label = Label(control_frame, text="Select X-axis:")
        x_label.pack(side='left', padx=5)
        x_menu = ttk.Combobox(control_frame, textvariable=x_var, values=axis_options)
        x_menu.pack(side='left', padx=5)

        # Dropdown for Y-axis selection
        y_var = StringVar(value="Power (P)")
        y_label = Label(control_frame, text="Select Y-axis:")
        y_label.pack(side='left', padx=5)
        y_menu = ttk.Combobox(control_frame, textvariable=y_var, values=axis_options)
        y_menu.pack(side='left', padx=5)

        # Save original limits for reset
        original_x_min = min(D_t_values)
        original_x_max = max(D_t_values)

        # Entry fields for setting X-axis min and max
        x_min_label = Label(control_frame, text="X-axis Min:")
        x_min_label.pack(side='left', padx=5)
        x_min_entry = Entry(control_frame, width=5)
        x_min_entry.pack(side='left', padx=5)

        x_max_label = Label(control_frame, text="X-axis Max:")
        x_max_label.pack(side='left', padx=5)
        x_max_entry = Entry(control_frame, width=5)
        x_max_entry.pack(side='left', padx=5)

        # Update plot function to apply new selections and X-axis limits
        def update_plot(*args):
            current_x_label = x_var.get()
            current_y_label = y_var.get()

            # Clear the axes and replot with the selected data
            ax.clear()

            # Choose data for X-axis and set units
            if current_x_label == "Bandwidth (B)":
                x_values = B_values
                ax.set_xlabel("B (MHz)")
            elif current_x_label == "System Temperature (T_sys)":
                x_values = T_sys_values
                ax.set_xlabel("T_sys (K)")
            elif current_x_label == "Distance (z)":
                x_values = z_values
                ax.set_xlabel("z (m)")
            elif current_x_label == "Transmitter Diameter (D_t)":
                x_values = D_t_values
                ax.set_xlabel("D_t (m)")
            elif current_x_label == "Power (P)":
                x_values = power_values_D_t
                ax.set_xlabel("P (W)")
            elif current_x_label == "Diameter (D)":
                x_values = D_values
                ax.set_xlabel("D (m)")
            elif current_x_label == "Frequency (nu)":
                x_values = nu_values
                ax.set_xlabel("nu (GHz)")

            # Choose data for Y-axis and set units
            if current_y_label == "Power (P)":
                y_values = power_values_D_t
                ax.scatter(x_values, y_values, color='blue', s=10)
                ax.set_ylabel("P (W)")
            elif current_y_label == "Bandwidth (B)":
                y_values = B_values
                ax.scatter(x_values, y_values, color='red', s=10)
                ax.set_ylabel("B (MHz)")
            elif current_y_label == "System Temperature (T_sys)":
                y_values = T_sys_values
                ax.scatter(x_values, y_values, color='green', s=10)
                ax.set_ylabel("T_sys (K)")
            elif current_y_label == "Distance (z)":
                y_values = z_values
                ax.scatter(x_values, y_values, color='purple', s=10)
                ax.set_ylabel("z (m)")
            elif current_y_label == "Transmitter Diameter (D_t)":
                y_values = D_t_values
                ax.scatter(x_values, y_values, color='orange', s=10)
                ax.set_ylabel("D_t (m)")
            elif current_y_label == "Diameter (D)":
                y_values = D_values
                ax.scatter(x_values, y_values, color='cyan', s=10)
                ax.set_ylabel("D (m)")
            elif current_y_label == "Frequency (nu)":
                y_values = nu_values
                ax.scatter(x_values, y_values, color='yellow', s=10)
                ax.set_ylabel("nu (GHz)")

            # Set limits for X-axis based on user input
            try:
                x_min = float(x_min_entry.get())
                x_max = float(x_max_entry.get())
                ax.set_xlim(x_min, x_max)
            except ValueError:
                ax.set_xlim(min(x_values), max(x_values))

            ax.set_title(f"Step 6: {current_y_label} vs {current_x_label}")
            ax.legend()
            canvas.draw()

        # Button to trigger plot update
        update_button = Button(control_frame, text="Update X-axis Limits", command=update_plot)
        update_button.pack(side='left', padx=10)

        # Reset button to restore original limits
        def reset_limits():
            x_min_entry.delete(0, 'end')
            x_max_entry.delete(0, 'end')
            ax.set_xlim(original_x_min, original_x_max)
            canvas.draw()

        reset_button = Button(control_frame, text="Reset X-axis Limits", command=reset_limits)
        reset_button.pack(side='left', padx=10)

        # Trigger plot update when X or Y axis selection changes
        x_var.trace("w", update_plot)
        y_var.trace("w", update_plot)
        x_min_entry.bind("<Return>", lambda event: update_plot())
        x_max_entry.bind("<Return>", lambda event: update_plot())

        # Save button to save the graph as an image
        save_button = Button(control_frame, text="Save Graph", command=lambda: save_graph_as_image(fig))
        save_button.pack(pady=10)

    except ValueError:
        messagebox.showwarning("Input Error", "Please enter valid numbers for Step 6 inputs.")

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


window = Tk()

icon_image = PhotoImage(file=relative_to_assets("HS_logo.png"))
window.iconphoto(False, icon_image)

window.geometry("1200x976")
window.configure(bg="#808080")
window.title("HoloSim")

canvas = Canvas(
    window,
    bg="#808080",
    height=976,
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
image_D_t = load_and_display_image(canvas, "D_t.png", 1000.0, 511.0)
image_T_sys = load_and_display_image(canvas, "T_sys.png", 997.0, 550.0)
image_T_R_z = load_and_display_image(canvas, "T_R_z.png", 1000.0, 529.0)
image_B = load_and_display_image(canvas, "B.png", 930.0, 571.0)

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
    canvas.create_text(850.0, 502.0, anchor="nw", text="Transmitter diameter      : ", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(1140, 502.0, anchor="nw", text="m", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(850.0, 520.0, anchor="nw", text="Distance between", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(1140, 520.0, anchor="nw", text="m", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(1030.0, 520.0, anchor="nw", text=":", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(850.0, 540.0, anchor="nw", text="System temperature        : ", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(1140, 540.0, anchor="nw", text="K", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(850.0, 562.0, anchor="nw", text="Bandwidth     :", fill="#000000", font=("Inter SemiBold", 14))
    canvas.create_text(1140, 562.0, anchor="nw", text="MHz", fill="#000000", font=("Inter SemiBold", 14))
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

entry_image_8 = load_and_display_image(canvas, "entry_8.png", 1090.0, 530.5)
entry_8 = create_entry(canvas, "entry_8.png", 1090.0, 530.5, 77.0, 17.0)

entry_image_9 = load_and_display_image(canvas, "entry_9.png", 1090.0, 551.5)
entry_9 = create_entry(canvas, "entry_9.png", 1090.0, 551.5, 77.0, 17.0)

entry_image_10 = load_and_display_image(canvas, "entry_10.png", 1090.0, 572.5)
entry_10 = create_entry(canvas, "entry_10.png", 1090.0, 572.5, 77.0, 17.0)

entry_image_12 = load_and_display_image(canvas, "entry_12.png", 305.0, 489.0)
entry_12 = create_entry(canvas, "entry_12.png", 305.0, 489.0, 77.0, 17.0)

entry_image_16 = load_and_display_image(canvas, "entry_16.png", 690.0, 490.0)
entry_16 = create_entry(canvas, "entry_16.png", 690.0, 490.0, 77.0, 17.0)

entry_image_17 = load_and_display_image(canvas, "entry_17.png", 1090.0, 510.0)
entry_17 = create_entry(canvas, "entry_17.png", 1090.0, 510.0, 77.0, 17.0)

entry_image_19 = load_and_display_image(canvas, "entry_19.png", 1110.0, 224.0)
entry_19 = create_entry(canvas, "entry_19.png", 1110.0, 224.0, 78.0, 16.0)


button_image_1 = PhotoImage(file=relative_to_assets("button1.png"))
button_1 = Button(image=button_image_1, borderwidth=0, highlightthickness=0, command=calculate_all_steps, relief="flat")
button_1.place(x=1040, y=730, width=148, height=46)

button_image_2 = PhotoImage(file=relative_to_assets("button_2.png")) # Connect the updated handler
button_2 = Button(image=button_image_2, borderwidth=0, highlightthickness=0, command=calculate_and_display_step1, relief="flat")
button_2.place(x=30.0, y=300.0, width=135, height=50)

button_image_3 = PhotoImage(file=relative_to_assets("button_3.png")) # Connect Button 3 to the calculate_and_display_step4 function
button_3 = Button(image=button_image_3, borderwidth=0, highlightthickness=0, command=calculate_and_display_step4, relief="flat")
button_3.place(x=30, y=630, width=135, height=50)

button_image_4 = PhotoImage(file=relative_to_assets("button_4.png"))
button_4 = Button(image=button_image_4, borderwidth=0, highlightthickness=0, command=calculate_and_display_step5, relief="flat")
button_4.place(x=450, y=630, width=135, height=50)

button_image_5 = PhotoImage(file=relative_to_assets("button_5.png"))
button_5 = Button(image=button_image_5, borderwidth=0, highlightthickness=0, command=calculate_and_display_step2,  relief="flat")
button_5.place(x=450.0, y=300.0, width=135, height=50)

button_image_6 = PhotoImage(file=relative_to_assets("button_6.png"))
button_6 = Button(image=button_image_6,borderwidth=0,highlightthickness=0,command=calculate_and_display_step3, relief="flat") # Connect Button 6 to the calculate_and_display_step3 function
button_6.place(x=860.0, y=300, width=135, height=50)

button_image_7 = PhotoImage(file=relative_to_assets("button_7.png"))
button_7 = Button(image=button_image_7, borderwidth=0, highlightthickness=0, command=calculate_and_display_step6,  relief="flat") # Connect to Step 6 function
button_7.place(x=860.0, y=630.0, width=135, height=50)


button_image_8 = PhotoImage(file=relative_to_assets("button_8.png"))
button_8 = PhotoImage(file=relative_to_assets("button_8.png"))
#button_8 = Button(image=button_image_8, borderwidth=0, highlightthickness=0, command=lambda: [print("Step 1 Graph button pressed"), plot_step1(float(entry_2.get()), float(entry_5.get()))],relief="flat")
#button_8.place(x=170.0, y=300.0, width=135, height=50)

button_image_9 = PhotoImage(file=relative_to_assets("button_9.png"))
#button_9 = Button(image=button_image_9, borderwidth=0, highlightthickness=0, command=lambda: [print("Step 2 Graph button pressed"), plot_step2(float(entry_1.get()), calculate_step2(float(entry_1.get())))],  relief="flat")
#button_9.place(x=590.0, y=300.0, width=135, height=50)

button_image_10 = PhotoImage(file=relative_to_assets("button_10.png"))
#button_10 = Button(image=button_image_10,borderwidth=0,highlightthickness=0,command=lambda: [print("Step 3 Graph button pressed"), plot_step3(*calculate_step3(calculate_step2(float(entry_1.get().strip()))), calculate_step2(float(entry_1.get().strip())))],relief="flat")
#button_10.place(x=1000.0, y=300.0, width=135, height=50)

button_image_11 = PhotoImage(file=relative_to_assets("button_11.png"))
#button_11 = Button(image=button_image_11, borderwidth=0, highlightthickness=0, command=lambda: [print("Step 4 Graph button pressed"), plot_step4(*calculate_step4(float(entry_1.get())))],relief="flat")
#button_11.place(x=170, y=630, width=135, height=50)

button_image_12 = PhotoImage(file=relative_to_assets("button_12.png"))
#button_12 = Button(image=button_image_12, borderwidth=0, highlightthickness=0, command=lambda: [print("Step 5 Graph button pressed"), plot_step5(*calculate_step5(calculate_step2(float(entry_1.get()))))],  relief="flat")
#button_12.place(x=590, y=630, width=135, height=50)

button_image_13 = PhotoImage(file=relative_to_assets("button_13.png"))
#button_13 = Button(image=button_image_13, borderwidth=0, highlightthickness=0, command=lambda: [print("Step 6 Graph button pressed"), plot_step6(calculate_step6(),float(entry_17.get()),float(entry_18.get()), float(entry_8.get()),float(entry_9.get()),float(entry_10.get()))],  relief="flat")
#button_13.place(x=1000.0, y=630.0, width=135, height=50)

Massive_Help_Button_image = PhotoImage(file=relative_to_assets("Massive_Help_Button.png"))
Massive_Help_Button = Button(image=Massive_Help_Button_image, borderwidth=0, highlightthickness=0, command=open_help_window, relief="flat")
Massive_Help_Button.place(x=1050.0, y=20.0, width=141, height=36)

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



# Example button setup to trigger the graph plotting for each step
button_plot_step1 = Button(image=button_image_8, borderwidth=0, highlightthickness=0,
                           command=plot_step1_graph, relief="flat")
button_plot_step1.place(x=170.0, y=300.0, width=135, height=50)

button_plot_step2 = Button(image=button_image_9, borderwidth=0, highlightthickness=0,
                           command=plot_step2_graph, relief="flat")
button_plot_step2.place(x=590.0, y=300.0, width=135, height=50)

button_plot_step3 = Button(image=button_image_10, borderwidth=0, highlightthickness=0,
                           command=plot_step3_graph, relief="flat")
button_plot_step3.place(x=1000.0, y=300.0, width=135, height=50)

button_plot_step4 = Button(image=button_image_11, borderwidth=0, highlightthickness=0,
                           command=plot_step4_graph, relief="flat")
button_plot_step4.place(x=170, y=630, width=135, height=50)

button_plot_step5 = Button(image=button_image_12, borderwidth=0, highlightthickness=0,
                           command=plot_step5_graph, relief="flat")
button_plot_step5.place(x=590, y=630, width=135, height=50)


button_plot_step6 = Button(image=button_image_13, borderwidth=0, highlightthickness=0,
                           command=plot_step6_graph, relief="flat")
button_plot_step6.place(x=1000.0, y=630.0, width=135, height=50)



window.mainloop()