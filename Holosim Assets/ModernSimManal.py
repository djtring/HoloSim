import tkinter as tk
from tkinter import ttk, filedialog
from pathlib import Path
from tkinter import Tk, Canvas, Entry, Button, PhotoImage, messagebox, Toplevel, Label
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

#from scipy.integrate import dblquad 

OUTPUT_PATH = Path(__file__).parent
ASSETS_PATH = OUTPUT_PATH / Path(r"/Users/manalsiddiqui/Desktop/Short Gui/build/assets/frame0")

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
        D_a = float(entry_18.get().strip())  # Receiver (aperture) diameter
        z = float(entry_8.get().strip())  # Distance between transmitter and receiver
        T_sys = float(entry_9.get().strip())  # System temperature
        B = float(entry_10.get().strip())  # Bandwidth in MHz
        
        # SNR from previous step (or set a default value)
        SNR = 28.23  # This should match your value from Step 5
        
        # Constants from the paper
        k = 1.38 * 10**-23  # Boltzmann constant in J/K
        G = 33  # Gain in dB

        # Noise floor
        noise_floor = 10 * np.log10(k * B * 10**6 * T_sys / (1 * 10**-3))
        print(f"Noise floor: {noise_floor}")

        # Wavelength (lam)
        lam = 3 * 10**8 / (float(entry_5.get()) * 10**9)  # Frequency from Step 1 (GHz)
        print(f"Wavelength (lambda): {lam}")

        # Beam waist at the transmitter
        w_0 = D_t / 2  # Beam waist
        z_R = np.pi * w_0**2 / lam  # Rayleigh range
        w_z = w_0 * np.sqrt(1 + (z / z_R)**2)  # Beam radius at the receiver location
        print(f"Beam waist (w_0): {w_0}, Rayleigh range (z_R): {z_R}, Beam radius (w_z): {w_z}")

        # Define integrand for the overlap integral
        def integrand(r, theta):
            return r * np.exp(-r**2 / w_z**2)

        # Calculate the overlap integral over the aperture
        overlap, _ = 1 #dblquad(integrand, 0, D_a / 2, lambda theta: 0, lambda theta: 2 * np.pi)
        print(f"Overlap integral: {overlap}")

        # Calculate beam-coupling efficiency
        eta = (overlap / (np.pi * w_z**2 / 2))**2
        eta_dB = 10 * np.log10(eta)
        print(f"Beam-coupling efficiency (eta): {eta}, Beam-coupling efficiency (eta_dB): {eta_dB}")

        # Calculate transmitter output power in dBm
        P_dB = noise_floor + SNR + eta_dB + G
        print(f"Transmitter output power (P_dB): {P_dB}")

        # Convert to watts
        P = (1 * 10**-3) * 10**(P_dB / 10)
        print(f"Transmitter output power: {P} W")

        return P

    except ValueError:
        messagebox.showerror("Input Error", "Please enter valid numbers for Step 6 inputs.")
        return None


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
    ax.scatter(data_dict[x_label_text], data_dict[y_label_text], color='blue', label=f'Step {step} Data')
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
        ax.scatter(data_dict[current_x_label], data_dict[current_y_label], color='blue', label=f'{current_x_label} vs {current_y_label}')
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

# Function to save the graph as an image
def save_graph_as_image(fig):
    # Ask the user where they want to save the image
    file_path = filedialog.asksaveasfilename(defaultextension='.png', filetypes=[("PNG files", "*.png"), ("JPEG files", "*.jpg"), ("All files", "*.*")])
    if file_path:
        fig.savefig(file_path)
        messagebox.showinfo("Save Successful", f"Graph saved as {file_path}")

# Modify plot_step to pass input field labels as axis names
def plot_step(result, step, x_label_text, y_label_text):
    if result is not None:
        if isinstance(result, tuple):  # Handle multiple results like in step 1
            open_graph_window(step, list(range(1, len(result) + 1)), result, x_label_text, y_label_text)
        else:  # Handle single result
            open_graph_window(step, [1], [result], x_label_text, y_label_text)

# Update the rest of the plotting functions similarly
def plot_step1(D, nu):
    if D and nu:
        # Calculate values for Step 1
        R_F, R_min, R_react = calculate_step1(D, nu)

        # Prepare the data dictionary
        data_dict = {
            "Frequency (ν)": [nu],
            "Diameter (D)": [D],
            "Far-Field Cutoff (R_F)": [R_F],
            "Minimum Distance (R_min)": [R_min],
            "Reactive Near-Field (R_react)": [R_react]
        }

        # Call the graph window with the data dictionary
        open_graph_window(1, data_dict, "Frequency (ν)", "Diameter (D)")

# Modify plot_step2 function to display the exact graph for Step 2 when the graph button is clicked
def plot_step2(a, delta_d):
    if a and delta_d:
        # Prepare the data dictionary for Step 2
        data_dict = {
            "Aperture (a)": [a],
            "Spatial Resolution (Δd)": [delta_d],
            "Diameter (D)": [float(entry_2.get())],   # D from Step 1
            "Frequency (ν)": [float(entry_5.get())],  # Frequency (ν) from Step 1
            "Far-field Cutoff (R_F)": [calculate_step1(float(entry_2.get()), float(entry_5.get()))[0]],  # R_F from Step 1
            "Minimum Distance (R_min)": [calculate_step1(float(entry_2.get()), float(entry_5.get()))[1]],  # R_min from Step 1
            "Reactive Near-field (R_react)": [calculate_step1(float(entry_2.get()), float(entry_5.get()))[2]]  # R_react from Step 1
        }

        # Filter only the required keys for Step 2
        allowed_keys = ["Aperture (a)", "Spatial Resolution (Δd)", "Diameter (D)", "Frequency (ν)", "Far-field Cutoff (R_F)", "Minimum Distance (R_min)", "Reactive Near-field (R_react)"]
        filtered_data_dict = {k: data_dict[k] for k in allowed_keys}

        # Call open_graph_window for Step 2
        open_graph_window(2, filtered_data_dict, "Aperture (a)", "Spatial Resolution (Δd)")

def plot_step3(t_int, t_map, delta_d):
    try:
        # Retrieve the actual input values from Step 3 fields
        f_1 = float(entry_3.get().strip())  # Primary Beam Taper Factor
        f_apo = float(entry_6.get().strip())  # Apodization Smoothing Factor
        f_osr = float(entry_7.get().strip())  # Oversampling Factor
        dtheta = float(entry_19.get().strip())  # Rotation Rate of the Antenna

        # Step 3 relevant data for X and Y axes:
        data_dict = {
            "f₁ primary beam taper factor": [f_1],
            "f_apo apodization smoothing factor": [f_apo],
            "f_osr oversampling factor": [f_osr],
            "δθ rotation rate of the antenna": [dtheta],
            "t_int integration time": [t_int],
            "t_map mapping time": [t_map]
        }

        # Set initial values for X and Y labels
        current_x_label = "f₁ primary beam taper factor"
        current_y_label = "t_int integration time"

        # Open the graph window with the data dictionary and default axis labels
        open_graph_window(3, data_dict, current_x_label, current_y_label)

    except ValueError:
        messagebox.showerror("Input Error", "Please enter valid numbers for the fields in Step 3.")

def plot_step4(theta_ext, theta_b, theta_sr, theta_ss):
    if theta_ext is not None:
        # Prepare the data dictionary with only the required values for Step 4
        data_dict = {
            "theta_ext": [theta_ext],
            "theta_b": [theta_b],
            "theta_sr": [theta_sr],
            "theta_ss": [theta_ss]
        }

        # Set initial values for X and Y labels
        current_x_label = "theta_ext"
        current_y_label = "theta_b"

        # Open the graph window with the filtered data dictionary and default axis labels
        open_graph_window(4, data_dict, current_x_label, current_y_label)

def plot_step5(theta_point, N_row, SNR):
    # Retrieve the delta_z value from the entry box
    delta_z = float(entry_16.get().strip())
    
    # Prepare the data dictionary, including delta_z (Surface Deformation)
    data_dict = {
        "theta_point": [theta_point],
        "N_row": [N_row],
        "SNR": [SNR],
        "δz surface deformation": [delta_z]  # Add this value to the data dictionary
    }
    
    # Open the graph window with the new data_dict, and set default x and y labels
    open_graph_window(5, data_dict, "theta_point", "SNR")

def plot_step6(P, D_t, D_a, z, T_sys, B):
    if P is not None:
        # Prepare the data dictionary with Step 6 values
        data_dict = {
            "Transmitter diameter (D_t)": [D_t],
            "Aperture diameter (D_a)": [D_a],
            "Distance between T & R (z)": [z],
            "System temperature (T_sys)": [T_sys],
            "Bandwidth (B)": [B],
            "Transmitter output power (P)": [P]
        }

        # Set initial values for X and Y labels
        current_x_label = "Transmitter diameter (D_t)"
        current_y_label = "Transmitter output power (P)"

        # Open the graph window with the filtered data dictionary and default axis labels
        open_graph_window(6, data_dict, current_x_label, current_y_label)

# Step 1 Calculate button action with refined input validation
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

        # Display the new results and store the item IDs in the list
        result_text_items.append(canvas.create_text(70.0, 190.0, anchor="nw", text=f"{R_F:.2f} m", fill="#000000", font=("Inter Medium", 13)))
        result_text_items.append(canvas.create_text(80.0, 210.0, anchor="nw", text=f"{R_min:.2f} m", fill="#000000", font=("Inter Medium", 13)))
        result_text_items.append(canvas.create_text(90.0, 230.0, anchor="nw", text=f"{R_react:.2f} m", fill="#000000", font=("Inter Medium", 13)))

    except ValueError:
        # Show a warning message if inputs are invalid
        messagebox.showwarning(
            "Input Error",
            "Please enter valid numerical values for D and ν (Frequency)."
        )

def open_help_window_1():
    # Create a new window for the help content
    new_window = Toplevel()
    new_window.title("Step 1 Help Window")
    new_window.geometry("900x600")  # Adjust the size as needed
    new_window.configure(bg="#1f1f1f")  # Set the background color

    new_window.resizable(False, False)

    # Load and display the help image
    help_image = PhotoImage(file=relative_to_assets("/Users/manalsiddiqui/Desktop/Short Gui/build/assets/frame0/Step1_Help.png"))
    help_label = Label(new_window, image=help_image)
    help_label.image = help_image  # Keep a reference to avoid garbage collection
    help_label.pack(expand=True)  # Center the image

    new_window.mainloop()


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
        step2_text_items.append(canvas.create_text(525.0, 190.0, anchor="nw", text=f"{delta_d:.2f} cm", fill="#000000", font=("Inter Medium", 13)))

    except ValueError:
        # Show a warning message if the input is invalid
        messagebox.showwarning(
            "Input Error",
            "Please enter a valid numerical value for a."
        )


def open_help_window_2():
    new_window = Toplevel()  
    new_window.title("Step 2 Help Window") 
    new_window.geometry("900x600")  
    new_window.configure(bg="#1f1f1f")

    new_window.resizable(False, False)

    help_image = PhotoImage(file=relative_to_assets("/Users/manalsiddiqui/Desktop/Short Gui/build/assets/frame0/Step2_help.png"))  
    help_label = Label(new_window, image=help_image)
    help_label.image = help_image  
    help_label.pack(expand=True)  

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
        step3_text_items.append(canvas.create_text(775.0, 235.0, anchor="nw", text=f"{t_int:.2f} s", fill="#000000", font=("Inter SemiBold", 13)))
        step3_text_items.append(canvas.create_text(905.0, 235.0, anchor="nw", text=f"{t_map:.2f} hr", fill="#000000", font=("Inter SemiBold", 13)))

    except ValueError:
        # Show a warning message if the input is invalid
        messagebox.showwarning(
            "Input Error",
            "Please enter valid numerical values for Step 3."
        )

def open_help_window_3(): 
    new_window = Toplevel()
    new_window.title("Step 3 Help Window")
    new_window.geometry("900x600")
    new_window.configure(bg= "#1f1f1f")

    new_window.resizable(False,False)

    help_image = PhotoImage(file=relative_to_assets("/Users/manalsiddiqui/Desktop/Short Gui/build/assets/frame0/step3_help.png"))
    help_label = Label(new_window, image=help_image)
    help_label.image = help_image 
    help_label.pack(expand=True)

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

        # Display the new results and store the item IDs in the list
        step4_text_items.append(canvas.create_text(160.0, 410.0, anchor="nw", text=f"{theta_ext:.2f} deg", fill="#000000", font=("Inter Medium", 13)))
        step4_text_items.append(canvas.create_text(160.0, 430.0, anchor="nw", text=f"{theta_b:.2f} arcsec", fill="#000000", font=("Inter Medium", 13)))
        step4_text_items.append(canvas.create_text(160.0, 450.0, anchor="nw", text=f"{theta_sr:.2f} arcsec", fill="#000000", font=("Inter Medium", 13)))
        step4_text_items.append(canvas.create_text(160.0, 470.0, anchor="nw", text=f"{theta_ss:.2f} arcsec", fill="#000000", font=("Inter Medium", 13)))

    except ValueError:
        # Show a warning message if inputs are invalid
        messagebox.showwarning(
            "Input Error",
            "Please enter valid numerical values for frequency (ν), D, and f_oss in Step 4."
        )


def open_help_window_4(): 
    new_window = Toplevel()
    new_window.title("Step 4 Help Window")
    new_window.geometry("900x600")
    new_window.configure(bg= "#1f1f1f")

    new_window.resizable(False,False)

    help_image = PhotoImage(file=relative_to_assets("/Users/manalsiddiqui/Desktop/Short Gui/build/assets/frame0/Step4_Help.png"))
    help_label = Label(new_window, image=help_image)
    help_label.image = help_image 
    help_label.pack(expand=True)


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
        step5_text_items.append(canvas.create_text(500.0, 420.0, anchor="nw", text=f"{theta_point:.2f} deg", fill="#000000", font=("Inter Medium", 13)))
        step5_text_items.append(canvas.create_text(500.0, 440.0, anchor="nw", text=f"{N_row:.2f}", fill="#000000", font=("Inter Medium", 13)))
        step5_text_items.append(canvas.create_text(500.0, 460.0, anchor="nw", text=f"{SNR:.2f} dB", fill="#000000", font=("Inter Medium", 13)))

    except ValueError:
        # Show a warning message if inputs are invalid
        messagebox.showwarning(
            "Input Error",
            "Please enter a valid numerical value for delta_z in Step 5."
        )

def open_help_window_5(): 
    new_window = Toplevel()
    new_window.title("Step 5 Help Window")
    new_window.geometry("900x600")
    new_window.configure(bg= "#1f1f1f")

    new_window.resizable(False,False)

    help_image = PhotoImage(file=relative_to_assets("/Users/manalsiddiqui/Desktop/Short Gui/build/assets/frame0/Step5_Help.png"))
    help_label = Label(new_window, image=help_image)
    help_label.image = help_image 
    help_label.pack(expand=True)

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
                900.0, 439.0, anchor="nw", text=f"{result:.4e} W", fill="#000000", font=("Inter Medium", 13)
            )
        )

    except ValueError:
        messagebox.showwarning(
            "Input Error",
            "Please enter valid numerical values for Step 6 inputs."
        )

def open_help_window_6(): 
    new_window = Toplevel()
    new_window.title("Step 6 Help Window")
    new_window.geometry("900x600")
    new_window.configure(bg= "#1f1f1f")

    new_window.resizable(False,False)

    help_image = PhotoImage(file=relative_to_assets("x"))
    help_label = Label(new_window, image=help_image)
    help_label.image = help_image 
    help_label.pack(expand=True)

window = Tk()

icon_image = PhotoImage(file="/Users/manalsiddiqui/Desktop/Short Gui/build/assets/frame0/HoloSim_logo.png")
window.iconphoto(False, icon_image)

window.geometry("1000x700")
window.configure(bg="#808080")
window.title("HoloSim")

canvas = Canvas(
    window,
    bg="#808080",
    height=700,
    width=1000,
    bd=0,
    highlightthickness=0,
    relief="ridge"
)

canvas.place(x=0, y=0)
canvas.create_text(
    13.0,
    9.0,
    anchor="nw",
    text="HoloSim.",
    fill="#FFFFFF",
    font=("Inter", 36, "bold")
)


image_image_1 = PhotoImage(file=relative_to_assets("image_1.png"))
image_1 = canvas.create_image(164.0, 164.0, image=image_image_1)

image_image_2 = PhotoImage(file=relative_to_assets("image_2.png"))
image_2 = canvas.create_image(501.0, 164.0, image=image_image_2)

entry_image_1 = PhotoImage(file=relative_to_assets("entry_1.png"))
entry_bg_1 = canvas.create_image(614.0, 105.5, image=entry_image_1)
entry_1 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_1.place(x=587.5, y=97.0, width=53.0, height=15.0)

image_image_3 = PhotoImage(file=relative_to_assets("image_3.png"))
image_3 = canvas.create_image(838.0, 164.0, image=image_image_3)

image_image_4 = PhotoImage(file=relative_to_assets("image_4.png"))
image_4 = canvas.create_image(164.0, 401.0, image=image_image_4)

image_image_5 = PhotoImage(file=relative_to_assets("image_5.png"))
image_5 = canvas.create_image(499.0, 401.0, image=image_image_5)

image_image_6 = PhotoImage(file=relative_to_assets("image_6.png"))
image_6 = canvas.create_image(838.0, 398.0, image=image_image_6)

canvas.create_rectangle(0.0, 525.0, 1003.0, 700.0, fill="#B6B9C9", outline="")

canvas.create_text(24.0, 65.0, anchor="nw", text="Step 1: Transmitter Placement", fill="#000000", font=("Inter Medium", 15 * -1))
canvas.create_text(24.0, 190.0, anchor="nw", text="R_F = ", fill="#000000", font=("Inter Medium", 13))
canvas.create_text(24.0, 210.0, anchor="nw", text="R_min = ", fill="#000000", font=("Inter Medium", 13))
canvas.create_text(24.0, 230.0, anchor="nw", text="R_react = ", fill="#000000", font=("Inter Medium", 13))
canvas.create_text(367.0, 65.0, anchor="nw", text="Step 2: Spatial Resolution Δd", fill="#000000", font=("Inter Medium", 15 * -1))
canvas.create_text(367.0, 190.0, anchor="nw", text="Spatial Resolution (Δd) < ", fill="#000000", font=("Inter Medium", 13))
canvas.create_text(699.0, 65.0, anchor="nw", text="Step 3: Grid Point Integration Time", fill="#000000", font=("Inter Medium", 15 * -1))
canvas.create_text(24.0, 302.0, anchor="nw", text="Step 4: Map Sampling Intervals θₛᵣ & θₛₛ", fill="#000000", font=("Inter Medium", 15 * -1))
canvas.create_text(70.0, 410.0, anchor="nw", text= "theta_ext = ", fill = "#000000", font=("inter Medium", 13 ))
canvas.create_text(70.0, 430.0, anchor="nw", text= "theta_b = ", fill = "#000000", font=("inter Medium", 13 ))
canvas.create_text(70.0, 450.0, anchor="nw", text= "theta_sr = ", fill = "#000000", font=("inter Medium", 13 ))
canvas.create_text(70.0, 470.0, anchor="nw", text= "theta_ss = ", fill = "#000000", font=("inter Medium", 13 ))
canvas.create_text(367.0, 307.0, anchor="nw", text="Step 5: Pointing Accuracy & SNR", fill="#000000", font=("Inter Medium", 15 * -1))
canvas.create_text(698.0, 294.0, anchor="nw", text="Step 6: Transmitter output power P", fill="#000000", font=("Inter Medium", 15 * -1))
canvas.create_text(380.0, 342.0, anchor="nw", text="δz surface deformation:", fill="#000000", font=("Inter Medium", 13 * -1)) 
canvas.create_text(21.0, 97.0, anchor="nw", text="D (Diameter in meters):\n\n", fill="#000000", font=("Inter SemiBold", 13 * -1))
canvas.create_text(357.0, 96.0, anchor="nw", text="a (Aperture in square meters):", fill="#000000", font=("Inter SemiBold", 13 * -1))

entry_image_2 = PhotoImage(file=relative_to_assets("entry_2.png"))
entry_bg_2 = canvas.create_image(241.0, 105.0, image=entry_image_2)
entry_2 = Entry(
    bd=0,
    bg="#D9D9D9",
    fg="#000716",
    highlightthickness=0
)
entry_2.place(x=202.0, y=96.0, width=78.0, height=16.0)

entry_image_3 = PhotoImage(file=relative_to_assets("entry_3.png"))
entry_bg_3 = canvas.create_image(944.0, 96.5, image=entry_image_3)
entry_3 = Entry(
    bd=0,
    bg="#D9D9D9",
    fg="#000716",
    highlightthickness=0
)
entry_3.place(x=905.5, y=87.0, width=77.0, height=17.0) # first entry for step 3 

entry_image_5 = PhotoImage(file=relative_to_assets("entry_5.png")) #step 1 entry 1 
entry_bg_5 = canvas.create_image(241.0, 129.5, image=entry_image_5)
entry_5 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_5.place(x=202.5, y=120.0, width=77.0, height=17.0)

entry_image_6 = PhotoImage(file=relative_to_assets("entry_6.png")) # middle step 3 
entry_bg_6 = canvas.create_image(944.0, 121.0, image=entry_image_6)
entry_6 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_6.place(x=905.0, y=112.0, width=78.0, height=16.0)
entry_6.insert(0, "1.3") # Sets the default value to 1.3 on startup

entry_image_19 = PhotoImage(file=relative_to_assets("entry_19.png"))
entry_bg_19 = canvas.create_image(944.0, 170.0, image=entry_image_19)
entry_19 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_19.place(x=905.5, y=162.0, width=78.0, height=16.0)  # last step 3

entry_image_7 = PhotoImage(file=relative_to_assets("entry_7.png")) # bottom step 3 
entry_bg_7 = canvas.create_image(944.0, 144.5, image=entry_image_7)
entry_7 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_7.place(x=905.5, y=136.0, width=77.0, height=17.0)

entry_image_12 = PhotoImage(file=relative_to_assets("entry_12.png"))
entry_bg_12 = canvas.create_image(270.0, 335.5, image=entry_image_12)
entry_12 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_12.place(x=230.5, y=328.0, width=77.0, height=17.0) # step 4 entry (Use this one)

entry_image_16 = PhotoImage(file=relative_to_assets("entry_16.png"))
entry_bg_16 = canvas.create_image(580.0, 350.0, image=entry_image_16)
entry_16 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_16.place(x=540.5, y=341.0, width=77.0, height=17.0) # entry for step 5

# Entry 17
entry_image_17 = PhotoImage(file=relative_to_assets("entry_17.png"))
entry_bg_17 = canvas.create_image(948.0, 326.5, image=entry_image_17)  # Shifted 10 pixels up
entry_17 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_17.place(x=909.5, y=317.0, width=77.0, height=17.0)  # Shifted 10 pixels up (D_t entry)

# Entry 18
entry_image_18 = PhotoImage(file=relative_to_assets("entry_18.png"))
entry_bg_18 = canvas.create_image(948.0, 351.5, image=entry_image_18)  # Shifted 10 pixels up
entry_18 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_18.place(x=909.5, y=342.0, width=77.0, height=17.0)  # Shifted 10 pixels up (D_a entry)

# Entry 8
entry_image_8 = PhotoImage(file=relative_to_assets("entry_8.png"))
entry_bg_8 = canvas.create_image(948.0, 376.5, image=entry_image_8)  # Shifted 10 pixels up
entry_8 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_8.place(x=909.5, y=367.0, width=77.0, height=17.0)  # Shifted 10 pixels up (T & R (z)) entry

# Entry 9
entry_image_9 = PhotoImage(file=relative_to_assets("entry_9.png"))
entry_bg_9 = canvas.create_image(948.0, 401.5, image=entry_image_9)  # Shifted 10 pixels up
entry_9 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_9.place(x=909.5, y=392.0, width=77.0, height=17.0)  # Shifted 10 pixels up (T_sys entry)

# Entry 10
entry_image_10 = PhotoImage(file=relative_to_assets("entry_10.png"))
entry_bg_10 = canvas.create_image(948.0, 426.5, image=entry_image_10)  # Shifted 10 pixels up
entry_10 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_10.place(x=909.5, y=417.0, width=77.0, height=17.0)  # Shifted 10 pixels up (B entry)


button_image_1 = PhotoImage(file=relative_to_assets("button_1.png"))
button_1 = Button(image=button_image_1, borderwidth=0, highlightthickness=0, command=lambda: print("button_1 clicked"), relief="flat", bg="#B6B9C9", activebackground="#B6B9C9")
button_1.place(x=500.2979125976562, y=600.8173217773438, width=250.48828125, height=70.5707778930664)


#Calc all text here
canvas.create_text(50, 550, anchor="nw", text="R_f (Meters) = ", fill = "#000000", font=("Inter Medium", 15))
canvas.create_text(50, 570, anchor="nw", text="R_f (Meters) = ", fill = "#000000", font=("Inter Medium", 15))
canvas.create_text(50, 590, anchor="nw", text="R_react (Meters) = ", fill = "#000000", font=("Inter Medium", 15))

canvas.create_text(50, 610, anchor="nw", text="Spatial Reslution (Δd) < ", fill="#000000", font=("Inter Medium", 15))

canvas.create_text(50, 630, anchor="nw", text="t_int (Seconds) = ", fill="#000000", font=("Inter Medium", 15))
canvas.create_text(50, 650, anchor="nw", text="t_map (Hours) = ", fill="#000000", font=("Inter Medium", 15))

canvas.create_text(270, 550, anchor="nw", text="theta_ext (Degrees)=", fill="#000000", font=("Inter Medium", 15))
canvas.create_text(270, 570, anchor="nw", text="theta_b (arcsec) = ", fill="#000000", font=("Inter Medium", 15))
canvas.create_text(270, 590, anchor="nw", text="theta_sr (arcsec) = ", fill="#000000", font=("Inter Medium", 15))
canvas.create_text(270, 610, anchor="nw", text="theta_ss (arcsec) = ", fill="#000000", font=("Inter Medium", 15))

canvas.create_text(270, 630, anchor="nw", text="theta_point (Degrees) <  ", fill="#000000", font=("Inter Medium", 15))
canvas.create_text(270, 650, anchor="nw", text="N_row = ", fill="#000000", font=("Inter Medium", 15))
canvas.create_text(500, 550, anchor="nw", text="SNR (dB) >  ", fill="#000000", font=("Inter Medium", 15))

canvas.create_text(500, 570, anchor="nw", text="Transmitter output power (W) = ", fill="#000000", font=("Inter Medium", 15))

button_image_2 = PhotoImage(file=relative_to_assets("button_2.png"))
button_2 = Button(
    image=button_image_2,
    borderwidth=0,
    highlightthickness=0,
    command=calculate_and_display_step1,  # Connect the updated handler
    relief="flat"
)
button_2.place(
    x=54.0,
    y=150.0,
    width=96.98621368408203,
    height=33.50228500366211
)

button_image_3 = PhotoImage(file=relative_to_assets("button_3.png"))
button_3 = Button(
    image=button_image_3,
    borderwidth=0,
    highlightthickness=0,
    command=calculate_and_display_step4,  # Connect Button 3 to the calculate_and_display_step4 function
    relief="flat"
)
button_3.place(x=67.35586547851562, y=360.26483154296875, width=71.3999252319336, height=26.21917724609375)

button_image_4 = PhotoImage(file=relative_to_assets("button_4.png"))
button_4 = Button(image=button_image_4, borderwidth=0, highlightthickness=0, command=calculate_and_display_step5, relief="flat")
button_4.place(x=396.96966552734375, y=380.26483154296875, width=71.39994812011719, height=26.21917724609375)

button_image_5 = PhotoImage(file=relative_to_assets("button_5.png"))
button_5 = Button(
    image=button_image_5, 
    borderwidth=0, 
    highlightthickness=0, 
    command=calculate_and_display_step2,  
    relief="flat"
)
button_5.place(x=404.0, y=126.0, width=96.25379180908203, height=33.50228500366211)

button_image_6 = PhotoImage(file=relative_to_assets("button_6.png"))
button_6 = Button(
    image=button_image_6,
    borderwidth=0,
    highlightthickness=0,
    command=calculate_and_display_step3,  # Connect Button 6 to the calculate_and_display_step3 function
    relief="flat"
)
button_6.place(x=730.8827514648438, y=185.96347045898438, width=95.98621368408203, height=33.50228500366211)

button_image_7 = PhotoImage(file=relative_to_assets("button_7.png"))
button_7 = Button(
    image=button_image_7, 
    borderwidth=0, 
    highlightthickness=0, 
    command=calculate_and_display_step6,  # Connect to Step 6 function
    relief="flat"
)
button_7.place(x=723.0, y=458.0, width=95.25379180908203, height=33.50228500366211)


button_image_8 = PhotoImage(file=relative_to_assets("button_8.png"))
button_8 = PhotoImage(file=relative_to_assets("button_8.png"))
button_8 = Button(
    image=button_image_8, 
    borderwidth=0, 
    highlightthickness=0, 
    command=lambda: [print("Step 1 Graph button pressed"), plot_step1(float(entry_2.get()), float(entry_5.get()))],
    relief="flat"
)
button_8.place(x=164.0, y=150.0, width=93.15172576904297, height=33.50228500366211)

button_image_9 = PhotoImage(file=relative_to_assets("button_9.png"))
button_9 = Button(
    image=button_image_9, 
    borderwidth=0, 
    highlightthickness=0, 
    command=lambda: [print("Step 2 Graph button pressed"), plot_step2(float(entry_1.get()), calculate_step2(float(entry_1.get())))],  
    relief="flat"
)
button_9.place(x=506.0, y=126.0, width=93.15172576904297, height=33.50228500366211)

button_image_10 = PhotoImage(file=relative_to_assets("button_10.png"))
button_10 = Button(
    image=button_image_10,
    borderwidth=0,
    highlightthickness=0,
    command=lambda: [print("Step 3 Graph button pressed"), plot_step3(*calculate_step3(calculate_step2(float(entry_1.get().strip()))), calculate_step2(float(entry_1.get().strip())))],
    relief="flat"
)
button_10.place(x=848.397216796875, y=185.96347045898438, width=93.15172576904297, height=33.50228500366211)

button_image_11 = PhotoImage(file=relative_to_assets("button_11.png"))
button_11 = Button(
    image=button_image_11, 
    borderwidth=0, 
    highlightthickness=0, 
    command=lambda: [print("Step 4 Graph button pressed"), plot_step4(*calculate_step4(float(entry_1.get())))],
    relief="flat"
)
button_11.place(x=146.17654418945312, y=360.26483154296875, width=69.76347351074219, height=26.21917724609375)

button_image_12 = PhotoImage(file=relative_to_assets("button_12.png"))
button_12 = Button(
    image=button_image_12, 
    borderwidth=0, 
    highlightthickness=0, 
    command=lambda: [print("Step 5 Graph button pressed"), plot_step5(*calculate_step5(calculate_step2(float(entry_1.get()))))],  
    relief="flat"
)
button_12.place(x=487.2551574707031, y=380.26483154296875, width=69.76347351074219, height=26.21917724609375)

button_image_13 = PhotoImage(file=relative_to_assets("button_13.png"))
button_13 = Button(
    image=button_image_13, 
    borderwidth=0, 
    highlightthickness=0, 
    command=lambda: [print("Step 6 Graph button pressed"), 
                     plot_step6(calculate_step6(),float(entry_17.get()),float(entry_18.get()), float(entry_8.get()),float(entry_9.get()),float(entry_10.get()))],  
    relief="flat"
)



canvas.create_text(690.0, 90.0, anchor="nw", text="f_1 primary beam taper factor: ", fill="#000000", font=("Inter SemiBold", 12))
canvas.create_text(690.0, 112.0, anchor="nw", text="f_apo apodization smoothing factor: ", fill="#000000", font=("Inter SemiBold", 11))
canvas.create_text(688.0, 135.0, anchor="nw", text="f_osr oversampling factor:", fill="#000000", font=("Inter SemiBold", 12))
canvas.create_text(688.0, 160.0, anchor="nw", text="δθ rotation rate of the antenna:", fill="#000000", font=("Inter SemiBold", 12))
canvas.create_text(730.0, 235.0, anchor="nw", text="t_int = ", fill="#000000", font=("Inter SemiBold", 13))
canvas.create_text(850.0, 235.0, anchor="nw", text="t_map = ", fill="#000000", font=("Inter SemiBold", 13))
canvas.create_text(17.0, 330.0, anchor="nw", text="fₒₛₛ oversampling factor along row:", fill="#000000", font=("Inter SemiBold", 12 * -1))
canvas.create_text(21.0, 116.0, anchor="nw", text="ν (Frequency in GHz):", fill="#000000", font=("Inter SemiBold", 13 * -1))

canvas.create_text(690.0, 319.0, anchor="nw", text="Transmitter diameter (D_t): ", fill="#000000", font=("Inter SemiBold", 13 * -1))
canvas.create_text(690.0, 342.0, anchor="nw", text="Aperture diameter (D_a): ", fill="#000000", font=("Inter SemiBold", 13 * -1))
canvas.create_text(690.0, 365.0, anchor="nw", text="Distance between T & R (z): ", fill="#000000", font=("Inter SemiBold", 13 * -1))
canvas.create_text(690.0, 392.0, anchor="nw", text="System temperature (Tₛᵧₛ): ", fill="#000000", font=("Inter SemiBold", 13 * -1))
canvas.create_text(690.0, 415.0, anchor="nw", text="Bandwidth (B): ", fill="#000000", font=("Inter SemiBold", 13 * -1))
canvas.create_text(700.0, 439.0, anchor="nw", text="Transmitter output power = ", fill="#000000", font=("Inter SemiBold", 13 * -1))

canvas.create_text(400.0, 420.0, anchor="nw", text="theta_point < ", fill="#000000", font=("Inter Medium", 13))
canvas.create_text(400.0, 440.0, anchor="nw", text="N_row = ", fill="#000000", font=("Inter Medium", 13))
canvas.create_text(400.0, 460.0, anchor="nw", text="SNR > ", fill="#000000", font=("Inter Medium", 13))

button_image_14 = PhotoImage(file=relative_to_assets("button_14.png"))
button_14 = Button(image=button_image_14,borderwidth=0,highlightthickness=0,command=open_help_window_1, relief="flat")
button_14.place(x=295.21929931640625, y=238.2511444091797, width=18.63034439086914, height=20.3926944732666)

button_image_15 = PhotoImage(file=relative_to_assets("button_15.png"))
button_15 = Button(image=button_image_15, borderwidth=0, highlightthickness=0, command=open_help_window_6, relief="flat")
button_15.place(x=967.3448486328125, y=474.2237548828125, width=18.63034439086914, height=20.3926944732666)

button_image_16 = PhotoImage(file=relative_to_assets("button_16.png"))
button_16 = Button(image=button_image_14, borderwidth = 0,highlightthickness=0, command=open_help_window_3, relief= "flat")
button_16.place(x=967.3448486328125, y=239.707763671875, width=18.63034439086914, height=20.3926944732666)

button_image_17 = PhotoImage(file=relative_to_assets("button_17.png"))
button_17 = Button(image=button_image_17, borderwidth=0, highlightthickness=0, command=open_help_window_5, relief="flat")
button_17.place(x=630.5654907226562, y=474.2237548828125, width=18.63034439086914, height=20.3926944732666)

button_image_18 = PhotoImage(file=relative_to_assets("button_18.png"))
button_18 = Button(image=button_image_18, borderwidth=0, highlightthickness=0, command=open_help_window_4, relief="flat")
button_18.place(x=295.21929931640625, y=472.76715087890625, width=18.63034439086914, height=20.3926944732666)

button_image_19 = PhotoImage(file=relative_to_assets("button_19.png"))
button_19 = Button(image=button_image_19, borderwidth=0, highlightthickness=0, command=open_help_window_2, relief="flat")
button_19.place(x=630.5654907226562, y=239.707763671875, width=18.63034439086914, height=20.3926944732666)
window.resizable(False, False)


window.mainloop()