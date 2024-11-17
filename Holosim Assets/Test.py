from pathlib import Path
from tkinter import Tk, Canvas, Entry, Button, PhotoImage, messagebox
import numpy as np

# Function for Step 1: Far-field and reactive near-field cutoffs
def calculate_step1(D, nu):
    c = 3 * 10**8  # Speed of light in m/s
    lam = c/(nu*10**9)  # Wavelength in meters
    R_F = (2*D**2)/lam  # Far-field cutoff in meters
    R_min = 5*D  # Minimum transmitter distance in meters
    R_react = 0.62*np.sqrt(D**3/lam)  # Reactive near-field cutoff in meters
    return R_F, R_min, R_react

# Function for Step 2: Spatial resolution
def calculate_step2(a):
    delta_d = a*np.sqrt(2)/2 *10**2  # Spatial resolution in cm
    return delta_d

# Function for Step 3: Integration time and map time
def calculate_step3(f_1, f_apo, f_osr, dtheta, nu, D, delta_d):
    t_int = 6.2*10**4*f_osr*f_1*f_apo**2 / (dtheta*nu*D)  # Integration time in seconds
    t_map = 171768 * f_osr * f_1 * f_apo**2 * D / (dtheta * nu * delta_d**2)  # Total map time in hours
    return t_int, t_map

# Function for Step 4: Angular extent, beamwidth, and sampling intervals
def calculate_step4(f_1, f_apo, f_osr, f_oss, delta_d, nu, D):
    theta_ext = 1717.7*f_1*f_apo/(delta_d*nu)  # Angular extent in degrees
    theta_b = 61836.6*f_1/(nu*D)  # Beamwidth in arcsec
    theta_sr = theta_b/f_osr  # Sampling interval between rows in arcsec
    theta_ss = theta_b/f_oss  # Sampling interval along scan in arcsec
    return theta_ext, theta_b, theta_sr, theta_ss

# Function for Step 5: On-axis power measurement (example placeholder)
def calculate_step5(D, delta_d, nu, delta_z, M_0):
    power = M_0 * (D ** 2) / (delta_d * delta_z * nu)  # Placeholder calculation
    return power

# Function for Step 6: Surface deformation and frequency effects (example placeholder)
def calculate_step6(delta_z, nu):
    result = delta_z * nu * 10**-6  # Placeholder calculation
    return result

OUTPUT_PATH = Path(__file__).parent
ASSETS_PATH = OUTPUT_PATH / Path(r"/Users/mustafa/Desktop/Short Gui/build/assets/frame0")

def relative_to_assets(path: str) -> Path:
    return ASSETS_PATH / Path(path)

window = Tk()
window.geometry("1000x700")
window.configure(bg="#808080")
window.title("Holosim")

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
    text="Holosim.",
    fill="#000000",
    font=("Inter ExtraBold", 36 * -1)
)

# Image placeholders
image_image_1 = PhotoImage(file=relative_to_assets("image_1.png"))
image_1 = canvas.create_image(164.0, 164.0, image=image_image_1)

image_image_2 = PhotoImage(file=relative_to_assets("image_2.png"))
image_2 = canvas.create_image(501.0, 164.0, image=image_image_2)

image_image_3 = PhotoImage(file=relative_to_assets("image_3.png"))
image_3 = canvas.create_image(838.0, 164.0, image=image_image_3)

image_image_4 = PhotoImage(file=relative_to_assets("image_4.png"))
image_4 = canvas.create_image(164.0, 401.0, image=image_image_4)

image_image_5 = PhotoImage(file=relative_to_assets("image_5.png"))
image_5 = canvas.create_image(499.0, 401.0, image=image_image_5)

image_image_6 = PhotoImage(file=relative_to_assets("image_6.png"))
image_6 = canvas.create_image(838.0, 398.0, image=image_image_6)

canvas.create_rectangle(0.0, 525.0, 1003.0, 700.0, fill="#B6B9C9", outline="")

# Step 1
canvas.create_text(24.0, 65.0, anchor="nw", text="Step 1:", fill="#000000", font=("Inter Medium", 15 * -1))

entry_image_1 = PhotoImage(file=relative_to_assets("entry_1.png"))
entry_bg_1 = canvas.create_image(614.0, 105.5, image=entry_image_1)
entry_1 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_1.place(x=587.5, y=97.0, width=53.0, height=15.0)

entry_image_2 = PhotoImage(file=relative_to_assets("entry_2.png"))
entry_bg_2 = canvas.create_image(614.0, 137.5, image=entry_image_2)
entry_2 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_2.place(x=587.5, y=130.0, width=53.0, height=15.0)

def calculate_and_display_step1():
    D = float(entry_1.get())
    nu = float(entry_2.get())
    R_F, R_min, R_react = calculate_step1(D, nu)
    messagebox.showinfo("Step 1 Result", f"R_F = {R_F:.2f} m\nR_min = {R_min:.2f} m\nR_react = {R_react:.2f} m")

button_1 = Button(window, image=PhotoImage(file=relative_to_assets("button_1.png")), command=calculate_and_display_step1, borderwidth=0, highlightthickness=0, relief="flat")
button_1.place(x=587.5, y=160.0, width=53.0, height=15.0)

# Step 2
canvas.create_text(363.0, 65.0, anchor="nw", text="Step 2:", fill="#000000", font=("Inter Medium", 15 * -1))

entry_image_3 = PhotoImage(file=relative_to_assets("entry_3.png"))
entry_bg_3 = canvas.create_image(614.0, 230.0, image=entry_image_3)
entry_3 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_3.place(x=587.5, y=230.0, width=53.0, height=15.0)

def calculate_and_display_step2():
    a = float(entry_3.get())
    delta_d = calculate_step2(a)
    messagebox.showinfo("Step 2 Result", f"δ_d = {delta_d:.2f} cm")

button_2 = Button(window, image=PhotoImage(file=relative_to_assets("button_2.png")), command=calculate_and_display_step2, borderwidth=0, highlightthickness=0, relief="flat")
button_2.place(x=587.5, y=260.0, width=53.0, height=15.0)

# Step 3
canvas.create_text(710.0, 65.0, anchor="nw", text="Step 3:", fill="#000000", font=("Inter Medium", 15 * -1))

entry_image_4 = PhotoImage(file=relative_to_assets("entry_4.png"))
entry_bg_4 = canvas.create_image(614.0, 330.0, image=entry_image_4)
entry_4 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_4.place(x=587.5, y=330.0, width=53.0, height=15.0)

entry_image_5 = PhotoImage(file=relative_to_assets("entry_5.png"))
entry_bg_5 = canvas.create_image(614.0, 360.0, image=entry_image_5)
entry_5 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_5.place(x=587.5, y=360.0, width=53.0, height=15.0)

entry_image_6 = PhotoImage(file=relative_to_assets("entry_6.png"))
entry_bg_6 = canvas.create_image(614.0, 390.0, image=entry_image_6)
entry_6 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_6.place(x=587.5, y=390.0, width=53.0, height=15.0)

def calculate_and_display_step3():
    D = float(entry_4.get())
    nu = float(entry_5.get())
    delta_d = float(entry_6.get())
    f_1, f_apo, f_osr, dtheta = 1.13, 1.3, 2.2, 300  # Example fixed values
    t_int, t_map = calculate_step3(f_1, f_apo, f_osr, dtheta, nu, D, delta_d)
    messagebox.showinfo("Step 3 Result", f"t_int = {t_int:.2f} s\nt_map = {t_map:.2f} hours")

button_3 = Button(window, image=PhotoImage(file=relative_to_assets("button_3.png")), command=calculate_and_display_step3, borderwidth=0, highlightthickness=0, relief="flat")
button_3.place(x=587.5, y=420.0, width=53.0, height=15.0)

# Step 4
canvas.create_text(24.0, 298.0, anchor="nw", text="Step 4:", fill="#000000", font=("Inter Medium", 15 * -1))

entry_image_7 = PhotoImage(file=relative_to_assets("entry_7.png"))
entry_bg_7 = canvas.create_image(614.0, 490.0, image=entry_image_7)
entry_7 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_7.place(x=587.5, y=490.0, width=53.0, height=15.0)

entry_image_8 = PhotoImage(file=relative_to_assets("entry_8.png"))
entry_bg_8 = canvas.create_image(614.0, 520.0, image=entry_image_8)
entry_8 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_8.place(x=587.5, y=520.0, width=53.0, height=15.0)

entry_image_9 = PhotoImage(file=relative_to_assets("entry_9.png"))
entry_bg_9 = canvas.create_image(614.0, 550.0, image=entry_image_9)
entry_9 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_9.place(x=587.5, y=550.0, width=53.0, height=15.0)

def calculate_and_display_step4():
    D = float(entry_7.get())
    nu = float(entry_8.get())
    delta_d = float(entry_9.get())
    f_1, f_apo, f_osr, f_oss = 1.13, 1.3, 2.2, 15  # Example fixed values
    theta_ext, theta_b, theta_sr, theta_ss = calculate_step4(f_1, f_apo, f_osr, f_oss, delta_d, nu, D)
    messagebox.showinfo("Step 4 Result", f"θ_ext = {theta_ext:.2f} deg\nθ_b = {theta_b:.2f} arcsec\nθ_sr = {theta_sr:.2f} arcsec\nθ_ss = {theta_ss:.2f} arcsec")

button_4 = Button(window, image=PhotoImage(file=relative_to_assets("button_4.png")), command=calculate_and_display_step4, borderwidth=0, highlightthickness=0, relief="flat")
button_4.place(x=587.5, y=580.0, width=53.0, height=15.0)

# Step 5
canvas.create_text(363.0, 298.0, anchor="nw", text="Step 5:", fill="#000000", font=("Inter Medium", 15 * -1))

entry_image_10 = PhotoImage(file=relative_to_assets("entry_10.png"))
entry_bg_10 = canvas.create_image(614.0, 97.0, image=entry_image_10)
entry_10 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_10.place(x=700.0, y=97.0, width=53.0, height=15.0)

entry_image_11 = PhotoImage(file=relative_to_assets("entry_11.png"))
entry_bg_11 = canvas.create_image(614.0, 130.0, image=entry_image_11)
entry_11 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_11.place(x=700.0, y=130.0, width=53.0, height=15.0)

entry_image_12 = PhotoImage(file=relative_to_assets("entry_12.png"))
entry_bg_12 = canvas.create_image(614.0, 160.0, image=entry_image_12)
entry_12 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_12.place(x=700.0, y=160.0, width=53.0, height=15.0)

entry_image_13 = PhotoImage(file=relative_to_assets("entry_13.png"))
entry_bg_13 = canvas.create_image(614.0, 190.0, image=entry_image_13)
entry_13 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_13.place(x=700.0, y=190.0, width=53.0, height=15.0)

entry_image_14 = PhotoImage(file=relative_to_assets("entry_14.png"))
entry_bg_14 = canvas.create_image(614.0, 220.0, image=entry_image_14)
entry_14 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_14.place(x=700.0, y=220.0, width=53.0, height=15.0)

def calculate_and_display_step5():
    D = float(entry_10.get())
    delta_d = float(entry_11.get())
    nu = float(entry_12.get())
    delta_z = float(entry_13.get())
    M_0 = float(entry_14.get())
    power = calculate_step5(D, delta_d, nu, delta_z, M_0)
    messagebox.showinfo("Step 5 Result", f"Power = {power:.2f} W")  # Example placeholder output

button_5 = Button(window, image=PhotoImage(file=relative_to_assets("button_5.png")), command=calculate_and_display_step5, borderwidth=0, highlightthickness=0, relief="flat")
button_5.place(x=700.0, y=280.0, width=53.0, height=15.0)

# Step 6
canvas.create_text(710.0, 298.0, anchor="nw", text="Step 6:", fill="#000000", font=("Inter Medium", 15 * -1))

entry_image_15 = PhotoImage(file=relative_to_assets("entry_15.png"))
entry_bg_15 = canvas.create_image(614.0, 330.0, image=entry_image_15)
entry_15 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_15.place(x=700.0, y=330.0, width=53.0, height=15.0)

entry_image_16 = PhotoImage(file=relative_to_assets("entry_16.png"))
entry_bg_16 = canvas.create_image(614.0, 360.0, image=entry_image_16)
entry_16 = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
entry_16.place(x=700.0, y=360.0, width=53.0, height=15.0)

def calculate_and_display_step6():
    delta_z = float(entry_15.get())
    nu = float(entry_16.get())
    result = calculate_step6(delta_z, nu)
    messagebox.showinfo("Step 6 Result", f"Result = {result:.6f}")  # Example placeholder output

button_6 = Button(window, image=PhotoImage(file=relative_to_assets("button_6.png")), command=calculate_and_display_step6, borderwidth=0, highlightthickness=0, relief="flat")
button_6.place(x=700.0, y=420.0, width=53.0, height=15.0)

window.mainloop()
