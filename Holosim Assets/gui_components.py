# gui_components.py

import os
import tkinter as tk
from pathlib import Path
from tkinter import Tk, Canvas, Entry, Button, PhotoImage, messagebox, Toplevel, Label
from PIL import Image, ImageTk

class GUI:
    def __init__(self, window, relative_to_assets):
        self.window = window
        self.relative_to_assets = relative_to_assets

        self.window.geometry("1200x940")
        self.window.configure(bg="#808080")
        self.window.title("HoloSim")

        icon_image = PhotoImage(file=self.relative_to_assets("HS_logo.png"))
        self.window.iconphoto(False, icon_image)

        self.canvas = Canvas(
            self.window,
            bg="#808080",
            height=940,
            width=1200,
            bd=0,
            highlightthickness=0,
            relief="ridge"
        )
        self.canvas.place(x=0, y=0)

        # Title text
        self.canvas.create_text(
            13.0,
            15.0,
            anchor="nw",
            text="HoloSim.",
            fill="#FFFFFF",
            font=("Inter", 42, "bold")
        )

        # Now define helper methods similar to holosim.py but as class methods:
        # We'll call these methods to create entries/buttons/images.
        
        # We'll store references to entry and button images if needed to avoid garbage collection:
        self.images = []
        self.entry_images = []
        self.button_images = []
        
        # Load and display image helper method
    def load_and_display_image(self, image_path, x, y):
        image = PhotoImage(file=self.relative_to_assets(image_path))
        self.images.append(image)  # keep a reference
        self.canvas.create_image(x, y, image=image)
        return image

    def create_entry(self, image_path, x, y, width, height, default_text=""):
        entry_image = PhotoImage(file=self.relative_to_assets(image_path))
        self.entry_images.append(entry_image)
        self.canvas.create_image(x, y, image=entry_image)
        entry = Entry(bd=0, bg="#D9D9D9", fg="#000716", highlightthickness=0)
        entry.place(x=x - width / 2, y=y - height / 2, width=width, height=height)
        entry.insert(0, default_text)
        return entry

    def create_button(self, image_path, x, y, width, height, command):
        button_image = PhotoImage(file=self.relative_to_assets(image_path))
        self.button_images.append(button_image)
        btn = Button(image=button_image, borderwidth=0, highlightthickness=0, command=command, relief="flat")
        btn.place(x=x, y=y, width=width, height=height)
        return btn

    def build_gui(self, 
                  calculate_all_steps_cmd, calculate_and_display_step1_cmd, calculate_and_display_step2_cmd, 
                  calculate_and_display_step3_cmd, calculate_and_display_step4_cmd, calculate_and_display_step5_cmd, 
                  calculate_and_display_step6_cmd, open_help_window_cmd, open_Reference_Window_cmd):
       

        
        # images from holosim.py
        image_1 = self.load_and_display_image("Image_1.png", 200.0, 230.0)
        image_2 = self.load_and_display_image("Image_2.png", 600.0, 230.0)
        image_3 = self.load_and_display_image("Image_3.png", 1000.0, 230.0)
        image_4 = self.load_and_display_image("Image_4.png", 200.0, 560.0)
        image_5 = self.load_and_display_image("Image_5.png", 600.0, 560.0)
        image_6 = self.load_and_display_image("Image_6.png", 1000.0, 560.0)
        image_Line_1 = self.load_and_display_image("Line_1.png", 550.0, 920.0)
        image_Line_2 = self.load_and_display_image("Line_2.png", 600.0, 790.0)

        # Bottom rectangle
        rect = self.canvas.create_rectangle(0.0, 725.0, 1200.0, 976.0, fill="#B6B9C9", outline="")

        # Step 1 latex
        image_D = self.load_and_display_image("D.png", 185.0, 159.0)
        image_nu = self.load_and_display_image("nu.png", 205.0, 180.0)
        image_R_F = self.load_and_display_image("R_F.png", 100.0, 219.0)
        image_R_min = self.load_and_display_image("R_min.png", 105.0, 239.0)
        image_R_react = self.load_and_display_image("R_react.png", 107.0, 259.0)

        # Step 2 latex
        image_a = self.load_and_display_image("a.png", 654.0, 159.0)
        image_small_delt_d = self.load_and_display_image("small_delt_d.png", 550.0, 239.0)

        # Step 3 latex
        image_f_1 = self.load_and_display_image("f_1.png", 988.0, 158.0)
        image_f_apo = self.load_and_display_image("f_apo.png", 980.0, 178.0)
        image_f_osr = self.load_and_display_image("f_osr.png", 1025.0, 203.0)
        image_Theta_dot = self.load_and_display_image("Theta_dot.png", 998.0, 220.0)
        image_t_int = self.load_and_display_image("t_int.png", 968.0, 260.0)
        image_t_map = self.load_and_display_image("t_map.png", 970.0, 285.0)

        # Step 4 latex
        image_f_oss = self.load_and_display_image("f_oss.png", 234.0, 490.0)
        image_step4_outputs = self.load_and_display_image("step4_outputs.png", 100.0, 560.0)

        # Step 5 latex
        image_small_delta_z = self.load_and_display_image("small_delta_z.png", 621.0, 490.0)
        image_Theta_Point_less = self.load_and_display_image("Theta_Point_less.png", 540.0, 538.0)
        image_N_Row_Eq = self.load_and_display_image("N_Row_Eq.png", 535.0, 563.0)
        image_SNR_Less = self.load_and_display_image("SNR_Less.png", 534.0, 588.0)

        # Step 6 latex
        image_D_t = self.load_and_display_image("D_t.png", 1000.0, 491.0)
        image_T_sys = self.load_and_display_image("T_sys.png", 997.0, 530.0)
        image_Z_Final = self.load_and_display_image("Z_Final.png", 975.0, 510.0)
        image_B = self.load_and_display_image("B.png", 990.0, 551.0)

        # inputs latex
        image_Input_D = self.load_and_display_image("Input_D.png", 30.0, 820.0)
        image_Input_nu = self.load_and_display_image("Input_nu.png", 30.0, 850.0)
        image_Input_a = self.load_and_display_image("Input_a.png", 30.0, 880.0)
        image_Input_f_1 = self.load_and_display_image("Input_f_1.png", 30.0, 910.0)
        image_Input_f_apo = self.load_and_display_image("Input_f_apo.png", 180.0, 820.0)
        image_Input_f_osr = self.load_and_display_image("Input_f_osr.png", 180.0, 850.0)
        image_Input_theta_dot = self.load_and_display_image("Input_theta_dot.png", 170.0, 880.0)
        image_Input_f_oss = self.load_and_display_image("Input_f_oss.png", 180.0, 910.0)
        image_Input_delta_z = self.load_and_display_image("Input_delta_z.png", 330.0, 820.0)
        image_Input_D_t = self.load_and_display_image("Input_D_t.png", 330.0, 850.0)
        image_Input_z = self.load_and_display_image("Input_z.png", 324.0, 880.0)
        image_Input_T_sys = self.load_and_display_image("Input_T_sys.png", 330.0, 910.0)
        image_Input_B = self.load_and_display_image("Input_B.png", 445.0, 820.0)

        self.canvas.create_text(10.0, 730.0, anchor="nw", text="Final Parameters", fill="#000000", font=("Inter", 20, "bold"))
        self.canvas.create_text(15.0, 765.0, anchor="nw", text="Inputs:", fill="#000000", font=("Inter Medium", 17))
        self.canvas.create_text(570.0, 765.0, anchor="nw", text="Outputs:", fill="#000000", font=("Inter Medium", 17))

        # outputs latex
        image_R_min_Equals = self.load_and_display_image("R_min_Equals.png", 607.0, 820.0)
        image_Delta_lower_d_less = self.load_and_display_image("Delta_lower_d_less.png", 595.0, 850.0)
        image_t_int_equals = self.load_and_display_image("t_int_equals.png", 596.0, 880.0)
        image_t_map_equals = self.load_and_display_image("t_map_equals.png", 601.0, 910.0)
        image_Theta_ext_Equals = self.load_and_display_image("Theta_ext_Equals.png", 830.0, 820.0)
        image_Theta_b_equals = self.load_and_display_image("Theta_b_equals.png", 823.0, 850.0)
        image_Theta_sr_equals = self.load_and_display_image("Theta_sr_equals.png", 825.0, 880.0)
        image_Theta_ss_equals = self.load_and_display_image("Theta_ss_equals.png", 825.0, 910.0)
        image_Theta_point_less = self.load_and_display_image("Theta_point_less.png", 1040.0, 820.0)
        image_N_row_equals = self.load_and_display_image("N_row_equals.png", 1035.0, 850.0)
        image_SNR = self.load_and_display_image("SNR.png", 1031.0, 880.0)

        self.canvas.create_text(1007.0, 895.0, anchor="nw", text="P >", fill="#000000", font=("Inter Medium", 20, "italic"))

        # Create text elements for steps
        self.create_text_elements()

        # Entries (use the same calls as holosim.py)
        self.entry_1 = self.create_entry("entry_1.png", 715.0, 160.0, 85.0, 15.0)
        self.entry_2 = self.create_entry("entry_2.png", 265.0, 158.0, 78.0, 16.0)
        self.entry_3 = self.create_entry("entry_3.png", 1095.0, 156.5, 77.0, 17.0)
        self.entry_5 = self.create_entry("entry_5.png", 265.0, 180.0, 77.0, 17.0)
        self.entry_6 = self.create_entry("entry_6.png", 1095.0, 180.0, 78.0, 16.0, default_text="1.3")
        self.entry_7 = self.create_entry("entry_7.png", 1095.0, 202.0, 77.0, 17.0)
        self.entry_8 = self.create_entry("entry_8.png", 1090.0, 510.5, 77.0, 17.0)
        self.entry_9 = self.create_entry("entry_9.png", 1090.0, 531.5, 77.0, 17.0)
        self.entry_10 = self.create_entry("entry_10.png", 1090.0, 552.5, 77.0, 17.0)
        self.entry_12 = self.create_entry("entry_12.png", 305.0, 489.0, 77.0, 17.0)
        self.entry_16 = self.create_entry("entry_16.png", 690.0, 490.0, 77.0, 17.0)
        self.entry_17 = self.create_entry("entry_17.png", 1090.0, 490.0, 77.0, 17.0)
        self.entry_19 = self.create_entry("entry_19.png", 1095.0, 224.0, 78.0, 16.0)

        # Buttons (pass commands that were given as parameters)
        button_image_1 = PhotoImage(file=self.relative_to_assets("button1.png"))
        self.button_1 = Button(image=button_image_1, borderwidth=0, highlightthickness=0, command=calculate_all_steps_cmd, relief="flat")
        self.button_1.place(x=1040, y=730, width=148, height=46)
        self.button_images.append(button_image_1)

        button_image_2 = PhotoImage(file=self.relative_to_assets("button_2.png"))
        self.button_2 = Button(image=button_image_2, borderwidth=0, highlightthickness=0, command=calculate_and_display_step1_cmd, relief="flat")
        self.button_2.place(x=110.0, y=300.0, width=135, height=50)
        self.button_images.append(button_image_2)

        button_image_3 = PhotoImage(file=self.relative_to_assets("button_3.png"))
        self.button_3 = Button(image=button_image_3, borderwidth=0, highlightthickness=0, command=calculate_and_display_step4_cmd, relief="flat")
        self.button_3.place(x=110, y=630, width=135, height=50)
        self.button_images.append(button_image_3)

        button_image_4 = PhotoImage(file=self.relative_to_assets("button_4.png"))
        self.button_4 = Button(image=button_image_4, borderwidth=0, highlightthickness=0, command=calculate_and_display_step5_cmd, relief="flat")
        self.button_4.place(x=530, y=630, width=135, height=50)
        self.button_images.append(button_image_4)

        button_image_5 = PhotoImage(file=self.relative_to_assets("button_5.png"))
        self.button_5 = Button(image=button_image_5, borderwidth=0, highlightthickness=0, command=calculate_and_display_step2_cmd,  relief="flat")
        self.button_5.place(x=530.0, y=300.0, width=135, height=50)
        self.button_images.append(button_image_5)

        button_image_6 = PhotoImage(file=self.relative_to_assets("button_6.png"))
        self.button_6 = Button(image=button_image_6, borderwidth=0, highlightthickness=0, command=calculate_and_display_step3_cmd, relief="flat")
        self.button_6.place(x=940.0, y=300, width=135, height=50)
        self.button_images.append(button_image_6)

        button_image_7 = PhotoImage(file=self.relative_to_assets("button_7.png"))
        self.button_7 = Button(image=button_image_7, borderwidth=0, highlightthickness=0, command=calculate_and_display_step6_cmd, relief="flat")
        self.button_7.place(x=940.0, y=630.0, width=135, height=50)
        self.button_images.append(button_image_7)

        Guide_Button_image = PhotoImage(file=self.relative_to_assets("Massive_Help_Button.png"))
        self.Guide_Button = Button(image=Guide_Button_image, borderwidth=0, highlightthickness=0, command=open_help_window_cmd, relief="flat")
        self.Guide_Button.place(x=1050.0, y=20.0, width=140, height=36)
        self.button_images.append(Guide_Button_image)

        button_image_14 = PhotoImage(file=self.relative_to_assets("button_14.png"))
        self.button_14 = Button(image=button_image_14,borderwidth=0,highlightthickness=0,command=open_Reference_Window_cmd, relief="flat")
        self.button_14.place(x=350, y=350, width=30, height=32)
        self.button_images.append(button_image_14)

        button_image_15 = PhotoImage(file=self.relative_to_assets("button_15.png"))
        self.button_15 = Button(image=button_image_15, borderwidth=0, highlightthickness=0, command=open_Reference_Window_cmd, relief="flat")
        self.button_15.place(x=1150.3448486328125, y=680.2237548828125, width=30, height=32)
        self.button_images.append(button_image_15)

        button_image_16 = PhotoImage(file=self.relative_to_assets("button_16.png"))
        self.button_16 = Button(image=button_image_14, borderwidth=0, highlightthickness=0, command=open_Reference_Window_cmd, relief="flat")
        self.button_16.place(x=1150.3448486328125, y=350, width=30, height=32)
        

        button_image_17 = PhotoImage(file=self.relative_to_assets("button_17.png"))
        self.button_17 = Button(image=button_image_17, borderwidth=0, highlightthickness=0, command=open_Reference_Window_cmd, relief="flat")
        self.button_17.place(x=740, y=680, width=30, height=32)
        self.button_images.append(button_image_17)

        button_image_18 = PhotoImage(file=self.relative_to_assets("button_18.png"))
        self.button_18 = Button(image=button_image_18, borderwidth=0, highlightthickness=0, command=open_Reference_Window_cmd, relief="flat")
        self.button_18.place(x=350, y=680, width=30, height=32)
        self.button_images.append(button_image_18)

        button_image_19 = PhotoImage(file=self.relative_to_assets("button_19.png"))
        self.button_19 = Button(image=button_image_19, borderwidth=0, highlightthickness=0, command=open_Reference_Window_cmd, relief="flat")
        self.button_19.place(x=740, y=350, width=30, height=32)
        self.button_images.append(button_image_19)

      
        self.window.resizable(False, False)

    def create_text_elements(self):
        self.canvas.create_text(55.0, 100.0, anchor="nw", text="Step 1: Transmitter Placement", fill="#000000", font=("Inter Medium", 20))
        self.canvas.create_text(50.0, 150.0, anchor="nw", text="Aperture Diameter  \n\n", fill="#000000", font=("Inter SemiBold", 14))
        self.canvas.create_text(315, 150, anchor="nw", text="m", fill="#000000", font=("Inter SemiBold", 14))
        self.canvas.create_text(50.0, 170.0, anchor="nw", text="Transmitter frequency  ", fill="#000000", font=("Inter SemiBold", 14))
        self.canvas.create_text(315, 170, anchor="nw", text="GHz", fill="#000000", font=("Inter SemiBold", 14))

        self.canvas.create_text(485.0, 100.0, anchor="nw", text="Step 2: Spatial Resolution   ", fill="#000000", font=("Inter Medium", 20))
        self.canvas.create_text(420.0, 150.0, anchor="nw", text="Distance between corner adjustors    :", fill="#000000", font=("Inter SemiBold", 14))
        self.canvas.create_text(770, 150, anchor="nw", text="m", fill="#000000", font=("Inter SemiBold", 14))
        self.canvas.create_text(570.0, 230.0, anchor="nw", text="< ", fill="#000000", font=("Inter Medium", 15))

        self.canvas.create_text(830.0, 85.0, anchor="nw", text="Step 3: Grid Point Integration Time and", fill="#000000", font=("Inter Medium", 20))
        self.canvas.create_text(930, 110, anchor="nw", text="Total Map Time", fill="#000000", font=("Inter Medium", 20))
        self.canvas.create_text(820.0, 150.0, anchor="nw", text="Primary beam taper factor      : ", fill="#000000", font=("Inter SemiBold", 13))
        self.canvas.create_text(820.0, 170.0, anchor="nw", text="Apodization smoothing          :", fill="#000000", font=("Inter SemiBold", 13))
        self.canvas.create_text(820.0, 193.0, anchor="nw", text="Oversampling factor btwn rows        :", fill="#000000", font=("Inter SemiBold", 13))
        self.canvas.create_text(820.0, 215.0, anchor="nw", text="Rotation rate of the antenna     :", fill="#000000", font=("Inter SemiBold", 13))
        self.canvas.create_text(1145, 215.0, anchor="nw", text="arc/sec", fill="#000000", font=("Inter SemiBold", 13))

        self.canvas.create_text(45.0, 410.0, anchor="nw", text="Step 4: Map Angular Extent and", fill="#000000", font=("Inter Medium", 20))
        self.canvas.create_text(105, 435, anchor="nw", text="Sampling Intervals", fill="#000000", font=("Inter Medium", 20))
        self.canvas.create_text(20.0, 480.0, anchor="nw", text="Oversampling factor along row        :", fill="#000000", font=("Inter SemiBold", 14))

        self.canvas.create_text(460.0, 410.0, anchor="nw", text="Step 5: Pointing Accuracy and", fill="#000000", font=("Inter Medium", 20))
        self.canvas.create_text(520, 435, anchor='nw', text="SNR Requirement", fill="#000000", font=("Inter Medium", 20))
        self.canvas.create_text(463.0, 480.0, anchor="nw", text="Surface deformation       :", fill="#000000", font=("Inter Medium", 15))
        self.canvas.create_text(740, 480, anchor="nw", text="μm", fill="#000000", font=("Inter Medium", 15))

        self.canvas.create_text(850.0, 430.0, anchor="nw", text="Step 6: Transmitter output power", fill="#000000", font=("Inter Medium", 20))
        self.canvas.create_text(850.0, 482.0, anchor="nw", text="Transmitter diameter      : ", fill="#000000", font=("Inter SemiBold", 14))
        self.canvas.create_text(1140, 482.0, anchor="nw", text="m", fill="#000000", font=("Inter SemiBold", 14))
        self.canvas.create_text(850.0, 500.0, anchor="nw", text="Distance between", fill="#000000", font=("Inter SemiBold", 14))
        self.canvas.create_text(1140, 500.0, anchor="nw", text="m", fill="#000000", font=("Inter SemiBold", 14))
        self.canvas.create_text(980.0, 500.0, anchor="nw", text=" :", fill="#000000", font=("Inter SemiBold", 14))
        self.canvas.create_text(850.0, 520.0, anchor="nw", text="System temperature        : ", fill="#000000", font=("Inter SemiBold", 14))
        self.canvas.create_text(1140, 520.0, anchor="nw", text="K", fill="#000000", font=("Inter SemiBold", 14))
        self.canvas.create_text(847.0, 542.0, anchor="nw", text=" Detector Bandwidth     :", fill="#000000", font=("Inter SemiBold", 14))
        self.canvas.create_text(1140, 542.0, anchor="nw", text="MHz", fill="#000000", font=("Inter SemiBold", 14))
        self.canvas.create_text(930.0, 600.0, anchor="nw", text="P", fill="#000000", font=("Inter SemiBold", 15, "italic"))
        self.canvas.create_text(940.0, 600.0, anchor="nw", text=" > ", fill="#000000", font=("Inter SemiBold", 15))

    def display_results(self, R_F, R_min, R_react, delta_d, t_int, t_map,
                        theta_ext, theta_b, theta_sr, theta_ss,
                        theta_point, N_row, SNR, transmitter_power):
        self.canvas.delete("results")

        self.canvas.create_text(645.0, 810.0, anchor="nw", text=f"{R_min:.2f} m", fill="#000000", font=("Inter Medium", 15), tags="results")
        self.canvas.create_text(645, 840, anchor="nw", text=f"{delta_d:.2f} cm", fill="#000000", font=("Inter Medium", 15), tags="results")
        self.canvas.create_text(645, 870, anchor="nw", text=f"{t_int:.2f} s", fill="#000000", font=("Inter Medium", 15), tags="results")
        self.canvas.create_text(645, 900, anchor="nw", text=f"{t_map:.2f} hr", fill="#000000", font=("Inter Medium", 15), tags="results")

        self.canvas.create_text(860, 810, anchor="nw", text=f"{theta_ext:.2f} deg", fill="#000000", font=("Inter Medium", 15), tags="results")
        self.canvas.create_text(860, 840, anchor="nw", text=f"{theta_b:.2f} arcsec", fill="#000000", font=("Inter Medium", 15), tags="results")
        self.canvas.create_text(860, 870, anchor="nw", text=f"{theta_sr:.2f} arcsec", fill="#000000", font=("Inter Medium", 15), tags="results")
        self.canvas.create_text(860, 900, anchor="nw", text=f"{theta_ss:.2f} arcsec", fill="#000000", font=("Inter Medium", 15), tags="results")

        self.canvas.create_text(1080, 810, anchor="nw", text=f"{theta_point:.2f} deg", fill="#000000", font=("Inter Medium", 15), tags="results")
        self.canvas.create_text(1080, 840, anchor="nw", text=f"{N_row:.2f}", fill="#000000", font=("Inter Medium", 15), tags="results")
        self.canvas.create_text(1080, 870, anchor="nw", text=f"{SNR:.2f} dB", fill="#000000", font=("Inter Medium", 15), tags="results")
        self.canvas.create_text(1080, 900, anchor="nw", text=f"{transmitter_power:.4e} W", fill="#000000", font=("Inter Medium", 15), tags="results")

    def display_inputs(self):
        try:
            self.canvas.delete("inputs")

            inputs = [
                {"value": self.entry_2.get().strip(), "unit": "m", "x": 60.0, "y": 810.0},
                {"value": self.entry_5.get().strip(), "unit": "GHz", "x": 60.0, "y": 840.0},
                {"value": self.entry_1.get().strip(), "unit": "m", "x": 60.0, "y": 870.0},
                {"value": self.entry_3.get().strip(), "unit": "", "x": 60.0, "y": 900.0},
                {"value": self.entry_6.get().strip(), "unit": "", "x": 210.0, "y": 810.0},
                {"value": self.entry_7.get().strip(), "unit": "", "x": 210.0, "y": 840.0},
                {"value": self.entry_19.get().strip(), "unit": "arc/sec", "x": 210.0, "y": 875.0},
                {"value": self.entry_12.get().strip(), "unit": "", "x": 210.0, "y": 903.0},
                {"value": self.entry_16.get().strip(), "unit": "μm", "x": 360.0, "y": 810.0},
                {"value": self.entry_17.get().strip(), "unit": "m", "x": 360.0, "y": 840.0},
                {"value": self.entry_8.get().strip(), "unit": "m", "x": 360.0, "y": 870.0},
                {"value": self.entry_9.get().strip(), "unit": "K", "x": 360.0, "y": 900.0},
                {"value": self.entry_10.get().strip(), "unit": "MHz", "x": 475.0, "y": 810.0},
            ]

            print("Displaying inputs:")
            for item in inputs:
                value_str = item["value"]
                unit = item["unit"]
                x = item["x"]
                y = item["y"]

                if value_str == "":
                    print(f"Warning: Empty input detected at coordinates ({x}, {y})")

                text = f"{value_str} {unit}".strip()
                self.canvas.create_text(x, y, anchor="nw", text=text, fill="#000000", font=("Inter Medium", 15), tags="inputs")
                print(f"Displayed: '{text}' at ({x}, {y})")

        except Exception as e:
            print(f"Error displaying inputs: {e}")