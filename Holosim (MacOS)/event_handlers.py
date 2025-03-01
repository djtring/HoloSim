# This class acts like a bridge between the GUI and the calculator class.
# It handles all events triggered by the user. Each method retrieves the necessary input from the GUI,
# calls the corresponding method from the calculator class to compute results, 
# and then updates the GUI accordingly.

import os
import sys
from tkinter import messagebox
from calculator import Calculator

class EventHandelers:
    def __init__(self, gui, calculator):
        self.gui = gui
        self.calculator = calculator

    def open_help_window(self):
        base_path = os.path.dirname(__file__)
        pdf_path = os.path.join(base_path, 'Submillimeter_Wave_Holography_for_Large_Dish_Antennas.pdf')

        if os.name == 'posix':
            try:
                if sys.platform == 'darwin':  # macOS
                    os.system(f'open "{pdf_path}"')
                else:  # Linux
                    os.system(f'xdg-open "{pdf_path}"')
            except Exception as e:
                messagebox.showerror("Error", f"Could not open help PDF: {e}")
        else:
            messagebox.showerror("Error", "Unsupported operating system.")

    def open_Reference_Window(self):
        base_path = os.path.dirname(__file__)
        pdf_path = os.path.join(base_path, 'HoloSim_Quick_Reference.pdf')

        if os.name == 'posix':
            try:
                if sys.platform == 'darwin':
                    os.system(f'open "{pdf_path}"')
                else:
                    os.system(f'xdg-open "{pdf_path}"')
            except Exception as e:
                messagebox.showerror("Error", f"Could not open reference PDF: {e}")
        else:
            messagebox.showerror("Error", "Unsupported operating system.")

    def calculate_all_steps(self):
        try:
            D_input = self.gui.entry_2.get().strip()
            nu_input = self.gui.entry_5.get().strip()

            if not D_input or not nu_input:
                raise ValueError("Empty input for step 1")

            D = float(D_input)
            nu = float(nu_input)

            R_F, R_min, R_react = self.calculator.calculate_step1(D, nu)

            a_input = self.gui.entry_1.get().strip()
            if not a_input:
                raise ValueError("Empty input for step 2")
            a = float(a_input)

            delta_d = self.calculator.calculate_step2(a)

            f_1 = float(self.gui.entry_3.get().strip())
            f_apo = float(self.gui.entry_6.get().strip())
            f_osr = float(self.gui.entry_7.get().strip())
            dtheta = float(self.gui.entry_19.get().strip())

            t_int, t_map = self.calculator.calculate_step3(delta_d, D, nu, f_1, f_apo, f_osr, dtheta)

            f_oss_input = self.gui.entry_12.get().strip()
            if not f_oss_input:
                raise ValueError("Empty input for step 4")
            f_oss = float(f_oss_input)

            theta_ext, theta_b, theta_sr, theta_ss = self.calculator.calculate_step4(delta_d, nu, D, f_oss)

            delta_z_input = self.gui.entry_16.get().strip()
            if not delta_z_input:
                raise ValueError("Empty input for step 5")
            delta_z = float(delta_z_input)

            theta_point, N_row, SNR = self.calculator.calculate_step5(delta_d, delta_z, D, nu, f_apo, f_osr, f_oss)

            D_t_input = self.gui.entry_17.get().strip()
            z_input = self.gui.entry_8.get().strip()
            T_sys_input = self.gui.entry_9.get().strip()
            B_input = self.gui.entry_10.get().strip()

            if not D_t_input or not z_input or not T_sys_input or not B_input:
                raise ValueError("Empty input for step 6")

            D_t = float(D_t_input)
            z = float(z_input)
            T_sys = float(T_sys_input)
            B = float(B_input)

            transmitter_power = self.calculator.calculate_step6(SNR, D_t, D, z, T_sys, B, nu)

            # Clear old final results
            self.gui.canvas.delete("results")

            # Display final (all steps) results with "results" tag
            self.gui.display_results(R_F, R_min, R_react, delta_d, t_int, t_map,
                                     theta_ext, theta_b, theta_sr, theta_ss,
                                     theta_point, N_row, SNR, transmitter_power)

            
            self.gui.display_inputs()

        except ValueError as e:
            messagebox.showerror("Calculation Error", str(e))

    def calculate_and_display_step1(self):
        try:
            D_input = self.gui.entry_2.get().strip()
            nu_input = self.gui.entry_5.get().strip()

            if not D_input or not nu_input:
                raise ValueError("Empty input for step 1")

            D = float(D_input)
            nu = float(nu_input)

            R_F, R_min, R_react = self.calculator.calculate_step1(D, nu)

            # Clear previous Step 1 results before showing new ones
            self.gui.canvas.delete("step1_results")

            # Display step 1 results with a unique tag
            self.gui.canvas.create_text(145.0, 210.0, anchor="nw", text=f"{R_F:.2f} m", fill="#000000", font=("Inter Medium", 13), tags="step1_results")
            self.gui.canvas.create_text(145.0, 230.0, anchor="nw", text=f"{R_min:.2f} m", fill="#000000", font=("Inter Medium", 13), tags="step1_results")
            self.gui.canvas.create_text(145.0, 250.0, anchor="nw", text=f"{R_react:.2f} m", fill="#000000", font=("Inter Medium", 13), tags="step1_results")

            # No display_inputs() call here, as requested

        except ValueError as e:
            messagebox.showerror("Calculation Error", str(e))

    def calculate_and_display_step2(self):
        try:
            a_input = self.gui.entry_1.get().strip()

            if not a_input:
                raise ValueError("Empty input for step 2")

            a = float(a_input)
            delta_d = self.calculator.calculate_step2(a)

            # Clear previous Step 2 results
            self.gui.canvas.delete("step2_results")

            # Display step 2 results
            self.gui.canvas.create_text(590.0, 230.0, anchor="nw", text=f"{delta_d:.2f} cm", fill="#000000", font=("Inter Medium", 15), tags="step2_results")

            # No display_inputs()

        except ValueError as e:
            messagebox.showerror("Calculation Error", str(e))

    def calculate_and_display_step3(self):
        try:
            D_input = self.gui.entry_2.get().strip()
            nu_input = self.gui.entry_5.get().strip()

            if not D_input or not nu_input:
                raise ValueError("Empty input from previous step")

            D = float(D_input)
            nu = float(nu_input)

            a = float(self.gui.entry_1.get().strip())
            delta_d = self.calculator.calculate_step2(a)

            f_1 = float(self.gui.entry_3.get().strip())
            f_apo = float(self.gui.entry_6.get().strip())
            f_osr = float(self.gui.entry_7.get().strip())
            dtheta = float(self.gui.entry_19.get().strip())

            t_int, t_map = self.calculator.calculate_step3(delta_d, D, nu, f_1, f_apo, f_osr, dtheta)

            # Clear previous Step 3 results
            self.gui.canvas.delete("step3_results")

            self.gui.canvas.create_text(1010.0, 250.0, anchor="nw", text=f"{t_int:.2f} s", fill="#000000", font=("Inter SemiBold", 13), tags="step3_results")
            self.gui.canvas.create_text(1010.0, 275.0, anchor="nw", text=f"{t_map:.2f} hr", fill="#000000", font=("Inter SemiBold", 13), tags="step3_results")

            # No display_inputs()

        except ValueError as e:
            messagebox.showerror("Calculation Error", str(e))

    def calculate_and_display_step4(self):
        try:
            nu_input = self.gui.entry_5.get().strip()
            D_input = self.gui.entry_2.get().strip()
            f_oss_input = self.gui.entry_12.get().strip()

            if not nu_input or not D_input or not f_oss_input:
                raise ValueError("Empty input")

            nu = float(nu_input)
            D = float(D_input)
            f_oss = float(f_oss_input)

            a = float(self.gui.entry_1.get().strip())
            delta_d = self.calculator.calculate_step2(a)

            theta_ext, theta_b, theta_sr, theta_ss = self.calculator.calculate_step4(delta_d, nu, D, f_oss)

            # Clear previous Step 4 results
            self.gui.canvas.delete("step4_results")

            self.gui.canvas.create_text(130.0, 520.0, anchor="nw", text=f"{theta_ext:.2f} deg", fill="#000000", font=("Inter Medium", 14), tags="step4_results")
            self.gui.canvas.create_text(130.0, 543.0, anchor="nw", text=f"{theta_b:.2f} arcsec", fill="#000000", font=("Inter Medium", 14), tags="step4_results")
            self.gui.canvas.create_text(130.0, 568.0, anchor="nw", text=f"{theta_sr:.2f} arcsec", fill="#000000", font=("Inter Medium", 14), tags="step4_results")
            self.gui.canvas.create_text(130.0, 590.0, anchor="nw", text=f"{theta_ss:.2f} arcsec", fill="#000000", font=("Inter Medium", 14), tags="step4_results")

            # No display_inputs()

        except ValueError as e:
            messagebox.showerror("Calculation Error", str(e))

    def calculate_and_display_step5(self):
        try:
            delta_z_input = self.gui.entry_16.get().strip()
            if not delta_z_input:
                raise ValueError("Empty input")
            delta_z = float(delta_z_input)

            a = float(self.gui.entry_1.get().strip())
            delta_d = self.calculator.calculate_step2(a)

            D = float(self.gui.entry_2.get().strip())
            nu = float(self.gui.entry_5.get().strip())
            f_apo = 1.3
            f_osr = 2.2
            f_oss = float(self.gui.entry_12.get().strip())

            theta_point, N_row, SNR = self.calculator.calculate_step5(delta_d, delta_z, D, nu, f_apo, f_osr, f_oss)

            # Clear previous Step 5 results
            self.gui.canvas.delete("step5_results")

            self.gui.canvas.create_text(575.0, 528.0, anchor="nw", text=f"{theta_point:.2f} deg", fill="#000000", font=("Inter Medium", 14), tags="step5_results")
            self.gui.canvas.create_text(575.0, 555.0, anchor="nw", text=f"{N_row:.2f}", fill="#000000", font=("Inter Medium", 14), tags="step5_results")
            self.gui.canvas.create_text(575.0, 580.0, anchor="nw", text=f"{SNR:.2f} dB", fill="#000000", font=("Inter Medium", 14), tags="step5_results")

            # No display_inputs()

            return SNR

        except ValueError as e:
            messagebox.showerror("Calculation Error", str(e))
            return None

    def calculate_and_display_step6(self):
        try:
            SNR = self.calculate_and_display_step5()
            if SNR is None:
                return

            D = float(self.gui.entry_2.get().strip())
            nu = float(self.gui.entry_5.get().strip())

            D_t_input = self.gui.entry_17.get().strip()
            z_input = self.gui.entry_8.get().strip()
            T_sys_input = self.gui.entry_9.get().strip()
            B_input = self.gui.entry_10.get().strip()

            if not D_t_input or not z_input or not T_sys_input or not B_input:
                raise ValueError("Empty input for Step 6")

            D_t = float(D_t_input)
            z = float(z_input)
            T_sys = float(T_sys_input)
            B = float(B_input)

            result = self.calculator.calculate_step6(SNR, D_t, D, z, T_sys, B, nu)

            if result is None:
                return

            # Clear previous Step 6 results
            self.gui.canvas.delete("step6_results")

            self.gui.canvas.create_text(960.0, 600.0, anchor="nw", text=f"{result:.4e} W", fill="#000000", font=("Inter Medium", 14), tags="step6_results")

            # No display_inputs()

        except ValueError as e:
            messagebox.showwarning("Input Error", str(e))