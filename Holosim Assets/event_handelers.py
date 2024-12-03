from tkinter import messagebox
from calculator import Calculator


class EventHandelers: 
    def __init__(self, gui, calculator):
        self.gui = gui
        self.calculator = calculator

    def calculate_all_steps(self):
        try:
            D_input = self.gui.entry_2.get().strip() 
            nu_input = self.gui.entry_5.get().strip()

            if not D_input or not nu_input:
                raise ValueError("Error", "Empty unput for step 1")
            
            D = float(D_input)
            nu = float(nu_input)

            R_F, R_min, R_react = self.calculator.calculate_step1(D, nu)

            a_input = self.gui.entry_1.get().strip()
            
            if not a_input: 
                raise ValueError("Error", "Empty input for step 2")
            
            a = float(a_input)

            delta_d = self.calculator.calculate_step2(a)

            f_1 = float(self.gui.entry_3.get().strip())
            f_apo = float(self.gui.entry_6.get().strip())
            f_osr = float(self.gui.entry_7.get().strip())
            dtheta = float(self.gui.entry_17.get().strip())

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

            transmitter_power = self.calculator.calculate_step6(SNR, D_t, D, z, T_sys, B, nu)

            self.gui.canvanas.delete("results")

            self.gui.display_results(R_F, R_min, R_react, delta_d, t_int, t_map, theta_ext, theta_b, theta_sr, theta_ss, theta_point, N_row, SNR, transmitter_power)

            self.gui.display_inputs()
        
        except ValueError as e:
            messagebox.showerror("Calculation Error", str(e))
                



                