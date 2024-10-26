import tkinter as tk
import sys
import os
from pathlib import Path
from calculator import Calculator
from event_handlers import EventHandelers
from gui_components import GUI

class HoloSimApp:
    def __init__(self):
        # Initialize the main Tk window
        self.window = tk.Tk()

        # Determine the base path dynamically (Handles both script & .exe mode)
        if getattr(sys, '_MEIPASS', False):  
            self.assets_path = Path(sys._MEIPASS) / "assets"
        else:
            self.assets_path = Path(__file__).resolve().parent / "assets"

        def relative_to_assets(path: str) -> Path:
            """ Returns the full path to assets, working for both .py and .exe modes """
            return self.assets_path / path

        # Create GUI, Calculator, and EventHandelers instances
        self.gui = GUI(self.window, relative_to_assets)
        self.calculator = Calculator()
        self.handlers = EventHandelers(self.gui, self.calculator)

        # Build the GUI and connect buttons to the event handlers

        self.gui.build_gui(
            calculate_all_steps_cmd=self.handlers.calculate_all_steps,
            calculate_and_display_step1_cmd=self.handlers.calculate_and_display_step1,
            calculate_and_display_step2_cmd=self.handlers.calculate_and_display_step2,
            calculate_and_display_step3_cmd=self.handlers.calculate_and_display_step3,
            calculate_and_display_step4_cmd=self.handlers.calculate_and_display_step4,
            calculate_and_display_step5_cmd=self.handlers.calculate_and_display_step5,
            calculate_and_display_step6_cmd=self.handlers.calculate_and_display_step6,
            open_help_window_cmd=self.handlers.open_help_window,
            open_Reference_Window_cmd=self.handlers.open_Reference_Window
        )

    def run(self):
        # Start the Tkinter main loop
        self.window.mainloop()