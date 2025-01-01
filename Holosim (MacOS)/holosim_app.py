# holosim_app.py

import tkinter as tk
from pathlib import Path
from calculator import Calculator
from event_handlers import EventHandelers
from gui_components import GUI

class HoloSimApp:
    def __init__(self):
        # Initialize the main Tk window
        self.window = tk.Tk()

        # Determine the output path and assets path based on current file location
        self.output_path = Path(__file__).resolve().parent
        self.assets_path = self.output_path

        def relative_to_assets(path: str) -> Path:
            return self.assets_path / Path(path)

        # Create the GUI, Calculator, and EventHandelers instances
        self.gui = GUI(self.window, relative_to_assets)
        self.calculator = Calculator()
        self.handlers = EventHandelers(self.gui, self.calculator)

        # Build the GUI and connect buttons to the event handlers
        # These commands correspond to methods defined in EventHandelers
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