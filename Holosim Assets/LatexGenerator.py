import tkinter as tk
from tkinter import filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

# Function to save the plot as an image
def save_plot():
    # Ask the user where they want to save the file
    file_path = filedialog.asksaveasfilename(defaultextension=".png", 
                                             filetypes=[("PNG files", "*.png"), 
                                                        ("All files", "*.*")])
    if file_path:
        # Save the figure to the specified file path with tight bounding box
        fig.savefig(file_path, transparent=True, bbox_inches='tight', pad_inches=0)

# Create the main application window and set its size
root = tk.Tk()
root.title("Matplotlib LaTeX in Tkinter")
root.geometry("300x200")  # Set a smaller window size

# Create a figure and axis for plotting with smaller figure size
fig, ax = plt.subplots(figsize=(2, 4))  # Adjust figure size for stacked lines

# Make the figure background transparent
fig.patch.set_alpha(0)  # Set transparency for the figure background
ax.set_facecolor('none')  # Ensure the axis background is also transparent

latex_text_ext = r"$B:$"

latex_font_size = 7 # Adjust font size

# Set the LaTeX formatted equation
ax.text(0.5, 0.5, latex_text_ext, fontsize=latex_font_size, va='center', ha='center')

# Remove axes for a clean appearance
ax.axis('off')

# Create a canvas widget to display the plot in Tkinter
canvas = FigureCanvasTkAgg(fig, master=root)

# Remove the background of the canvas
canvas.get_tk_widget().config(bg='white', highlightthickness=0)
canvas.draw()

# Pack the canvas widget into the Tkinter window
canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

# Create a Save button
save_button = tk.Button(root, text="Save", command=save_plot)
save_button.pack(pady=10)

# Start the Tkinter event loop
root.mainloop()