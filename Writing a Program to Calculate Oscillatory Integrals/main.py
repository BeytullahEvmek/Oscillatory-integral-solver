# Library (just one for this page)
import tkinter as tk
# Importing the functions that display each page/module of the app
from airy_function import show_airy_function
from filon_type import show_filon_type
from stationary_phase import show_stationary_phase

#Removes all widgets inside content_frame so that it removes the old content in the new page.
def clear_frame():
    for widget in content_frame.winfo_children():
        widget.destroy()

#It makes it so that program can activate the other parts of the programs via buttons
def show_page(page_func):
    # clear_frame works already but just in case if it doesn't work after returning from another page
    clear_frame()
    page_func(content_frame, show_main_page)

#Closes the app
def exit_app(): #Closes the app
    root.quit()

# Initialize the main Tkinter window
root = tk.Tk()
root.title("Oscillatory Integral") #Name
root.configure(bg="#e6ffe6")  # Background color

# Gets resolution of the app
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()

root.state('zoomed')  # Maximizes the window

# Makes the frame and makes borders to look fancy, also colors outside green as well so it doesn't look weird
frame = tk.Frame(root, bg="#ccffcc", bd=10, relief="sunken")
frame.place(relx=0.5, rely=0.5, anchor="center", relwidth=0.98, relheight=0.98)

# Makes inside of the frame and colors it green
content_frame = tk.Frame(frame, bg="#ccffcc")
content_frame.place(relx=0.5, rely=0.5, anchor="center", relwidth=0.9, relheight=0.9)

#Main Page with the declaration of parent_frame
def show_main_page(parent_frame):

    clear_frame() #Calls clear_frame()

    button_font_size = int(screen_height * 0.04) #Size of the ,text inside buttons
    button_font = ("Arial", button_font_size) #Type of the fond

    # Splits the screen to rows and columns
    for i in range(4):
        parent_frame.grid_columnconfigure(i, weight=1, uniform="button_col")
    for i in range(4):
        parent_frame.grid_rowconfigure(i, weight=1)

    #Label
    tk.Label(parent_frame, text="Welcome", font=("Helvetica", 80, "bold"), bg="#ccffcc").grid(row=0, column=1, columnspan=2, pady=30, sticky="nsew")

    # Buttons, including exit
    btn1 = tk.Button(parent_frame, text="Airy Function\nApproximation", font=button_font, command=lambda: show_page(show_airy_function))
    btn2 = tk.Button(parent_frame, text="Stationary Phase\nApproximation", font=button_font, command=lambda: show_page(show_stationary_phase))
    btn3 = tk.Button(parent_frame, text="Filon Quadrature\nApproximation", font=button_font, command=lambda: show_page(show_filon_type))
    btn_exit = tk.Button(parent_frame, text="Exit", font=("Arial", 40), bg="#ffcccc", command=exit_app)

    # Position of the buttons (sticky element makes it so that it spread the button)
    btn1.grid(row=1, column=0, columnspan=2, sticky="nsew", padx=10, pady=10)
    btn2.grid(row=1, column=2, columnspan=2, sticky="nsew", padx=10, pady=10)
    btn3.grid(row=2, column=0, columnspan=2, sticky="nsew", padx=10, pady=10)
    btn_exit.grid(row=2, column=2, columnspan=2, sticky="nsew", padx=10, pady=10)


# Shows the main menu
show_main_page(content_frame)

# Start Tkinter event loop
root.mainloop()
