import tkinter
import tkinter.messagebox

def showMessage():
	tkinter.messagebox.showinfo("Welcome!", "Hello, there, Javier!")

def menu_action():
	sys.sterr.write("Executed menu!")

main = tkinter.Tk()

menubar = tkinter.Menu(main)

filemenu = tkinter.Menu(menubar)
filemenu.add_command(label = "New", command = menu_action)
filemenu.add_command(label = "Open", command = menu_action)
filemenu.add_command(label = "Save", command = menu_action)
filemenu.add_command(label = "Save as", command = menu_action)
filemenu.add_separator()
filemenu.add_command(label = "Quit", command = main.quit)
helpmenu = tkinter.Menu(menubar)
helpmenu.add_command(label = "Help", command = menu_action)
helpmenu.add_separator()
helpmenu.add_command(label = "About", command = menu_action)
menubar.add_cascade(label = "File", menu = filemenu)
menubar.add_cascade(label = "Help", menu = helpmenu)
main.config(menu = menubar)

b = tkinter.Button(main, text = "Greetings!", command = showMessage)

b.pack()

main.mainloop()

#tkinter.filedialog.askopenfilename( filetypes=[("Python files", "*.py")],
#title = “Open Python File” )