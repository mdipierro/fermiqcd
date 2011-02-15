from Tkinter import *
import Pmw

root = Tk()
root.title('MenuBar')
Pmw.initialise()

balloon = Pmw.Balloon(root)

menuBar = Pmw.MenuBar(root, hull_relief=RAISED,  hull_borderwidth=1,
                      balloon=balloon)
menuBar.pack(fill=X)

menuBar.addmenu('Buttons', 'Simple Commands')
for i in range(100):
	menuBar.addmenuitem('Buttons', 'command', str(i),
	        	    font=('StingerLight', 14), label=str(i))
root.mainloop()