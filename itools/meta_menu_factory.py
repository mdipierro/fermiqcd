from Tkinter import *
import tkSimpleDialog

class EntryDialog(tkSimpleDialog.Dialog):
    def __init__(self, root, title,label, variable):
        self.root=root
        self.title_=title
        self.label=label
        self.textvariable=variable
        tkSimpleDialog.Dialog.__init__(self,root)
    def body(self,root):
        self.title(self.title_)
        Label(root, text=self.label).grid(row=0)
        self.e1 = Entry(root,textvariable=self.textvariable)
        self.e1.grid(row=1)
        return self.e1 # initial focus

class ScrollDialog(tkSimpleDialog.Dialog):
    def __init__(self, root, title, label, variable, from_=0, to=100, resolution=1):
        self.root=root
        self.title_=title        
        self.label=label        
        self.textvariable=variable
        self.from_=from_
        self.to=to
        self.resolution=resolution
        tkSimpleDialog.Dialog.__init__(self,root)
    def body(self,root):
        self.title(self.title_)
        Label(root, text=self.label).grid(row=0)
        #self.e1 = Entry(root,textvariable=self.textvariable)
        #self.e1.grid(row=1)
        self.s1=Scale(root,from_=self.from_,to=self.to,resolution=self.resolution,variable=self.textvariable,orient=HORIZONTAL,length=200)
        self.s1.grid(row=1,sticky='news')
        return self.s1 # initial focus
    
class MenuFactory:
    def __init__(self,root,items):
        self.root=root
        self.items=items
        self.menu=Menu(root)
        for item in items:            
            self.build(self.menu,item[0],item[1])

    def build(self,menu,label,command):
        if type(command)==type([]):
            submenu=Menu(menu,tearoff=0)
            for item in command:
                if item[0]=='---':
                    submenu.add_separator()
                else:
                    self.build(submenu,item[0],item[1])
            menu.add_cascade(label=label,menu=submenu)
        elif type(command)==type(None):
            menu.add_command(label=label,command=command,state=DISABLED)
        elif type(command)==type(BooleanVar()) and type(command.get())==type(True):
            menu.add_checkbutton(label=label,variable=command)
        elif type(command)==type(StringVar()) and type(command.get())==type(""):
            menu.add_radiobutton(label=label,variable=command,value=label)
        else:
            menu.add_command(label=label,command=command)

    def connect(self):
        self.root.config(menu=self.menu)
        
# display the menu


def test_all():
    root=Tk()
    tc=lambda:0
    bv=BooleanVar()
    sv=StringVar()
    dv=DoubleVar()
    sv.set("B")
    menu=[ ["File", [["Subsection1", tc],
                     ["---"],
                     ["Subsection2", None],
                     ["Subsection2", tc],
                     ["Check", bv]]],
           ["Submenus", [["Subsection1",tc],
                         ["Subsection2",[["Subsubsection1",tc],
                                         ["Subsubsection2",[['A',tc],
                                                            ['B',tc]]]]]]],
           ["Radios", [["A", sv],
                       ["B", sv],
                       ["C", sv]]],
           ["Dialogs", [["Dialog", lambda: EntryDialog(root,'Title','What?',sv)],
                        ["Scroll", lambda: ScrollDialog(root,'Title','What?',dv)]]]]
    m=MenuFactory(root,menu)
    m.connect()
    root.mainloop()

if __name__=='__main__':
    test_all()
