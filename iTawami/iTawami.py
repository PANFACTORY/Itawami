#**********************************************************
#Title      :iTawami
#Author     :Tanabe Yuta
#Date       :2019/07/25
#Copyright  :(C)2019 TanabeYuta
#**********************************************************


import numpy as np
import math
import tkinter as tk
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib.backends.backend_tkagg  import FigureCanvasTkAgg
import scipy.linalg
from numba import jit


#********** Get deflection distribution (support all edges) **********
def GetW(_a, _b, _xi, _ita, _m, _n, _x, _y, _D):
    m = np.arange(1, _m+1)
    n = np.arange(1, _n+1)

    #----------Culculating coefficient----------
    c = 4.0*np.tile(np.sin(math.pi*_xi*m/_a), (_n, 1)).T*np.tile(np.sin(math.pi*_ita*n/_b), (_m, 1))/(_a*_b)
    d = np.square(np.square(np.tile(m, (_n, 1)).T)/(_a**2.0) + np.square(np.tile(n, (_m, 1)))/(_b**2.0)) 
    mn = c/d
    
    #----------Culculating deflection----------
    xm = np.sin(math.pi*np.tile(_x, (_m, 1)).T*np.tile(m, (_x.size, 1))/_a)
    ny = np.sin(math.pi*np.tile(n, (_y.size, 1)).T*np.tile(_y, (_n, 1))/_b)
    return -np.dot(xm, np.dot(mn, ny))/math.pi**4.0/_D


#********** Get deflection distribution (fix all edges) **********
@jit
def GetW2(_a, _b, _xi, _ita, _m, _n, _x, _y, _D):
    m = np.arange(1, _m+1)
    n = np.arange(1, _n+1)

    #----------Culculating coefficient----------
    a = np.zeros([_m*_n, _m*_n])
    b = np.zeros(_m*_n)
    for i in range(1, _m+1):
        for j in range(1, _n+1):
            ii = (i-1)*_n+(j-1)
            a[ii][ii] += 3.0*(i/_a)**4.0+3.0*(j/_b)**4.0+2.0*(i/_a)**2.0*(j/_b)**2.0

            for k in range(1, _n+1):
                jj = (i-1)*_n+(k-1)
                if(k != j):
                    a[ii][jj] += 2.0*(i/_a)**4.0
            for k in range(1, _m+1):
                jj = (k-1)*_n+(j-1)
                if(k != i):
                    a[ii][jj] += 2.0*(j/_b)**4.0

            b[ii] = (1.0 - math.cos(2.0*math.pi*i*_xi/_a))*(1.0 - math.cos(2.0*math.pi*j*_ita/_b))
    a *= 4.0*_D*math.pi**4.0*_a*_b

    A = scipy.linalg.lu_factor(a)
    mn = np.reshape(scipy.linalg.lu_solve(A, b), (_m, _n))
    
    #----------Culculate deflection----------
    xm = 1.0 - np.cos(2.0*math.pi*np.tile(_x, (_m, 1)).T*np.tile(m, (_x.size, 1))/_a)
    ny = 1.0 - np.cos(2.0*math.pi*np.tile(n, (_y.size, 1)).T*np.tile(_y, (_n, 1))/_b)
    return -np.dot(xm, np.dot(mn, ny))


#********** Control GUI **********
class Window:
    #----------Initializer----------
    def __init__(self):
        #.....Parameters for simuation.....
        self.a = 2.0
        self.b = 4.0
        t = 0.1
        self.xi = 0.5*self.a
        self.ita = 0.5*self.b
        E = 210000.0
        V = 0.3
        self.D = E*t**3.0/(12.0*(1.0-V**2.0))
        self.wmax = 0.0

        #.....Parameters for glaphics.....
        self.angle1 = 0.0
        self.angle2 = 0.0
        self.mode = "s"
        self.capture = ""

        #.....Setting GUI.....
        self.root = tk.Tk()
        self.root.title("iTawami")

        #.....Binding Keys.....
        self.root.bind('<Up>', self.UpKey)
        self.root.bind('<Down>', self.DownKey)
        self.root.bind('<Right>', self.RightKey)
        self.root.bind('<Left>', self.LeftKey)
        self.root.bind('<Key-s>', self.SKey)        #Draw all edges supported plate
        self.root.bind('<Key-f>', self.FKey)        #Draw all edges fixed plate
        self.root.bind('<Key-a>', self.AKey)        #Modify value of a
        self.root.bind('<Key-b>', self.BKey)        #Modify value of b
        self.root.bind('<KeyPress>', self.Key)      

        #.....Initial plot.....
        self.Init_graph()
        self.PlotPlate()
        self.Draw()

        self.root.mainloop()


    #----------Initialize graph----------
    def Init_graph(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_xlabel("x") 
        self.ax.set_ylabel("y") 
        self.ax.set_zlabel("w") 
        self.Canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.Canvas.get_tk_widget().grid(row=0, column=0, rowspan=10)

    
    #----------Make and plot plate(support all edges)----------
    def PlotPlate(self):
        x = np.arange(0, self.a, self.a/100.0)
        y = np.arange(0, self.b, self.b/100.0)
        X, Y = np.meshgrid(x, y)
        Z = GetW(self.a, self.b, self.xi, self.ita, 100, 100, x, y, self.D)
        self.wmax = Z.min()
        self.ax.plot_surface(X,Y,Z.T)
        self.ax.set_title("All edges supported") 


    #----------Make and plot plate(fix all edges)----------
    def PlotPlate2(self):
        x = np.arange(0, self.a, self.a/100.0)
        y = np.arange(0, self.b, self.b/100.0)
        X, Y = np.meshgrid(x, y)
        Z = GetW2(self.a, self.b, self.xi, self.ita, 100, 100, x, y, self.D)
        self.wmax = Z.min()
        self.ax.plot_surface(X,Y,Z.T)
        self.ax.set_title("All edges fixed") 
    

    #----------Draw plots----------
    def Draw(self):
        self.ax.view_init(self.angle1, self.angle2)
        self.ax.text2D(0, 0.1, "a=" + '{:.2g}'.format(self.a), transform=self.ax.transAxes)
        self.ax.text2D(0, 0.05, "b=" + '{:.2g}'.format(self.b), transform=self.ax.transAxes)
        self.ax.text2D(0, 0, "wmax=" + '{:.4g}'.format(self.wmax), transform=self.ax.transAxes)
        self.Canvas.draw()


    def UpKey(self, event):
        self.angle1 += 10
        self.Draw()

    def DownKey(self, event):
        self.angle1 -= 10
        self.Draw()

    def RightKey(self, event):
        self.angle2 += 10
        self.Draw()

    def LeftKey(self, event):
        self.angle2 -= 10
        self.Draw()

    def SKey(self, event):
        self.mode = "s"
        self.Init_graph()
        self.PlotPlate()
        self.Draw()

    def FKey(self, event):
        self.mode = "f"
        self.Init_graph()
        self.PlotPlate2()
        self.Draw()

    def AKey(self, event):
        try:
            self.a = float(self.capture)  
        except:
            self.capture = ""
        self.xi = 0.5*self.a
        self.Init_graph()
        if(self.mode == "s"):
            self.PlotPlate()
        else:
            self.PlotPlate2()
        self.Draw()
        self.capture = ""

    def BKey(self, event):
        try:
            self.b = float(self.capture)  
        except:
            self.capture = ""
        self.ita = 0.5*self.b
        self.Init_graph()
        if(self.mode == "s"):
            self.PlotPlate()
        else:
            self.PlotPlate2()
        self.Draw()
        self.capture = ""

    def Key(self, event):
        if(event.char.isdecimal() or event.char == "."):
            self.capture += event.char        


#********** Main process **********
if __name__ == "__main__":
    #----------入力値の受け取り----------
    w = Window()
