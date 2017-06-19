import numpy as np
import pandas as pd
from math import sin,cos,tan,asin,acos,atan
from pint import UnitRegistry
u = UnitRegistry()

def RotateZ(body, omega):
    omega *= -1
    rmat = [[cos(omega), -sin(omega), 0],
            [sin(omega), +cos(omega), 0],
            [0,          0,          +1]]
    return Rotate(body,rmat)
def RotateY(body, omega):
    omega *= -1
    rmat = [[cos(omega), 0,          -sin(omega)],
            [0,          1,          0          ],
            [sin(omega), 0,          +cos(omega)]]
    return Rotate(body,rmat)
def RotateX(body, omega):
    omega *= -1
    rmat = [[1,          0,           0],
            [0,          cos(omega), -sin(omega)],
            [0,          sin(omega), +cos(omega)]]
    return Rotate(body,rmat)
def Offset(body, offset):
    return Modify(body, lambda p: np.add(p, offset)) # keeps units
def Rotate(body, rmat):
    return Modify(body, lambda p: (p * np.mat(rmat)).tolist()[0]) # keeps units
def Modify(body, f):
    for point in body:
        point[0], point[1], point[2] = f(point)

def Unit(head, tail):
    vector = np.mat(head)-np.mat(tail)
    vector = vector/np.linalg.norm(vector)
    return vector.tolist()[0]
def Length(head,tail):
    return sum([(h-t)**2 for h,t in zip(head,tail)])**.5


def calculate():
    
    # Part dimensions
    # ---------------
    
    boom_jib =[18788, 0, -742]*u.mm
    boom_cyljib =   [12588, 0, -1267]*u.mm
    boom_cylboom =  [5800, 0, -1450]*u.mm
    boom_cog =      [8600, 0, -500]*u.mm
    boom_pedestal1 =  [0, 1280, 0]*u.mm
    boom_pedestal2 =  [0, -1280, 0]*u.mm
    boom = [boom_jib, boom_cyljib, boom_cylboom, boom_cog, boom_pedestal1, boom_pedestal2]

    jib_cyljib =    [1631, 0, 1236]*u.mm
    jib_load =      [14957, 0, 1150]*u.mm
    jib_cog =       [5200, 0, 600]*u.mm
    jib_boom1 =     [0, 697, 0]*u.mm
    jib_boom2 =     [0, -697, 0]*u.mm
    jib = [jib_cyljib, jib_load, jib_cog, jib_boom1, jib_boom2]

    pedestal_boom =     [0,0,0]*u.mm
    pedestal_cylboom =  [1300,0,-2500]*u.mm
    pedestal = [pedestal_boom, pedestal_cylboom]

    boomjib = boom+jib
    pbj = pedestal+boom+jib

    # Component assembly
    # ------------------
    
    RotateY     (jib,       180*u.deg+beta)
    Offset      (jib,       boom_jib)
    RotateY     (boomjib,   alfa)
    RotateZ     (pbj,       gamma)

    # Define constraints
    # ------------------

    x,y,z = [[1.,0.,0.]*u.t,
             [0.,1.,0.]*u.t,
             [0.,0.,1.]*u.t]
    cb = Unit(boom_cylboom,pedestal_cylboom)*u.t
    cj = Unit(jib_cyljib,boom_cyljib)*u.t

    RotateZ([x,y,z,cb],gamma)
    d1 = boom_pedestal1
    d2 = boom_pedestal2
    d3 = boom_cylboom

    constraintsFD = [[x,d1],
                   [y,d1],
                   [z,d1],
                   [x,d2],
                   [z,d2],
                   [cb,d3]]
    constraintsFM = np.transpose(
        np.vstack(np.hstack((F, np.cross(F,D))) for F,D in constraintsFD)
        ) # losing units here: +cross, ~stacks
    
    # Define loads
    # ------------
    
    loadsFD = [[acc*m_jib,     jib_cog],
             [acc*m_boom,    boom_cog],
             [acc*m_cyl*2,   np.add(boom_cyljib,jib_cyljib)/2],
             [acc*m_load,    jib_load]]
    loadsFM = sum(np.hstack((F,np.cross(D,F))) for F,D in loadsFD)
    
    # Solve for constraints
    # ---------------------

    reaction = np.linalg.solve(constraintsFM,loadsFM)

    # Post process
    # ------------

##    print(loadsFM.astype(int))
##    print(constraintsFM.astype(int))
##    print(constraintsFM)
##    print(reaction.astype(int))
    
    return {"alfa":alfa.to(u.deg).magnitude,
            "beta":beta.to(u.deg).magnitude,
            "gamma":gamma.to(u.deg).magnitude,
            "m_load":m_load.to(u.t).magnitude,
            "reachx":Length(pedestal_boom[0:2],jib_load[0:2]).to(u.m).magnitude,
            "reachz":jib_load[2].to(u.m).magnitude,
            "cblen":Length(boom_cylboom,pedestal_cylboom).to(u.m).magnitude,
            "cjlen":Length(jib_cyljib,boom_cyljib).to(u.m).magnitude,
            **dict(zip(["FX","FY","FZ","MX","MY","MZ"],loadsFM)),
            **dict(zip(["RX1","RY1","RZ1","RX2","RZ2","RC3"],reaction))}

def main():
    global alfa, beta, gamma, m_jib, m_boom, m_cyl, m_load, acc, df, data

    alfas =  list(range(0,  75+1, 15))*u.deg
    betas =  list(range(0, 120+1, 30))*u.deg
    gammas = list(range(0, 360+1, 30))*u.deg
   
    m_jib = 6*u.t
    m_boom = 16*u.t
    m_cyl = 2.6*u.t
    m_loads = [5,10]*u.t

    acc = np.array([0.21,  0.48,    -1.17])

    data = []

    for gamma in gammas:
        gamma = 180*u.deg
        for beta in betas:
            beta = 0*u.deg
            for alfa in alfas:
                alfa = 0*u.deg
                for m_load in m_loads:
                    data.append(calculate())
                    break
                break
            break
        break
                    
    df = pd.DataFrame(data)
    cols = ['gamma','alfa','beta','reachx','RC3','cblen']
    # print(df[df.reachx > 24][df.m_load == 5][cols])
    print(df[cols])

       

main()
