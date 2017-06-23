import matplotlib.pyplot as plt
import logging
import numpy as np
import pandas as pd
from math import sin,cos,tan,asin,acos,atan
from pint import UnitRegistry
u = UnitRegistry()

def Scale(body,factor):
    '''Not used'''
    rmat = [[factor,0,0],
            [0,factor,0],
            [0,0,factor]]
    return Rotate(body,rmat) 
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
    for point in body.values():
        point[0], point[1], point[2] = f(point)

def Unit(head, tail):
    vector = np.mat(head)-np.mat(tail)
    vector = vector/np.linalg.norm(vector)
    return vector.tolist()[0]
def Length(head,tail):
    return sum([(h-t)**2 for h,t in zip(head,tail)])**.5
def FD2FM(FD):
    return [np.hstack((F,np.cross(D,F))) for F,D in FD]

def plot_geometry(geom,xaxis,yaxis,width,height):
    '''Depricated'''
    h = '+-'+'--'*width+'-+\n'
    d = h
    vs = geom.values()
    xs = [v[xaxis] for v in vs]
    ys = [v[yaxis] for v in vs]
    xmax = max(     xs)    .m
    xmin =     min( xs)    .m
    ymax = max(         ys).m
    ymin =     min(     ys).m
    f = min(height/(ymax-ymin),width/(xmax-xmin))
    for y in range(height,0-1,-1):
        d += '|'
        for x in range(0,width+1,1):
            for k in geom:
                q1 = int((geom[k][xaxis].m-xmin)*f) == x 
                q2 = int((geom[k][yaxis].m-ymin)*f) == y
                if q1 and q2:
                    d += k.split('_')[0][0]+k.split('_')[1][0]
                    break
            else:
                d += '  '
        d += '|\n'
    d += h
    print(d)
def geometry():
    
    # Part dimensions
    # ---------------
    
    boom = {"boom_jib" :        [18788, 0, -742]*u.mm,
            "boom_cyljib" :     [12588, 0, -1267]*u.mm,
            "boom_cylboom" :    [5800, 0, -1450]*u.mm,
            "boom_cog" :        [8600, 0, -500]*u.mm,
            "boom_pedestal1" :  [0, 1280, 0]*u.mm,
            "boom_pedestal2" :  [0, -1280, 0]*u.mm}

    jib = {"jib_cyljib" :       [1631, 0, 1236]*u.mm,
           "jib_load" :         [14957, 0, 1150]*u.mm,
           "jib_cog" :          [5200, 0, 600]*u.mm,
           "jib_boom1" :        [0, 697, 0]*u.mm,
           "jib_boom2" :        [0, -697, 0]*u.mm}

    pedestal = {"pedestal_boom" :    [0,0,0]*u.mm,
                "pedestal_cylboom" :  [1300,0,-2500]*u.mm}

    boomjib = {**boom,**jib}
    pbj =     {**pedestal,**boomjib}

    # Component assembly
    # ------------------
    
    RotateY     (jib,       180*u.deg+beta)
    Offset      (jib,       pbj['boom_jib'])
    RotateY     (boomjib,   alfa)
    RotateZ     (pbj,       gamma)
    
    return pbj

def calculate(geom):
    
    # Define constraints
    # ------------------

    CS = {"cs_x": [1.,0.,0.]*u.t,
          "cs_y": [0.,1.,0.]*u.t,
          "cs_z": [0.,0.,1.]*u.t}
    RotateZ(CS,gamma)
    x,y,z = CS['cs_x'], CS['cs_y'], CS['cs_z']
    
    
    # Constraint 1
    d1 = geom["boom_pedestal1"]
    d2 = geom["boom_pedestal2"]
    d3 = geom["boom_cylboom"]
    c = Unit(geom['boom_cylboom'],geom["pedestal_cylboom"])*u.t

    constraintsFD1 = [[x,d1],
                      [y,d1],
                      [z,d1],
                      [x,d2],
                      [z,d2],
                      [c,d3]]
    constraintsFM1 = np.transpose(np.vstack(FD2FM(constraintsFD1))) # losing units here: +cross, ~stacks

    # Constraint 2
    d1 = geom['jib_boom1']
    d2 = geom['jib_boom2']
    d3 = geom['jib_cyljib']
    c = Unit(geom['jib_cyljib'],geom['boom_cyljib'])*u.t

    constraintsFD2 = [[x,d1],
                      [y,d1],
                      [z,d1],
                      [x,d2],
                      [z,d2],
                      [c,d3]]
    constraintsFM2 = np.transpose(np.vstack(FD2FM(constraintsFD2))) # losing units here: +cross, ~stacks

    # Define loads
    # ------------

    # Loads 1
    loadsFD1 = [[acc*m_jib,     geom['jib_cog']],
             [acc*m_boom,    geom['boom_cog']],
             [acc*m_cyl*2,   np.add(geom['boom_cyljib'],geom['jib_cyljib'])/2],
             [acc*m_load,    geom['jib_load']]]
    loadsFM1 = sum(FD2FM(loadsFD1))

    # Loads 2
    loadsFD2 = [[acc*m_jib,     geom['jib_cog']],
             [acc*m_load,    geom['jib_load']]]
    loadsFM2 = sum(FD2FM(loadsFD2))
    
    # Solve for constraints
    # ---------------------

    reaction1 = np.linalg.solve(constraintsFM1,loadsFM1)
    reaction2 = np.linalg.solve(constraintsFM2,loadsFM2)

    # Post process
    # ------------
    for a in [loadsFM2, constraintsFM2, reaction2,loadsFM1, constraintsFM1, reaction1]:
        logging.debug(np.array_str(a,max_line_width=120))

    return {"alfa":alfa.to(u.deg).magnitude,
            "beta":beta.to(u.deg).magnitude,
            "gamma":gamma.to(u.deg).magnitude,
            "m_load":m_load.to(u.t).magnitude,
            "reachx":Length(geom['pedestal_boom'][0:2],geom['jib_load'][0:2]).to(u.m).magnitude,
            "reachz":geom['jib_load'][2].to(u.m).magnitude,
            "cblen":Length(geom['boom_cylboom'],geom['pedestal_cylboom']).to(u.m).magnitude,
            "cjlen":Length(geom['jib_cyljib'],geom['boom_cyljib']).to(u.m).magnitude,
            **dict(zip(["FXb","FYb","FZb","MXb","MYb","MZb"],loadsFM1)),
            **dict(zip(["RX1b","RY1b","RZ1b","RX2b","RZ2b","RC3b"],reaction1)),
            **dict(zip(["FXj","FYj","FZj","MXj","MYj","MZj"],loadsFM2)),
            **dict(zip(["RX1j","RY1j","RZ1j","RX2j","RZ2j","RC3j"],reaction2))}

def iterator():
    global m_load
##    m_load = 0*u.t
    g = geometry()
    while True:
        r = calculate(g)
        q1 = r["RC3b"] <   89*2
        q2 = r["RC3b"] > -124*2
        q3 = r["RC3j"] <   89*2
        q4 = r["RC3j"] > -124*2
        if q1 and q2 and q3 and q4:
            logging.info(f"passed: {alfa} {beta} {m_load}")
            m_load += 1*u.t
            good_r = r
        else:
            logging.info(f"final: {alfa} {beta} {m_load}")
            return good_r

def main(key=None):
    global alfa, beta, gamma, accs, m_jib, m_boom, m_cyl, m_load, acc, df, data

    alfas =  list(range(0,  75+1, 5))*u.deg # boom angle
    betas =  list(range(0, 120+1, 5))*u.deg # jib angle
    gammas = list(range(0, 360+1, 45))*u.deg # slew angle
   
    m_jib = 6*u.t
    m_boom = 16*u.t
    m_cyl = 2.6*u.t
    m_loads = [5,10]*u.t

    acc = np.array([-.48, 0.0, -1.17])
    data = []

    if key=='all':
        for gamma in gammas:
            for beta in betas:
                for alfa in alfas:
                    for m_load in m_loads:
                        data.append(calculate(geometry()))
        df = pd.DataFrame(data)
        df.to_pickle('fbd')
    elif key=='iterator':
        for gamma in gammas:
            for beta in betas:
                for alfa in alfas:
                    m_load = 0*u.t
                    data.append(iterator())
        df = pd.DataFrame(data)
        df.to_pickle('fbdi')
    else:
        gamma = 0*u.deg
        beta = 90*u.deg
        alfa = 35*u.deg
        m_load = 5*u.t
        data.append(calculate(geometry()))
        df = pd.DataFrame(data)
   
       
# logging.basicConfig(level=logging.INFO)
if __name__ == '__main__':
    main(key='iterator')
# plot_geometry(geometry(),xaxis=0,yaxis=2,width=39,height=30)
##pd.DataFrame(geometry()).T.plot.scatter(0,2, subplots=True)
##pd.DataFrame(geometry()).T.plot.scatter(1,2, subplots=True)
##pd.DataFrame(geometry()).T.plot.scatter(0,1, subplots=True)
plt.show()

