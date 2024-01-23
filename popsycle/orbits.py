import os, sys, math
import numpy as np
from astropy.table import Table, Column, MaskedColumn
import matplotlib.pyplot as plt
import astropy.constants as c
from astropy import units as u
from copy import deepcopy
from random import choice


class Orbit(object):
    """
    Exerpted from orbits.py from gcworks
    
    Methods
    -------
    kep2xyz :
	Generates (r, v, a) in AU, AU/yr, and AU/yr^2 respectively from keplarian parameters.
    
    eccen_anomaly :
        Calculates the eccentricity anomaly.
    """
    def kep2xyz(self, epochs):
        """
        Generates (r, v, a) in AU, AU/yr, and AU/yr^2 respectively from keplarian parameters. 
        
        epoch: numpy array
           Array of times (in years) that you would like to compute the parameters for
        
        mass: float or array-like
           Mass of primary object
        

        An example call is:

        >>> orb = orbits.Orbit()
        orb.w = omega # [degrees]
        orb.o = bigOm # [degrees]
        orb.i = incl # [degrees]
        orb.e = e_mag # float between 0 and 1
        orb.p = p # [years]
        orb.t0 = t0 # [years] This is initial
        orb.mass = mass # [Msun]

        (r, v, a) = orb.kep2xyz(array([refTime]))
        
        """

        GM = self.mass * c.G.to("cm3/(Msun s2)")

        epoch_num = len(epochs)

        # meanMotion in radians per year
        meanMotion = 2.0 * math.pi / self.p
        
        # Semi-major axis in AU
        axis = (self.p**2 * self.mass)**(1.0/3.0)

        ecc_sqrt = np.sqrt(1.0 - self.e**2)        
        
        #switch angular quantities to radians
        w = math.radians(self.w)
        o = math.radians(self.o)
        i = math.radians(self.i)
        
        # Mean anomaly
        mean_anomaly = meanMotion * (epochs - self.t0)

        #----------
        # Now for each epoch we compute the x and y positions
        #----------

        # Eccentric anomaly
        E = self.eccen_anomaly(mean_anomaly, self.e)
            
        cos_E = np.cos(E)
        sin_E = np.sin(E)
        
        Edot = meanMotion / (1.0 - (self.e * cos_E))

        X = cos_E - self.e
        Y = ecc_sqrt * sin_E

        #----------
        # Calculate Thiele-Innes Constants
        #----------
        cos_bigOm = np.cos(o)
        sin_bigOm = np.sin(o)
        cos_i = np.cos(i)
        sin_i = np.sin(i)
        
        cos_om = np.cos(w)
        sin_om = np.sin(w)

        self.conA = axis * (cos_om * cos_bigOm  - sin_om * sin_bigOm * cos_i)
        self.conB = axis * (cos_om * sin_bigOm  + sin_om * cos_bigOm * cos_i)
        self.conC = axis * (sin_om * sin_i)
        self.conF = axis * (-sin_om * cos_bigOm - cos_om * sin_bigOm * cos_i)
        self.conG = axis * (-sin_om * sin_bigOm + cos_om * cos_bigOm * cos_i)
        self.conH = axis * (cos_om * sin_i)
        
        # initialize zero arrays for r, v, and a
        r = np.zeros((epoch_num, 3), dtype='float64')
        v = np.zeros((epoch_num, 3), dtype='float64')
        a = np.zeros((epoch_num, 3), dtype='float64')

        r[:,0] = (self.conB * X) + (self.conG * Y)
        r[:,1] = (self.conA * X) + (self.conF * Y)
        r[:,2] = (self.conC * X) + (self.conH * Y)
        
        v[:,0] = Edot * ((-self.conB * sin_E) + (self.conG * ecc_sqrt * cos_E))
        v[:,1] = Edot * ((-self.conA * sin_E) + (self.conF * ecc_sqrt * cos_E))
        v[:,2] = Edot * ((-self.conC * sin_E) + (self.conH * ecc_sqrt * cos_E))
        
    
        # Calculate accleration
        for ii in range(epoch_num):
            rmag_cm = (np.sqrt( (r[ii,:]**2).sum() )*(u.au)).to("cm").value
            a[ii,:] = -GM * (r[ii,:]*(u.au)).to("cm").value / rmag_cm**3
        
        # from cm/s^2 to AU/yr^2
        a = (a*(u.cm/u.s**2)).to("au/yr2").value

        return (r, v, a)


    def eccen_anomaly(self, m, ecc, thresh=1e-10):
        """
        Calculates the eccentricity anomaly
        
	Parameters
	----------        
	m: numpy array
            Mean anomalies
            
        ecc: float between 0-1
            The eccentricity of the orbit
        """
        # set default values

        if (ecc < 0. or ecc >= 1.):
            raise Exception('Eccentricity must be >= 0 and < 1')
            
        #
        # Range reduction of m to -pi < m <= pi
        #
        mx = m.copy()

        ## ... m > pi
        zz = (np.where(mx > math.pi))[0]
        mx[zz] = mx[zz] % (2.0 * math.pi)
        zz = (np.where(mx > math.pi))[0]
        mx[zz] = mx[zz] - (2.0 * math.pi)
        
        # ... m < -pi
        zz = (np.where(mx <= -math.pi))[0]
        mx[zz] = mx[zz] % (2.0 * math.pi)
        zz = (np.where(mx <= -math.pi))[0]
        mx[zz] = mx[zz] + (2.0 * math.pi)

        #
        # Bail out for circular orbits...
        #
        if (ecc == 0.0):
            return mx

        aux   = (4.0 * ecc) + 0.50
        alpha = (1.0 - ecc) / aux

        beta = mx/(2.0*aux)
        aux = np.sqrt(beta**2 + alpha**3)
   
        z=beta+aux
        zz=(np.where(z <= 0.0))[0]
        z[zz]=beta[zz]-aux[zz]

        test=abs(z)**0.3333333333333333

        z =  test.copy()
        zz = (np.where(z < 0.0))[0]
        z[zz] = -z[zz]

        s0=z-alpha/z
        s1 = s0-(0.0780 * s0**5) / (1.0 + ecc)
        e0 = mx + ecc*((3.0 * s1) - (4.0 * s1**3))

        se0=np.sin(e0)
        ce0=np.cos(e0)

        f  = e0 - (ecc*se0) - mx
        f1 = 1.0 - (ecc*ce0)
        f2 = ecc*se0
        f3 = ecc*ce0
        f4 = -1.0 * f2
        u1 = -1.0 * f/f1
        u2 = -1.0 * f/(f1 + 0.50*f2*u1)
        u3 = -1.0 * f/(f1 + 0.50*f2*u2
                 + 0.166666666666670*f3*u2*u2)
        u4 = -1.0 * f/(f1 + 0.50*f2*u3
                 + 0.166666666666670*f3*u3*u3
                 + 0.0416666666666670*f4*u3**3)

        eccanom=e0+u4

        zz = (np.where(eccanom >= 2.00*math.pi))[0]
        eccanom[zz]=eccanom[zz]-2.00*math.pi
        zz = (np.where(eccanom < 0.0))[0]
        eccanom[zz]=eccanom[zz]+2.00*math.pi

        # Now get more precise solution using Newton Raphson method
        # for those times when the Kepler equation is not yet solved
        # to better than 1e-10
        # (modification J. Wilms)

        mmm = mx.copy()
        ndx = (np.where(mmm < 0.))[0]
        mmm[ndx] += (2.0 * math.pi)
        diff = eccanom - ecc*np.sin(eccanom) - mmm

        ndx = (np.where(abs(diff) > 1e-10))[0]
        for i in ndx:
            # E-e sinE-M
            fe = eccanom[i]-ecc*np.sin(eccanom[i])-mmm[i]
            # f' = 1-e*cosE
            fs = 1.0 - ecc*np.cos(eccanom[i])
            oldval=eccanom[i]
            eccanom[i]=oldval-fe/fs

            loopCount = 0
            while (abs(oldval-eccanom[i]) > thresh):
                # E-e sinE-M
                fe = eccanom[i]-ecc*np.sin(eccanom[i])-mmm[i]
                # f' = 1-e*cosE
                fs = 1.0 - ecc*np.cos(eccanom[i])
                oldval=eccanom[i]
                eccanom[i]=oldval-fe/fs
                loopCount += 1
                
                if (loopCount > 10**6):
                    msg = 'eccen_anomaly: Could not converge for e = %f' % ecc
                    raise EccAnomalyError(msg)

            while (eccanom[i] >=  math.pi):
                eccanom[i] = eccanom[i] - (2.0 * math.pi)
                
            while (eccanom[i] < -math.pi ):
                eccanom[i] = eccanom[i] + (2.0 * math.pi)

        return eccanom
    
class EccAnomalyError(Exception):
    def __init__(self, message):
        self.message = message

def a_to_P(mass, a):
    """
    Goes from semimajor axis in AU to period in years
    
    Parameters
    ---------- 
    mass: float or array-like
        Primary object mass in Msun.
    
    a: float or array-like
        Semimajor axis in AU.
        
    Returns
    ------- 
    period: float or array-like
        Orbital period in years.
    """
    
    G_units = c.G.to("AU3/(M_sun*year2)").value
    period = (a**3*4*(np.pi**2)/G_units/mass)**(1/2)
    return period

def add_positions(ss):
    """
    Adds x and y positions randomly in a box of length and width 40000 AU for each system.
    
    Parameters
    ----------   
    ss: astropy table
        Star system table without positions
        
    Returns
    -------
    ss_temp: astropy table
        Star system table with positions added
    """
    ss_temp = deepcopy(ss)
    
    ss_temp.add_column( Column(np.zeros(len(ss_temp), dtype=float), name='x', description='AU') )
    ss_temp.add_column( Column(np.zeros(len(ss_temp), dtype=float), name='y', description='AU') )
    
    sign_x = np.array([choice([-1,1]) for i in range(len(ss_temp))])
    sign_y = np.array([choice([-1,1]) for i in range(len(ss_temp))])
    ss_temp['x'] = sign_x*20000*np.random.rand(len(ss_temp))
    ss_temp['y'] = sign_y*20000*np.random.rand(len(ss_temp))
    
    return ss_temp

def add_mult_positions(companions, ss_pos, logAge, add_velocities = False):
    """
    Adds x and y positions of multiple companions by transforming keplerian parameters to xyz in AU
    using code origially from gcworks and random initial times. Then adding them to the random posiiton of the primary object.
    
    Parameters
    ----------    
    ss_pos: astropy table
        Star system table with positions added with add_positions()
    
    companion: astropy table
        Companion table without positions
    
    logAge: float or int
        Log of age of cluster with age in years
        
    Returns
    -------
    companion_temp: astropy table
        Companion table with positions added
        
    """
    companions_temp = deepcopy(companions)
    
    companions_temp.add_column( Column(np.zeros(len(companions_temp), dtype=float), name='x', description='AU') )
    companions_temp.add_column( Column(np.zeros(len(companions_temp), dtype=float), name='y', description='AU') )
    
    if add_velocities == True:
        companions_temp.add_column( Column(np.zeros(len(companions_temp), dtype=float), name='v_x', description='AU/yr') )
        companions_temp.add_column( Column(np.zeros(len(companions_temp), dtype=float), name='v_y', description='AU/yr') )
        
    orb = Orbit()
    for i in companions_temp:
        orb.w = i['omega'] #degrees
        orb.o = i['Omega'] #degrees
        orb.i = i['i'] #degrees
        orb.e = i['e'] #between 0 and 1
        orb.p = a_to_P(ss_pos[i['system_idx']]['systemMass'],10**i['log_a']) #year
        orb.t0 = (10**logAge)*np.random.rand() #year
        orb.mass = ss_pos[i['system_idx']]['mass']
        
        (r, v, a) = orb.kep2xyz(np.array([10**logAge]))
               
        x = r[0][0]
        y = r[0][1]
        
        #putting positions relative to primary object
        i['x'] = ss_pos[i['system_idx']]['x'] + x #AU
        i['y'] = ss_pos[i['system_idx']]['y'] + y #AU
        
        if add_velocities == True:
            i['v_x'] = v[0][0] #AU/yr
            i['v_y'] = v[0][1] #AU/yr
    
    
    return companions_temp

def distance_to_center_of_mass(ss_pos, companions_pos):
    """
    Adds extra column to star system and companions table with x and y distance to the center of mass in AU.
    Assumes hierarchical triples (two closest stars orbit their center of mass and triple orbits them) and no quads+.
    
    Parameters
    ----------    
    ss_pos: astropy table
        Star system table with positions added with add_positions()
    
    companion_pos: astropy table
        Companion table with positions added with add_mult_positions()
        
    Returns
    ------- 
    ss_pos_temp: astropy table
        Star system table with distance to center of mass in AU added
    
    companion_pos_temp: astropy table
        Companion table with distance to center of mass in AU added
    """
    
    ss_pos_temp = deepcopy(ss_pos)
    companions_pos_temp = deepcopy(companions_pos)
    
    ss_pos_temp.add_column( Column(np.zeros(len(ss_pos_temp), dtype=float), name='com_x') )
    ss_pos_temp.add_column( Column(np.zeros(len(ss_pos_temp), dtype=float), name='com_y') )
    companions_pos_temp.add_column( Column(np.zeros(len(companions_pos_temp), dtype=float), name='com_x') )
    companions_pos_temp.add_column( Column(np.zeros(len(companions_pos_temp), dtype=float), name='com_y') )
    
    companions_pos_temp.sort(['system_idx','log_a'])
    
    for i in range(len(ss_pos_temp)):
        if ss_pos_temp[i]['isMultiple'] == True:
            companion_indicies = np.where(companions_pos_temp['system_idx'] == i)[0]
            
            primary_mass = ss_pos_temp[i]['mass']
            companion_masses = companions_pos_temp[companion_indicies]['mass']
            
            primary_x = ss_pos_temp[i]['x']
            companion_x = companions_pos_temp[companion_indicies]['x']
            com_x = (primary_mass*primary_x + companion_masses[0]*companion_x[0])/(primary_mass + companion_masses[0])
            
            primary_y = ss_pos_temp[i]['y']
            companion_y = companions_pos_temp[companion_indicies]['y']
            com_y = (primary_mass*primary_y + companion_masses[0]*companion_y[0])/(primary_mass + companion_masses[0])
            
            ss_pos_temp[i]['com_x'] = com_x - primary_x
            companions_pos_temp[companion_indicies[0]]['com_x'] = com_x - companion_x[0]
            ss_pos_temp[i]['com_y'] = com_y - primary_y
            companions_pos_temp[companion_indicies[0]]['com_y'] = com_y - companion_y[0]
            
            # Assumes hierarchical triples 
            if len(companion_indicies) == 2:
                center_mass = primary_mass + companion_masses[0]
                com_x_out = (center_mass*com_x + companion_masses[1]*companion_x[1])/(center_mass + companion_masses[1])
                com_y_out = (center_mass*com_y + companion_masses[1]*companion_y[1])/(center_mass + companion_masses[1])
                
                companions_pos_temp[companion_indicies[1]]['com_x'] = com_x_out - companion_x[1]
                companions_pos_temp[companion_indicies[1]]['com_y'] = com_y_out - companion_y[1]
                
    return ss_pos_temp, companions_pos_temp
                
            
def plot_projected_cluster(ss_pos, companions_pos):
    """
    Plots projected cluster with lines between companions and primary stars
    
    Parameters
    ----------   
    ss_pos: astropy table
        Star system table with positions added with add_positions()
    
    companion_pos: astropy table
        Companion table with positions added with add_mult_positions()
    """
    plt.figure(figsize=(10,10))
    plt.plot(ss_pos['x'], ss_pos['y'],linestyle='none',marker='o' )
    plt.plot(companions_pos['x'], companions_pos['y'],linestyle='none',marker='.' )
    
    #makes lines between companion and primary star
    for i in companions_pos:
        plt.plot([i['x'], ss_pos[i['system_idx']]['x']],[i['y'], ss_pos[i['system_idx']]['y']],color='grey',linewidth=1)
        
    plt.xlabel("x (AU)")
    plt.ylabel("y (AU)")
    plt.show()
    
    return

def plot_companion_orbit(ss, companions_pos, logAge, t0 = None, system = None):
    """
    Plots the orbit of one system assuming the primary object is at (0,0). By default random companion and initial time.
    
    Parameters
    ----------       
    ss: astropy table
        Star system table (does not matter if it has positions or not)
    
    companion_pos: astropy table
        Companion table with positions added with add_mult_positions()
    
    logAge: float or int
        Log of age of cluster with age in years
        
    t0: float or int, optional
        Initial time of creation of the system in years.
        Default random.
    
    system: int, optional
        Index of desired companion in companion_pos table.
        Default random.
    """
    
    if system == None:
        system = np.random.randint(len(companions_pos))
    if t0 == None:
        t0 = (10**logAge)*np.random.rand() #year

    companion = companions_pos[system]
    
    orb = Orbit()
    orb.w = companion['omega'] #degrees
    orb.o = companion['Omega'] #degrees
    orb.i = companion['i'] #degrees
    orb.e = companion['e'] #between 0 and 1
    orb.p = a_to_P(ss[companion['system_idx']]['mass'],10**companion['log_a']) #year
    orb.t0 = t0 #year
        
    (r, v, a) = orb.kep2xyz(np.linspace(1, 10**logAge), mass=ss[companion['system_idx']]['mass'])
        
    for i in r:
        plt.plot(i[0],i[1], marker='.', color = 'blue')
        if i[0] == r[-1][0]:
            plt.plot(i[0],i[1], marker='*', color = 'gold', markersize=15)
            
    plt.xlabel("x (AU)")
    plt.ylabel("y (AU)")
    plt.show()
    
    return
