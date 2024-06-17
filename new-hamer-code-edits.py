import numpy as np
import h5py
from scipy.spatial import KDTree
import time, re, json
from numpyencoder import NumpyEncoder
from mpmath import mp, mpf, nstr
from mpmath import isnan as mp_isnan
from astropy.table import Table

##########
# Conversions
##########
masyr_to_degday = 1.0 * (1.0e-3 / 3600.0) * (1.0 / 365.25)
kms_to_kpcs = 1.0 * (3.086 * 10**16) ** -1
kms_to_kpcday = 1.0 * (3.086 * 10**16) ** -1 * 86400.0
au_to_kpc = 4.848 * 10**-9

# (15 min / 1825 days)^2 where 15 min is Roman's fiducial cadence
# and 1825 days is Roman's fiducial observational baseline
## change: instead of using Roman time, we will use obs_time as specified by calc_events (DAYS)
transit15minSq = mpf('3.25786e-11')

# set precision
mp.dps = 30

# region coordinate transformations

def GLat_exact(px, py, pz):
    return np.arcsin(pz / np.sqrt(px**2 + py**2 + pz**2))

def GLon_exact(px, py, pz):
    return np.arctan(py / px)

def spherical_exact(px, py, pz):
    return(
        np.sqrt(px**2 + py**2 + pz**2),
        GLat_exact(px, py, pz),
        GLon_exact(px, py, pz)
    )

def cartesian_exact(rad, glat, glon):
    return(
        rad * np.cos(glat * np.pi / 180) * np.cos(glon * np.pi / 180),        
        rad * np.cos(glat * np.pi / 180) * np.sin(glon * np.pi / 180),
        rad * np.sin(glat * np.pi / 180)
    )

#endregion

# this holds all the information about some circle with a linear trajectory
# r is radius, p1 is initial center of circle, v is 2D velocity vector, t is duration
# returns (circle1, circle2, rectangle points)
def RRectPath(r, p1, v, t):
	# p1 (p2) is center of circle at the beginning (end)
    p2 = p1 + v * t

    # v is 2D vector of proper motion, if v = 0, then the shape is just a 
    # stationary circle and zero-size rectangle
    if (v[0] == 0 and v[1] == 0):
        return ((p1, r), (p1, r), (p1, p1, p1, p1))

    # unit vector perpendicular to trajectory
    uPerp = np.array([-v[1], v[0]]) / np.sqrt(np.dot(v,v))

    # boundary points of rectangular portion
    rectPoints = (p1 + uPerp * r, p2 + uPerp * r, p2 - uPerp * r, p1 - uPerp * r)

    return ((p1, r), (p2, r), rectPoints)

# after writing many functions, had to start using precision numbers.
# this is a decorator to map all the arguments of a function f(x1,x2...) --> f(mpf(x1),mpf(x2)...)
def make_precise(func):
    def inner(*args):
        return(func(*map(mpf, args)))

    return inner

# the squared difference between  the two solutions of the quadratic equation
# formulas here were solved analytically in mathematica and copy/pasted
@make_precise
def tSolSqDiff(r1, r2, x1, y1, vx1, vy1, x2, y2, vx2, vy2):
    # -r1 -r2 + sqrt(...) < 0 condition
    if (np.power(r1 + r2,2) <= (np.power(vy1*x1 - vy2*x1 - vy1*x2 + vy2*x2 - vx1*y1 + vx2*y1 + vx1*y2 - 
        vx2*y2,2)/
        (np.power(vx1,2) - 2*vx1*vx2 + np.power(vx2,2) + np.power(vy1,2) - 2*vy1*vy2 + 
        np.power(vy2,2)))):
        # no solution
        return(np.nan)
    

    return(4*(-4*r1*r2*vx1*vx2 - 4*r1*r2*vy1*vy2 - 4*vy1*vy2*x1*x2 + 2*vx1*vy1*x1*y1 - 2*vx2*vy1*x1*y1 - 2*vx1*vy2*x1*y1 + 
     2*vx2*vy2*x1*y1 - 2*vx1*vy1*x2*y1 + 2*vx2*vy1*x2*y1 + 2*vx1*vy2*x2*y1 - 2*vx2*vy2*x2*y1 - 2*vx1*vy1*x1*y2 + 
     2*vx2*vy1*x1*y2 + 2*vx1*vy2*x1*y2 - 2*vx2*vy2*x1*y2 + 2*vx1*vy1*x2*y2 - 2*vx2*vy1*x2*y2 - 2*vx1*vy2*x2*y2 + 
     2*vx2*vy2*x2*y2 - 4*vx1*vx2*y1*y2 - 2*vx1*vx2*np.power(r1,2) - 2*vy1*vy2*np.power(r1,2) - 2*vx1*vx2*np.power(r2,2) - 
     2*vy1*vy2*np.power(r2,2) + 2*r1*r2*np.power(vx1,2) + 2*y1*y2*np.power(vx1,2) + np.power(r1,2)*np.power(vx1,2) + 
     np.power(r2,2)*np.power(vx1,2) + 2*r1*r2*np.power(vx2,2) + 2*y1*y2*np.power(vx2,2) + 
     np.power(r1,2)*np.power(vx2,2) + np.power(r2,2)*np.power(vx2,2) + 2*r1*r2*np.power(vy1,2) + 
     2*x1*x2*np.power(vy1,2) + np.power(r1,2)*np.power(vy1,2) + np.power(r2,2)*np.power(vy1,2) + 
     2*r1*r2*np.power(vy2,2) + 2*x1*x2*np.power(vy2,2) + np.power(r1,2)*np.power(vy2,2) + 
     np.power(r2,2)*np.power(vy2,2) + 2*vy1*vy2*np.power(x1,2) - np.power(vy1,2)*np.power(x1,2) - 
     np.power(vy2,2)*np.power(x1,2) + 2*vy1*vy2*np.power(x2,2) - np.power(vy1,2)*np.power(x2,2) - 
     np.power(vy2,2)*np.power(x2,2) + 2*vx1*vx2*np.power(y1,2) - np.power(vx1,2)*np.power(y1,2) - 
     np.power(vx2,2)*np.power(y1,2) + 2*vx1*vx2*np.power(y2,2) - np.power(vx1,2)*np.power(y2,2) - 
     np.power(vx2,2)*np.power(y2,2))*np.power(-2*vx1*vx2 - 2*vy1*vy2 + np.power(vx1,2) + np.power(vx2,2) + 
     np.power(vy1,2) + np.power(vy2,2),-2))

# solve for the value of t when the circles are exactly touching
# solutions correspond to the beginning of event (first overlap of two circles)
# and end of event (when two circles stop overlapping)
# tSol used in rrQuadSolve to find the beginning and end t for an event
@make_precise
def tSol(solutionIndex, r1, r2, x1, y1, vx1, vy1, x2, y2, vx2, vy2):

    if (solutionIndex == 1):
        sqrtSign = -1
    elif(solutionIndex == 2):
        sqrtSign = 1
    else:
        raise Exception('solutionIndex must 1 or 2, corresponding to the two different solutions of the quadratic equation')

    # check condition for solution existing

    cond1 = r1 + r2 + np.sqrt(np.power(vy1*x1 - vy2*x1 - vy1*x2 + vy2*x2 - vx1*y1 + vx2*y1 + vx1*y2 - 
        vx2*y2,2)/
        (np.power(vx1,2) - 2*vx1*vx2 + np.power(vx2,2) + np.power(vy1,2) - 2*vy1*vy2 + 
        np.power(vy2,2)))

    cond2 = -r1 - r2 + np.sqrt(np.power(vy1*x1 - vy2*x1 - vy1*x2 + vy2*x2 - vx1*y1 + vx2*y1 + vx1*y2 - 
        vx2*y2,2)/
        (np.power(vx1,2) - 2*vx1*vx2 + np.power(vx2,2) + np.power(vy1,2) - 2*vy1*vy2 + 
        np.power(vy2,2)))

    if (cond1 > 0 and cond2 > 0):
        return(np.nan)
        

    return(-2*vx1*x1 + 2*vx2*x1 + 2*vx1*x2 - 
		2*vx2*x2 - 2*vy1*y1 + 
		2*vy2*y1 + 2*vy1*y2 - 
		2*vy2*y2 + 
		sqrtSign * np.sqrt(np.power(2*vx1*x1 - 
		2*vx2*x1 - 2*vx1*x2 + 
		2*vx2*x2 + 2*vy1*y1 - 
		2*vy2*y1 - 2*vy1*y2 + 
		2*vy2*y2,2) - 
		4*(np.power(vx1,2) - 
		2*vx1*vx2 + 
		np.power(vx2,2) + 
		np.power(vy1,2) - 
		2*vy1*vy2 + np.power(vy2,2))*
		(-np.power(r1,2) - 2*r1*r2 - 
		np.power(r2,2) + 
		np.power(x1,2) - 2*x1*x2 + 
		np.power(x2,2) + 
		np.power(y1,2) - 2*y1*y2 + 
		np.power(y2,2))))/(
		2.*(np.power(vx1,2) - 2*vx1*vx2 + 
		np.power(vx2,2) + 
		np.power(vy1,2) - 2*vy1*vy2 + 
		np.power(vy2,2))
	)

# solve for the times t1,t2 where two circles overlap given two RRectPaths
# only done after comparing the square differences
def rrQuadSolve(rr1, rr2):

	d1 = rr1[1][0] - rr1[0][0]
	d2 = rr2[1][0] - rr2[0][0]

	t1 = tSol(1, rr1[0][1], rr2[0][1], *rr1[0][0], *d1, *rr2[0][0], *d2)

	t2 = tSol(2, rr1[0][1], rr2[0][1],
		*rr1[0][0], *d1,
		*rr2[0][0], *d2
	)
	return((t1, t2))


# xxx: does this implement a t = -T/2 start?
def rrQuadDiff(rr1, rr2):

	d1 = rr1[1][0] - rr1[0][0]
	d2 = rr2[1][0] - rr2[0][0]
        
	return(tSolSqDiff(
        rr1[0][1], rr2[0][1],
		*rr1[0][0], *d1,
		*rr2[0][0], *d2
	))

# pretend all stars have the radius of the sun
# returns angular size in radians
## change: use actual radius of star -- ask how to calculate this from other attributes of the star
solar_radius_kpc = 2.26E-11 ## this is already changed in synthetic.py, still need to change it here...
def star_size(rad): ## rad is radial distance away
    return solar_radius_kpc / rad

# mass in solar masses, distances in kpc (all galaxia defaults)
# returns angular in radians
def einstein_radius(d_lens, d_source, mass):
    # https://en.wikipedia.org/wiki/Einstein_radius

    d_LS = np.abs(d_source - d_lens)

    theta_E_arcsec = np.sqrt(mass/np.power(10, 11.09)) / np.sqrt(d_lens * d_source / (d_LS * 1E6))

    if(np.isnan(theta_E_arcsec)):
        print(d_LS)
        raise Exception('somethings gone wrong, probably precision lost')

    return(np.pi/180 * theta_E_arcsec/3600)

## take a patch as made in popsycle, make the kdtree from all objects in the patch, iterate over the whole patch
def processPatch_quad_starTree(patch, duration):

    sources = patch[:]
    lenses = patch[:]
    # diagnostic: num sources, num lenses
    print('num sources: %d, lenses: %d' % (len(sources), len(lenses))) ## should be the same number
    
    # location of final point in spherical coords
    def end_movement_spherical_noCartesian(pVec, vVec):
        return spherical_exact(
            cartesian_exact(pVec['rad'], pVec['glat'], pVec['glon'])[0] + vVec['vx'] * kms_to_kpcday * duration,
            cartesian_exact(pVec['rad'], pVec['glat'], pVec['glon'])[1] + vVec['vy'] * kms_to_kpcday * duration,
            cartesian_exact(pVec['rad'], pVec['glat'], pVec['glon'])[2] + vVec['vz'] * kms_to_kpcday * duration
        )
    
    def mid_movement_spherical_noCartesian(pVec, vVec):
        return spherical_exact(
            cartesian_exact(pVec['rad'], pVec['glat'], pVec['glon'])[0] + vVec['vx'] * kms_to_kpcday * duration/2,
            cartesian_exact(pVec['rad'], pVec['glat'], pVec['glon'])[1] + vVec['vy'] * kms_to_kpcday * duration/2,
            cartesian_exact(pVec['rad'], pVec['glat'], pVec['glon'])[2] + vVec['vz'] * kms_to_kpcday * duration/2
        )

    # end position of sources and lenses respectively
    endPosSph_sources = np.asarray(end_movement_spherical_noCartesian(sources[['rad','glat','glon']], sources[['vx', 'vy', 'vz']])).T
    endPosSph_lenses = np.asarray(end_movement_spherical_noCartesian(lenses[['rad','glat','glon']], lenses[['vx', 'vy', 'vz']])).T

    midPosSph_sources = np.asarray(mid_movement_spherical_noCartesian(sources[['rad','glat','glon']], sources[['vx', 'vy', 'vz']])).T
    midPosSph_lenses = np.asarray(mid_movement_spherical_noCartesian(lenses[['rad','glat','glon']], lenses[['vx', 'vy', 'vz']])).T

    def to_radians(rgg):
        return (rgg['rad'], rgg['glat']*np.pi/180, rgg['glon']*np.pi/180)

    # starting position of sources and lenses respectively
    startPosSph_sources = np.asarray(to_radians(sources[['rad','glat','glon']])).T
    startPosSph_lenses = np.asarray(to_radians(lenses[['rad','glat','glon']])).T
    

    ## array of total displacements for lenses, all displacement is related to lenses
    total_disp_arr = np.sqrt((startPosSph_lenses - endPosSph_lenses)[:, 1]**2 + (startPosSph_lenses - endPosSph_lenses)[:, 2]**2)
    max_disp = total_disp_arr.max()
    disp_95 = np.percentile(total_disp_arr, 95)
    radius_cut_95 = 2*disp_95 ## double because lens and source could each have about this much motion 
    radius_cut_max = 2*max_disp
    ## could have more, but that case is currently being neglected, as this will catch almost all cases (especially if taken from midpoint positions)
    print('==== displacements ====')    
    print('max: %s' % max_disp)
    print('95th percentile: %s' % disp_95)
    print('============================')


    ## change: I believe maxSphRadius is equivalent to radius_cut, so we can use the input from calc_events instead
    ## note: radius_cut is input with arcseconds, and maxSphRadius is currently in radians, so there needs to be a unit conversion
    maxSphRadius = 1e-5 # to do: find dynamically instead of by eye
    kdt = KDTree(midPosSph_sources[:, 1:3]) ## now uses midpoint instead of start to capture more relevant potential events
    print('number of objects (sources) in the kdtree: %d' % len(midPosSph_sources[:, 1:3]))

    startTime = time.time()
    totalNearbyObjects = 0 ## counts the number of nearby (potential) sources near a lens
    totalIsolatedObjects = 0
    totalLensingEvents = 0

    resultData = {'transitData': []}

    # assuming there are fewer lenses (ffps), we iterate over ffps but query sources from kdtree
    # to do: choose whether to loop over sources or lenses automatically based on which is smaller
    ## note: this doesnt matter because we will have the same # of potential sources and lenses, as each obj can be a lens or source

    lens_id = []
    sorc_id = []
    ## list for now cus appending is better, can always make it an array later
    less_95_case = 0 ## for testing purposes
    greater_95_case = 0
    for i, lens in enumerate(lenses):
        # if (i%1000 == 0):
        #     print('i at %d' % i)
        if total_disp_arr[i] <= disp_95:
            ## if the ith lens has a displacement lower than disp_95, we use the 95th percentile threshold value
            results = kdt.query_ball_point((midPosSph_lenses[i][1], midPosSph_lenses[i][2]), 2*disp_95) ## lens' glat, glon (coords)
            results = [index for index in results if index != i] ## dont count the same object as a source if its already the lens
            less_95_case += 1
        else: ##in the >95% case
            results = kdt.query_ball_point((midPosSph_lenses[i][1], midPosSph_lenses[i][2]), min(maxSphRadius, radius_cut_max))
            results = [index for index in results if index != i]
            greater_95_case += 1
        """if i in range(6):
            print('i is ', i, 'results are', results)
            print('lens id is', lens['obj_id'])
            print(lens)
            print(sources[i])
            print('==================================')"""
        # loop through nearby sources to see if trajectory contours overlap
        # if so, find duration of event
        for res in results:
            # res is the id in the original list of the star
            ## res is actually the INDEX of the object in the original patch list
            if(lens['rad'] < startPosSph_sources[res][0]): ## does this line mean that duplicates cant be found?
                # lens is nearer than the source
                totalNearbyObjects += 1
                lens_id.append(i) ## this is the index of the object in the original patch list
                sorc_id.append(res) ## this is the index of the object in the original patch list
                ## lens[0] for actual obj_id of lens
                ## this copies behavior of _calc_event_cands_radius, but only when fed a big patch

                ## ============ END of cands_radius ====================
                # starting coordinates
                lensCoord_start = startPosSph_lenses[i, 1:3]
                sourceCoord_start = startPosSph_sources[res][1:3]

                lensCoord_end  = endPosSph_lenses[i, 1:3]
                sourceCoord_end = endPosSph_sources[res][1:3]

                d_lens =  lensCoord_end - lensCoord_start
                d_source = sourceCoord_end - sourceCoord_start

                # this is here since the time in popsycle might go from t=-1/2 to t=1/2
                # so for purposes of comparing our output to theirs, make this the same.
                lensCoord_start = lensCoord_start - d_lens/2
                sourceCoord_start = sourceCoord_start - d_source/2

                # trajectory contours of lens and source respectively
                rr_lens = RRectPath(
                    2*einstein_radius(lens['rad'], startPosSph_sources[res][0], lens['mass']), ## change: should be system mass
                    np.array(lensCoord_start),
                    np.array(d_lens),
                1)
                ## 2 einstein radii will be an input that can be changed later, just use 2 as placeholder for now
                ## change: einstein_radius in synthetic.py is in MAS, nicks code assumes its in RADIANS
                rr_source = RRectPath(
                    star_size(startPosSph_sources[res][0]),
                    np.array(sourceCoord_start),
                    np.array(d_source),
                1)

                deltaTSq = rrQuadDiff(rr_lens, rr_source)

                ## change:(?) see if this try catch can be improved, but if not, its fine for now
                # some better code would remove the need for this try catch
                # but it works fine, so i'm leaving it in
                try:
                    if(np.isnan(deltaTSq)):
                        # the two circles never overlap at the same time
                        continue
                except:
                    # chill
                    pass

                if (deltaTSq < transit15minSq):
                    # event is too short duration to be seen with 15 min cadence
                    continue
                
                # possible collision, now solve for the exact times where the circles are touching.
                (t1, t2) = rrQuadSolve(rr_lens, rr_source)

                if((t1 >= 0 and t1 <= 1) or (t2 >= 0 and t2 <= 1)):

                    print('comparing delta t:   ',t2 - t1, np.sqrt(deltaTSq))

                    totalLensingEvents +=1
                    '''
                    print('======== lensing event ========')
                    print('start/end time of transit: %s, %s' % (t1, t2))
                    #print('source id: %s' % stars[res]['obj_id'])
                    #print('lens id: %s' % ffp['obj_id'])
                    #print('ffp index: %s' % i)
                    #print('star index: %s' % res)
                    print('einstein radius: %s' % einstein_radius(ffp['rad'], sphCoords[res][0], ffp['mass']))
                    print('star radius: %s' % star_size(sphCoords[res][0]))
                    print('star coord: %s' % starCoord)
                    print('dStar: %s' % d_star)
                    '''
                    
                    # start and end time of event
                    times_np = np.array(mp.matrix([t1,t2]).tolist(), dtype=np.float32)

                    # output data to write to FITS file
                    ## need to check how this fits table compares to the popsycle one
                    resultData['transitData'].append({
                        'id_lens': i,
                        'id_source': res,
                        'source obj_id': sources[res]['obj_id'],
                        'lens obj_id': lens['obj_id'],
                        't1': times_np[0],
                        't2': times_np[1],
                        'm_lens': lens['mass'],
                        'm_source': sources[res]['mass']
                    })

                    if(resultData['transitData'][-1]['id_source'] > len(sources)):
                        print('=!==!==!==!==!==!==!==!==!==!=')
                        print(resultData['transitData'][-1])
                        raise Exception('somehow source id larger than # of sources')
                    
        # there are no sources nearby this lens
        if len(results) == 0:
            totalIsolatedObjects += 1

    endTime = time.time()
    print('search radius is %s rad (adjust manually if necessary)' % maxSphRadius)
    print('completed search for transit events in %s s' % (endTime-startTime))
    print('total close sources found: %d. total lonely lenses found: %d' % (totalNearbyObjects, totalIsolatedObjects))
    ## totalNearbyObjects may still be higher than total # of objects because sources may be close to multiple other
    print('total lensing events: %s' % totalLensingEvents)
    print('less_95_case = %s' % less_95_case)
    print('greater_95_case = %s' % greater_95_case)
    ## change: print statements should be modified to match current implementation (late priority)
    print('length of sorc_id is', len(sorc_id))
    print(sorc_id[0:20])
    print('length of lens_id is', len(lens_id))
    print(lens_id[0:20])

    resultData['sourceCoordsStart'] = startPosSph_sources[ [tr['id_source'] for tr in resultData['transitData']] ]
    resultData['sourceCoordsEnd'] = endPosSph_sources[ [tr['id_source'] for tr in resultData['transitData']] ]

    resultData['lensCoordsStart'] = startPosSph_lenses[ [tr['id_lens'] for tr in resultData['transitData']] ]
    resultData['lensCoordsEnd'] = endPosSph_lenses[ [tr['id_lens'] for tr in resultData['transitData']] ]

    return(resultData)

## note: likely dont need this function at all because # of lenses = # of sources
## but if we do use it, same notes as for previous function apply
# identical to the above function, but lenses outnumber sources, so lenses are in the kd-tree
# and sources are looped over
def processPatch_quad_ffpTree(patch, duration):

    stars = patch[np.where(patch['rem_id'] == 0)]
    ffps  = patch[np.where(patch['rem_id'] == 104)] # 105 for ffps, our hard-coded PopSyCLE ids

    print('num stars: %d, ffps: %d' % (len(stars), len(ffps)))
    
    def end_movement_spherical_noCartesian(pVec, vVec):
        return spherical_exact(
            cartesian_exact(pVec['rad'], pVec['glat'], pVec['glon'])[0] + vVec['vx'] * kms_to_kpcday * duration,
            cartesian_exact(pVec['rad'], pVec['glat'], pVec['glon'])[1] + vVec['vy'] * kms_to_kpcday * duration,
            cartesian_exact(pVec['rad'], pVec['glat'], pVec['glon'])[2] + vVec['vz'] * kms_to_kpcday * duration
        )

    endPosSph = np.asarray(end_movement_spherical_noCartesian(stars[['rad','glat','glon']], stars[['vx', 'vy', 'vz']])).T
    endPosSph_ffps = np.asarray(end_movement_spherical_noCartesian(ffps[['rad','glat','glon']], ffps[['vx', 'vy', 'vz']])).T

    def to_radians(rgg):
        return (rgg['rad'], rgg['glat']*np.pi/180, rgg['glon']*np.pi/180)

    sphCoords = np.asarray(to_radians(stars[['rad','glat','glon']])).T
    sphCoords_ffps = np.asarray(to_radians(ffps[['rad','glat','glon']])).T
    
    print('==== max displacement ====')    
    print('stars: %s' % np.sqrt(np.max((sphCoords - endPosSph)[:, 1]**2 + (sphCoords - endPosSph)[:, 2]**2)))

    print('ffps: %s' % np.sqrt(np.max((sphCoords_ffps - endPosSph_ffps)[:, 1]**2 + (sphCoords_ffps - endPosSph_ffps)[:, 2]**2)))
    print('use these numbers to set the spherical search radius in the kdtree')
    print('============================')
    



    maxSphRadius = 1e-6
    kdt = KDTree(sphCoords_ffps[:, 1:3])
    print('number of guys in the kdtree: %d' % len(sphCoords_ffps[:, 1:3]))

    startTime = time.time()
    totalNearbyObjects = 0
    totalIsolatedObjects = 0
    totalLensingEvents = 0

    resultData = {'transitData': []}

    for i, star in enumerate(stars):
        # if (i%1000 == 0):
        #     print('i at %d' % i)
        '''
        print('star coords from the list')
        print((star['glat']*np.pi/180, star['glon']*np.pi/180))

        print('star coords from sphCoords')
        print(sphCoords[i,1:3])

        raise Exception('done')
        '''
        results = kdt.query_ball_point(sphCoords[i,1:3], maxSphRadius)

        
        
        for res in results:
            # res is the id in the original list of the ffps
            if(star['rad'] > sphCoords_ffps[res][0]):
                # lens is nearer than the source
                totalNearbyObjects += 1

                # starting coordiantes
                ffpCoord  = sphCoords_ffps[res, 1:3]
                ffpCoord_end  = endPosSph_ffps[res, 1:3]

                starCoord = sphCoords[i, 1:3]                
                starCoord_end = endPosSph[i, 1:3]

                d_ffp =  ffpCoord_end - ffpCoord
                d_star = starCoord_end - starCoord

                ffpCoord = ffpCoord - d_ffp/2
                starCoord = starCoord - d_star/2


                rr_ffp = RRectPath(
                    einstein_radius(sphCoords_ffps[res, 0], sphCoords[i, 0], ffps[res]['mass']),
                    np.array(ffpCoord),
                    np.array(d_ffp),
                1)

                rr_star = RRectPath(
                    star_size(sphCoords[i][0]),
                    np.array(starCoord),
                    np.array(d_star),
                1)

                deltaTSq = rrQuadDiff(rr_ffp, rr_star)

                if(mp_isnan(deltaTSq)):
                        # no collision
                        continue

                try:
                    if(np.isnan(deltaTSq)):
                        # no collision
                        continue
                except:
                    # chill
                    pass

                if (deltaTSq < transit15minSq):
                    # too short
                    continue
                
                (t1, t2) = rrQuadSolve(rr_ffp, rr_star)

                try:
                    if((t1 >= 0 and t1 <= 1) or (t2 >= 0 and t2 <= 1)):
                        print('hit')
                except TypeError as e:
                    print('this shouldn\'t happen')
                    print(e, t1, t2, deltaTSq)
                    continue
                

                if(t1 >= 0 and t2 <= 1):

                    print('comparing delta t:   ',t2 - t1, np.sqrt(deltaTSq))

                    totalLensingEvents +=1
                    '''
                    print('======== lensing event ========')
                    print('start/end time of transit: %s, %s' % (t1, t2))
                    #print('source id: %s' % stars[res]['obj_id'])
                    #print('lens id: %s' % ffp['obj_id'])
                    #print('ffp index: %s' % i)
                    #print('star index: %s' % res)
                    print('einstein radius: %s' % einstein_radius(ffp['rad'], sphCoords[res][0], ffp['mass']))
                    print('star radius: %s' % star_size(sphCoords[res][0]))
                    print('star coord: %s' % starCoord)
                    print('dStar: %s' % d_star)
                    '''
                    
                    times_np = np.array(mp.matrix([t1,t2]).tolist(), dtype=np.float32)

                    resultData['transitData'].append({
                        'id_ffp': res,
                        'id_star': i,
                        'star obj_id': star['obj_id'],
                        'ffp obj_id': ffps[res]['obj_id'],
                        't1': times_np[0],
                        't2': times_np[1],
                        'm_ffp': ffps[res]['mass'],
                        'm_star': star['mass']
                    })

                    if(resultData['transitData'][-1]['id_star'] > len(stars)):
                        print('=!==!==!==!==!==!==!==!==!==!=')
                        print(resultData['transitData'][-1])
                        raise Exception('somehow star id larger than # of stars')

        if len(results) == 0:
            totalIsolatedObjects += 1

    endTime = time.time()
    print('search radius is %s rad' % maxSphRadius)
    print('did search for transits in %s s' % (endTime-startTime))
    print('total close stars found: %d. total lonely ffps found: %d' % (totalNearbyObjects, totalIsolatedObjects))
    print('total lensing events: %s' % totalLensingEvents)

    resultData['starCoordsStart'] = sphCoords[ [tr['id_star'] for tr in resultData['transitData']] ]
    resultData['starCoordsEnd'] = endPosSph[ [tr['id_star'] for tr in resultData['transitData']] ]

    resultData['ffpCoordsStart'] = sphCoords_ffps[ [tr['id_ffp'] for tr in resultData['transitData']] ]
    resultData['ffpCoordsEnd'] = endPosSph_ffps[ [tr['id_ffp'] for tr in resultData['transitData']] ]

    return(resultData)

## do we want output to be a .json? dont we want a .fits table? ask about this.
## pretty sure processPatch_quad_ffp/starTree return fits files, but computeEvents turns them into a json file
## change: even so, these fit files are not complete with what we need
# takes in an hdf5 file from PopSyCLE with catalog of all objects and performs analysis
# writes output as json
def computeEvents(hdf5_filename):
    hf = h5py.File(hdf5_filename, "r")

    binRE = re.compile(r'l(\d{1,2})b(\d{1,2})')
    patchKeys = list(filter(binRE.match, hf.keys()))
    print('patches in this file are %s' % patchKeys)

    results = {}
    for k in patchKeys:
        ## with the old implementation, it only ran through the first small patch
        print('processing patch %s' % k)
        patchResult = processPatch_quad_starTree(hf[k], 1000) ## computeEvents is currently hard coded to duration of 1000
        results[k] = patchResult
    
    jsonFilename = '%s-mk2.json' % hdf5_filename[:-3]
    print('writing output to file: %s' % jsonFilename)

    with open(jsonFilename, 'w') as f:
        json.dump(results, f, cls=NumpyEncoder)

# actual run step
# to do: make filename an input
## computeEvents('pbh-binl9b9_earth_01_fdm_5_long.h5')
## commented out for ease of testing

# useful diagnostic function to compare to what PopSyCLE finds as events in
# a given patch. we find the same events exactly but faster
def convertFits(fn, relevantIdx=104):
    ## fn is a FITS table
    ## written such that it ONLY works for ffps. Could be modified later to be general, but I don't see a reason to (yet)
    ## change: if we want to use this function, need to modify relevantIdx and make it work for not only ffp lenses
    duration = 1825
    
    t = Table.read(fn)

    goodIdx = np.where(t['rem_id_L'] == relevantIdx)[0]
    goodRows = t[goodIdx]
    
    def end_movement_spherical_noCartesian(pVec, vVec, suffix):
        def s(str):
            return(str+suffix)

        return spherical_exact(
            cartesian_exact(pVec[s('rad')], pVec[s('glat')], pVec[s('glon')])[0] + vVec[s('vx')] * kms_to_kpcday * duration,
            cartesian_exact(pVec[s('rad')], pVec[s('glat')], pVec[s('glon')])[1] + vVec[s('vy')] * kms_to_kpcday * duration,
            cartesian_exact(pVec[s('rad')], pVec[s('glat')], pVec[s('glon')])[2] + vVec[s('vz')] * kms_to_kpcday * duration
        )

    endPosSph = np.asarray(end_movement_spherical_noCartesian(goodRows[['rad_S','glat_S','glon_S']], goodRows[['vx_S', 'vy_S', 'vz_S']], '_S')).T
    endPosSph_ffps = np.asarray(end_movement_spherical_noCartesian(goodRows[['rad_L','glat_L','glon_L']], goodRows[['vx_L', 'vy_L', 'vz_L']], '_L')).T

    def to_radians(rgg, suffix):
        return (rgg['rad%s' % suffix], rgg['glat%s' % suffix]*np.pi/180, rgg['glon%s' % suffix]*np.pi/180)

    sphCoords = np.asarray(to_radians(goodRows[['rad_S','glat_S','glon_S']], '_S')).T
    sphCoords_ffps = np.asarray(to_radians(goodRows[['rad_L','glat_L','glon_L']], '_L')).T

    totalEvents = 0
    resultData = {
        'transitData': [],
        'good lens obj_ids': [],
        'pop t_E values': [],        
        'my t_E values': []
        }

    for i in range(len(goodRows)):
        # starting coordiantes
        ffpCoord  = sphCoords_ffps[i, 1:3]
        ffpCoord_end  = endPosSph_ffps[i, 1:3]

        starCoord = sphCoords[i, 1:3]                
        starCoord_end = endPosSph[i, 1:3]

        d_ffp =  ffpCoord_end - ffpCoord
        d_star = starCoord_end - starCoord

        ffpCoord = ffpCoord - d_ffp/2
        starCoord = starCoord - d_star/2

        #d_lens, d_source, mass

        rr_ffp = RRectPath(
            einstein_radius(sphCoords_ffps[i, 0], sphCoords[i, 0], goodRows[i]['mass_L']),
            np.array(ffpCoord),
            np.array(d_ffp),
        1)

        rr_star = RRectPath(
            star_size(sphCoords[i][0]),
            np.array(starCoord),
            np.array(d_star),
        1)

        deltaTSq = rrQuadDiff(rr_ffp, rr_star)

        resultData['transitData'].append({
            'star obj_id': goodRows[i]['obj_id_S'],
            'ffp obj_id': goodRows[i]['obj_id_L'],
            'm_ffp': goodRows[i]['mass_L'],
            'm_star': goodRows[i]['mass_S'],
            'deltaTSq': nstr(deltaTSq, n=13)
        })

        resultData['starCoordsStart'] = sphCoords
        resultData['starCoordsEnd'] = endPosSph

        resultData['ffpCoordsStart'] = sphCoords_ffps
        resultData['ffpCoordsEnd'] = endPosSph_ffps


        if(mp_isnan(deltaTSq)):
            # no collision

            #print('nan -- radial distances L,S %s %s' % (sphCoords_ffps[i, 0],sphCoords[i, 0]))
            #print('nan -- einsten radius %s' % einstein_radius(sphCoords_ffps[i, 0], sphCoords[i, 0], goodRows[i]['mass_L']))
            #print('------------')
            continue

        vRel = (d_ffp - d_star) / duration
        speedRel = np.sqrt(vRel[0]**2 + vRel[1]**2)
        t_E = einstein_radius(sphCoords_ffps[i, 0], sphCoords[i, 0], goodRows[i]['mass_L']) / speedRel

        resultData['good lens obj_ids'].append(goodRows[i]['obj_id_L'])
        resultData['pop t_E values'].append(goodRows[i]['t_E'])
        resultData['my t_E values'].append(t_E)

        if (deltaTSq < transit15minSq):
            # too short
            continue
        
        (t1, t2) = rrQuadSolve(rr_ffp, rr_star)

        if((t1 >= 0 and t1 <= 1) or (t2 >= 0 and t2 <= 1)):
            print(1825*(t2-t1))
            totalEvents += 1

        
    print('total event count: %d' % totalEvents)

    jsonFilename = 'fits-data%s.json' % fn[:-5]
    print('writing output to file: %s' % jsonFilename)
    
    with open(jsonFilename, 'w') as f:
        json.dump(resultData, f, cls=NumpyEncoder)

#convertFits('30_mass_small.fits')
