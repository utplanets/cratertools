from sklearn.neighbors import NearestNeighbors
import numpy as np
import pandas as pd

import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

_thresh_longlat2 = 0.25
_thresh_rad = 0.25


def metrics(table):
    """Calculate the metrics for a set of crater matches.

    Arguments: 
        table : Pandas table with the columns
                N_match
                N_detect
                N_csv
    Returns:
        table : modified table with the additional columns
                frac_new_csv : fraction of unmatched craters relative to the matched + diff
                frac_new_detect : fraction of detected creaters that are new (FPR)
                precision : The ratio of true positives to all positives
                recall : The ratio of true positives to all craters
                f1 : The harmonic mean of precision and recall.
    """

    tp = table["N_match"]
    fp = table["N_detect"]-table["N_match"]
    fn = table["N_csv"] - table["N_match"]
    g = table["N_csv"]
    diff = table["N_detect"] - table["N_match"]
    table["frac_new_csv"] = diff/(table["N_csv"] +diff)
    table["frac_new_detect"] = diff/(table["N_detect"])
    p = tp/(tp+fp)
    r = tp/(fn+tp)
    table["precision"] = p
    table["recall"] = r
    table["f1"] = 2*(r*p)/(r + p)
    return table

def errmet(data, arad=3391):
    ckd = 180. / (np.pi*arad)
    xp, xg = data["A_Long"], data["B_Long"]
    yp, yg = data["A_Lat"], data["B_Lat"]

    rp, rg = data["A_Diameter (km)"], data["B_Diameter (km)"]
    lm = 0.5*(yp+yg)
    dlon = np.abs(xp-xg)*np.cos(np.deg2rad(lm))/(rg*ckd)
    dlat = np.abs(yp-yg)/(rg*ckd)
    dr = np.abs(rp-rg)/rg
    return np.percentile(dlon,[25,50,75]),np.percentile(dlat,[25,50,75]),np.percentile(dr,[25,50,75])

def global_metric(matched_data, matched_count, craters_B, arad=3391):
    q=(np.array(matched_count) == 0)

    tp = np.sum(np.array(matched_count) > 0)
    fp = np.sum(np.array(matched_count) == 0)
    fn = len(craters_B) - tp

    mm = pd.DataFrame(metrics(dict(N_match=tp,
                                   N_detect=tp+fp,
                                   N_csv=tp+fn)),
                      index=['data'])
    mm['precision']*=100
    mm['recall']*=100

    met = errmet(matched_data, arad=arad)
    m25, med, m75 = [np.array(a) for a in zip(*met)]

    mm["err_lo"], mm["err_la"], mm["err_r"] = med

    m25 = pd.DataFrame(m25,
                       index=["err_lo","err_la","err_r"],
                       columns=["data 25"])
    m75 = pd.DataFrame(m75,
                       index=["err_lo","err_la","err_r"],
                       columns=["data 75"])

    return pd.concat([mm.T,m25,m75],axis=1)

def filter_craters(craters, indices, matches=1,
                   radius=3391, 
                   thresh_rad=_thresh_rad,
                   thresh_longlat2=_thresh_longlat2):

    k2d = 180. / (np.pi * radius)       # km to deg                                                                             
    duplicate_craters = dict()
    reverse_craters = dict()
    from tqdm import tqdm
    used = set()
    for iloc, index in enumerate(indices[:]):
        loc, remainder = iloc, index
        all_lo, all_la, all_rad = craters[index].T
        lo, la, rad = craters[iloc].T
        Long, Lat, Rad = all_lo[:], all_la[:], all_rad[:]
        
        minr = np.where(rad < Rad, rad, Rad)
        
        la_m = 0.5 * (la + all_la)
        dL = (((Long - lo) * np.cos(np.deg2rad(la_m)) / (0.5 * minr * k2d))**2
              + ((Lat - la) / (0.5 * minr * k2d))**2)
        dR = np.abs(Rad - rad) / minr
        w = np.where((dR < thresh_rad) & (dL < thresh_longlat2))
        if len(w[0]) >= matches:
            #proposed_used = index[w]
            #new_used = set(proposed_used) - used
            #used+=new_used
            save = []
            for i in index[w]:
                if i not in reverse_craters:
                    save.append(i)
                    reverse_craters[i] = loc
            if len(save):
                duplicate_craters[loc] = save
    rc = np.zeros(len(duplicate_craters)*3).reshape((len(duplicate_craters),3))

    for i,k in enumerate(duplicate_craters.keys()):
        rc[i] = craters[duplicate_craters[k]][0]#.mean(axis=0)

    rc = pd.DataFrame(rc,
                      columns=['x','y','r'])
    
    return rc, reverse_craters, duplicate_craters

def rep_filter_unique_craters(incraters,
                              x='x', y='y', r='r', matches=1,
                              max_neighbours=8,max_iter=8,
                              radius=3391, 
                              thresh_rad=_thresh_rad,
                              thresh_longlat2=_thresh_longlat2):

    craters = incraters[~incraters.duplicated(subset=[x,y,r])]
    iteration = 0
    rcs = 0

    while iteration < max_iter:
        kn_craters = np.array([np.deg2rad(craters[x])*radius*np.cos(np.deg2rad(craters[y])),
                  np.deg2rad(craters[y])*radius,
                  360*craters[r]/(radius*2*np.pi*np.cos(np.deg2rad(craters[y])))
                 ]).T

        neighbours = max_neighbours if max_neighbours < len(craters) else len(craters)

        nbrs = NearestNeighbors(n_neighbors=neighbours, algorithm='ball_tree').fit(kn_craters)

        data,indices = nbrs.kneighbors(kn_craters)

        cc=craters[[x,y,r]].values

        craters, un, dup = filter_craters(cc,
                                          indices,
                                          matches=matches,
                                          radius=radius,
                                          thresh_rad=thresh_rad,
                                          thresh_longlat2=thresh_longlat2)
        craters.columns = [x,y,r]

        if len(craters) == rcs or iteration == max_iter:
            break
        rcs = len(craters)
        iteration += 1
  
    return craters, un, dup, data, indices

def match_craters(craters_A, craters_B, indices, thresh_longlat2, thresh_rad, radius=3391):
    """Match craters between two datasets in three dimensions.

    Iterate through each pair and compare against a threshold for detection.
    Arguments:
        craters_A : pandas table of crater locations
        craters_B : pandas table of crater locations
        thresh_longlat2 : spatial threshold in terms of fractions of radius
        thresh_rad : radius threshold in terms of fractions of radius
    Keywords:
        radius : planetary radius in km
    Returns:
        So many things
    """

    k2d = 180. / (np.pi * radius)       # km to deg                         

    # set up the output
    duplicate_craters = dict()
    reverse_craters = dict()
    unique = []

    for loc, index in enumerate(indices[:]):
        # extract the target crater from the A list, and possible matches from
        # the B list according to the indices
        lo, la, rad = craters_A[loc, :]
        all_lo, all_la, all_rad = craters_B[index].T

        Long, Lat, Rad = all_lo[:], all_la[:], all_rad[:]

        # The minimum radius in out comparison for each pair
        minr = np.where(rad < Rad, rad, Rad)
        la_m = np.mean(all_la) 

        # The difference in (longitude/radius)**2+(long/radius)**2
        # as a fraction of the minimum radius
        dL = (((Long - lo) * np.cos(np.pi * la_m / 180.) / (0.5*minr * k2d))**2
              + ((Lat - la) / (0.5*minr * k2d))**2)
        # The difference in radius/diameter relative to the minimum radius
        dR = np.abs(Rad - rad) / minr
        # Locations where the matches are possible
        w = np.where((dR < thresh_rad) & (dL < thresh_longlat2))
        #If there is a match for this crater, store the location in the
        # 'duplicate_craters' dictionary (to remember the location)
        if len(w[0]) > 0:
            save = index[w]
            duplicate_craters[loc] = save
        # Append to the 'unique' list the number of matches that were possible
        # (>1 is a sign of duplicates or clusters)
        unique.append(len(w[0]))
    # Loop over
    # Create a new table with 7 columns, enough rows for the matching craters
    ld = len(duplicate_craters)
    rc = np.zeros(ld*7).reshape((ld, 7))
    idx = np.zeros(ld, dtype=int)
    # Iterate through the duplicate_craters, for each match
    # copy the location data from craters_B and craters_A, and the number
    # of matches into the array created above
    for i, k in enumerate(duplicate_craters.keys()):
        for j in duplicate_craters[k]:
            rc[i] = np.hstack([craters_B[duplicate_craters[k]][0],
                              craters_A[k],
                              len(duplicate_craters[k])])
        idx[i] = k  # store the key in the 'index' array

    # Create a dataframe out of the array, using the crater indices
    # as the index column
    rc = pd.DataFrame(rc, columns=['Bx', 'By', 'Br',
                                   'Ax', 'Ay', 'Ar',
                                   'l'], index=idx)
    # return the dataframe,
    # an empty dict (this time),
    # the dictionary of duplicate craters
    # a list of the unique craters 
    return rc, reverse_craters, duplicate_craters, unique

def kn_match_craters(craters_A, craters_B, x='x', y='y', r='r',  max_neighbours=10, radius=3391, 
                     thresh_rad=_thresh_rad, thresh_longlat2=_thresh_longlat2):
    """Match craters between two data sets in three dimensions.

    Uses a K Nearest Neighbours algorithm to cluster small groups, then matches the closest
    crater in lat,lon,diameter dimensions.

    Arguments:
        craters_A : pandas table of crater locations
        craters_B : pandas table of crater locations
    Keywords:
        x : longitude dimension name
        y : latitude dimension name
        r : radius dimension name
        max_neightbours: number of neighbours to consider in the KNN algorithm
        radius : planetary radius in km
    Returns:
        So many things


    """
    logger = logging.getLogger()
    logger.info("kn_match_craters: max_neighbours={}, radius={}, thresh_rad={}, thresh_longlat={}".format(max_neighbours, radius, thresh_rad, thresh_longlat2))
    ca = np.cos(np.deg2rad(craters_A[y]))
    kn_A = np.array([np.deg2rad(craters_A[x])*radius*ca,
                     np.deg2rad(craters_A[y])*radius,
                     360*craters_A[r]/(radius*2*np.pi*ca)
                     ]).T

    cb = np.cos(np.deg2rad(craters_B[y]))
    kn_B = np.array([np.deg2rad(craters_B[x])*radius*cb,
                     np.deg2rad(craters_B[y])*radius,
                     360*craters_B[r]/(radius*2*np.pi*cb)
                     ]).T

    # Sanity check for number of neighbours
    neighbours = min([len(craters_A), len(craters_B), max_neighbours])

    cA = craters_A[[x, y, r]].values
    cB = craters_B[[x, y, r]].values

    if len(craters_A) > max_neighbours and len(craters_B) > max_neighbours:
        nbrs = NearestNeighbors(n_neighbors=neighbours,
                                algorithm='ball_tree').fit(kn_B)
        data, indices = nbrs.kneighbors(kn_A)
    else:
        indices = [np.arange(len(cB))]*len(cA)
        data = None
    rc, un, dup, unique = match_craters(cA, cB, indices,
                                        thresh_longlat2,
                                        thresh_rad,
                                        radius=radius)
    # rc.columns = ["robbins_"+x, "robbins_"+y, "robbins_"+r, x, y, r, 'l']
    rc.columns = ["B_"+x, "B_"+y, "B_"+r, "A_"+x, "A_"+y, "A_"+r, 'l']
    # return the dataframe 'rc' sorted by index
    # an empty dict (this time), <- really?
    # the dictionary of duplicate craters
    # a list of the unique craters
    # The KNN data and indices arrays
    return rc.sort_index(), un, dup, unique, data, indices
