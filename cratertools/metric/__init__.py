from sklearn.neighbors import NearestNeighbors

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
        all_lo, all_la, all_rad = crater_B[index].T

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

def kn_match_craters(craters_A, craters_B, x='x', y='y', r='r',  max_neighbours=10, radius=3391):
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
    neighbours = max([len(craters_A), len(craters_B), max_neighbours])

    cA = craters_A[[x, y, r]].values
    cb = craters_B[[x, y, r]].values

    if len(craters_A) > max_neighbours and len(craters_B) > max_neighbours:
        nbrs = NearestNeighbors(n_neighbors=neighbours,
                                algorithm='ball_tree').fit(kn_B)
        data, indices = nbrs.kneighbors(kn_A)
    else:
        indices = [np.arange(len(cB))]*len(cA)
        data = None

    rc, un, dup, unique = match_craters(cA, cB, indices,
                                        _thresh_longlat2,
                                        _thresh_rad,
                                        radius=radius)
    # rc.columns = ["robbins_"+x, "robbins_"+y, "robbins_"+r, x, y, r, 'l']
    rc.columns = ["B_"+x, "B_"+y, "B_"+r, "A_"+x, "A_"+y, "A_"+r, 'l']
    # return the dataframe 'rc' sorted by index
    # an empty dict (this time), <- really?
    # the dictionary of duplicate craters
    # a list of the unique craters
    # The KNN data and indices arrays
    return rc.sort_index(), un, dup, unique, data, indices
