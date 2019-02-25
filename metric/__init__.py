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

