import pandas as pd
import numpy as np
import RNA
import math
import rna_toolkit
from pseudo_decomposition import Peudo_Decom


def extract_features_pseudoknot_free(seq, bp):
    feature_dict = {}
    seq_temp = seq

    feature_dict['gc_perentage'] = float('nan')

    feature_dict['ent_3'] = float('nan')
    feature_dict['ent_4'] = float('nan')
    feature_dict['ent_5'] = float('nan')
    feature_dict['ent_6'] = float('nan')
    feature_dict['ent_7'] = float('nan')
    feature_dict['ent_8'] = float('nan')

    feature_dict['bfe_per']  = float('nan')
    feature_dict['kfe_per']  = float('nan')

    if rna_toolkit.entropy_max(len(seq_temp), 3) > 0:
        ent_3 = rna_toolkit.entropy(seq_temp, 3) / rna_toolkit.entropy_max(len(seq_temp), 3)

    if rna_toolkit.entropy_max(len(seq_temp), 4) > 0:
        ent_4 = rna_toolkit.entropy(seq_temp, 4) / rna_toolkit.entropy_max(len(seq_temp), 4)

    if rna_toolkit.entropy_max(len(seq_temp), 5) > 0:
        ent_5 = rna_toolkit.entropy(seq_temp, 5) / rna_toolkit.entropy_max(len(seq_temp), 5)

    if rna_toolkit.entropy_max(len(seq_temp), 6) > 0:
        ent_6 = rna_toolkit.entropy(seq_temp, 6) / rna_toolkit.entropy_max(len(seq_temp), 6)

    if rna_toolkit.entropy_max(len(seq_temp), 7) > 0:
        ent_7 = rna_toolkit.entropy(seq_temp, 7) / rna_toolkit.entropy_max(len(seq_temp), 7)

    if rna_toolkit.entropy_max(len(seq_temp), 8) > 0:
        ent_8 = rna_toolkit.entropy(seq_temp, 8) / rna_toolkit.entropy_max(len(seq_temp), 8)

    gc_perentage = (seq_temp.count('G') + seq_temp.count('C')) / float(len(seq_temp))
    


    a = RNA.fold_compound(seq_temp)
    a.pf()
    (s, mfe) = a.mfe()

    bfe,kfe = Peudo_Decom(bp,seq,'a')
    if mfe != 0:
        bfe_per = abs(bfe - mfe) / abs(mfe)
        feature_dict['bfe_per']  = bfe_per

    if bfe != 0:
        kfe_per = abs(bfe - kfe) / abs(bfe)
        feature_dict['kfe_per']  = kfe_per


    feature_dict['ent_3'] = ent_3
    feature_dict['ent_4'] = ent_4
    feature_dict['ent_5'] = ent_5
    feature_dict['ent_6'] = ent_6
    feature_dict['ent_7'] = ent_7
    feature_dict['ent_8'] = ent_8

    feature_dict['gc_perentage'] = gc_perentage
    
    df = pd.DataFrame(feature_dict)
    return df


def extract_features_pseudoknotted(seq, bp):
    feature_dict = {}
    seq_temp = seq

    feature_dict['gc_perentage'] = float('nan')

    feature_dict['ent_3'] = float('nan')
    feature_dict['ent_4'] = float('nan')
    feature_dict['ent_5'] = float('nan')
    feature_dict['ent_6'] = float('nan')
    feature_dict['ent_7'] = float('nan')
    feature_dict['ent_8'] = float('nan')

    feature_dict['bfe_per']  = float('nan')
    feature_dict['kfe_per']  = float('nan')

    if rna_toolkit.entropy_max(len(seq_temp), 3) > 0:
        ent_3 = rna_toolkit.entropy(seq_temp, 3) / rna_toolkit.entropy_max(len(seq_temp), 3)

    if rna_toolkit.entropy_max(len(seq_temp), 4) > 0:
        ent_4 = rna_toolkit.entropy(seq_temp, 4) / rna_toolkit.entropy_max(len(seq_temp), 4)

    if rna_toolkit.entropy_max(len(seq_temp), 5) > 0:
        ent_5 = rna_toolkit.entropy(seq_temp, 5) / rna_toolkit.entropy_max(len(seq_temp), 5)

    if rna_toolkit.entropy_max(len(seq_temp), 6) > 0:
        ent_6 = rna_toolkit.entropy(seq_temp, 6) / rna_toolkit.entropy_max(len(seq_temp), 6)

    if rna_toolkit.entropy_max(len(seq_temp), 7) > 0:
        ent_7 = rna_toolkit.entropy(seq_temp, 7) / rna_toolkit.entropy_max(len(seq_temp), 7)

    if rna_toolkit.entropy_max(len(seq_temp), 8) > 0:
        ent_8 = rna_toolkit.entropy(seq_temp, 8) / rna_toolkit.entropy_max(len(seq_temp), 8)

    gc_perentage = (seq_temp.count('G') + seq_temp.count('C')) / float(len(seq_temp))
    


    a = RNA.fold_compound(seq_temp)
    a.pf()
    (s, mfe) = a.mfe()

    bfe,kfe = Peudo_Decom(bp,seq,'a')
    if mfe != 0:
        bfe_per = abs(bfe - mfe) / abs(mfe)
        feature_dict['bfe_per']  = bfe_per

    if bfe != 0:
        kfe_per = abs(bfe - kfe) / abs(bfe)
        feature_dict['kfe_per']  = kfe_per


    feature_dict['ent_3'] = ent_3
    feature_dict['ent_4'] = ent_4
    feature_dict['ent_5'] = ent_5
    feature_dict['ent_6'] = ent_6
    feature_dict['ent_7'] = ent_7
    feature_dict['ent_8'] = ent_8

    feature_dict['gc_perentage'] = gc_perentage
    
    df = pd.DataFrame(feature_dict)
    return df
