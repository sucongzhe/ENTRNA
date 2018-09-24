import numpy as np
import pickle
import pandas as pd
from sklearn import preprocessing
from scipy.spatial import distance
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.model_selection import KFold
import itertools
from sklearn.metrics import accuracy_score
from scipy.spatial import distance
import extract_features
from sklearn.metrics import confusion_matrix
from itertools import combinations



def bulid_model(real_RNA_loc,folder_simulation_result):
    try:
        feature_set = ['ent_3','ent_8','bfe_per']
        df_real = pd.read_csv(real_RNA_loc, header = 0)
        df_v1 = df_real[df_real['bfe_per'] < 2]
        file_name_list = df_v1['file_name'].tolist()
        df_v1['length'] = df_v1['seq'].apply(lambda x:len(x))
        result_training_accuracy_list = []
        result_test_accuracy_list = []

        mat_all = df_v1[feature_set].as_matrix()

        train = list(range(len(file_name_list)))
        mat_train_real = mat_all[train]
        file_name_train_list = [file_name_list[i] for i in train]
        x_negative_list = []
        for ii in range(len(file_name_train_list)):
            x_true = mat_train_real[ii]
            df_temp = pd.read_csv(folder_simulation_result+file_name_train_list[ii].replace('.ct','.csv'), header = 0)
            mat_temp_v1 = df_temp[feature_set].as_matrix()

            distance_list = []
            for row_i in mat_temp_v1:
                if np.isnan(row_i).any() == True:
                    distance_list.append(0)
                else:
                    dst = distance.euclidean(row_i,x_true)
                    distance_list.append(dst)
            x_negative_list.append(mat_temp_v1[np.argmax(distance_list)])
        x_train = np.concatenate((mat_train_real,x_negative_list), axis=0)
        y_train = np.concatenate((np.ones(len(train)),np.zeros(len(train))), axis=0)
        clf = LogisticRegression()        
        clf.fit(x_train,y_train)
        y_pred = clf.predict(x_train)
        train_acc = accuracy_score(y_train,y_pred)
        tn, fp, fn, tp = confusion_matrix(y_train, y_pred).ravel()
        return clf
    except ValueError:
        print 'failed!'


    
def pred_foldability(df_pred,clf):
    mat_x = df_pred[['ent_3','ent_8','bfe_per']].as_matrix()
    for row_i in mat_x:
        if np.isnan(row_i).any() == True:
            return float('nan')
            
        else:
            return clf.predict_proba([row_i,])[0][1]



def entrna_main(seq,sec_str,real_RNA_loc = "RNASTRAND_pseudoknotted_feature.csv", folder_simulation_result = "~/RNASTRAND_extract_feature/"):
    df_pred = extract_features.extract_features(seq,sec_str)
    clf = bulid_model(real_RNA_loc,folder_simulation_result)
    foldability = pred_foldability(df_pred,clf)
    return foldability

