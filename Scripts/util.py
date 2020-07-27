import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import RidgeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.model_selection import KFold, StratifiedKFold
from imblearn.over_sampling import SMOTE
from sklearn.svm import SVC
from collections import Counter
from sklearn.cluster import KMeans
from sklearn.metrics import accuracy_score
from scipy.spatial import distance_matrix
import seaborn as sns
import lightgbm as lgb
from sklearn.covariance import ShrunkCovariance, GraphicalLasso, LedoitWolf, OAS, MinCovDet
from sklearn.mixture import GaussianMixture
#from KB_Clustering import *
from sklearn.metrics import confusion_matrix
from sklearn.utils.multiclass import unique_labels

from sklearn.metrics import recall_score
import umap.umap_ as umap
import warnings
warnings.filterwarnings('ignore')


#%%
def KMeansClustering(train,k,q):
    train = train.values[:,:q]
    kmeans = KMeans(n_clusters=k, random_state=0).fit(train)
    return kmeans.labels_

def ClusteringAcc(label1, label2):
    acc1 = accuracy_score(label1, label2)
    for i in range(len(label1)):
        label1[i] -= 1
        if label1[i] == -1:
            label1[i] = 1
    acc2 = accuracy_score(label1, label2)
    if acc1>acc2:
        return acc1
    else:
        return acc2

def KBClustering(train, q):
    Data = train.values[:,:q]
    D = distance_matrix(Data,Data)
    
    B = get_cont_mat(D)
    LD = get_local_depth(B) 
    CC, E = get_clusters(B)   
    return CC, E, D
            
            

def kFoldRF(train, nfold):
    oof = np.zeros(len(train))
    X = train.drop('target', axis=1)
    y = train['target']

    kf = KFold(n_splits=nfold,shuffle=True, random_state=42)
    for trn_idx, val_idx in kf.split(X,y):
        clf = RandomForestClassifier()
        clf.fit(X.iloc[trn_idx], y.iloc[trn_idx])
        oof[val_idx] = clf.predict(X.iloc[val_idx])

    acc = accuracy_score(oof, train['target'])
    print('RF oof acc score: {}'.format(accuracy_score(oof, train['target'])))
    
    return oof, acc


def kFoldRidge(train, nfold):
    oof = np.zeros(len(train))
    X = train.drop('target', axis=1)
    y = train['target']

    kf = KFold(n_splits=nfold,shuffle=True,  random_state=42)
    for trn_idx, val_idx in kf.split(X,y):
        clf = RidgeClassifier()
        clf.fit(X.iloc[trn_idx], y.iloc[trn_idx])
        oof[val_idx] = clf.predict(X.iloc[val_idx])

    acc = accuracy_score(oof, train['target'])
    print('Ridge oof acc score: {}'.format(accuracy_score(oof, train['target'])))
    
    return oof, acc


def kFoldKNN(train, nfold):
    oof = np.zeros(len(train))
    X = train.drop('target', axis=1)
    y = train['target']

    kf = KFold(n_splits=nfold, shuffle=True,  random_state=42)
    for trn_idx, val_idx in kf.split(X,y):
        clf = KNeighborsClassifier()
        clf.fit(X.iloc[trn_idx], y.iloc[trn_idx])
        oof[val_idx] = clf.predict(X.iloc[val_idx])

    acc = accuracy_score(oof, train['target'])
    print('KNN oof acc score: {}'.format(accuracy_score(oof, train['target'])))
    
    return oof, acc

def get_cov_estimator(cov_type):
    if cov_type == 'LW':
        model = LedoitWolf()
    elif cov_type == 'OAS':
        model = OAS()
    elif cov_type == 'MCD':
        model = MinCovDet()
    elif cov_type[:2] == 'SC':
        shrinkage = float(cov_type.split('_')[1])
        model = ShrunkCovariance(shrinkage=shrinkage)
    elif cov_type[:2] == 'GL':
        alpha = float(cov_type.split('_')[1])
        model = GraphicalLasso(alpha=alpha)
    return model


def get_mean_cov(x, y, n_clusters_per_class=1, cov_type='LW'):
    model = get_cov_estimator(cov_type)
    ones = (y == 1).astype(bool)
    x2 = x[ones]
    model.fit(x2)
    p1 = model.precision_
    m1 = model.location_

    onesb = (y == 0).astype(bool)
    x2b = x[onesb]
    model.fit(x2b)
    p2 = model.precision_
    m2 = model.location_

    ms = np.stack([m1] * n_clusters_per_class + [m2] * n_clusters_per_class)
    ps = np.stack([p1] * n_clusters_per_class + [p2] * n_clusters_per_class)
    return ms, ps

def kFoldGMM(train, nfold):
    oof = np.zeros(len(train))
    skf = KFold(n_splits=n_splits,
                          random_state=42,
                          shuffle=True)
    for train_idx, valid_idx in skf.split(train, train['target']):
        ms, ps = get_mean_cov(train[train_idx],
                              train.loc[train_idx]['target'].values,
                              n_clusters_per_class=n_clusters_per_class,
                              cov_type=cov_type)

        gm = GaussianMixture(n_components=2 * n_clusters_per_class,
                             init_params=init_params,
                             covariance_type='full',
                             max_iter=100,
                             n_init=n_init,
                             precisions_init=ps,
                             random_state=42)
        gm.fit(train[train_idx])
        oof[valid_idx] += gm.predict_proba(train[valid_idx])[:, :2].sum(1)
    print('GMM oof acc score: {}'.format(accuracy_score(oof, train['target'])))
    return oof

def kFoldQDA(train, nfold):
    oof = np.zeros(len(train))
    X = train.drop('target', axis=1)
    y = train['target']

    kf = KFold(n_splits=nfold,shuffle=True, random_state=42)
    for trn_idx, val_idx in kf.split(X,y):
        clf = QuadraticDiscriminantAnalysis()
        clf.fit(X.iloc[trn_idx], y.iloc[trn_idx])
        oof[val_idx] = clf.predict(X.iloc[val_idx])

    acc = accuracy_score(oof, train['target'])
    print('QDA oof acc score: {}'.format(accuracy_score(oof, train['target'])))
    
    return oof, acc

def kFoldSVC(train, nfold):
    oof = np.zeros(len(train))
    X = train.drop('target', axis=1)
    y = train['target']

    kf = KFold(n_splits=nfold,shuffle=True, random_state=42)
    for trn_idx, val_idx in kf.split(X,y):
        clf = SVC()
        clf.fit(X.iloc[trn_idx], y.iloc[trn_idx])
        oof[val_idx] = clf.predict(X.iloc[val_idx])

    acc = accuracy_score(oof, train['target'])
    print('SVC oof acc score: {}'.format(accuracy_score(oof, train['target'])))
    
    return oof, acc



def display_importances(feature_importance_df_):
    feature_importance_df_ = feature_importance_df_.reset_index()
    cols = feature_importance_df_[["feature", "importance"]].groupby("feature").mean().sort_values(by="importance", ascending=False)[:40].index
    best_features = feature_importance_df_.loc[feature_importance_df_.feature.isin(cols)]

    plt.figure(figsize=(20, 10))
    sns.barplot(x="importance", y="feature", data=best_features, orient = 'h', order = best_features.feature)
    plt.title('LightGBM Features (avg over folds)')
    plt.tight_layout()
    plt.savefig('FeatureImportance.png')
    return best_features

def kfold_lightgbm(train_df, num_folds, feat=None, stratified = False, multi = False, numClass=6):

    # Cross validation model
    if stratified:
        folds = StratifiedKFold(n_splits= num_folds, shuffle=True, random_state=326)
    else:
        folds = KFold(n_splits= num_folds, shuffle=True, random_state=326)

    # Create arrays and dataframes to store results
    oof_preds = np.zeros(train_df.shape[0])
    feature_importance_df = pd.DataFrame()
    
    if feat is not None:
        feats = [f for f in feat if f not in ['target']]
    else:
        feats = [f for f in train_df.columns if f not in ['target']]

    # k-fold
    for n_fold, (train_idx, valid_idx) in enumerate(folds.split(train_df[feats], train_df['target'])):
        train_x, train_y = train_df[feats].iloc[train_idx], train_df['target'].iloc[train_idx]
        valid_x, valid_y = train_df[feats].iloc[valid_idx], train_df['target'].iloc[valid_idx]

        # set data structure
        lgb_train = lgb.Dataset(train_x,
                                label=train_y,
                                free_raw_data=False)
        lgb_test = lgb.Dataset(valid_x,
                               label=valid_y,
                               free_raw_data=False)

        # params optimized by optuna     
        if multi:
            params = {'num_leaves': 16,
             'objective':'multiclass',
             'num_class': numClass,
             'max_depth': -1,
             'learning_rate': 0.01,
             "boosting": "gbdt",
             "metric": 'multi_logloss',
             "verbosity": -1,
             "random_state": 2019}
        else:
            params = {'num_leaves': 32,
             'objective':'binary',
             'max_depth': -1,
             'learning_rate': 0.01,
             "boosting": "gbdt",
             "metric": 'binary_logloss',
             "verbosity": -1,
             "random_state": 2019}
            
        
        reg = lgb.train(
                        params,
                        lgb_train,
                        valid_sets=[lgb_train, lgb_test],
                        valid_names=['train', 'test'],
                        num_boost_round=10000,
                        early_stopping_rounds= 200,
                        verbose_eval=-1
                        )
        
        if not multi:
            oof_preds[valid_idx] = reg.predict(valid_x, num_iteration=reg.best_iteration)
        else:
            sub = reg.predict(valid_x, num_iteration=reg.best_iteration)
            tmp = np.zeros(len(sub))
            for i in range(len(tmp)):
                tmp[i] = np.argmax(sub[i])
            oof_preds[valid_idx] = tmp
            
            
        fold_importance_df = pd.DataFrame()
        fold_importance_df["feature"] = feats
        fold_importance_df["importance"] = np.log1p(reg.feature_importance(importance_type='gain', iteration=reg.best_iteration))
        fold_importance_df["fold"] = n_fold + 1
        
        feature_importance_df = pd.concat([feature_importance_df, fold_importance_df], axis=0)
        #print('Fold %2d accuracy : %.6f' % (n_fold + 1, (valid_y, oof_preds[valid_idx])))
        del reg, train_x, train_y, valid_x, valid_y

    # display importances

    feature_importance_df = feature_importance_df.groupby('feature').agg({'importance':['mean']})
    feature_importance_df.columns = ['importance']
    feature_importance_df = feature_importance_df.sort_values(by='importance', ascending=False)
    display_importances(feature_importance_df)
    
    if not multi:
        oof_preds[oof_preds<0.5] = 0
        oof_preds[oof_preds>=0.5] = 1
    acc = accuracy_score(oof_preds, train_df['target'])
    print('LGBM oof acc score: {}'.format(accuracy_score(oof_preds, train_df['target'])))
    return feature_importance_df, oof_preds, acc


def ScaleColumns(df, scafac):
    assert df.shape[1] == len(scafac)
    for i in range(len(scafac)):
        col = df.columns[i]
        df[col] = df[col]*scafac[i]
        
    return df


def TrainSeries(df, nfold, smote=False):
    CV_scores = []
    if smote:
        sm = SMOTE(random_state=42)
        X_resamp_tr, y_resamp_tr = sm.fit_resample(df.drop(['target'], axis=1), df['target'])
        X_resamp_tr = pd.DataFrame(X_resamp_tr)
        y_resamp_tr = pd.DataFrame({"target": y_resamp_tr})
        df = X_resamp_tr
        df['target'] = y_resamp_tr['target']
        
    print('Dummy classifier of target mode: {}'.format(accuracy_score([1]*len(df), df['target'])))
    CV_scores += [accuracy_score([1]*len(df), df['target'])]
    oof1, acc1 = kFoldRF(df, nfold)
    #oof2, acc2 = kFoldQDA(df, nfold)
    oof3, acc3 = kFoldKNN(df, nfold)
    oof4, acc4 = kFoldRidge(df, nfold)
    oof5, acc5 = kFoldSVC(df, nfold)
    #oof7, acc7 = kFoldGMM(df, 10)
    featImp, oof6, acc6 = kfold_lightgbm(df,num_folds=nfold, stratified=True, multi=True, numClass=61)
    CV_scores += [acc1, acc6, acc3, acc4, acc5]
    
    return CV_scores


def plot_confusion_matrix(y_true, y_pred, classes,
                          normalize=False,
                          title=None,
                          cm = None,
                          cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if cm is None:
        if not title:
            if normalize:
                title = 'Normalized confusion matrix'
            else:
                title = 'Confusion matrix, without normalization'
    
        # Compute confusion matrix
        cm = confusion_matrix(y_true, y_pred)
        # Only use the labels that appear in the data
        classes = classes[unique_labels(y_true, y_pred)]
        if normalize:
            cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
            print("Normalized confusion matrix")
        else:
            print('Confusion matrix, without normalization')
    
        #print(cm)

    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           xticklabels=classes, yticklabels=classes,
           title=title,
           ylabel='True label',
           xlabel='Predicted label')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()
    return ax,cm
    
