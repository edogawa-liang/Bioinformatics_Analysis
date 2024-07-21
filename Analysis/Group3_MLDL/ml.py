from sklearn.metrics import  roc_auc_score, accuracy_score, precision_score, recall_score, f1_score, confusion_matrix
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
import lightgbm as lgb
from catboost import CatBoostClassifier
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle



# ML model
# 一次一個部位 一個模型
def ml_model(X_train, X_test, y_train, y_test, model, cancer, smote=False):
    '''
    X_train, X_test, y_train, y_test: numpy array
    model: str, 'DecisionTree', 'RandomForest', 'LightGBM', 'CatBoost'
    cancer: str, 'Bone', 'Brain', 'Kidney', 'Liver', 'Lung'
    smote: bool, default False
    '''
    if smote:
        i = 0
    else:
        i = ['Bone', 'Brain', 'Kidney', 'Liver', 'Lung'].index(cancer)
    
    # training
    if model == 'DecisionTree':
        clf = DecisionTreeClassifier(random_state=0)
        clf.fit(X_train, y_train[:, i])

    elif model == 'RandomForest':
        clf = RandomForestClassifier(random_state=0)
        clf.fit(X_train, y_train[:, i])
    
    elif model == 'LightGBM':
        clf = lgb.LGBMClassifier(random_state=0, metric='binary_logloss', max_depth=3, verbose=-1)
        clf.fit(X_train, y_train[:, i], eval_set=[(X_train, y_train[:, i])], eval_metric='binary_logloss')
    
    elif model == 'CatBoost':
        clf = CatBoostClassifier(random_state=0, iterations=150, verbose=0)
        clf.fit(X_train, y_train[:, i], plot=False, eval_set=[(X_train, y_train[:, i])])

    # 儲存模型
    if smote:
        with open(f'Group3_MLDL/model/{cancer}_{model}_smote.pkl', 'wb') as file:
            pickle.dump(clf, file)
    else:
        with open(f'Group3_MLDL/model/{cancer}_{model}.pkl', 'wb') as file:
            pickle.dump(clf, file)

    
    # 預測
    y_pred = clf.predict(X_test)
    y_prob = clf.predict_proba(X_test)[:, 1]

    acc = accuracy_score(y_test[:, i], y_pred)
    precision = precision_score(y_test[:, i], y_pred, zero_division=0)
    recall = recall_score(y_test[:, i], y_pred)
    f1 = f1_score(y_test[:, i], y_pred)
    auc = roc_auc_score(y_test[:, i], y_prob)

    return acc, precision, recall, f1, auc


# 多個部位 一個模型
# 一次跑完所有轉移部位
def check_transfer(X_train, X_test, y_train, y_test, model):
    """
    X_train, X_test, y_train, y_test: numpy array
    model: str, 'DecisionTree', 'RandomForest', 'LightGBM', 'CatBoost'
    plot: bool, default False (適用於LightGBM, CatBoost畫出training loss)
    """
    metrics = {}
    cancers = ['Bone', 'Brain', 'Kidney', 'Liver', 'Lung']

    for i in range(y_train.shape[1]):
        acc, precision, recall, f1, auc = ml_model(X_train, X_test, y_train, y_test, model, cancer=cancers[i])
        metrics[f'{cancers[i]}'] = {
        'Accuracy': acc,
        'Precision': precision,
        'Recall': recall,
        'F1': f1,
        'AUC': auc
        }
    
    df_transfer = pd.DataFrame.from_dict(metrics, orient='index').round(4) 
    df_transfer.to_csv(f'Group3_MLDL/result/{model}_transfer.csv')
    
    return df_transfer


# smote 後再重訓練模型
# 一次一個部位, 多個模型
def afterSMOTE_models(models, cancer_type, X_train_smote, X_test, y_train_smote, y_test):
    '''
    models: list, ['DecisionTree', 'RandomForest', 'LightGBM', 'CatBoost']
    cancer_type: str, 'Bone', 'Brain', 'Kidney', 'Liver', 'Lung'
    X_train_smote, X_test, y_train_smote, y_test: numpy array
    '''
    metrics = ['Accuracy', 'Precision', 'Recall', 'F1', 'AUC']

    print(123)
    smo = {'Model': []}  # 初始化包含 'Model' 鍵的字典
    for metric in metrics:
        smo[metric] = []  # 為每個度量指標添加一個空列表

    for model in models:
        acc, precision, recall, f1, auc = ml_model(X_train_smote, X_test, y_train_smote.reshape(-1, 1), y_test, model, cancer=cancer_type, smote=True) # SMOTE只針對一個部位: i=0
        smo['Model'].append(model)
        smo['Accuracy'].append(acc)
        smo['Precision'].append(precision)
        smo['Recall'].append(recall)
        smo['F1'].append(f1)
        smo['AUC'].append(auc)

    df = pd.DataFrame(smo).round(4) 
    df.set_index('Model', inplace=True)

    print(f"{cancer_type} model after SMOTE:")
    df.to_csv(f'Group3_MLDL/result/{cancer_type}_smote.csv')
    
    return df


# feature importance
def featureplot(importance, k, X_colname, title):
    '''
    importance: 
    k: number of features
    X_colname: colname of X (list)
    '''

    indices_of_largest = np.argsort(np.abs(importance), )[-k:][::-1]
    k_importance = [X_colname[i] for i in indices_of_largest]
    
    plt.figure(figsize=(6, 6))
    plt.barh(k_importance, np.abs(importance)[indices_of_largest])
    plt.xlabel('Importance', size=12)
    plt.title(title, size=13)
    plt.yticks(range(k), labels=k_importance)
    plt.savefig(f'Group3_MLDL/plot/feat_{title}.png', dpi=300, bbox_inches='tight')
    
    # return k_importance, np.abs(importance)[indices_of_largest]


# 選擇AUC最高的模型
def choose_model(df_dt, df_rf, df_lgbm, df_catb, df):
    df_dt['Source'] = 'DecisionTree'
    df_rf['Source'] = 'RandomForest'
    df_lgbm['Source'] = 'LightGBM'
    df_catb['Source'] = 'CatBoost'
    combined_df = pd.concat([df_dt.T, df_rf.T, df_lgbm.T, df_catb.T], axis=1)
    df['Source'] = ['DecisionTree_Smote', 'RandomForest_Smote', 'LightGBM_Smote', 'CatBoost_Smote']
    df = df.T
    df.columns = ["Brain"]*4
    final_df = pd.concat([combined_df, df], axis=1)
    cancers = ['Bone', 'Brain', 'Kidney', 'Liver', 'Lung']
    results = []

    for i in range(5):
        d = final_df.loc[:,final_df.columns==cancers[i]]
        dd = pd.to_numeric(d.loc['AUC'])
        result = d.iloc[:, dd.tolist().index(max(dd.tolist()))]
        results.append(result)

    df = pd.concat(results, axis=1)
    df.to_csv('Group3_MLDL/result/choose_model.csv')
    return df