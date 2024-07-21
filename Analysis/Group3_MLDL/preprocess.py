import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from collections import Counter
from sklearn.preprocessing import MinMaxScaler
import seaborn as sns
import matplotlib.pyplot as plt
import os

# 在當前目錄新增資料夾: model, result, plot
if not os.path.exists("Group3_MLDL/model"):
    os.makedirs("Group3_MLDL/model")
if not os.path.exists("Group3_MLDL/result"):
    os.makedirs("Group3_MLDL/result")
if not os.path.exists("Group3_MLDL/plot"):
    os.makedirs("Group3_MLDL/plot")


def preprocess_data(genes_path, meta_path, resp_path, stratify=True):
    # input data
    data_genes = pd.read_csv(genes_path)
    data_meta = pd.read_csv(meta_path)
    data_resp = pd.read_csv(resp_path)

    # 移除duplicate
    data_genes = data_genes.drop_duplicates()
    data_meta = data_meta.drop_duplicates()
    data_resp = data_resp.drop_duplicates() 

    # 將三個資料集根據ID合併
    # data1 = pd.merge(data_genes, data_meta, how='inner', left_on='Unnamed: 0', right_on='ID')
    # data2 = pd.merge(data1, data_resp, how='inner', left_on='Unnamed: 0', right_on='ID')
    data = pd.concat([data_genes, data_meta, data_resp], axis=1)

    # 移除ID欄位
    # data = data2.drop(['Unnamed: 0', 'ID_x', 'ID_y'], axis=1)
    data = data.drop(['ID'], axis=1)

    # 看primary_disease的種類
    # 將出現次數少於10次的primary_disease改名為'Others'
    data['primary_disease'] = data['primary_disease'].apply(lambda x: 'Others' if data['primary_disease'].value_counts()[x] < 10 else x)
    data['primary_disease'].value_counts()

    # 將primary_disease轉換成one-hot encoding 建立新欄位
    cancers = data['primary_disease'].value_counts().index.tolist()
    for cancer in cancers:
        data[cancer] = data['primary_disease'].apply(lambda x: 1 if cancer in x else 0)
    data.drop('primary_disease', axis=1, inplace=True)

    # 將類別變數轉乘One-hot encoding
    conditions = {'Metastasis': 1, 'NonMetastasis': 0}
    data['bone'] = data['bone'].map(conditions)
    data['brain'] = data['brain'].map(conditions)
    data['kidney'] = data['kidney'].map(conditions)
    data['liver'] = data['liver'].map(conditions)
    data['lung'] = data['lung'].map(conditions)

    # 將y獨立出來, 有多個y
    data_y = data[['bone', 'brain', 'kidney', 'liver', 'lung']]
    data_X = data.drop(['bone', 'brain', 'kidney', 'liver', 'lung'], axis=1)

    # 若每行>0.5的數量小於2, 則移除
    data_X = data_X.loc[:, (data_X > 0.5).sum(axis=0) > 2]
    columns = data_X.columns

    # min-max normalization
    # scaler = MinMaxScaler()
    # data_X = scaler.fit_transform(data_X)

    # 將x, y 轉為numpy array
    data_X = data_X.to_numpy()
    data_y = data_y.to_numpy()

    if stratify:
        # 拆分訓練集和測試集(分層抽樣)
        X_train, X_test, y_train, y_test = split(data_X, data_y)
    else:
        # 隨機抽樣
        X_train, X_test, y_train, y_test = train_test_split(data_X, data_y, test_size=0.1, random_state=123)


    return X_train, X_test, y_train, y_test, columns


# 分層抽樣
def split(data_X, data_y):

    # 計算組合並給新組合編號
    group = Counter(tuple(row) for row in data_y)
    group_idx = {comb: idx for idx, comb in enumerate(group)}

    # 新組合編號
    data_y2 = np.array([group_idx[tuple(row)] for row in data_y])
    unique_classes = np.unique(data_y2) # 共25個組合
    
    # 只有一個樣本的類別直接給train
    single_sample_classes = [cls for cls in unique_classes if np.sum(data_y2 == cls) == 1]
    remove_idx = [idx for idx, cls in enumerate(data_y2) if cls in single_sample_classes] #只有一個樣本的類別的index

    # 將這個index的data_X data_y挑出來 給train
    X1_train = data_X[remove_idx, :]
    y1_train = data_y[remove_idx]

    # 移除這些index
    data_X = np.delete(data_X, remove_idx, axis=0)
    data_y = np.delete(data_y, remove_idx, axis=0)
    data_y2 = np.delete(data_y2, remove_idx)

    # 剩下的做分層抽樣
    X_train, X_test, y_train, y_test = train_test_split(data_X, data_y, test_size=0.11, stratify=data_y2, random_state=123)

    # 將剛剛挑出來的data併回train
    X_train = np.concatenate((X_train, X1_train), axis=0)
    y_train = np.concatenate((y_train, y1_train), axis=0)

    return X_train, X_test, y_train, y_test


# 繪製分層抽樣後的圖(也可繪製隨機抽樣的圖)
def stratify_plot(y_train, y_test):
    cancer = ['Bone', 'Brain', 'Kidney', 'Liver', 'Lung']
    plt.figure(figsize=(12, 8))
    for i in range(5):
        df_train = pd.DataFrame(y_train[:, i], columns=['Category'])
        df_train['Type'] = 'Train'
        df_test = pd.DataFrame(y_test[:, i], columns=['Category'])
        df_test['Type'] = 'Test'
        df = pd.concat([df_train, df_test])
        ax = plt.subplot(2, 3, i+1)
        sns.countplot(data=df, x='Type', hue='Category', ax=ax)
        plt.title(f'{cancer[i]}')

        for p in ax.patches:
            if p.get_height() > 0:
                ax.annotate(f'{int(p.get_height())}', (p.get_x() + p.get_width() / 2., p.get_height()),
                            ha='center', va='center', fontsize=10, color='black', xytext=(0, 5),
                            textcoords='offset points')
    plt.tight_layout()  # 調整每個子圖的間距以防重疊


# 繪製分層抽樣後的圖(水平堆疊的)
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def stratify_stack_plot(y_train, y_test):
    cancer = ['Bone', 'Brain', 'Kidney', 'Liver', 'Lung']
    plt.figure(figsize=(12, 8))
    for i in range(5):
        df_train = pd.DataFrame(y_train[:, i], columns=['Category'])
        df_train['Type'] = 'Train'
        df_test = pd.DataFrame(y_test[:, i], columns=['Category'])
        df_test['Type'] = 'Test'
        
        df = pd.concat([df_train, df_test])
        df['Type'] = pd.Categorical(df['Type'], categories=['Train', 'Test'], ordered=True)

        proportion_df = df.groupby(['Type', 'Category']).size().reset_index(name='Count')
        total_by_type = proportion_df.groupby('Type')['Count'].sum()
        proportion_df['Percentage'] = proportion_df.apply(lambda x: 100 * x['Count'] / total_by_type[x['Type']], axis=1)

        ax = plt.subplot(2, 3, i+1)
        categories = proportion_df['Category'].drop_duplicates().sort_values()
        bottom = 0 

        for category in categories:
            category_data = proportion_df[proportion_df['Category'] == category]
            bars = plt.bar(category_data['Type'], category_data['Percentage'], bottom=bottom, label=f'Category {category}')
            bottom += category_data['Percentage'].values  

            for bar in bars:
                y_val = bar.get_y() + bar.get_height() / 2  
                x_val = bar.get_x() + bar.get_width() / 2  
                plt.annotate(f'{bar.get_height():.1f}%', (x_val, y_val), ha='center', va='center', fontsize=10, color='white')

        ax.axhline(50, color='red', linewidth=1.5, linestyle='--')
        plt.title(f'{cancer[i]}')
        plt.ylabel('Percentage')
        plt.legend()

    plt.tight_layout()  
    plt.show()