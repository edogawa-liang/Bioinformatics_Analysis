import torch.nn as nn
import torch
import numpy as np
import random
from sklearn.metrics import  roc_auc_score, accuracy_score, precision_score, recall_score, f1_score, confusion_matrix
import matplotlib.pyplot as plt


# 設定隨機種子
seed = 123
torch.manual_seed(seed)
np.random.seed(seed)
random.seed(seed)

class MultiHeadNN(nn.Module):
    def __init__(self, input_dim, shared_dim, output_dim):
        super(MultiHeadNN, self).__init__()
        # 共享層
        self.shared_layer = nn.Sequential(
            nn.Linear(input_dim, shared_dim),
            nn.ReLU(),
            # nn.Dropout(0.5),
            # nn.Linear(shared_dim, shared_dim),
            # nn.ReLU()
        )
        # 獨立頭
        self.output1 = nn.Linear(shared_dim, output_dim)
        self.output2 = nn.Linear(shared_dim, output_dim)
        self.output3 = nn.Linear(shared_dim, output_dim)
        self.output4 = nn.Linear(shared_dim, output_dim)
        self.output5 = nn.Linear(shared_dim, output_dim)

    def forward(self, x):
        shared = self.shared_layer(x)
        out1 = self.output1(shared)
        out2 = self.output2(shared)
        out3 = self.output3(shared)
        out4 = self.output4(shared)
        out5 = self.output5(shared)
        
        return out1, out2, out3, out4, out5
    

# 訓練模型
def train(model, X_train_tensor, y_train_tensor, criterion, optimizer, batch_size=4):
    model.train()
    permutation = torch.randperm(X_train_tensor.size()[0])
    
    total_loss = 0    
    for i in range(0, X_train_tensor.size()[0], batch_size):
        indices = permutation[i:i+batch_size]
        batch_X, batch_y = X_train_tensor[indices], y_train_tensor[indices]

        optimizer.zero_grad()
        outputs = model(batch_X)
        loss = sum(criterion(output, batch_y[:, idx].unsqueeze(1)) for idx, output in enumerate(outputs)) # 加總5個頭的損失
        loss.backward()
        optimizer.step()
        total_loss += loss.item()

    return total_loss / X_train_tensor.size()[0] # 每個樣本的平均損失


# 測試模型
def test(model, X_test_tensor):
    model.eval()
    with torch.no_grad():
        outputs = model(X_test_tensor)
        sigmoid_outputs = [torch.sigmoid(output).cpu().numpy() for output in outputs]  # Convert logits to probabilities
    
    return sigmoid_outputs


# 評估模型
def evaluate(outputs, labels, cancers, thr=0.5):
    epoch_metrics = {}
    for idx, output in enumerate(outputs):

        pred = output > thr
        label = labels[:, idx]
        
        auc = roc_auc_score(label, output)
        acc = accuracy_score(label, pred)
        pr = precision_score(label, pred, zero_division=0)
        re = recall_score(label, pred, zero_division=0)
        f1 = f1_score(label,pred, zero_division=0)
        cm = confusion_matrix(label, pred)
        
        epoch_metrics[f'{cancers[idx]}'] = {
            'AUC': auc,
            'ACC': acc,
            'Precision': pr,
            'Recall': re,
            'F1': f1,
            'Confusion Matrix': cm
        }         

    return epoch_metrics

# 繪製圖表
def evaluate_plot(performance_history, use_metric, cancers):
    evaluate_data = {cancer: [] for cancer in cancers}

    for epoch_data in performance_history:
        for head, metrics in epoch_data.items():
            evaluate_data[head].append(metrics[use_metric])

    # plt.figure(figsize=(6, 4))
    for head, values in evaluate_data.items():
        plt.plot(values, label=head)

    plt.title(f'Testing {use_metric}')
    plt.xlabel('Epoch')
    plt.ylabel(f'{use_metric}')
    # plt.xticks(range(len(performance_history)))
    plt.legend(loc="upper right")   
    # plt.show()