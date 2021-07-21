import matplotlib.pyplot as plt
from sklearn.datasets import make_regression as mr
from sklearn.model_selection import train_test_split as tts
from regression import linear_r
import numpy as np

X,y = mr(n_samples=100, n_features=1, noise=20, random_state=4)
X_train, X_test, y_train, y_test = tts(X,y, test_size=0.2, train_size=0.8)


def mse(y_true, y_pred):
    return np.mean((y_true-y_pred) **2)


reg = linear_r(lr=0.1)
reg.fit(X_train, y_train)
predicted = reg.predict(X_test)
mse_val = mse(y_test, predicted)
print(mse_val)
plt.scatter(X, y, edgecolor='k', marker='.')
plt.scatter(X_test, predicted, c='r', edgecolors='k')
plt.show()