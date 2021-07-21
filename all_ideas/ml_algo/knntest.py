import os

from sklearn import datasets
from sklearn.model_selection import train_test_split
import numpy as np
import matplotlib.pyplot as plt
# from os import chdir
# chdir('/Users/anmol_gorakshakar/python/machine_learning/')

from knn import knn

iris = datasets.load_iris()
X, y = iris.data, iris.target
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, train_size=0.8)

clf = knn(k=2)
clf.fit(X_train, y_train)
predict = clf.predict(X_test)
acc = np.sum(predict == y_test)/len(y_test)
print(acc)

plt.scatter(predict, y_test, edgecolor='k', c=y_test)
plt.show()