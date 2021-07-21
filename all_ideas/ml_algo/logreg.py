# y = wx+b
import numpy as np


class Logreg:
    def __init__(self, lr=0.001, n_iters=1000):
        self.lr = lr
        self.n_iters = n_iters
        self.w = None
        self.b = None

    def fit(self, X, y):
        n_samples, n_features = X.shape
        self.w = np.zeros(n_features)
        self.b = 0
        for _ in range(self.n_iters):
            # approximate y using linear regression model
            linear_model = np.dot(X, self.w) + self.b
            # apply sigmoid function
            predicted_y = self._sigmod(linear_model)

            # compute derivatives
            dw = (1/n_samples) * np.dot(X.T, (predicted_y-y))
            db = (1/n_samples) * sum(predicted_y-y)

            # change parameters
            self.w -= self.lr * dw
            self.b -= self.lr * db
        return

    def predict(self, X):
        linear_model = np.dot(X, self.w) + self.b
        predicted_y = self._sigmod(linear_model)
        predicted_class = [1 if x > 0.5 else 0 for x in predicted_y]
        return np.array(predicted_class)

    def _sigmod(self, x):
        return 1 / (1 + np.exp(-x))

