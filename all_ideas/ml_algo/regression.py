# y = wx+b
import numpy as np


class linear_r:
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
            y_predicted = np.dot(X, self.w) + self.b
            # Apply Gradient Decent
            dw = (1/n_samples) * np.dot(X.T, (y_predicted - y))
            db = (1/n_samples) * np.sum(y_predicted - y)

            self.w -= self.lr * dw
            self.b -= self.lr * db
        pass

    def predict(self, X):
        y_predicted = np.dot(X, self.w) + self.b
        return y_predicted