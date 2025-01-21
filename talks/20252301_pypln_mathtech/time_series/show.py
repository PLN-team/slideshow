import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parameters
n_points = 1000  # Number of time points
mean = [2, 2]   # Means of X1 and X2
std_dev = [0.001, 0.001]  # Standard deviations of X1 and X2
correlation = 0.2  # Correlation coefficient
period = 0.02
generate_poiss = True

# Covariance matrix
cov_matrix = [[std_dev[0]**2, correlation * std_dev[0] * std_dev[1]],
              [correlation * std_dev[0] * std_dev[1], std_dev[1]**2]]
T = 500
# Generate Gaussian time series
time = np.linspace(0, T, n_points)  # Time values
X = np.random.multivariate_normal(mean, cov_matrix, size=n_points)

# Add non-linear trends
X1 = X[:, 0] + np.sin(period * time)  # Sinusoidal trend for X1
X2 = X[:, 1] + np.cos(period*time)            # Quadratic trend for X2
if generate_poiss is True:
    X1 = np.random.poisson(lam = np.exp(X1))
    X2 = np.random.poisson(lam = np.exp(X2))

# Plot in 3D
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')

# 3D plot with non-linear trends
ax.scatter(time, X1, X2, label="Gaussian Time Series with Non-Linear Trend", color='b', alpha=0.7)
ax.set_xlabel('Time (t)')
ax.set_ylabel('X1(t)')
ax.set_zlabel('X2(t)')
ax.set_title('3D Representation of Gaussian Time Series with Non-Linear Trend')

plt.legend()
plt.show()
