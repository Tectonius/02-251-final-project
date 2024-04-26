import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


# with open('scored_pairs.pkl', 'rb') as f:
#     scored_pairs = pickle.load(f)


# Example input data (replace this with your actual data)
scored_pairs = [
    (90, 85),
    (92, 82),
    (88, 80),
    (95, 87),
    (87, 81)
]

# Extracting scores for optimized and default hyperparameters
optimized_scores = [pair[0] for pair in scored_pairs]
default_scores = [pair[1] for pair in scored_pairs]

# Plotting a box plot to compare distributions
plt.figure(figsize=(8, 6))
sns.boxplot(data=[optimized_scores, default_scores], palette="Set3")
plt.xticks([0, 1], ['Optimized', 'Default'])
plt.title("Comparison of ClustalW Performance")
plt.ylabel("Alignment Score")
plt.show()

# Calculating summary statistics
optimized_mean = np.mean(optimized_scores)
default_mean = np.mean(default_scores)
optimized_std = np.std(optimized_scores)
default_std = np.std(default_scores)

print("Summary Statistics:")
print("Optimized Hyperparameters - Mean:", optimized_mean, " Std Dev:", optimized_std)
print("Default Hyperparameters - Mean:", default_mean, " Std Dev:", default_std)

# Performing statistical test (e.g., paired t-test) to check for significant difference
from scipy.stats import ttest_rel

t_statistic, p_value = ttest_rel(optimized_scores, default_scores)
print("\nPaired t-test results:")
print("t-statistic:", t_statistic)
print("p-value:", p_value)
if p_value < 0.05:
    print("The difference in scores is statistically significant.")
else:
    print("The difference in scores is not statistically significant.")
