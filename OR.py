from scipy.stats import fisher_exact
import numpy as np

# Define the contingency table (2x2)
contingency_table = np.array([[20, 30], [10, 40]])

# Perform Fisher's exact test
odds_ratio, p_value = fisher_exact(contingency_table)

print("Fisher's exact test:")
print("Odds Ratio:", odds_ratio)
print("p-value:", p_value)

# Calculate the odds ratio manually
a, b, c, d = contingency_table.flatten()
odds_ratio_manual = (a * d) / (b * c)

print("\nOdds Ratio (Manual Calculation):", odds_ratio_manual)

# Calculate the confidence intervals using the Mantel-Haenszel method
from statsmodels.stats.contingency_tables import Table2x2

table = Table2x2(contingency_table)
conf_int = table.oddsratio_confint(method='normal')

print("\nConfidence Intervals:")
print("Lower Bound:", conf_int[0])
print("Upper Bound:", conf_int[1])
#%%
