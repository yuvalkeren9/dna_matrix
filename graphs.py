import numpy as np
import matplotlib.pyplot as plt

# Assume n = 200
dna_strand_len = 200

def f(X, n):
    return ((n * X) / (1 - (3/4)**X)) - n*X

def g(X, n):
    return 1 - ((n*X) / (f(X,n) + n*X))

# Generate x values (natural numbers only)
x_values = np.arange(1, 30)  # x from 1 to 20 (since x is a natural number)
y_values = f(x_values, dna_strand_len)

# Plot the function expected number of extra cycles graph
plt.figure(figsize=(8, 5))
plt.plot(x_values, y_values, 'bo-', label=rf'$f(x) = \frac{{{dna_strand_len} \cdot x}}{{1 - (3/4)^x}}$')
# Labels and title
plt.xlabel("x (Natural numbers)")
plt.ylabel("f(x)")
plt.title("Plot of expected number of extra cycles")
plt.legend()
plt.grid(True)
even_ticks = x_values[x_values % 4 == 0]  # Filter to get only even numbers
plt.xticks(even_ticks)

# Show the plot
plt.show()

#Plot the percentage graph
y_values_percentage = g(x_values, dna_strand_len)
plt.figure(figsize=(8, 5))
plt.plot(x_values, y_values_percentage, 'bo-', label=rf'$f(x) = \frac{{{dna_strand_len} \cdot x}}{{1 - (3/4)^x}}$')
# Labels and title
plt.xlabel("x (Natural numbers)")
plt.ylabel("f(x)")
plt.title("Plot of the extra cycles with respect to the least cycles needed")
plt.legend()
plt.grid(True)
even_ticks = x_values[x_values % 4 == 0]  # Filter to get only even numbers
plt.xticks(even_ticks)

# Show the plot
plt.show()
