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
plt.plot(x_values, y_values, 'go-', label=rf'$f(x) = \frac{{{dna_strand_len} \cdot X}}{{1 - \left( \frac{{3}}{{4}} \right)^x}} - {dna_strand_len} \cdot x$')
# Labels and title
plt.xlabel("X")
plt.ylabel("f(x)")
plt.title("Expected Number of Overhead Cycles, n=200")
plt.legend()
plt.grid(True)
even_ticks = x_values[x_values % 4 == 0]  # Filter to get only even numbers
plt.xticks(even_ticks)

# Show the plot
plt.show()

#Plot the percentage graph
y_values_percentage = g(x_values, dna_strand_len)
plt.figure(figsize=(8, 5))
plt.plot(x_values, y_values_percentage, 'bo-', label=rf'$g(x) = \frac{{{dna_strand_len} \cdot x}}{{f(x)}}$')
# Labels and title
plt.xlabel("X")
plt.ylabel("g(x)")
plt.title("Percentage of Overhead Cycles with respect to Number of Bases")
plt.legend()
plt.grid(True)
even_ticks = x_values[x_values % 4 == 0]  # Filter to get only even numbers
plt.xticks(even_ticks)

# Show the plot
plt.show()

#multiple together

n_values = [5, 50, 100, 200]

# Plot f(x) for different n values
plt.figure(figsize=(8, 5))
for n in n_values:
    y_values = f(x_values, n)
    label = r'$f(x) \text n={}$'.format(n)
    plt.plot(x_values, y_values, label=f"n={n}")

# Labels and title
plt.xlabel("X")
plt.ylabel("f(x)")
plt.title("Expected Number of Overhead Cycles for Different n Values")
plt.grid(True)
plt.legend()


# Show the plot
plt.show()

