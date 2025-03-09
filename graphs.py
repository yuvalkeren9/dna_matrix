import numpy as np
import matplotlib.pyplot as plt


###########################     PROBABILITY GRAPHS START ###############################

# Assume n = 200
dna_strand_len = 200

def plot_simple_graph(x_values, y_values, x_label, y_label, title):
    """
    Plot a graph

    :param x_values: The x values.
    :param y_values: The y values.
    :param x_label: The x label.
    :param y_label: The y label.
    :param title: The title to use for the graph.
    """
    # Plot the function expected number of extra cycles graph
    plt.figure(figsize=(8, 5))
    plt.plot(x_values, y_values, 'go-')
    # Labels and title
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.grid(True)

    # Show the plot
    plt.show()


def calc_exp_num_of_overhead_cycles(X, n):
    """
    Calculate the expected number of overhead cycles

    :param X: The number of strands.
    :param n: The length of a DNA strand.
    """
    return ((n * X) / (1 - (3/4)**X)) - n*X


def calc_perc_with_resp_to_bases_num(X, n):
    """
    Calculate the percentage of overhead cycles with respect to the number of bases.

    :param X: The number of strands.
    :param n: The length of a DNA strand.
    """
    return 1 - ((n*X) / (calc_exp_num_of_overhead_cycles(X,n) + n*X))

def plot_graph(func_to_use, num_of_strands, dna_strand_len):
    """
    Plot a graph with a given function and its parameters

    :param func_to_use: The function to apply.
    :param num_of_strands: The number of strands.
    :param dna_strand_len: The length of a DNA strand.
    """
    # Generate x values (natural numbers only)
    x_values = np.arange(1, 30)  # x from 1 to 20 (since x is a natural number)
    y_values = calc_exp_num_of_overhead_cycles(x_values, dna_strand_len)

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

if __name__ == '__main__':
    plot_graph(1,1, 200)
    #Plot the percentage graph
    x_values = np.arange(1, 30)  # x from 1 to 20 (since x is a natural number)
    y_values_percentage = calc_perc_with_resp_to_bases_num(x_values, dna_strand_len)
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
        y_values = calc_exp_num_of_overhead_cycles(x_values, n)
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

###########################     PROBABILITY GRAPHS END ###############################



