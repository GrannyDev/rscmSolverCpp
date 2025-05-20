import matplotlib
# Use the PGF backend for LaTeX compatibility
matplotlib.use("pgf")
import matplotlib.pyplot as plt

# Configure PGF parameters for LaTeX
# Note: pgf.preamble expects a single string, not a list
plt.rcParams.update({
    "text.usetex": True,
    "pgf.rcfonts": False,                # disable Matplotlib’s auto‐fonting
    "font.size": 10,                   # font size
})
# List to store cost values
costs = []

# Read the log file and extract the cost from each line
with open('log.txt', 'r') as file:
    for line in file:
        parts = line.strip().split(':')
        if len(parts) == 2:
            try:
                cost = int(parts[1].strip())
                costs.append(cost)
            except ValueError:
                continue

# Print statistics if available
if costs:
    print(f"Lowest cost: {min(costs)}")
    print(f"Highest cost: {max(costs)}")
    print(f"Average cost: {sum(costs) / len(costs):.2f}")
print(f"Number of costs: {len(costs)}")

# Check for data before plotting
if not costs:
    print("No valid cost data found.")
else:
    # Create a histogram plot of the costs with reduced height
    plt.figure(figsize=(8, 3))  # Adjust the figure size to reduce height
    plt.hist(costs, bins=15, edgecolor='#5e81ac', color='#b48ead')
    plt.xlabel("Area Estimation")
    plt.ylabel("Frequency")
    # Save the figure in PGF format for inclusion in LaTeX
    plt.savefig("histogram.pgf")
    # print("Plot saved as 'histogram.pgf'.")
    # Show the plot
    plt.show()
