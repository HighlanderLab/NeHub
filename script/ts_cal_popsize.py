import numpy as np

# Choose time intervals (e.g., [0, 100], [100, 500], [500, 1000] generations)
time_intervals = [(0, 100), (100, 500), (500, 1000)]

# Dictionary to store the inferred population sizes for each interval
pop_size_by_interval = {interval: [] for interval in time_intervals}

# Iterate through the tree sequence
for tree in ts.trees():
    for interval in time_intervals:
        start, end = interval
        coalescence_times = []

        # Get coalescence times within the specified interval
        for node in tree.nodes():
            coalescence_time = tree.time(node)
            if start <= coalescence_time < end:  # Check if within the time interval
                coalescence_times.append(coalescence_time)

        # Calculate the population size for this interval
        if coalescence_times:
            mean_coalescence_time = np.mean(coalescence_times)
            Ne = mean_coalescence_time / 2
            pop_size_by_interval[interval].append(Ne)

# Average the population sizes for each interval
for interval, pop_sizes in pop_size_by_interval.items():
    if pop_sizes:
        avg_pop_size = np.mean(pop_sizes)
        print(f"Inferred population size between {interval[0]} and {interval[1]} generations: {avg_pop_size}")
