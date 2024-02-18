import matplotlib.pyplot as plt
import networkx as nx

# This function plots the histogram for a list of values
def plot_hist_values(values, xlabel, ylabel, title, nbins=20, size=(6, 4)):
    plt.figure(figsize=size)
    plt.hist(values, bins=nbins)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()

def color_graph_by_attribute(G, pos, attrib, title, weights, cmap=plt.cm.Oranges):
    # Map attrib values to colors using a colormap
    # Here, we use the 'viridis' colormap, but you can choose others like 'plasma', 'inferno', etc.
    colors = cmap(attrib)
    # Plot the graph
    plt.figure(figsize=(10, 10))  # Set the figure size|
    nx.draw_networkx(G, pos, node_color=attrib, cmap=cmap, with_labels=False, node_size=7, width=weights, alpha=1)
    # plt.colorbar(plt.cm.ScalarMappable(cmap=plt.cm.Oranges), label='Eigenvector Centrality', norm=norm)
    plt.title(title)
    plt.show()