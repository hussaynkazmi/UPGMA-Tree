#UPGMA TREE FROM DISTANCE MATRIX
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform

# === Define your distance matrix ===
labels = ["A", "B", "C", "D"]
dist_matrix = np.array([
    [0, 4, 6, 8],
    [4, 0, 6, 8],
    [6, 6, 0, 11],
    [8, 8, 11, 0]
])

# Convert to condensed form (since scipy needs upper triangular vector)
condensed = squareform(dist_matrix)

# Perform UPGMA clustering (average linkage)
Z = linkage(condensed, method='average')

# === Convert linkage distances to UPGMA node heights ===
# In UPGMA the node (cluster) age is half the distance at which clusters are
# joined, so convert the linkage distances (Z[:,2]) to heights = distance/2
Z_upgma = Z.copy()
Z_upgma[:, 2] = Z[:, 2] / 2.0

# === Plot dendrogram ===
plt.figure(figsize=(8,6))
dn = dendrogram(
    Z_upgma,
    labels=labels,
    color_threshold=0.0,
    above_threshold_color='black',
    leaf_font_size=12,
)

# Annotate individual branch lengths on the dendrogram
# For each merge in the linkage matrix, compute and show the branch lengths
# to its children (difference between parent and child heights)
try:
    icoord = dn['icoord']
    dcoord = dn['dcoord']
    
    # Get node heights from linkage matrix for computing branch lengths
    n_samples = len(labels)
    heights = {i: 0.0 for i in range(n_samples)}  # leaf heights = 0
    for i, row in enumerate(Z_upgma):
        node_id = n_samples + i
        heights[node_id] = row[2]  # height of this merged node
        
    # Annotate each branch segment
    for i, d in enumerate(dcoord):
        for j in range(4):
            if j == 1 or j == 2:  # skip middle segments
                continue
            if d[j] != d[j+1]:  # vertical branch segment
                x = icoord[i][j]
                y0, y1 = sorted([d[j], d[j+1]])
                ymid = (y0 + y1) / 2
                # Branch length is the height difference
                branch_len = y1 - y0
                if branch_len > 0:  # don't label zero-length branches
                    plt.text(x, ymid, f"{branch_len:.2f}", va='center',
                            ha='right' if j == 0 else 'left',
                            fontsize=8, backgroundcolor='white')

except Exception:
    # If annotation fails for any reason, continue without crashing the script
    pass
plt.title("UPGMA (Rooted) Phylogenetic Tree", fontsize=14, fontweight='bold')
plt.xlabel("Sequences")
plt.ylabel("Evolutionary Distance")
plt.grid(axis='y', linestyle='--', alpha=0.6)
# Save figure to file to avoid blocking (useful when running non-interactively)
plt.savefig('UPGMA_tree.png', bbox_inches='tight', dpi=200)
plt.show()

# Print the linkage matrices for verification
print("Original linkage matrix (Z) [merge indices, merge indices, distance, sample_count]:\n", Z)
print('\nUPGMA-adjusted linkage matrix (node heights = distance/2):\n', Z_upgma)


def linkage_to_newick(Z_linkage, labels_list):
    """Convert a linkage matrix (with node heights in Z_linkage[:,2]) to a Newick string.

    Assumes SciPy linkage ordering where newly created cluster id = n + k for merge k.
    """
    n = len(labels_list)
    # children mapping: node_id -> (left, right)
    children = {}
    heights = {i: 0.0 for i in range(n)}

    for k, row in enumerate(Z_linkage):
        left = int(row[0])
        right = int(row[1])
        node_id = n + k
        children[node_id] = (left, right)
        heights[node_id] = float(row[2])

    root = n + Z_linkage.shape[0] - 1

    def build(node, parent_height=None):
        if node < n:
            blen = parent_height - heights[node] if parent_height is not None else 0.0
            return f"{labels_list[node]}:{blen:.6f}"
        left, right = children[node]
        left_s = build(left, heights[node])
        right_s = build(right, heights[node])
        if parent_height is None:
            return f"({left_s},{right_s})"
        else:
            blen = parent_height - heights[node]
            return f"({left_s},{right_s}):{blen:.6f}"

    newick = build(root, None) + ";"
    return newick


try:
    newick = linkage_to_newick(Z_upgma, labels)
    print('\nUPGMA Newick string:\n', newick)
except Exception as e:
    print('Could not generate Newick:', e)
