# 🌳 UPGMA Tree Construction using Python

A Python implementation of the **UPGMA (Unweighted Pair Group Method with Arithmetic Mean)** algorithm for constructing a **rooted phylogenetic tree** from a distance matrix.

The script takes a pairwise distance matrix, performs hierarchical clustering using the UPGMA algorithm, visualizes the resulting phylogenetic tree, annotates branch lengths, and exports the tree in **Newick format**.

---

## 📌 Features

- Construct rooted phylogenetic trees using the UPGMA algorithm
- Accepts a user-defined distance matrix
- Generates publication-quality dendrograms
- Displays branch lengths
- Exports the tree in Newick format
- Saves the tree as a high-resolution PNG image
- Prints linkage matrices for verification

---

## 📂 Repository Structure

```
UPGMA-Tree/
│
├── UPGMA TREE.py          # Main Python script
├── UPGMA_tree.png         # Example output tree
├── README.md
├── requirements.txt
├── LICENSE
└── .gitignore
```

---

## 📦 Requirements

Python 3.9+

Install the required packages using:

```bash
pip install -r requirements.txt
```

or

```bash
pip install numpy scipy matplotlib
```

---

## 🚀 Usage

Run the script:

```bash
python "UPGMA TREE.py"
```

The script will:

- Read the predefined distance matrix
- Perform UPGMA clustering
- Construct a rooted phylogenetic tree
- Display branch lengths
- Save the tree as:

```
UPGMA_tree.png
```

- Print the linkage matrix
- Print the generated Newick string

---

## 🧬 Example Distance Matrix

```text
      A  B  C  D
A     0  4  6  8
B     4  0  6  8
C     6  6  0 11
D     8  8 11  0
```

---

## 📈 Example Output

The program generates a rooted phylogenetic tree similar to:

```
          ┌── A
      ┌───┤
      │   └── B
──────┤
      │
      └──── C
           │
           └──── D
```

The tree is also exported as a high-resolution PNG image.

---

## 🌲 Newick Output

Example:

```
(((A:1.000000,B:1.000000):0.500000,C:1.500000):0.750000,D:2.250000);
```

This format can be imported into many phylogenetic visualization tools such as:

- FigTree
- iTOL
- MEGA
- Dendroscope

---

## 🛠 Technologies Used

- Python
- NumPy
- SciPy
- Matplotlib

---

## 📖 Algorithm

UPGMA (Unweighted Pair Group Method with Arithmetic Mean) is a hierarchical clustering algorithm widely used in computational biology and phylogenetics.

Workflow:

1. Compute the pairwise distance matrix.
2. Merge the two closest taxa.
3. Update cluster distances using arithmetic means.
4. Repeat until a single rooted tree is obtained.

---

## 🎯 Applications

- Molecular Evolution
- Phylogenetics
- Comparative Genomics
- Evolutionary Biology
- Bioinformatics Education
- Hierarchical Clustering

---

## 📜 License

This project is licensed under the MIT License.

---

## 👨‍💻 Author

**Hussayn Kazmi**

GitHub: https://github.com/hussaynkazmi

---

⭐ If you found this project useful, consider giving it a star!