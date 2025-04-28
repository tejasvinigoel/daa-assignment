# Implementing Efficient Algorithms for Densest Subgraph Discovery
This repository contains implementations of two algorithms for clique densest subgraph in graphs, along with datasets and their analysis.

| Name | ID | Individual Contributions |
|------|-------|-------|
| Tanisha Agarwal | 2022A7PS0078H | Implementation of Algorithm 1 |
| Thrisha Arrabolu | 2022A7PS0127H | Implementation of Algorithm 4 |
| Snigdha Barui | 2022A7PS0215H | Implementation of Algorithm 4 & Documentation and analysis |
| Tejasvini Goel | 2022A7PS1672H | Implementation of Algorithm 1 & Website Development |
| B.Vaishnavi | 2022A7PS1357H | Analysis of Algorithms & Documentation and analysis |

## Report link:
https://docs.google.com/document/d/1wVR-fJXoiuKwuarOJ8yLDRcS12AT7vNz_sp4gcVc0ss/edit?tab=t.0

## Website link:
https://tejasvinigoel.github.io/daa-assignment2/

## Edited Datasets link:


## Code Execution:
- Ensure the datasets are in the same directory as the source code.
- Once all files are downloaded in the local system, on terminal:

## ⚙️ Compilation Instructions(on windows)

Use the following commands to compile the algorithms:

```bash
# Compile algo_1 algorithm
g++ algo_1.cpp -o algo_1

# Compile algo_4 algorithm
g++ algo_4.cpp -o algo_4



# Run algorithm 1
./algo_1 input.txt

# Run algorithm 4
./algo_4 input.txt
 
```
---
## Observations:
### This assignment has analyzed and compared three prominent algorithms for maximal clique enumeration: Bron–Kerbosch with pivot, Bron–Kerbosch with degeneracy ordering, and the arboricity-based method. Our findings reveal that each algorithm has its strengths for different graph types:

- Bron–Kerbosch with pivot serves as a reliable baseline for small to medium-sized graphs.
- Degeneracy ordering excels on large sparse graphs, particularly social and web networks, offering near-optimal performance.
- The arboricity-based method is most efficient for massive sparse and planar graphs, such as road networks.

All three algorithms produced consistent results, validating their correctness. The arboricity-based method demonstrated superior memory efficiency, using 40% less than the Bron–Kerbosch variants.


  
  
 
