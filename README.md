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
https://docs.google.com/document/d/1DX92KSZ35noW9xD8V28d1JUWK_9UgQhGyu0-wapfgn4/edit?usp=sharing

## Website link:
https://cozy-sable-bedbeb.netlify.app/

## Edited Datasets link:
https://drive.google.com/drive/folders/1A-_tZQI6Uqf_COlaW6JMCRppisf5UVz1?usp=sharing

## Code Execution:
- Ensure the datasets are in the same directory as the source code.
- Once all files are downloaded in the local system, on terminal:

## ⚙️ Compilation Instructions(on windows)

Use the following commands to compile the algorithms:

```bash
# Compile ELS algorithm
g++ -03 -o els.cpp 

# Compile Tomita algorithm
g++ -03 -o tomita.cpp 

# Compile Chiba-Nishizeki algorithm 
g++ -03 -o chiba.cpp


# Run ELS algorithm
./els.exe

# Run Tomita algorithm
./tomita.exe

# Run Chiba-Nishizeki algorithm 
./chiba.exe
 
```
---
## Observations:
### This assignment has analyzed and compared three prominent algorithms for maximal clique enumeration: Bron–Kerbosch with pivot, Bron–Kerbosch with degeneracy ordering, and the arboricity-based method. Our findings reveal that each algorithm has its strengths for different graph types:

- Bron–Kerbosch with pivot serves as a reliable baseline for small to medium-sized graphs.
- Degeneracy ordering excels on large sparse graphs, particularly social and web networks, offering near-optimal performance.
- The arboricity-based method is most efficient for massive sparse and planar graphs, such as road networks.

All three algorithms produced consistent results, validating their correctness. The arboricity-based method demonstrated superior memory efficiency, using 40% less than the Bron–Kerbosch variants.


  
  
 
