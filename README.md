# Implementing Maximal Clique Enumeration Algorithms
This repository contains implementations of several algorithms for maximal clique enumeration in graphs, along with datasets and their analysis.

| Name | ID | Individual Contributions |
|------|-------|-------|
| Tanisha Agarwal | 2022A7PS0078H | Implementation of Chiba algorithm & Implementation of ELS algorithm |
| Thrisha Arrabolu | 2022A7PS0127H | Implementation of Tomita algorithm & Documentation and analysis |
| Snigdha Barui | 2022A7PS0215H | Implementation of ELS algorithm & Documentation and analysis |
| Tejasvini Goel | 2022A7PS1672H | Implementation of Chiba algorithm & Website Development |
| B.Vaishnavi | 2022A7PS1357H | Implementation of Tomita algorithm & Documentation and analysis |

## Report link:
https://docs.google.com/document/d/1DX92KSZ35noW9xD8V28d1JUWK_9UgQhGyu0-wapfgn4/edit?usp=sharing

## Website link:
https://cozy-sable-bedbeb.netlify.app/

## Edited Datasets link:
https://drive.google.com/drive/folders/1A-_tZQI6Uqf_COlaW6JMCRppisf5UVz1?usp=sharing

## Code Execution:
- Ensure the datasets are in the same directory as the source code.
- Once all files are downloaded in the local system, on terminal:

Compilation:
bash
# Compile ELS algorithm
g++ els.cpp -o els.exe

# Compile Tomita algorithm
g++ tomita.cpp -o tomita.exe

# Compile Chiba-Nishizeki algorithm 
g++ chiba.cpp -o chiba.exe

To run:
bash

#  Run the executable
# run ELS algorithm
./els.exe
# run Tomita algorithm
./tomita.exe
# run Chiba-Nishizeki algorithm 
./chiba.exe


### This assignment has analyzed and compared three prominent algorithms for maximal clique enumeration: Bron–Kerbosch with pivot, Bron–Kerbosch with degeneracy ordering, and the arboricity-based method. Our findings reveal that each algorithm has its strengths for different graph types:

- Bron–Kerbosch with pivot serves as a reliable baseline for small to medium-sized graphs.
- Degeneracy ordering excels on large sparse graphs, particularly social and web networks, offering near-optimal performance.
- The arboricity-based method is most efficient for massive sparse and planar graphs, such as road networks.

All three algorithms produced consistent results, validating their correctness. The arboricity-based method demonstrated superior memory efficiency, using 40% less than the Bron–Kerbosch variants.


  
  
 
