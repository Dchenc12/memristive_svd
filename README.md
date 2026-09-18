# Memristive singular value decomposition (MSVD)

## Overview

This repository contains source data and simulation code for the article "*Memristive Singular Value Decomposition*", which is published on *Nature Communications* with DOI: https://doi.org/10.1038/s41467-026-76272-2. 

The code demonstrates the MSVD workflow both with and without the selective representation enhanced architecture (SREA), providing comparative analysis of their solutions.

## Requirements

- Matlab 2022/ Matlab 2023
- toolboxes: Fixed-Point Designer, Statistics and Machine Learning Toolbox

## Getting Started

To run the code, please execute the script `msvd.m` with Matlab.

`iris_data.mat` contains the data of [IRIS](https://archive.ics.uci.edu/dataset/53/iris) and used for decomposition.

Output:

```
V ground truth is: 
    0.3616   -0.6565    0.5810    0.3173
   -0.0823   -0.7297   -0.5964   -0.3241
    0.8566    0.1758   -0.0725   -0.4797
    0.3588    0.0747   -0.5491    0.7511

S ground truth is: 
   25.0899    6.0079    3.4205    1.8785

Solving with SREA
V solved with SREA is: 
    0.3623   -0.6486    0.6283    0.6068
   -0.0803   -0.7418   -0.5681   -0.3221
    0.8586    0.1597   -0.0504    0.2334
    0.3537    0.0600   -0.5291    0.6881

S solved with SREA is: 
   24.9202    5.9107    3.4038    2.2869

Solving without SREA
V solved without SREA is: 
    0.3576   -0.1783   -0.5786    0.7456
   -0.0810   -0.7125   -0.2149    0.0566
    0.8600    0.6395   -0.5848   -0.5967
    0.3549    0.2271   -0.5264   -0.2913

S solved without SREA is: 
   24.7901    5.4033    4.5884    2.2926
```



## License

This code repository is covered under the **Apache 2.0** License.

