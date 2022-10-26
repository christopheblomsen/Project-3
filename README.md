# Project 3

The hand in for Project 3 in FYS4150.

## Code
Can be built from using the `makefile` in the `code` directory. 

```
make
```

This will compile and link the `main` and `ex9` programs which is what is mainly used. To run said programs 

```
./main
./ex9
```

Warning, `ex9` takes around 3 hours to run. So do not use it unless the data files are corrupted.

Lastly we have the `give_me_omega.cpp`

```
g++ -I include src/* give_me_omega.cpp -o give_me_omega -larmadillo -O2
```

Which can be run as
```
./give_me_omega
```

## Data files
All the data files can be found in `code/data`. Where we have saved all the neccessary datafiles needed for the figuers. 
In this directory you can also find the `python` scripts to plot said figures.

## Report directory
The report source code can be found in the `report` directory. To compile the pdf

```
pdflatex main
biber main
pdflatex main
pdflatex main
```

`refs.org` is the refrence file needed for `biber`
