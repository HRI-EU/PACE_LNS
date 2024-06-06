# PACE Solver
Heuristic solver for the [PACE 2024 competition](https://pacechallenge.org/2024/) on the one-sided crossing
minimization problem. The solver uses a large neighborhood search approach and is implemented in C/C++.

To build the solver, run ```make``` in its top directory. This should produce an executable ```main```.

The solver reads the problem from STDIN. It terminates if the execution time exceeds 285 seconds or if a solution with 
zero crossings is found, and outputs the best solution found to STDOUT. To solve a problem specified in a file problem.gr,
the solver can be called as follows:

```cat problem.gr | ./main``` 



