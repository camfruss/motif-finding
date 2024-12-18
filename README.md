# motif-finding

TODO
- refactor w/ GibbsSampler class
- SerialGS
- add docs
- refactor has\_converged, max\_iters
- experiment class
- timing harness
- MPI impl
- benchmark against actual algo
- writeup / plots

For serial portion of the report:
- can change data type of PWM to float, double, long double
- run correctness tests
- try out differennt convergence criteria (how long to reach X% correct)
	- fixed iterations
	- multisampling, 10 rounds w/ shorter limit
	- early stopping w/ consensus with upper limit
	- likelihood peaked
- measure variance

- do basic profiling of where the code spends time
- timing as a function of number of motifs included in sequences
- as a function of motif length
- as function of inaccuracies w/in motif
- timing with varying pseudocounts

For MPI:
- convergence
- scaling studies
