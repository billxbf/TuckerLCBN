# Tucker LCBN
The Tucker decomposition is a generalization of the matrix singular value decomposition to tensors, which are multidimensional arrays. We seek to apply the decomposition as a dimension reduction technique in order to analyze large functional magnetic resonance imaging (f-MRI) datasets of human brains.

Neuroscientists are particularly interested in correlation among different areas in the brain, but computing and storing pairwise correlations between all pairs of brain areas can be infeasible, especially when the data set includes multiple participants and multiple trials. The current practice is to downsample the data in order to reduce the number of brain areas in the data, but this process loses information. We show that the dimension reduction via Tucker decomposition can be computed without explicitly computing and storing all correlations, making data analysis with the original granularity feasible and efficient. We demonstrate the advantage of using the full granularity to answer scientific questions about the data, including classifying participants across multiple trials.

## Paperwork
[Efficient Computation of Tucker Decomposition of Correlation-Based Tensors](https://github.com/billxbf/TuckerLCBN/blob/master/LaTeX/paper.pdf)


### Notes
F-MRI data we use is private to protect privacy of participants. 
Special thanks to Dr. Grey Ballard for project advising and parallel implementation (TuckerMPI).
