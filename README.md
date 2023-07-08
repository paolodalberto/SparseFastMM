# SparseFastMM
We are playing with sparse matrix sparse matrix multiplication in C and other languages. 

We implement the matrix operation C = A*B: 

* If the matrices are adjiacent matrices: * is sum and + is min
* If the matrices are just sparse matrices: the * is multiplication and + is sum

The sparse matrices are represented in COO format: (row,col, value) and the reason is because the COO is the most flexible although it could be the least efficient. The value can be a element or it can be a BNxBM  compact array. 
In particular, we show that the block allows to exploits: spatial and temporal locality. The spatial locality affect also the vectorization of the computation making more appealing. The temporal locality allows to ammortize the expensive data communication across more operations.
In passing, the operation A*B implies that there is an order how A and B are stored and accessed. Understandbly, we could assume that the oder should not matter using a COO layout.  In practice, to read a  row of A and a column on B for a single dot product, we do not want to read the whole A and the Whole B. 

Two basic ideas emerge when implementing A*B: 
* We need to sort the matrices A and B so that A is frinedly in the rows and B is friendly in the columns
  * Sorting is a modern application and new algorithms are actually found today.
  * Sorting a matrix of N elements: will take N*Log(N), which is often more time consuming that the multiplication itself.     
* The spMspM is a merge sort: and the complexity could be summarized as O(N)*Degree.
  * the preparation of the data and some trickeries makes possbile to achieve such a tight bound
  * You cannot do better than writing the output.

This new commit is basically a study about sorting, parallel sorting, block sparsities, any sparsity and jsut a reminder how C computations are done: Every instruction, variable, every procedure has to be defined and tested. 
I (Paolo) I have found the parallellisation of sorting by a recursive algorithm and core assignment instructive: the computation split like a recursive algorithm and each branch is actually its own thread. Never did it before. 

