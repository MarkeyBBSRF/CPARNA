#######################################################################
COMPILING THE C++ FILES:
To compile the C++ file, execute the following command:
g++ -o mh.o  mh.cpp  util.cpp `gsl-config --cflags --libs`

RUNNING THE SOFTWARE:
1. To run the code on the sample data set "datadn.txt," execute the following command:

python evolve.py 'logProd.npy' 'datadn.txt' 'trees' 'top_k_trees' 'clonal_frequencies' 'llh_trace' 2000 5000 1234 0.907 1 -2 0.2 

# 'datadn.txt': file name of the input data
# 'trees': output folder name where the MCMC trees/samples are saved 
# 'top_k_trees': output file name to save top-k trees in text format. 
# 'clonalFrequencies': output file to save clonal frequencies
# 'llhtrace': output file name to save log likelihood trace

The last seven arguments are the number of MCMC samples, number of Metropolis-Hastings iterations, and a random seed number for initializing the MCMC sampler, normlization factor, purity of the sample, mean of log fold change of expression rate, standard deviation of log fold change of expression rate.




