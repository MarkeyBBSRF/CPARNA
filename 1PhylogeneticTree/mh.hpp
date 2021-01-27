
#include<vector>
#include<math.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>


#include "util.hpp"


void sample_cons_params(struct node nodes[],struct config conf,gsl_rng *rand,int tp);
double multi_param_post(struct node nodes[], struct datum data[], int old,struct config conf);
double param_post(struct node nodes[], struct datum data[], int old,struct config conf,int tp);
void update_params(struct node nodes[], struct config conf);
void get_pi(struct node nodes[], double pi[], struct config conf, int old, int tp);

void load_data(char fname[], struct datum data[], struct config conf);
void load_tree(char fname[], struct node nodes[], struct config conf);
void write_params(char fname[], struct node nodes[], struct config conf);

void mh_loop(struct node nodes[], struct datum data[], struct config conf);

struct config{
	int MH_ITR;
	float MH_STD;
	int NDATA; // no. of data points
    int NDELTA; // no. of data points
	int NNODES; // no. of nodes in the tree
	int TREE_HEIGHT; 	
	int NTPS; // no. of samples / time points
    std::vector<double> NRATE;
    std::vector<double> NORM;
    std::vector<double> PURITY;

};



struct datum{
	
	int id;
	vector<int> a,d;

    double log_ll(double phi, double rate, int tp, double normf,double purityf ){
        int NDELTA = a.size();    ///
        double ll[NDELTA];
        for(int i=0;i<NDELTA;i++){
            ll[i] = log_complete_ll(phi,rate,tp, normf,purityf);
        }
        return logsumexp(ll,NDELTA);
    }
    
    double log_complete_ll(double phi, double rate, int tp, double normf,double purityf){
        return log_poisson_likelihood(a[tp], d[tp], rate, phi,normf,purityf);
    }
    
};

struct node{
	int id;
	vector<double> param,pi;
	vector<double> param1,pi1; // dummy	
	int ndata;
	vector<int> dids;	
	int nchild;
	vector<int> cids; // children ids	
	int ht;	
};
