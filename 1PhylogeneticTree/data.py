from numpy import *
import scipy.stats as stat
from scipy.special import gammaln
import util2 as u
import math,time
import numpy
import numpy.random

import traceback

class Datum(object):
	def __init__(self, name,id,a, d,prod):
		self.name = name; self.id=id;
		self.a=a
		self.d=d
		self.prod=prod
			
	def _get_log_mix_wts(self,delta):
		log_den = gammaln(sum(delta)+1)
		log_mix_wts = []
		for i, d_i in enumerate(delta):
			log_num = gammaln(d_i+1)
			for j, d_j in enumerate(delta):
				if i!=j: log_num += gammaln(d_j)
			log_mix_wts.append(log_num-log_den)
		return log_mix_wts
		
	# for multiple samples
	def _log_likelihood(self,phi,rate,norm,purity):
		ntps = len(phi) # multi sample
		return sum([self.__log_likelihood__(phi[tp],rate,tp,norm[tp],purity[tp]) for tp in arange(ntps)])

	def __log_likelihood__(self,phi,rate,tp,norm,purity):
		ll = []
		ll.append(self.__log_complete_likelihood__(phi,rate,tp,norm,purity))
		return u.logsumexp(ll)
		
	# for multiple samples
	def _log_complete_likelihood(self, phi,rate,norm,purity):
		ntps = len(self.a)
		return sum([self.__log_complete_likelihood__(phi,rate, tp,norm,purity) for tp in arange(ntps)])
		
	def __log_complete_likelihood__(self,phi,rate,tp,norm,purity):
		#logmu = numpy.random.randn()+2
        #mu = 10 #exp(logmu)
        #print('rate=',rate)
		return u.log_poisson_likelihood(self.a[tp], self.d[tp], phi,rate,norm,purity,self.prod[self.d[tp]-1])
		#return self.a[tp] * math.log(mu) + (self.d[tp] - self.a[tp]) * math.log(1 - mu) + self._log_bin_norm_const[tp]
		
