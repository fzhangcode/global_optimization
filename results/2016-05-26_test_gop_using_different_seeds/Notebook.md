2016-05-26 Test gop using different seeds
=============================
Purpose
------------
Test GOP using different seeds on data sets with a full factorial arrangement to see if the optimal objective function value is close to 0 and the optimal x and optimal theta are global solutions.

Conclusions
-----------------  
The optimal objective function ||y-xOpt*thetaOpt||^2_2 is 0. 

Background
-----------------


Materials and Equipment
------------------------------

Experimental Protocol
---------------------------
    `python test_gop_seeds.py`
    `python table.py`

x constraints: (it might be faster is using 2*np.sum(abs(x_star)))

	cons = np.sum(abs(x_star))

generate different seeds by timing:

	SEED = int(time.time())
    np.random.seed(SEED)

initial values:

	np.random.random_sample((M,K)) 


parallel processing in "get unique regions using cell numeration" and "solve master problems":

	NUM_CORES = 64

the number of seeds to test on:

	MAXSEED = 20

Results
-----------
Check 'results.csv'


Finish running tests on the data sets of M=1, K=2, N=2.   
tolerance:

	e = 0.01 


It is running the data sets of M=1, K=3, n=2.  
tolerance:

	e = 0.1 



Archived Samples
-------------------------


Archived Computer Data
------------------------------

Prepared by: _______Fan Zhang____________ Date: _____2016.05.26______________


Witnessed by: ________________________