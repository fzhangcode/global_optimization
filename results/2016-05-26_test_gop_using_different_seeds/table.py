__author__ = 'fzhang'

import h5py
import numpy as np

# Generate a table to characterize the results from test_gop_seeds.

if __name__ == "__main__":
        
        results = []    
        MAXSEED = 20
        num_seed = 0
        
        table=open('table.csv', 'a') 

        for i in range(MAXSEED):
            h5Filename = "M1_K2_N2_SEED%d.hdf5" %i
            results = []

            f = h5py.File(h5Filename, 'r')            
            #SEED = f['SEED'][...]
            initialx = f['initialx'][...]
            xOpt = f['xOpt'][...]
            x_true = f['x_true'][...]
            thetaOpt = f['thetaOpt'][...]
            theta_true = f['theta_true'][...]
            objOpt = f['objOpt'][...]
            MLBD = f['MLBD'][...]
            SUBD = f['SUBD'][...]
            #time = f['time'][...]

            results.append([float(initialx[0][0]), float(initialx[0][1]), float(xOpt[-2][0][0]), float(xOpt[-2][0][1]), float(x_true[0][0]), float(x_true[0][1]),\
            float(thetaOpt[-1][0][0]), float(thetaOpt[-1][0][1]), float(thetaOpt[-1][1][0]), float(thetaOpt[-1][1][1]), float(theta_true[0][0]), float(theta_true[0][1]),\
            float(theta_true[1][0]), float(theta_true[1][1]), objOpt, MLBD[-1], SUBD[-1] ])

            items = "[initialx[0], initialx[1], xOpt[0], xOpt[1], x_true[0], x_true[1], thetaOpt[0][0], thetaOpt[0][1], thetaOpt[1][0], thetaOpt[1][1],\
            theta_true[0][0], theta_true[0][1], theta_true[1][0], theta_true[1][1], objOpt, MLBD, SUBD]"
            
            if i == 0 :
                np.savetxt(table, results, fmt ="%.3f", delimiter=",", header = str(items))           
            else:
                np.savetxt(table, results, fmt ="%.3f", delimiter=",")