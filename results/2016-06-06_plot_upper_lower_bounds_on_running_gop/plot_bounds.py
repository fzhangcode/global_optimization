import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import time

if __name__ == "__main__":
    
            
    # GOP 
    with h5.File('M1_K2_N10_SEED0.hdf5','r') as f:
        SUBD_star = f['SUBD'][...]
        MLBD_star = f['MLBD'][...]
        time_acc = f['time_acc'][...]
                

    # plot GOP bounds by timing
    fig, ax = plt.subplots(figsize=(15, 8))
    plt.figure(1)

    fs =15

    #MAXITER_gop = 29 # 30 second
    #MAXITER_gop = 46 # 1 min
    #MAXITER_gop = 84 # 2 mins
    #MAXITER_gop = 408 # 10 mins
    MAXITER_gop = np.shape(SUBD_star)[0] # 87.33 hours.

    t = np.arange(1, MAXITER_gop - 1, 1) 
    t = np.arange(1, MAXITER_gop - 1, 1) 
    # plot by time
    plt.plot(time_acc[t], SUBD_star[t], 'r-', linewidth = 6, markersize=10, alpha=0.7, label = 'GOP Upper Bound')
    plt.plot(time_acc[t], MLBD_star[t], 'bo-', linewidth = 6, markersize=10, alpha=0.7,  label = 'GOP Lower Bound')    
    xlb = plt.xlabel('Time (s)', fontsize = fs)
    
    # plot by iterations
    #plt.plot(t, SUBD_star[t], 'r-', linewidth = 8, markersize=7.5, label = 'GOP Upper Bound')
    #plt.plot(t, MLBD_star[t], 'b--', linewidth = 8, markersize=7.5, label = 'GOP Lower Bound')    
    #xlb = plt.xlabel('Iterations', fontsize = 24)
             
    #plt.xticks(np.arange(0, MAXITER_gop, 20))
    plt.ylim(-12, 1)
    plt.legend(fontsize = fs, loc = 'lower right')           
    plt.setp(plt.gca().get_xticklabels(), fontsize=fs)
    plt.setp(plt.gca().get_yticklabels(), fontsize=fs)
    plt.tight_layout()
    #plt.show()	
    plt.savefig('figure_4.png')