import numpy as np
#import numpy.random as rn
from scipy.stats import binom, entropy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='whitegrid')
#
import os
path = '/Users/nanorosas/Packages/ctw'
os.chdir(path)
from lz76 import LZ76
import jpype as jp
jv_path = '/Library/Java/JavaVirtualMachines/jdk-15.0.1.jdk/Contents/MacOS/libjli.dylib'
if not jp.isJVMStarted():
    jp.startJVM(jv_path, '-ea', '-Xmx2048m', '-Djava.class.path=vmm.jar:trove.jar')
    #jp.startJVM(jv_path, '-ea', '-Xmx8192m', '-Djava.class.path=vmm.jar:trove.jar')    
jstr = jp.JPackage('java.lang').String
assert(jp.isJVMStarted())
javify = lambda py_str, ab_dict: jstr(bytearray([ab_dict[s] for s in py_str]))
#
#
#
def ctw_entropy(X, vmm_order=5):
    # Calculate entropy rate of binary signal X via CTW
    alphabet = set(X)
    ab_size = len(alphabet)
    ab_dict = {cc: i for cc,i in zip(sorted(alphabet), range(ab_size))}
    #
    ## Initialise and train model
    vmm = jp.JPackage('vmm.algs').DCTWPredictor()
    vmm.init(ab_size, vmm_order)
    vmm.learn(javify(X, ab_dict))
    #
    ## Test the model
    res = vmm.logEval(javify( X, ab_dict)) / len(X)
    return res
#
#
def lz_entropy(X):
    return LZ76(X) * np.log2(len(X))/len(X)
#
#
def gen_gm(size):
    # Generate data according to the Golden Mean process (Markov order 1, entropy rate=2/3)
    output = np.zeros(size)
    output[0] = binom.rvs(n=1,p=0.5,size=1)
    for i in range(size-1):
        temp = output[i-1]
        if temp==0:
            output[i+1] = 1
        else:
            output[i+1] = binom.rvs(n=1,p=0.5,size=1)
    return output.astype(int)
#
#
def gen_even(size):
    # Generate data according to the Golden Mean process (Markov order 1, entropy rate=2/3)
    output = np.zeros(size)
    output[0] = binom.rvs(n=1,p=0.5,size=1)
    counter = 0
    for i in range(size-1):
        #
        if output[i]==1:
            counter = counter+1
        elif output[i]==0:
            counter = 0
        #
        if counter%2==0:
            output[i+1] = binom.rvs(n=1,p=0.5,size=1)
        else:
            output[i+1] = 1
        #
    return output.astype(int)
#
#
def gen_rrxor(size):
    # XOR process: the value at times t+2 = XOR(values at time t, value at t+1)
    output = binom.rvs(n=1,p=0.5,size=size)
    for i in np.arange(0,size-3,3):
        output[i+2] = (output[i] + output[i+1])%2
    return output
#
#
def gen_nxor(size, p=0.3):
    # noisy XOR process, with noise strength given by p
    output = np.zeros(size)
    output[0] = binom.rvs(n=1,p=0.5,size=1)
    for i in range(size-2):
        output[i+2] = XOR( XOR(output[i+1], output[i]), binom.rvs(n=1,p=p,size=1) )
    return output.astype(int)
#
def XOR(a,b):
    return (a+b)%2






##########################
# Test parameters
#
folder_res = '/Users/nanorosas/Projects/vmm/res/'
#
# type of input data
mode = 'even'
mode = 'gm'
mode = 'iid'
#mode = 'nxor'
#mode = 'rrxor' 
#
# CTW order
vmm_order = 10
#
# Signal lengths to be tested
n_vec = np.arange(100,3000,300)
# Number of iterations per length
nb_samples = 40
#
#
#
##########################
# Calculations
#
if mode == 'even':
    lz76_vals = np.vstack([[  lz_entropy( gen_even(size=n0) ) for n0 in n_vec] for _ in range(nb_samples)])
    ctw_vals  = np.vstack([[ ctw_entropy( gen_even(size=n0), vmm_order= vmm_order) for n0 in n_vec] for _ in range(nb_samples)])    
    true_h = 2/3
#
#
elif mode == 'gm':
    lz76_vals = np.vstack([[  lz_entropy( gen_gm(size=n0) ) for n0 in n_vec] for _ in range(nb_samples)])
    ctw_vals  = np.vstack([[ ctw_entropy( gen_gm(size=n0), vmm_order= vmm_order) for n0 in n_vec] for _ in range(nb_samples)])
    true_h = 2/3
#
#
elif mode == 'iid':
    p = 0.5
    lz76_vals = np.vstack([[  lz_entropy( binom.rvs(1, p, size=n0) ) for n0 in n_vec] for _ in range(nb_samples)])
    ctw_vals  = np.vstack([[ ctw_entropy( binom.rvs(1, p, size=n0), vmm_order= vmm_order) for n0 in n_vec] for _ in range(nb_samples)])    
    P = [ binom.pmf(x, 1, p) for x in range(2) ]
    true_h = entropy(P, base=2)
#
#
elif mode == 'rrxor':
    lz76_vals = np.vstack([[  lz_entropy( gen_rrxor(size=n0) ) for n0 in n_vec] for _ in range(nb_samples)])
    ctw_vals  = np.vstack([[ ctw_entropy( gen_rrxor(size=n0), vmm_order= vmm_order) for n0 in n_vec] for _ in range(nb_samples)])
    true_h = 0.5940493509345 # OF COURSE THIS IS NOT THE RIGHT RESULT...
#
#
elif mode == 'nxor':
    lz76_vals = np.vstack([[  lz_entropy( gen_nxor(size=n0) ) for n0 in n_vec] for _ in range(nb_samples)])
    ctw_vals  = np.vstack([[ ctw_entropy( gen_nxor(size=n0), vmm_order= vmm_order) for n0 in n_vec] for _ in range(nb_samples)])    
    true_h = 0.5940493509345 # NEITHER THIS
#
#
#
#
#
#
#################################
# Arranging results
#
df_lz76 = pd.DataFrame( data=lz76_vals, columns=n_vec).T
df_lz76['Method'] = 'LZ76'
df_ctw  = pd.DataFrame(  data=ctw_vals, columns=n_vec).T
df_ctw['Method'] = 'CTW'
df = pd.concat([df_ctw,df_lz76],axis=0)
df.index.name = 'Window size'
df = df.reset_index()
df_melt = df.melt( id_vars=['Window size','Method'], var_name = 'Iteration', value_name='Estimate')
#
#
##################################
# Save results
#
df.to_csv( folder_res + 'mode_' + mode + '_ctworder_' + str(vmm_order) + '_it_' + str(nb_samples) + '.csv')
df_melt.to_csv( folder_res + 'melt_mode_' + mode + '_ctworder_' + str(vmm_order) + '_it_' + str(nb_samples) + '.csv')
#
#
##################################
# Plotting
#
H = true_h*np.ones(len(n_vec))
xx = np.linspace(0,len(n_vec)-1,len(n_vec))
sns.boxplot( x = 'Window size', y= 'Estimate', hue='Method', palette="muted", data=df_melt)
plt.plot(xx,H,'k--')
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()


