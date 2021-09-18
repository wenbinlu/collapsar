import os

fdir = '/Users/wenbinlu/Documents/Research/collapsar/'
fname_list = ['profile12']

for fname in fname_list:
    output = os.popen('python collapsar.py ' + fdir + ' ' + fname).read().strip('\n').split('\t')
    Mbh, abh, Ewind, Eacc, Eacc_ADAF = [float(x) for x in output]
    print('%.3e\t%.3e\t%.3e\t%.3e\t%.3e' % (Mbh, abh, Ewind, Eacc, Eacc_ADAF))
    
    
