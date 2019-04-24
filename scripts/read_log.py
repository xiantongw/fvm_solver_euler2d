#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 01:00:24 2019

@author: xtwang
"""

import matplotlib.pyplot as plt

def read_log(logfilename):
    f = open(logfilename, 'r')
    line = f.readline()
    niter = []
    residual = []
    while line:
        line_content = line.split(' ')
        if line_content[0] == 'niter:':
            niter.append(int(line_content[1]))
            residual.append(float(line_content[3]))
        line = f.readline()
    f.close()
    return niter, residual

niter0, residual0 = read_log('bump0_fs_1st.log')
niter1, residual1 = read_log('bump1_fs_1st.log')
niter2, residual2 = read_log('bump2_fs_1st.log')
niter3, residual3 = read_log('bump3_fs_1st.log')
niter4, residual4 = read_log('bump4_fs_1st.log')
plt.figure()
plt.yscale('log')
plt.plot(niter0, residual0)
plt.plot(niter1, residual1)
plt.plot(niter2, residual2)
plt.plot(niter3, residual3)
plt.plot(niter4, residual4)
plt.xlabel('Number of Iterations')
plt.ylabel('Residual Norm')
plt.title('Free-stream and Free-stream preservation test')
plt.savefig('free_stream_1st.eps')

niter0, residual0 = read_log('bump0_fs.log')
niter1, residual1 = read_log('bump1_fs.log')
niter2, residual2 = read_log('bump2_fs.log')
niter3, residual3 = read_log('bump3_fs.log')
niter4, residual4 = read_log('bump4_fs.log')
plt.figure()
plt.yscale('log')
plt.plot(niter0, residual0)
plt.plot(niter1, residual1)
plt.plot(niter2, residual2)
plt.plot(niter3, residual3)
plt.plot(niter4, residual4)
plt.xlabel('Number of Iterations')
plt.ylabel('Residual Norm')
plt.title('Free-stream and Free-stream preservation test')
plt.savefig('free_stream_2nd.eps')


niter0, residual0 = read_log('bump0_ss_1st.log')
niter1, residual1 = read_log('bump1_ss_1st.log')
niter2, residual2 = read_log('bump2_ss_1st.log')
niter3, residual3 = read_log('bump3_ss_1st.log')
niter4, residual4 = read_log('bump4_ss_1st.log')
plt.figure()
plt.yscale('log')
plt.plot(niter0, residual0, label="bump0")
plt.plot(niter1, residual1, label="bump1")
plt.plot(niter2, residual2, label="bump2")
plt.plot(niter3, residual3, label="bump3")
plt.plot(niter4, residual4, label="bump4")
plt.xlabel("Number of Iterations")
plt.ylabel("Residual Norm")
plt.legend()
plt.title("Residual Norm History - First Order Scheme")
plt.savefig('1st_converge.eps')

niter0, residual0 = read_log('bump0_ss_2nd.log')
niter1, residual1 = read_log('bump1_ss_2nd.log')
niter2, residual2 = read_log('bump2_ss_2nd.log')
niter3, residual3 = read_log('bump3_ss_2nd.log')
niter4, residual4 = read_log('bump4_ss_2nd.log')
plt.figure()
plt.yscale('log')
plt.plot(niter0, residual0, label="bump0")
plt.plot(niter1, residual1, label="bump1")
plt.plot(niter2, residual2, label="bump2")
plt.plot(niter3, residual3, label="bump3")
plt.plot(niter4, residual4, label="bump4")
plt.xlabel("Number of Iterations")
plt.ylabel("Residual Norm")
plt.legend()
plt.title("Residual Norm History - Second Order Scheme")
plt.savefig('2nd_converge.eps')
