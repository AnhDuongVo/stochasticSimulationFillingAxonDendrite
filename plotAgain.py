# this programm extract the data from the txt files (number of molecules and time) and plot it again

import numpy as np
from scipy.special import factorial
import matplotlib.pyplot as plt
from scipy import stats,polyval, polyfit, linspace
from copy import deepcopy
import time

receptors_all = []

for diffusion in [[0,0]]: # set the diffusion rate 
    for lambdas in [[30,70]]: # set the lambda_alpha and lambda_beta
        # get data of number of molecules
        r = "["
        data_text = "R"+str(lambdas[0])+"_"+str(lambdas[1])+"_"+str(diffusion[0])+"_"+str(diffusion[1])
        data = open(data_text)
        for line in data:
            r += line.rstrip()
        r = r[:-1] 
        r += "]"

        data.close()
        store_molecules = np.array(eval(r))
        print(diffusion, lambdas)
        receptors_all.append([diffusion, lambdas, r])
        
        # get data of time difference
        t = []
        current_time = 0
        data_text = "time_difference_"+str(lambdas[0])+"_"+str(lambdas[1])+"_"+str(diffusion[0])+"_"+str(diffusion[1])
        data = open(data_text)
        for line in data:
            t.append(np.ones(50)*float(line))
            current_time += float(line)
        data.close()
        store_time_difference = np.array(t)
        print(t[:3])
        
        print(diffusion, lambdas)
        receptors_all.append([diffusion, lambdas, r])
        
        # calculate standard deviation
        average = sum(store_time_difference*store_molecules)
        average /= current_time 
        sd = np.sqrt(sum((store_molecules - average)**2*(store_time_difference/current_time)))

        x = range(1,51)
        y = average
        print(y)
        e = sd

        # linear regression
        m,b = polyfit(x, y, 1) 

        # plot average receptors - errorbars
        f = plt.figure(figsize=(4,3))
        font = {'family' : 'serif',
                'weight' : 'normal',
                'size'   : 12}
        plt.rc('font', **font)
        plt.rc('font', serif='Times New Roman') 
        plt.rc('text', usetex=True)

        X = np.linspace(1.0,50,500,endpoint=True)
        Y = 10*np.exp(-(1/lambdas[0]-1/lambdas[1])*X)
        for i in range(50):
            if i ==0:
                plt.plot(X,Y,'-',color='red',label=r'$f(x)=10 \cdot exp(-(\frac{1}{\lambda_\alpha}-\frac{1}{\lambda_\beta})x)$')
            else:
                plt.plot(X,Y,'-',color='red')
        plt.errorbar(x, y, e, linestyle='None',color="dimgrey", marker='^')
        label_text = r'$regression(x)='+str(np.round(m,2))+'x+'+str(np.round(b,2))+'$'
        plt.plot(x, m*x+b,'-',color='black',label=label_text) 
        plt.yticks([2,4,6,8,10,12,14,16])
        plt.ylabel('Average number of proteins', fontsize=12)
        plt.xlabel('Compartments', fontsize=12)
        title_text = r'$\lambda_\alpha = '+str(lambdas[0])+r', \lambda_\beta = '+str(lambdas[1])+', D_+ = '+str(diffusion[0])+', D_- = '+str(diffusion[1])+'$'
        plt.title(title_text, fontsize=12)
        plt.legend(loc=1, fontsize=7)
        plt.show()
        file_text = "regression"+str(lambdas[0])+"_"+str(lambdas[1])+"_"+str(diffusion[0])+"_"+str(diffusion[1])+".pdf"
        f.savefig(file_text, bbox_inches='tight')

        # plot average receptors - without errorbars
        f = plt.figure(figsize=(4,3))
        font = {'family' : 'serif',
                'weight' : 'normal',
                'size'   : 12}
        plt.rc('font', **font)
        plt.rc('font', serif='Times New Roman') 
        plt.rc('text', usetex=True)

        X = np.linspace(1.0,50,500,endpoint=True)
        Y = 10*np.exp(-(1/lambdas[0]-1/lambdas[1])*X)
        for i in range(50):
            if i ==0:
                plt.plot(X,Y,'-',color='red',label=r'$f(x)=10 \cdot exp(-(\frac{1}{\lambda_\alpha}-\frac{1}{\lambda_\beta})x)$')
            else:
                plt.plot(X,Y,'-',color='red')
        label_text = r'$regression(x)='+str(np.round(m,2))+'x+'+str(np.round(b,2))+'$'
        plt.plot(x, m*x+b,'-',color='black',label=label_text) 
        plt.yticks([2,4,6,8,10,12,14,16])
        plt.ylabel('Average number of proteins', fontsize=12)
        plt.xlabel('Compartments', fontsize=12)
        title_text = r'$\lambda_\alpha = '+str(lambdas[0])+r', \lambda_\beta = '+str(lambdas[1])+', D_+ = '+str(diffusion[0])+', D_- = '+str(diffusion[1])+'$'
        plt.title(title_text, fontsize=12)
        plt.legend(loc=1, fontsize=7)
        plt.show()
        file_text = "regression2"+str(lambdas[0])+"_"+str(lambdas[1])+"_"+str(diffusion[0])+"_"+str(diffusion[1])+".pdf"
        f.savefig(file_text, bbox_inches='tight')
