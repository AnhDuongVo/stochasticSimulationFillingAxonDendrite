# stochastic simulation of filling axons/dendrites with proteins

# this program plots the number of proteins of each compartment depending on diffusion, internalization, externalization rate

import numpy as np
from scipy.special import factorial
import matplotlib.pyplot as plt
from scipy import stats,polyval, polyfit, linspace
from copy import deepcopy
import time

###################################################
######## functions ################################
###################################################

def next_values(a0,a,r1,r2):
    """returns values for the next reaction like time difference and reaction according to Gillespie"""
    
    # calculate next time
    new_time_difference = (1/a0)*np.log(1/r1)
    
    # choose next reaction R_mu under the condition that
    # sum(a[i], i=0, mu-1) < r2*a0 <= sum(a[i], i=0, mu)
    mu = 0
    N = r2*a0 - a[mu] 
    
    while N > 0:
        mu += 1
        N = N - a[mu]
    
    return(new_time_difference, mu)
    
def calculate_hi(n,m):
    """calculates the hi with the help of binomial coeff and factorial() 
    since hi is defined as total number of distinct 
    combinations of Ri reactant molecules"""
    
    b=[0]*(n+1)
    b[0]=1
    for i in range(1,n+1):
        b[i]=1
        j=i-1
        while j>0:
            b[j]+=b[j-1]
            j-=1
    hi = b[m]*factorial(m)
    return(hi)

###################################################
######## Gillespie algorithm ######################
###################################################
    
def gillespie_algo(init, rates, sub_stoch, prod_stoch, tmax, n_max):
    """generates a statistically correct trajectory of a stochastic equation

    input:
    s_i = array([s1,...,sN]) number of slots
    init = array([w1,...,wN,e1,...,eN,p]) number of molecules of each species
    rates = array([c1,..cM]) rates of each reaction
    sub_stoch, prod_stoch = stochiometry of substrates and products in matrix form
    tmax = maximum time
    n_max = estimated maximum number of reactions]

    output:
    store_time = array([[t1],[t2],[t3],...]) current time of each intervall
    store_number_molecules = array([[number molecules reaction 0],[number molecules reaction 0],..])
    store_filling_fraction_av = array([F1,...,FN]) average of filling fraction of each synapse
    coefficient_variation = array([CV_1,...,CV_N]) average of coefficient of variation of each synapse
    """
    
    # ****************************   
    # step 0: initialisation
    # ****************************

    # generate a array of two random numbers for step 2
    r1 = np.random.random_sample(n_max)
    r2 = np.random.random_sample(n_max)

    # initialise constant parameters
    stoch = sub_stoch + prod_stoch 
    number_reactions = np.shape(stoch)[0] # number of reactions
    number_species = np.shape(stoch)[1] # number of species

    # initialise current parameters
    current_time = 0
    current_species = init # current number of molecules of each species
    n_counter = 0 # number of already occured reactions
    
    # initialise variables to store time and molecule numbers
    store_time = np.zeros(n_max)
    store_time[n_counter] = current_time
    store_number_molecules = np.zeros((n_max,number_species))
    store_number_molecules[n_counter,:] = current_species 
    store_time_difference = np.zeros((n_max,number_species))   
    store_reactions = [0,0,0,0]
    
    while (current_time < tmax) and (n_counter < (n_max-3)):
        
        # ****************************   
        # step 1: calculate ai and a0
        # ****************************   

        a = np.ones((number_reactions,1))

        for i in range(number_reactions):
            hi = 1  # h1 is defined as the number of distinct 
                    # combinations of Ri reactant molecules 
            for j in range(number_species):
                # check whether the reactant is involved in this reaction
                if sub_stoch[i,j] == 0:
                    continue
                else:
                    # check the reactant has molecules available
                    if current_species[j] <= 0: 
                        hi = 0
                        continue
                    else:
                        hi *= calculate_hi(int(current_species[j]),np.absolute(sub_stoch[i,j]))
                        
            a[i] = hi*rates[i]
        
        a0 = sum(a)

        # ****************************   
        # step 2: calculate the next time difference and reaction
        # ****************************   
        new_time_difference,next_r = next_values(a0,a,r1[n_counter],r2[n_counter])
        store_time_difference[n_counter,:] = np.zeros(number_species)
        store_time_difference[n_counter,:] += new_time_difference 

        # ****************************   
        # step 3: update the system
        # ****************************   

        # update time, number species, counter
        current_time += new_time_difference 
        print(current_time, n_counter)    
        current_species += np.transpose(stoch[next_r,:])
        n_counter += 1

        # store current system
        store_time[n_counter] = current_time
        store_number_molecules[n_counter,:] = current_species 
        
        # update number of reactions
        if next_r in range(50):
            store_reactions[1] += 1
        if next_r in range(50,100):
            store_reactions[0] += 1
        if next_r in range(100,149):
            store_reactions[3] += 1
        if next_r in range(149,198):
            store_reactions[2] += 1
        
    # prepare the final output
    store_time = store_time[:n_counter]
    store_number_molecules = store_number_molecules[:n_counter,:]
    store_time_difference = store_time_difference[:n_counter,:]
    
    return(store_time_difference, current_time, store_time, store_number_molecules, store_reactions)

###################################################
######## Run the main program #####################
###################################################

# set the stochiometry of the substrates and products
sub = []
for j in range(50):
    a = []
    for i in range(50):
        if i == j:
            a.append(-1)
            continue
        a.append(0)
    sub.append(a)
for j in range(50):
    a = []
    for i in range(50):
        a.append(0)
    sub.append(a)
for j in range(49):
    a = []
    for i in range(50):
        if i == j:
            a.append(-1)
            continue
        a.append(0)
    sub.append(a)
for j in range(49):
    a = []
    for i in range(50):
        if i == (j+1):
            a.append(-1)
            continue
        a.append(0)
    sub.append(a)

sub_stoch =np.array(sub)

prod = []
for j in range(50):
    a = []
    for i in range(50):
        a.append(0)
    prod.append(a)
for j in range(50):
    a = []
    for i in range(50):
        if i == j:
            a.append(1)
            continue
        a.append(0)
    prod.append(a)
for j in range(49):
    a = []
    for i in range(50):
        if i == (j+1):
            a.append(1)
            continue
        a.append(0)
    prod.append(a)
for j in range(49):
    a = []
    for i in range(50):
        if i == j:
            a.append(1)
            continue
        a.append(0)
    prod.append(a)

prod_stoch =np.array(prod)

# simulation

cv_r_data = open("cv_reactions","a+")
for diffusion in [[0,0],[1000,1000],[50000,50000],[100000,100000]]:  # set diffusion rates here
    for lambdas in [[30,70],[70,30],[50,50]]:                        # set lambda_alpha and lambda_beta here            
        print("Start time: "+time.strftime("%d.%m.%Y %H:%M:%S"))
        print(lambdas)
        
        # define the initial conditions
        tmax = 30         # set the end time 
        n_max = 1000    # estimate n_max for later arrays
        alpha_0 = 1000    # set alpha_0 and beta_0 here
        beta_0 = 100
        
        # prepare the rates
        rates = []
        
        lambda_alpha = lambdas[0]
        lambda_beta = lambdas[1]
        
        D_plus = diffusion[0]
        D_minus = diffusion[1]
        
        for x in range(50):
            rates.append(beta_0/(np.e**(x/lambda_beta)))

        for x in range(50):
            rates.append(alpha_0/(np.e**(x/lambda_alpha)))

        for x in range(49):
            rates.append(D_plus)

        for x in range(49):
            rates.append(D_minus)

        rates = np.array(rates)

        # prepare initial number of molecules of each species according alpha_0 and beta_0
        init1 = []
        for i in range(50):
            init1.append(alpha_0/beta_0)
        init1 = np.array(init1)

        # gillespie algorithm
        init = deepcopy(init1)
        results = gillespie_algo(init, rates, sub_stoch, prod_stoch, tmax, n_max)
        store_time_difference = results[0]
        current_time = results[1]
        store_time = results[2]
        store_molecules = results[3]
        store_reactions = results[4]

        # txt file for number of molecules and time
        name_sim = str(lambda_alpha)+"_"+str(lambda_beta)+"_"+str(D_plus)+"_"+str(D_minus)
        
        R_data = open("R"+name_sim,"a+")
        time_data = open("time_"+name_sim,"a+")
        time_difference_data = open("time_difference_"+name_sim,"a+")

        for i in range(len(store_molecules)):
            R_data.write(str(list(store_molecules[i]))+",\n")
            time_data.write(str(store_time[i])+"\n")
            time_difference_data.write(str(store_time_difference[i][1])+"\n")
        R_data.close()
        time_data.close()
        time_difference_data.close()
        
        # calculate standard deviation
        average = sum(store_time_difference*store_molecules)
        average /= current_time 
        sd = np.sqrt(sum((store_molecules - average)**2*(store_time_difference/current_time)))

        x = range(1,51)
        y = average
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
        
        # save cv of the average of the protein numbers of each compartment and number of reactions
        sd_compartments = np.std(average)
        av_compartments = sum(average)/50
        cv_r_data.write(name_sim+" CV: "+str(np.round(sd_compartments/av_compartments*100,3))+"\n")
        cv_r_data.write("alpha,beta,D+,D- reaction: "+str(store_reactions)+"\n")
        cv_r_data.write("\n")
        
        print("End run time: "+time.strftime("%d.%m.%Y %H:%M:%S"))
        
cv_r_data.close()
