# StochasticSimulationFillingAxonDendrite

stochasticSimulationProteinsInDendriteAxon.py
- simulates stochastically the filling process of proteins in axons and dendrites (model with 50 compartments)
- to set diffusion rates, lambda_alpha, lambda_beta go to line 236 and 237
- to set the number of reactions and maximum time go to line 242 and 243
- to set alpha_0 and beta_0 go to line 244 and 245
- output: a txt file with information about CV of the number of proteins in each compartment, a txt file with the data of time differences between two reactions, txt file with the data of the time, txt file with the data of the number of molecules after each time intervall, 2 plots of the average number of molecules of each compartment

plotAgain.py
- extracts the data from the txt files (number of molecules and time) and plots it again
- before using this file you need to run stochasticSimulationProteinsInDendriteAxon.py first to get the txt files
- to choose the simulation that you want to plot again please select the right diffusion rate and lambdas in line 12 and 13
- output: 2 plots of the average number of molecules of each compartment



