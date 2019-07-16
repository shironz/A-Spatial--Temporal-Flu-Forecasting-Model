This is a course project that we develop a spatical-temporal forecasting model to predict the flu in U.S. We use Bayesian method. For the requirements about this project, please see "Project1.pdf".


1, We first develop a simulation experiments to show that our method could recove the true parameters, the codes are in "simulation.r".

2, We then apply the codes to the real data. For the variables we use, please see "variables.txt". The cleaned data are in folder "data", in which "w.txt" is the adjacent matrix for U.S. states. All the applications with real data are contained in "flu-predict.r". We train our model using data from Octoboer 2010 to December 2016. We run our MCMC for 200K times.

3, Finally, after geting the MCMC chain, we use codes "outcome6.r" to summarize the estimates and also use our trained model to predict the number of flu from Janually 2017 to December 2017, and also compare the predicted flu with true flu.

4, We report our outcomes in "MATH8820_Project1.pdf".



