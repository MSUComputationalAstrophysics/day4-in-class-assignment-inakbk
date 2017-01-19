In class #4 

NOTES

Explicit method: To approximate n+1 we use info at step n

Implicit method: To approximate n+1 we use info of step n+1 which have been approximated with a less accurate model (and info from step n)


The difference in the relative change in energy between the Euler methods and the RK4 is [ 0.49  0.041  0.004] for the timesteps [0.1, 0.01, 0.001]

Question:
Assume you wish to maintain energy to a given level of accuracy - say 0.01% between t = 0 and t = 4π. How many time steps of each of your methods do you need to reach that level of accuracy? 

Answer: 

For Euler you need dt=0.005 and for RK4 you can have it as large as dt=0.7 (!). See figure 2

At dt=0.005 the accuracy of the RK4 does not improve with smaller timestep (see figure 1), it reaches machine precition. 

Question:

How many total floating-point operations is that for each method for the entire integration? 

Answer:





