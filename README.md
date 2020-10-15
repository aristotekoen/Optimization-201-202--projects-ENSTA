# OPT201-202-projects-ENSTA
A series of interdependent projects that we worked on during the course of 2 optimization classes taught by Jean-Charles Gilbert and implemented in MATLAB. The goal was to create an optimization solver capable of implementing many different algorithms such as Newton, quasi-Newton, Sequential quadratic programming, gradient descent and to use static or adaptive timesteps (Armijo rule). 

We were testing the algorithms on the problem of minimizing the potential energy of an articulated chain. We sometimes applied it on a chain fixed at both extremities with no floor below and then also used it by imposing a floor below the chain acting as constraints of the minimization problem. 

You will find a folder with the matlab scripts: 
  - sqp.m corresponds to the optimizer
  - ch.m corresponds to the problem of the articulated chain which allows us to modify the lengths, number of vertices, plot the chain etc. 
  - cholmod.m is a function to do a Cholesky factorization
  - plancher.m is a script that generates the floor under the chain for the constrained problem
  - chs.m is the implementation of the chain based on ch.m

You will then find the pdfs we submitted for each of the project phase. 
