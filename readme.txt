Distribution-free approach for endogenously selected samples: Code and Data

The text below  explains program files used for simulation studies and data analysis in Qian and Xie (2021).

PLEASE CITE AS:

Qian, Yi, and Hui Xie. "Simplifying bias correction for selective sampling: A unified distribution-free approach to handling endogenously selected samples."
 Marketing Science 41.2 (2022): 336-360.


There are five files as follows:

1. Rcode.R
R  Code to analyze the simulated and real-life application datasets. 

2. pfunc.R
Contains a set of R functions called by Rcode.R for handling endogenously selected samples using semiparametric odds ratio mdoel. 

3. Outsample.exe
Fortran application called by the R function fitsor() to estimate endogenously selected samples (aka Outcome-dependent samples)
using  semiparametric odds ratio model. THIS APPLICATIONS HAS BEEN RECOMPLIED ON MARCH 2023 TO BE RUN ON NON-NATIVE WINDOWS COMPUTERS.  

4. libiomp5md.dll
The DLL library used by the IMSL numeric library for the Fortran Application. This DLL enables the Fortran application to run on computers that do not have
Intel Fortran complier and IMSL numeric library installed on the local computers.  

5. Data File: visits.csv

	Description:
        incomec -  Household income;
        distc - Distance to Store
        distc2 - distc^2
        agec - Age of customer
        marriedc - Marriage status of the customer
        kidsc  -  Number of kids at home.
	totalvisit - total number of store visits within the first year. 
        The variables ending with the letter 'c' mean that these variables are standardized. Data are also masked to protect original data values. 


To run the program, download all files and save them under the same folder in Windows operating system.
Then source Rcode.R in R (and the Stata portion in the file separately in Stata) to run the programs. These codes were developed under R 3.6.3 and Stata 14. 

