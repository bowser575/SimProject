import numpy as np
import random as rand
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
#created a bernoulli class
 
class bernoulli():
    def pmf(x,p):
        """
        probability mass function        
        """
        f = p**x*(1-p)**(1-x)
        return f
    
    def mean(p):
        """
        expected value of bernoulli random variable
        """
        return p
    
    def var(p):
        """
        variance of bernoulli random variable
        """
        return p*(1-p)
    
    def std(p):
        """
        standart deviation of bernoulli random variable
        """
        return bernoulli.var(p)**(1/2)
    
    def rvs(p,size=1):
        """
        random variates
        """
        rvs = np.array([])
        for i in range(0,size):
            if np.random.rand() <= p:
                a=1
                rvs = np.append(rvs,a)
            else:
                a=0
                rvs = np.append(rvs,a)
        return rvs
def runTrial(p,numStudents): #one simulation
    students = np.zeros(numStudents)
    students[0] = 1 
    students_sick_tracker = np.zeros(21)
    students_sick_tracker[0] = 3
    dayCounter = 0 
    infected = True
    temp_prob =1-p # used to calculate prob for each day below 
    while(infected):
        for i in range(0,len(students)):
            n_zeros = np.count_nonzero(students)
            p = 1-((temp_prob) ** n_zeros)
            if(students[i] == 0 ):
                if np.random.rand() <= p:
                    students_sick_tracker[i] = 3
                    students[i] = 1
            else:
                students_sick_tracker[i] -= 1  #Subtract a day from time of being sick
                if students_sick_tracker[i] == 0:
                    students[i] = 0 # student is healthy
        dayCounter += 1 
        if(not np.any(students)): # Check if any student is still sick 
            infected = False
    return dayCounter

#rand.seed(1)
p=0.02 # probability of having an accident
bernoulli.mean(p) # return -> 0.2
bernoulli.var(p) # return -> 0.16
bernoulli.std(p) # return -> 0.4



counter = 0
dayList = np.array([])
while(counter < 1000):
    dayList = np.append(dayList, runTrial(p,21))
    counter+=1



print("mean ", np.mean(dayList))   
print("median ", np.median(dayList))   
print("std ", np.std(dayList))  
print("max ",np.max(dayList)) 
print("min ", np.min(dayList)) 


plt.hist(dayList)
plt.xlabel("Days For Pandemic to End")
plt.ylabel("Number of Trials")

plt.savefig("mygraph.png")