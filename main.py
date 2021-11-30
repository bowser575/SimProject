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
def runTrial(p,numStudents): #one simulation with reinfection 
    students = np.zeros(numStudents)
    students[0] = 1 
    students_sick_tracker = np.zeros(21)
    has_sick = np.zeros(21)
    has_sick[0] = True
    students_sick_tracker[0] = 3
    dayCounter = 0 
    infected = True
    temp_prob =1-p # used to calculate prob for each day below 
    while(infected):
        for i in range(0,len(students)):
            n_zeros = np.count_nonzero(students) #number of sick students
            p = 1-((temp_prob) ** n_zeros)
            if(students[i] == 0 ):
                if np.random.rand() <= p:
                    if has_sick[i] is not True:
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

def runTrialNoReinfection(p,numStudents): #one simulation w/o reinfection
    students = np.zeros(numStudents)
    students[0] = 1 
    students_sick_tracker = np.zeros(21)
    has_sick = np.zeros(21)
    has_sick[0]= True
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
                    if not has_sick[i]:
                        students_sick_tracker[i] = 3
                        students[i] = 1
                        has_sick[i] = True
            else:
                students_sick_tracker[i] -= 1  #Subtract a day from time of being sick
                if students_sick_tracker[i] == 0:
                    students[i] = 0 # student is healthy
        #matrix[index][dayCounter] = np.count_nonzero(students)
        dayCounter += 1 

        if(not np.any(students)): # Check if any student is still sick 
            infected = False
    return dayCounter#,matrix
def runTrialNoReinfectionOneDay(p,numStudents): #one simulation with shorter infection time 
    students = np.zeros(numStudents)
    students[20] = 1 #must start the student at the back of the array or else the loop breaks on the first day
    students_sick_tracker = np.zeros(21)
    has_sick = np.zeros(21)
    has_sick[20]= True
    students_sick_tracker[20] = 1
    dayCounter = 0 
    infected = True
    temp_prob =1-p # used to calculate prob for each day below 
    while(infected):
        for i in range(0,len(students)):
            n_zeros = np.count_nonzero(students)
            p = 1-((temp_prob) ** n_zeros)
            if(students[i] == 0 ):
                if np.random.rand() <= p:
                    if not has_sick[i]:
                        students_sick_tracker[i] = 1
                        students[i] = 1
                        has_sick[i] = True
            else:
                students_sick_tracker[i] -= 1  #Subtract a day from time of being sick
                if students_sick_tracker[i] == 0:
                    students[i] = 0 # student is healthy
        #matrix[index][dayCounter] = np.count_nonzero(students)
        dayCounter += 1 

        if(not np.any(students)): # Check if any student is still sick 
            infected = False
    return dayCounter
def runTrialNoReinfectionMatrix(p,numStudents,matrix,index): #one simulation with a matrix to keep track of avg number of students infected
    students = np.zeros(numStudents)
    students[0] = 3
    students_sick_tracker = np.zeros(21)
    has_sick = np.zeros(21)
    has_sick[0]= True
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
                    if not has_sick[i]:
                        students_sick_tracker[i] = 3
                        students[i] = 1
                        has_sick[i] = True
            else:
                students_sick_tracker[i] -= 1  #Subtract a day from time of being sick
                if students_sick_tracker[i] == 0:
                    students[i] = 0 # student is healthy
        matrix[index][dayCounter] = np.count_nonzero(students)
        dayCounter += 1 

        if(not np.any(students)): # Check if any student is still sick 
            infected = False
    return dayCounter,matrix

if __name__ == "__main__":
    #rand.seed(1) un comment if you want repeatable results
    p=0.02 # probability of getting sick
    matrix = np.zeros((1000,1000)) 
    counter = 0
    dayList = np.array([])
    while(counter < 1000):
        #temp,matrix  = runTrialNoReinfectionMatrix(p,21,matrix,counter) # matrix example
        temp  = runTrialNoReinfection(p,21)
        dayList = np.append(dayList, temp)
        counter+=1

    # maxDays = np.max(dayList)
    # newmatrix = matrix[:,:int(maxDays)]
    # avergae = newmatrix.mean(axis = 0)
    # np.savetxt("foo.csv", avergae, delimiter=",")

    plt.hist(dayList)
    plt.xlabel("Days For Pandemic to End")
    plt.ylabel("Number of Trials")
    plt.savefig("BaseCaseGraph.png")

    print("mean ", np.mean(dayList))   
    print("median ", np.median(dayList))   
    print("std ", np.std(dayList))  
    print("max ",np.max(dayList)) 
    print("min ", np.min(dayList)) 

    p = .01
    dayList = np.array([])
    bins = []

    while(p <.31):
        counter = 0
        print(p)
        tempList = np.array([])
        while(counter < 5000):
            tempList = np.append(dayList, runTrialNoReinfection(p,21))
            counter+=1
        dayList = np.append(dayList, np.mean(tempList))
        p+=.01
    yaxis = [.01,.02,.03,.04,.05,.06,.07,.08,.09,.1,.11,.12,.13,.14,.15,.16,.17,.18,.19,.2,.21,.22,.23,.24,.25,.26,.27,.28,.29,.30]




    plt.plot(yaxis, dayList)
    plt.xlabel('Prob')
    plt.ylabel('Avg # days')

    plt.savefig("BaseCaseExtended.png")