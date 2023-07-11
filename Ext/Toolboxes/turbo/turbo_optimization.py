import sys

from turbo import Turbo1  # import Turbo
import numpy as np
import json
import os
import GPy
import scipy.io
from pyDOE import lhs
from turbo import Turbo1
import time
import csv
from tqdm import tqdm # to track iterations
from torch.quasirandom import SobolEngine


# json file
def inputs_params():
    """
    In this function, we can get the information of the system from a json file, called "scenario.json".
    """

    parameters_file = os.getcwd() + "/scenario.json"
    with open('scenario.json', 'r') as outputfile:

         outputdata = json.load(outputfile)
         numberOfDimensions = outputdata['number_dimensions']
         inputParameters = outputdata['input_parameters']
         maxNumberOfIterations = outputdata['optimization_iterations']
         lowerBounds = [];
         upperBounds = [];

         for i in range(1, numberOfDimensions+1):
             readingInputs = {}
             currentInput = "x" + str(i)
             readingInputs["x" + str(i)] = inputParameters[currentInput]
             readingInputsRanges = readingInputs[currentInput]
             lowerBounds.append(readingInputsRanges[0]);
             upperBounds.append(readingInputsRanges[1]);
            
    return numberOfDimensions, maxNumberOfIterations, lowerBounds, upperBounds

# Check to see that objective output is calculated and saved in a csv file, waiting until the objective value is calculated from Breach 

def waiting_to_get_objectvalue_from_matlab(nameOfFileToWaitFor):

    while not os.path.exists (nameOfFileToWaitFor): 
       time.sleep(1)


class Levy:
    # Reading system information saved in a json file by calling function "inputs_params"
    numberOfDimensions, maxNumberOfIterations, lowerBounds, upperBounds = inputs_params()
    def __init__(self, dim = numberOfDimensions):
         self.dim = dim
         self.lb = lowerBounds
         self.ub = upperBounds

    def __call__(self, x):
        assert len(x) == self.dim
        assert x.ndim == 1
        assert np.all(x <= self.ub) and np.all(x >= self.lb)
        
        global numberOfIterations 
        nameOfCurrentInputFile = "turboToMatlab.csv"
       
        # Writing the input values in a csv file
        with open(nameOfCurrentInputFile, 'w') as writingInput:
            writingInputsToCSV = csv.writer(writingInput)
            for i in range(1, len(x)+1):
                convertingValuesToFloat = np.float(x[i-1])
                x[i-1] = convertingValuesToFloat
             
               
            writingInputsToCSV.writerow(x)
            writingInput.seek(0,2)                    
            size=writingInput.tell()               
            writingInput.truncate(size-2)
             
        #writingInput.close()

        # Tell MATLAB that Turbo has finished writing
        f = open('turboFinishedWriting.dummy', 'w')
        f.close()
        
        nameOfCurrentOutputFile = "matlabToTurbo.csv"
        nameOfFileToWaitFor = 'matlabFinishedWriting.dummy'
   
        # Waiting to get objective value from Breach by calling function "waiting_to_get_objectvalue_from_matlab" 
        waiting_to_get_objectvalue_from_matlab(nameOfFileToWaitFor) 

        # Opening and reading value function saved in a csv file 
        outData = open(nameOfCurrentOutputFile, "r")
        outObj = []
        for out in outData:
             outObj.append(float(out))
    
        outData.close()

        os.remove(nameOfCurrentOutputFile)
        os.remove(nameOfFileToWaitFor)
        numberOfIterations = numberOfIterations + 1 # Increasing numberOfIterations
      
        return outObj[0]




if __name__ == '__main__':
    global numberOfIterations   
    numberOfDimensions, maxNumberOfIterations, lowerBounds, upperBounds = inputs_params()
    numberOfIterations = 1 

    f = Levy(numberOfDimensions)

    turbo1 = Turbo1(
        f=f,  # Handle to objective function
        lb=np.array(lowerBounds),  # Numpy array specifying lower bounds
        ub=np.array(upperBounds),  # Numpy array specifying upper bounds
        n_init=2*numberOfDimensions,  # Number of initial bounds from an Latin hypercube design
        max_evals = maxNumberOfIterations,  # Maximum number of evaluations
        batch_size=10,  # How large batch size TuRBO uses
        verbose=True,  # Print information from each batch
        use_ard=True,  # Set to true if you want to use ARD for the GP kernel
        max_cholesky_size=2000,  # When we switch from Cholesky to Lanczos
        n_training_steps=50,  # Number of steps of ADAM to learn the hypers
        min_cuda=1024,  # Run on the CPU for small datasets
        device="cpu",  # "cpu" or "cuda"
        dtype="float64",  # float64 or float32
    )

    turbo1.optimize()
    sys.exit(0)