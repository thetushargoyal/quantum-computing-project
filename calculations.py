#Quantum Route Optimizer - calculations module. Supplement to bachelor thesis made as part of the Business Economics & IT programme at KEA, by Patrik Žori and Kieran Olivier Holm, for A.P. Møller - Mærsk, 2020

from __main__ import debug
import numpy as np
import time
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
from collections import defaultdict
from functools import partial

# def func_CreateDistanceMatrix(visitList,distanceData):  #This function will create an m by n matrix of distances between each of the countries to be visited, using distance data from the CERDI database
#     countryNum = len(visitList)
#     inputMatrix = np.zeros((countryNum,countryNum)) #Creates a symmetric 2D matrix filled with 0s, the dimensions being equal to the amount of selected countries
#     rowIndex=0
#     for x in visitList: #For every selected country...
#         addToMatrix = []
#         for y in visitList: #...pair with every selected country and...
#             if x == y:
#                 addToMatrix.append(0)   #...set the distance to 0 if country x and country y are the same or...
#             else:
#                 distance = distanceData[(x,y)]  #...use the distance from the database to connect the 2 countries.
#                 addToMatrix.append(distance)    #Append the latest distance to an array...
#         inputMatrix[rowIndex]=addToMatrix   #...then, once all distances for the first selected country have been retrieved, append the array as a row of the matrix.
#         rowIndex+=1
#     inputMatrix = np.vstack((np.zeros((1,countryNum)),inputMatrix)) #Adds a row of 0s at the top of the matrix to optimize encoding and improve result accuracy.
#     print(f'Distance matrix (an m by n matrix of distances between each of the countries to be visited):\n{inputMatrix}') if debug==True else ''
#     return inputMatrix


import numpy as np

def func_CreateDistanceMatrix(visitList, distanceData):
    countryNum = len(visitList)
    inputMatrix = np.zeros((countryNum, countryNum))

    for rowIndex, x in enumerate(visitList):
        for colIndex, y in enumerate(visitList):
            if x == y:
                inputMatrix[rowIndex, colIndex] = 0
            else:
                # Check if (x,y) pair exists, otherwise check (y,x)
                if (x, y) in distanceData:
                    inputMatrix[rowIndex, colIndex] = distanceData[(x, y)]
                elif (y, x) in distanceData:
                    inputMatrix[rowIndex, colIndex] = distanceData[(y, x)]
                else:
                    # Handle the missing pair
                    raise ValueError(f"No distance data for pair: ({x}, {y}) or ({y}, {x})")
                    
    # Add a row of zeros at the top of the matrix (if necessary for your algorithm)
    inputMatrix = np.vstack((np.zeros(countryNum), inputMatrix))
    
    # Debugging print statement
    if debug:
        print(f'Distance matrix:\n{inputMatrix}')
    
    return input

def func_createQUBO(distanceMatrix):    #Formats the distance matrix as a quadratic unconstrained binary optimization (QUBO) problem, needed to perform calculation on an annealer.
    maxDistance = np.max(np.array(distanceMatrix))
    normalizedMatrix = distanceMatrix / maxDistance #Normalizes the distance matrix, so that all values are between 0 and 1, proportionally to all values. Normalized data can be processed more efficiently.
    print(f'Normalized matrix (the distance matrix normalized to values between 0 and 1 for optimal calculation):\n{normalizedMatrix}') if debug==True else ''
    n = len(normalizedMatrix)
    x = partial(func_binarize, num_variables=n) #The partial() function will prepare a function for execution, but will not actually do so until the x is called elsewhere.
    qubo = defaultdict(float)
    if n < 10:  #Sets a multiplier constant, which will be used to apply weights to wanted/unwanted solutions, the constant depending on the size of the problem.
        const_mul = 8500
    else:
        const_mul = 17000
    #Adding time (rows of the distance matrix) constraints to make sure that no two countries are visited at the same time, row by row (lines 41 to 46 adapted from tsp-demo-unitary-fund). The first row is skipped:
    for row in range(1, n-1):
        for i in range(1, n-1):
            qubo[(x(row, i), x(row, i))] += -const_mul  #Using the func_binarize function, this will tie a combination of 2 countries to a single qubit, represented in binary form (for example: 0 0 0 1 0 0, meaning that this segment of the route is followed at time slot 4). While doing so, it adds the time constraint to qubits where the country is paired to itself to ensure that this option will have too high energy at the end of annealing to be considered a valid solution.
            for j in range(i+1, n-1):
                qubo[(x(row, i), x(row, j))] += 2 * const_mul
                qubo[(x(row, j), x(row, i))] += 2 * const_mul
    print(f'\nQUBO with time constraints:\n{qubo}') if debug==True else ''
    #Adding location (columns of the distance matrix) to make sure that the result does not suggest staying in the country the traveler is already located in, column by column (lines 49 to 54 adapted from tsp-demo-unitary-fund). The first row is skipped:
    for col in range(1, n-1):
        for i in range(1, n-1):
            qubo[(x(i, col), x(i, col))] += -const_mul  #Using the func_binarize function, this will tie a combination of 2 countries to a single qubit, represented in binary form. While doing so, it adds the location constraint to qubits to ensure that this option will have too high energy at the end of annealing to be considered a valid solution.
            for j in range(i+1, n-1):
                qubo[(x(i, col), x(j, col))] += 2 * const_mul
                qubo[(x(j, col), x(i, col))] += 2 * const_mul
    print(f'\nQUBO with time and locations constraints:\n{qubo}') if debug==True else ''
    #The objective function which maps the distance matrix onto the QUBO (lines 57 to 66 adapted from tsp-demo-unitary-fund):
    for i in range(1, n-1):
        dist_to_first =  normalizedMatrix[0, i]
        dist_to_last =  normalizedMatrix[-1, i]
        qubo[x(1, i), x(1, i)] += 1 * dist_to_first
        qubo[(x(n-2, i), x(n-2, i))] += 1 * dist_to_last
        for j in range(1, n-1):
            if i != j:
                for step in range(1, n-2):
                    dist = normalizedMatrix[i, j]
                    qubo[(x(step, i), x((step+1), j))] += 1 * dist
    print(f'\nFinal QUBO:\n{qubo}') if debug==True else ''
    return qubo,n

def func_solveTSPdwave(distanceMatrixLength,token,url,solver,qubo):  #Sends the QUBO to the annealer and retrieves the result
    if distanceMatrixLength > 7:    #The numReads variable specifies how many times the annealing will be performed on the QUBO, thus changing the amount of returned samples. The variable is scaled with the amount of countries to be visited, as more countries require more samples to get an accurate result.
        numReads = 4000
    elif distanceMatrixLength > 14:
        numReads = 7000
    else:
        numReads = 2000
    if '2000Q' in solver:   #Sets the chain strength, which scales with the number of countries. In addition, the Advantage system provides better results with a lower chain strength than the 2000Q system.
        chain = distanceMatrixLength*100
    else:
        chain = int(round((distanceMatrixLength*100)/1.5,0))
    timeStart = time.time() #Makes a timestamp of when the annealing process began
    try:
        binaryResult = EmbeddingComposite(DWaveSampler(token=token, endpoint=url, solver=solver)).sample_qubo(qubo, chain_strength=chain, num_reads=numReads)   #Sends the QUBO to the annealer
    except Exception as error:  #Attempts to catch an exception if one is returned by the D-Wave system
        print(f'An error has occured while trying to send the problem to the annealer: {error}')
        print('The calculation has been cancelled, returning to main menu.')
        return None,None,True
    print(f'\nBinarized optimal solution returned by the annealer:\n{binaryResult}') if debug==True else ''
    return binaryResult,timeStart,None

def func_decodeResult(binaryResult):   #Decodes the binary result returned by the annealer back to a human-readable route. Lines 92 to 102 adapted from quantum_tsp.
    distribution = {}
    prevEnergy = binaryResult.record[0].energy
    for record in binaryResult.record:  #Iterates hrough all solutions to generate a full solution distribution, and also finds the lowest-energy (optimal) result
        sample = record[0]
        binarySolution = [node for node in sample]
        solution = func_debinarize(binarySolution)
        distribution[tuple(solution)] = (record.energy, record.num_occurrences)
        if record.energy <= prevEnergy and not record.num_occurrences < binaryResult.record[0].num_occurrences/2:
            decodedResult = solution
        prevEnergy = record.energy
    timeEnd = time.time() #Makes a timestamp of when the annealing process ended
    return decodedResult,distribution,timeEnd

def func_binarize(i, j, num_variables): #Used to convert the input distance matrix into binary pairs of locations. Function adapted from tsp-demo-unitary-fund.
    return i * num_variables + j

def func_debinarize(binarySolution): #Decodes the binary result returned by the annealer into a list of location indexes. Function adapted from quantum_tsp.
    points_order = []
    number_of_points = int(np.sqrt(len(binarySolution)))
    for p in range(number_of_points):
        for j in range(number_of_points):
            if binarySolution[(number_of_points) * p + j] == 1:
                points_order.append(j)
    return points_order

def func_calculateCost(cost_matrix, solution):  #Used to calculate the cost of the optimal solution, or the cost of all solutions (depends on how the function was called in main.py). Function adapted from quantum_tsp.
    cost = 0
    for i in range(len(solution)):
        a = i%len(solution)
        b = (i+1)%len(solution)
        cost += cost_matrix[solution[a]][solution[b]]
    return cost

def func_solveTSPbruteforce(distanceMatrixLength,distanceMatrix):   #Used to calculate the optimal route using a classical brute-force approach. Function adapted from tsp-demo-unitary-fund.
    import itertools
    number_of_nodes = distanceMatrixLength-1
    initial_order = range(0, number_of_nodes)
    bruteStart = time.time()
    all_permutations = [list(x) for x in itertools.permutations(initial_order)]
    cost_matrix = distanceMatrix
    best_permutation = all_permutations[0]
    best_cost = func_calculateCost(cost_matrix, all_permutations[0])
    for permutation in all_permutations:
        current_cost = func_calculateCost(cost_matrix, permutation)
        if current_cost < best_cost:
            best_permutation = permutation
            best_cost = current_cost
    bruteTime = time.time() - bruteStart
    return best_permutation,best_cost,bruteTime
