## Optimising Chip Placement Using Quantum Annealing
This project aims to optimize the placement of nodes (components or blocks) on a chip canvas by leveraging the capabilities of quantum annealing. The primary objectives are to minimize wire length and congestion while ensuring that the chip density remains within acceptable limits. This is achieved by formulating the chip placement problem as a Quadratic Unconstrained Binary Optimization (QUBO) problem and solving it using D-Wave's Quantum Annealers.

## Problem Description
Chip placement is a critical and complex aspect of chip design, involving the positioning of macros (like SRAMs) and standard cells (logic gates) on a chip canvas. The optimization must consider power, performance, and area (PPA) metrics, as well as constraints on placement density and routing congestion. The complexity of the problem is due to the large size of netlist graphs, which can contain millions to billions of nodes, and the high computational cost of evaluating target metrics using traditional methods.

## Approach
### Steps
1. **Formulating the Problem as an Objective Function**:
   - The chip placement problem is represented as a mathematical expression (Objective function) that captures the constraints and goals.
   - The goal is to find the lowest values that satisfy this function, representing optimal solutions to the chip placement problem.
2. **Finding Good Solutions by Sampling**:
   - Sampling involves obtaining values from low-energy states of the objective function.
   - The optimal solution is expected to be among these samples.

### Quantum Annealing
Quantum annealing is used to find low-energy solutions for optimization problems by initializing qubits in a superposition state and gradually reducing quantum fluctuations. This allows qubits to transition from higher energy states to lower energy states, guided by a time-dependent Hamiltonian. Quantum tunneling enables qubits to overcome energy barriers, making quantum annealing effective for complex optimization problems.

### QUBO Representation

![Alt text](https://github.com/thetushargoyal/quantum-computing-project/blob/main/formulation.png)

## Code Implementation
The main functions of the implementation include:

1. **Scenario Setup**:
   Initializes the grid graph and defines points of interest (POIs), existing macro locations, and potential new macro locations.

2. **Distance Calculation**:
   Computes the Manhattan distance between points.

3. **Building the BQM**:
   Constructs the Binary Quadratic Model (BQM) for the problem.

4. **Solving the BQM**:
   Uses a sampler to solve the BQM and returns the new macro locations.

## Results
The implementation successfully identified viable new macro locations, demonstrating the potential of quantum annealing for chip placement optimization. The results include average distances to POIs and existing macros, and the distance between new macros.

## Conclusion
This project explored the application of quantum computing in chip placement optimization, a traditionally complex and computationally expensive problem. The use of quantum annealing provided promising results, and as quantum computing technology advances, it holds significant potential for further enhancing such optimization problems.
