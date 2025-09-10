#!/usr/bin/env python3
#
# Please look for "TODO" in the comments, which indicate where you
# need to write your code.
#
# Part 4: Solve the Coupled Simple Harmonic Oscillator Problem (1 point)
#
# * Objective:
#   Take the coupled harmonic oscillator problem we solved in class
#   and rewrite it using a well-structured Python class.
# * Details:
#   The description of the problem and the solution template can be
#   found in `hw1/p4.py`.
#
# From lecture `02w`, we solve systems of coupled harmonic oscillators
# semi-analytically by numerically solving eigenvalue problems.
# However, the code structure was not very clean, making the code hard
# to reuse.
# Although numerical analysis in general does not require
# object-oriented programming, it is sometime useful to package
# stateful caluation into classes.
# For this assignment, we will provide a template class.
# Your responsibility to implement the methods in the class.


import numpy as np
import scipy.linalg as la
from scipy.integrate import solve_ivp


class CoupledOscillators:
    """A class to model a system of coupled harmonic oscillators.

    Attributes:
        Omega (np.ndarray): array of angular frequencies of the normal modes.
        V     (np.ndarray): matrix of eigenvectors representing normal modes.
        M0    (np.ndarray): initial amplitudes of the normal modes.

    """

    def __init__(self, X0=[-0.5, 0, 0.5], m=1.0, k=1.0):
        """Initialize the coupled harmonic oscillator system.

        Args:
            X0 (list or np.ndarray): initial displacements of the oscillators.
            m  (float):              mass of each oscillator (assumed identical for all oscillators).
            k  (float):              spring constant (assumed identical for all springs).

        """
        N = len(X0)
        self.N = N
        self.mass = m
        self.X0 = X0
        self.k = k

        # Construct the stiffness matrix K
        stiffness = np.zeros((N,N))
        
        for i in range(N):
                stiffness[i][i] = 2. * k            # Diagonal elements equal to 2k
                if i>0:                             
                    stiffness[i-1][i] = -1. * k     # Off-diagonal elements equal to -k    
                if i< N - 1:
                    stiffness[i+1][i] = -1. * k

        # Solve the eigenvalue problem for K to find normal modes
        evals, modes = np.linalg.eig(stiffness/m)

        # Compute initial modal amplitudes M0 (normal mode decomposition)
        amps = np.linalg.solve(modes, X0)

        # Store angular frequencies and eigenvectors
        self.K = stiffness
        self.Omega = np.sqrt(np.abs(evals))
        self.Modes = modes
        self.M0 = amps

    def __call__(self, t):
        """Calculate the displacements of the oscillators at time t.

        Args:
            t (float): time at which to compute the displacements.

        Returns:
            np.ndarray: displacements of the oscillators at time t.

        """
        # Reconstruct the displacements from normal modes
        Xt = self.M0 * np.cos(self.Omega *t)
        return self.Modes @ Xt



if __name__ == "__main__":

    # Initialize the coupled oscillator system with default parameters
    co = CoupledOscillators()

    # Print displacements of the oscillators at each time step
    print("Time(s)  Displacements")
    print("----------------------")
    for t in np.linspace(0, 10, num=101):
        X = co(t)             # compute displacements at time t
        print(f"{t:.2f}", X)  # print values for reference
