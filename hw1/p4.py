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
        stiffness = np.zeros((N,N))
        mass_matrix = np.zeros((N,N))

        for i in range(N):
            for j in range(N):
                if i==j:
                    stiffness[i][j] = 2. * k
                    mass_matrix[i][j] = m
                else:
                    stiffness[i][j] = -1. * k
                    mass_matrix[i][j] = 0.
        
        self.K = stiffness
        self.X0 = X0

        # TODO: Construct the stiffness matrix K
        # TODO: Solve the eigenvalue problem for K to find normal modes
        # TODO: Store angular frequencies and eigenvectors
        # TODO: Compute initial modal amplitudes M0 (normal mode decomposition)

        mass_inv = np.linalg.inv(mass_matrix)
        evals, modes = np.linalg.eig(stiffness @ mass_inv)
        amps = modes.T @ (mass_matrix @ X0)


        self.Omega = np.sqrt(evals)
        self.mass = mass_matrix
        self.inv_mass = mass_inv
        self.Modes = modes
        self.M0 = amps

    def __call__(self, t):
        """Calculate the displacements of the oscillators at time t.

        Args:
            t (float): time at which to compute the displacements.

        Returns:
            np.ndarray: displacements of the oscillators at time t.

        """
        # TODO: Reconstruct the displacements from normal modes
        x_double_dot = np.zeros(self.N)
        for i in range(self.N):
            left = self.X0[i-1] if i > 0 else 0.0      # fixed boundary at left
            right = self.X0[i+1] if i < self.N-1 else 0.0   # fixed boundary at right
            x_double_dot[i] = (self.K @ self.inv_mass)[i] * (left - 2*self.X0[i] + right)

        soln = solve_ivp(x_double_dot, (0.,t), self.X0, t_eval = t)
        return soln.y


if __name__ == "__main__":

    # Initialize the coupled oscillator system with default parameters
    co = CoupledOscillators()

    # Print displacements of the oscillators at each time step
    print("Time(s)  Displacements")
    print("----------------------")
    for t in np.linspace(0, 10, num=101):
        X = co(t)             # compute displacements at time t
        print(f"{t:.2f}", X)  # print values for reference
