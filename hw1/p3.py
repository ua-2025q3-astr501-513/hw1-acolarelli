#!/usr/bin/env python3
#
# Please look for "TODO" in the comments, which indicate where you
# need to write your code.
#
# Part 3: Implement a Numerically Stable Quadratic Equation Solver (1 point)
#
# * Objective:
#   Implement a numerically stable quadratic equation solver that does
#   not catastrophic cancellation.
# * Details:
#   The description of the problem and the solution template can be
#   found in `hw1/p3.py`.
#
# From lecture `01w`, we learned about catastrophic cancellation---the
# significant loss of precision that occurs when subtracting two
# nearly equal numbers.
# This problem actually appeared in CK's research!
# While solving for the initial conditions of (unstable) spherical
# photon orbits around a black hole for the convergence test of
# [GRay2](https://ui.adsabs.harvard.edu/abs/2018ApJ...867...59C),
# catastrophic cancellation introduced errors so severe that photons
# would not remain on their spherical orbits for an radian.
# CK suspected a bug in the integrator and spent an entire month
# debugging the wrong part of the code.
# At the end, he realized the real problem.
# The standard quadratic formula we all learn in high school was
# simply not accurate enough for reliable numerical computation.
#
# Here, let's implement a numerically stable quadratic equation solver
# to overcome catastrophic cancellation.
# Please make sure that you take care of all the special cases.

def quadratic(a, b, c):
    """Numerically stable quadratic equation solver

    The standard quadratic formula

        x = (-b +- sqrt(b^2 - 4ac)) / (2a)

    is algebraically correct but can suffer from *catastrophic
    cancellation* when b^2 >> 4ac and the sign of b matches the
    chosen +-.
    In that case, subtracting two nearly equal numbers causes a large
    loss of precision.

    A more stable alternative is obtained by multiplying top and
    bottom by the conjugate, leading to two equivalent forms.
    To avoid cancellation, choose the version that keeps the
    subtraction well-separated:

        x1 = (-b - sign(b) * sqrt(b^2 - 4ac)) / (2a)
        x2 = (c / a) / x1

    This way, at least one root is always computed stably.

    Args:
        a, b, c: coefficients for the quadratic equation
                 a x^2 + b x + c = 0.

    Returns:
        x1, x2: the two roots of the quadratic equation.
                If there are two real roots, x1 < x2.
                If there is only one real root, x2 == None.
                If there is no real root, x1 == x2 == None.
    """

    # implement the stable quadratic equation solver here (I wrote this before I realized we could use numpy so it's not the most efficient, but it works)
    delta = b**2 - 4*a*c

    # First check the degree of the polynomial
    if a is None or a==0:
        if b is None or b==0:
            x1, x2 = None, None     # If of order zero (constant), there are no roots
        else:
            x1, x2 = (-c/b), None   # If of order one (linear), there is one root, which is just the y-intercept
    
    # Now check the discriminant
    elif delta>0:
        discrim = delta**0.5            # If discriminant is positive, we have two real roots
        if abs((4*a*c)/(b**2))<1e-4:    # Check if b^2 >> 4ac to avoid catastrophic cancellation
            if b>0:
                x1 = (-b - discrim)/(2*a)
                x2 = (c/a)/x1
            elif b<0:
                x2 = (-b + discrim)/(2*a)
                x1 = (c/a)/x2
        else:                           # Now the "normal" quadratic equation solution
            x1 = (-b - discrim)/(2*a)
            x2 (-b + discrim)/(2*a)

    elif delta==0:                      # If discriminant is zero, there is only one real root
        x1 = -b/(2*a)
        x2 = None

    elif delta<0:                       # If discriminant is negative, there are no real roots
        x1, x2 = None, None
    return x1, x2
