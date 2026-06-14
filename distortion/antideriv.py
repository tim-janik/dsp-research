#!/usr/bin/env python3

import sympy as sp
from sympy.utilities.codegen import codegen

# Define the symbolic variable
x = sp.Symbol('x')

# Define your exact Padé formula
f = (x * (27 + x**2)) / (27 + 9 * x**2)

# Compute 1st and 2nd antiderivatives
F1 = sp.integrate(f, x)
F2 = sp.integrate(F1, x)

print("--- Exact 1st Antiderivative ---")
sp.pprint(F1)

print("--- Exact 2nd Antiderivative ---")
sp.pprint(F2)

print("\n--- Optimized C++ Expression ---")
# This optimizes common subexpressions automatically
(name, c_code) = codegen(("evaluateF1", F1), "C", "waveshaper")[0]
print(c_code)
