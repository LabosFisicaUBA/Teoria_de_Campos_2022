# Online Python compiler (interpreter) to run Python online.
# Write Python 3 code in this online editor and run it.
# Get started with interactive Python!
# Supports Python Modules: builtins, math,pandas, scipy 
# matplotlib.pyplot, numpy, operator, processing, pygal, random, 
# re, string, time, turtle, urllib.request
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy as sp
from sympy import *
from sympy.tensor.tensor import TensorIndexType, TensorHead,TensorSymmetry
from sympy.tensor.tensor import tensor_indices

import warnings
from sympy.utilities.exceptions import SymPyDeprecationWarning
warnings.filterwarnings("ignore")
# message may be omitted to filter all SymPyDeprecationWarnings
#message=r”(?s).*<regex matching the warning message>”,category=SymPyDeprecationWarning,module=r”<regex matching your module>”


#Lorentz = TensorIndexType('Lorentz', dummy_name='L')      #Define el tensor
Lorentz = TensorIndexType('Lorentz')      #Define el tensor
asym2 = TensorSymmetry.fully_symmetric(-2)                #Tensor totalmente antisimetrico
A = TensorHead('A', [Lorentz, Lorentz], asym2)
mu, nu = tensor_indices('mu nu', Lorentz)                  #Define los indices del tensor
repl = {Lorentz: diag(1, -1, -1, -1)}                     #Define la metrica del tensor

pprint(repl)

P = TensorHead('P', [Lorentz], TensorSymmetry.no_symmetry(1))
pprint(P(-mu)*P(nu))
E, px, py, pz = symbols('E p_x p_y p_z', positive=True)
repl.update({P(mu): [E, px, py, pz]})
beta = symbols('beta', real=True)
gamma,c = symbols('gamma c', positive=True)

M_01 = TensorHead('M_01', [Lorentz, Lorentz], asym2)
M_02 = TensorHead('M_02', [Lorentz, Lorentz], asym2)
M_03 = TensorHead('M_03', [Lorentz, Lorentz], asym2)

etha_1, etha_2, etha_3 = symbols("etha_1:4",real=true)


repl.update({M_01(-mu,-nu): [ [cosh(etha_1), sinh(etha_1),0,0],[sinh(etha_1), cosh(etha_1), 0,0],[0,0,1,0],[0,0,0,1]]})
repl.update({M_02(-mu,-nu): [ [cosh(etha_2), 0 ,sinh(etha_2),0],[0,1,0,0],[sinh(etha_2),0,cosh(etha_2),0],[0,0,0,1]]})
repl.update({M_03(-mu,-nu): [ [cosh(etha_3), 0 ,0,sinh(etha_3)],[0,1,0,0],[0,0,1,0],[sinh(etha_3),0,0,cosh(etha_3)]]})

pprint(M_01)

pprint(' ')
pprint(M_01(mu,nu).replace_with_arrays(repl, [mu,nu]))
pprint(' ')

pprint(' ')
pprint(M_01(-mu,nu).replace_with_arrays(repl, [-mu,nu]))
pprint(' ')

#D = diag(1,-1,-1,-1)
#g = [list(D[0,:]),list(D[1,:]),list(D[2,:]),list(D[3,:])]
g = Array(diag(1,-1,-1,-1))
G = TensorHead('G', [Lorentz, Lorentz], TensorSymmetry.no_symmetry(2))
repl.update({G(-mu,-nu): g })

#G = Lorentz.metric

sigma, rho = tensor_indices('sigma rho', Lorentz)                  #Define los indices del tensor
#T = M_01(-mu,-sigma)*G(sigma,rho)*M_01(-rho,-nu)
T = M_01(-mu,-sigma)*G(sigma,rho)*M_01(-rho,-nu)
pprint(T)
pprint(type(T))
T1 = T(-mu,-nu).replace_with_arrays(repl, [-mu,-nu])
pprint(simplify(T1))
Ty = M_02(-mu,-sigma)*G(sigma,rho)*M_02(-rho,-nu)
Tz = M_03(-mu,-sigma)*G(sigma,rho)*M_03(-rho,-nu)
Ty = Ty(-mu,-nu).replace_with_arrays(repl, [-mu,-nu])
pprint(simplify(Ty))
Tz = Tz(-mu,-nu).replace_with_arrays(repl, [-mu,-nu])
pprint(simplify(Tz))
sigma_2, sigma_3 = tensor_indices('sigma_2 sigma_3', Lorentz)                  #Define los indices del tensor
rho_2, rho_3 = tensor_indices('rho_2 rho_3', Lorentz)                  #Define los indices del tensor
#Tr = M_03(-mu,-sigma_3)*M_02(sigma_3,sigma_2)*M_01(-sigma_2,-sigma)*G(sigma,rho)*M_01(-rho,-rho_2)*M_02(rho_2,rho_3)*M_03(-rho_3,-nu)
#pprint(simplify(Tr))
Tr = M_02(-mu,sigma_2)*M_01(-sigma_2,-sigma)*G(sigma,rho)*M_01(-rho,-rho_2)*M_02(rho_2,-nu)
Tr2 = Tr(-mu,-nu).replace_with_arrays(repl, [-mu,-nu])
pprint(simplify(Tr2))
