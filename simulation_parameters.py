import plotting_functions as ptf
import function_library as sf
import numpy as np

# ========================================================================================================
# SIMULATION PARAMETERS
# ========================================================================================================
sim_length = 1
time_step = 0.0001
relative_tolerance = 1e-6
absolute_tolerance = 1e-8

# ========================================================================================================
# USER SET PARAMETERS - NETWORK TOPOLOGY AND DYNAMICS
# ========================================================================================================
n = 5                # Number of nodes
N = 3                # Number of states in each node
P = [
    [1,0,0,0,0],
    [0,1,0,0,0],
    [0,0,1,0,0],     # Pinned node matrix
    [0,0,0,0,0],
    [0,0,0,0,0]
]
Gamma = [
    [1,0,0],
    [0,1,0],         # Inner Coupling Matrix
    [0,0,1]
]
A = [
    [-4,1,1,1,1],
    [1,-4,1,1,1],
    [1,1,-4,1,1],    # Outer Coupling Matrix | Rows must sum to 0 | Matrix must be symmetric
    [1,1,1,-4,0],
    [1,1,1,0,-3]
]
c = 40               # Coupling strength

sigma, r, b = 10.0, 28.0, 8.0/3.0   # Lorenz coefficients

alpha = 5000     # ETC scalar
beta = 0.4       # ETC decay rate

initial_conditions = np.random.uniform(-5, 5, (n*N) + 3)
# ========================================================================================================
# RESTRICTED PARAMETERS 
# ========================================================================================================

tau = 0.008         # Delay 

D = [
    [300,0,0,0,0],
    [0,300,0,0,0],
    [0,0,300,0,0],     # Control gain
    [0,0,0,0,0],
    [0,0,0,0,0]
]
# ========================================================================================================
# Computed Parameters
# ========================================================================================================
# Number of pinned nodes
l = 0
for i in range(0,len(P[0])):
    for j in range(0,len(P[0])):
        l += P[i][j]
# Lipchitz 
L = sf.lipschitz_constant(sigma, r, b, n_grid=300)
# OSL
rho = sf.osl_constant(sigma, r, b, n_grid=300)
# Maximum Delay
tau_max = sf.compute_tau_max(A, L, rho, Gamma, c, l)
# Minimum Control
minimum_control = sf.minimal_control_strength(A, rho, Gamma, c, l)

# ========================================================================================================
# PRINT THEOREM REQUIREMENTS
# ========================================================================================================
print('='*10)
print("Lipchitz:", L)
print('OSL:', rho)
print('Maximum delay:', tau_max)
print('Minimum control strenght:', minimum_control)
print('='*10)

