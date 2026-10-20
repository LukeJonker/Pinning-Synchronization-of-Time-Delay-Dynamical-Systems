# ========================================================================================================
# COMPUTATIONAL AND SUPPORTING FUNCTIONS 
# ========================================================================================================
import numpy as np

# ========================================================================================================
# GENERAL VERSITILE SUPPORT FUNCTIONS
# ========================================================================================================
'''
CHECK IF PINNED
Checked if the input node i is a pinned node or not.
'''
def check_pinned(i,P):
    return (P[i][i] == 1)

# ========================================================================================================
# EVENT TRIGGERED CONTROL FUNCTIONS
# ========================================================================================================
'''
GET MINIMUM CONTROL GAIN
Extracting the minimum control gain out of the control matrix
'''
def minimum_control_gain(D):
    num_nodes = len(D[0])

    # Locating the maximum value over the matrix
    min = 0
    for i in range(0, num_nodes):
        for j in range(0, num_nodes):
            if min <= D[i][j]:
                min = D[i][j]

    # Taking the minimum over all non zero matrix entries
    for i in range(0, num_nodes):
        for j in range(0, num_nodes):
            if min >= D[i][j] and D[i][j] != 0:
                min = D[i][j]
    return(min)

'''
GET PINNED NODE STATE VECTOR
This function takes in the state and the sample state, and returns the vectors of the state for the pinned nodes, setting all non pinned node states
to zero.
'''
def reduce_to_pinned_vector(state, sample_state, P):
    num_nodes = len(P[0])

    state_temp = state.copy()
    sample_state_temp = sample_state.copy()

    # Removing the target states from the state vector
    for i in range(0,3):
        state_temp[i] = 0
        sample_state_temp[i] = 0

    # For each node indexed as node target_state = y0,y1,y2 | node_0 = {y3,y4,y5} *(p[0][0]) etc...
    for i in range(0, num_nodes):
        state_temp[3*i + 3] *= P[i][i]
        state_temp[3*i + 4] *= P[i][i]
        state_temp[3*i + 5] *= P[i][i]
        sample_state_temp[3*i + 3] *= P[i][i]
        sample_state_temp[3*i + 4] *= P[i][i]
        sample_state_temp[3*i + 5] *= P[i][i]
    return(state_temp, sample_state_temp)

'''
[EPSILON] ETC sampling error e(t) - e(t_k)
'''
def get_etc_sampling_error(state, sample_state):
    target_state = state[:3]
    target_sample_state = sample_state[:3]

    state_error = []
    state_sample_error = []
    etc_sampling_error = []

    for i in range(0, len(state), 3):
        for j in range(0,3):
            state_error.append(state[i+j] - target_state[j])
            state_sample_error.append(sample_state[i+j] - target_sample_state[j])

    for i in range(0, len(state_error)):
        etc_sampling_error.append(state_error[i] - state_sample_error[i])
        
    return etc_sampling_error

'''
[ZETA] - COMPUTE ETC ERROR MEASURE
Computes the error measure based on the current state, the current etc error and the minimum control gain in the matrix [NEEDS TARGET STATE]
'''
def compute_etc_error_measure(state, sample_state, D, P):
    num_nodes = len(D[0])
    target_state = state[:3] 
    
    control_gain = minimum_control_gain(D)
    etc_sampling_error = np.linalg.norm(get_etc_sampling_error(state, sample_state), 1)

    state_error_sum = 0
    for i in range(num_nodes):
        if P[i][i] == 1: 
            idx_x, idx_y, idx_z = 3*i + 3, 3*i + 4, 3*i + 5
            state_error_sum += abs(state[idx_x] - target_state[0])
            state_error_sum += abs(state[idx_y] - target_state[1])
            state_error_sum += abs(state[idx_z] - target_state[2])
            
    measure = control_gain * etc_sampling_error * state_error_sum
    return measure

'''
ETC CONTROL UPDARTE CHECK
This funcition checks the current state against the last sampled state and if it breaches the triggering condition it updates the sampled state to
the current state, otherwise it keeps the sampled state.
'''
def ETC_control_update_check(state, sample_state, time, P, D, alpha, beta):
 
    if compute_etc_error_measure(state, sample_state, D, P) >= alpha*np.exp(-beta*time):
        # print('Updated control')
        return(state)
    else:
        return(sample_state)

# ========================================================================================================
# COMPUTING PARAMETERS
# ========================================================================================================
'''
Computing the total sychronization error at each time in the times subsequence
'''
def compute_sync_error(data, num_nodes, vars_per_node):
    time_steps = data.shape[0]
    
    # 1. Reshape the 2D array into a 3D array: [time, node, variable]
    # This neatly groups each node's x, y, z variables together
    reshaped_data = data.reshape(time_steps, num_nodes, vars_per_node)
    
    # 2. Pre-allocate an array to hold the total error at each time step
    error_over_time = np.zeros(time_steps)
    
    # 3. Loop through every unique pair of nodes
    for j in range(num_nodes):
        for k in range(j + 1, num_nodes):
            
            # Subtract the states of node k from node j across all time steps
            diff = reshaped_data[:, j, :] - reshaped_data[:, k, :]
            
            # Calculate the Euclidean distance (magnitude of the difference)
            dist = np.linalg.norm(diff, axis=1)
            
            # Add this pair's error to the total error at each time step
            error_over_time += dist
            
    return error_over_time



'''
Takes in the network topology and the computed OSL and Lipchitz coefficents, and computes what the theorem calls the maximum premissable delay 
that will still allow for synchronization.
'''
def compute_tau_max(A, L, rho, Gamma, c, l): 
    A = np.array(A)
    Gamma = np.array(Gamma)
    A22 = A[l:, l:]       
    A21 = A[l:, :l]        

    Omega22 = c * np.kron(A22, Gamma)
    Omega21 = c * np.kron(A21, Gamma)
    omega21 = np.linalg.norm(Omega21, 2)           
    omega22 = np.linalg.norm(Omega22, 2)           

    Omega22_sym = 0.5 * (Omega22 + Omega22.T)
    mu = np.max(np.linalg.eigvalsh(Omega22_sym))   

    numerator   = -mu - rho - omega21
    denominator = L + omega22 + omega21             

    tau_max = numerator / denominator
    return tau_max

'''
Computing the minimum control gain required across all pinned nodes to satisfy one of the theorems conditions
'''
def minimal_control_strength(A, rho, Gamma, c, l):
    A = np.array(A)
    Gamma = np.array(Gamma)
    A11 = A[:l:, :l]        
    A21 = A[l:, :l]        

    Omega11 = c * np.kron(A11, Gamma)
    Omega21 = c * np.kron(A21, Gamma)
    omega21 = np.linalg.norm(Omega21, 2)           
    omega11 = np.linalg.norm(Omega11, 2)    

    return(rho + omega11 + omega21)


'''
Functions that support and compute the Lipchitz coefficent and OSL coefficient based on the specified network topology
'''
def ellipsoid_params(sigma, r, b):
    c = (r + sigma) / 2
    a_x = c * np.sqrt(b / sigma)
    a_y = c * np.sqrt(b)
    a_z = c
    return c, a_x, a_y, a_z

def boundary_point(phi, theta, sigma, r, b):
    c, a_x, a_y, a_z = ellipsoid_params(sigma, r, b)
    x = a_x * np.sin(phi) * np.cos(theta)
    y = a_y * np.sin(phi) * np.sin(theta)
    z = c + a_z * np.cos(phi)
    return np.array([x, y, z])

def jacobian(state, sigma, r, b):
    x, y, z = state
    return np.array([
        [-sigma, sigma,  0.0],
        [ r - z,  -1.0,   -x],
        [    y,      x,   -b]
    ])

def lipschitz_constant(sigma, r, b, n_grid=300):
    phis   = np.linspace(0, np.pi,   n_grid)
    thetas = np.linspace(0, 2*np.pi, n_grid)
    best = -np.inf
    for phi in phis:
        for theta in thetas:
            state = boundary_point(phi, theta, sigma, r, b)
            J = jacobian(state, sigma, r, b)
            val = np.linalg.norm(J, 2)
            if val > best:
                best = val
    return best

def osl_constant(sigma, r, b, n_grid=300):
    phis   = np.linspace(0, np.pi,   n_grid)
    thetas = np.linspace(0, 2*np.pi, n_grid)
    best = -np.inf
    for phi in phis:
        for theta in thetas:
            state = boundary_point(phi, theta, sigma, r, b)
            J = jacobian(state, sigma, r, b)
            J_sym = 0.5 * (J + J.T)
            val = np.max(np.linalg.eigvalsh(J_sym))
            if val > best:
                best = val
    return best