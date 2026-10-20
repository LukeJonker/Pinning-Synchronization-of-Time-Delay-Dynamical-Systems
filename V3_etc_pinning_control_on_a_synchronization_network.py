# ========================================================================================================
# MAIN
# ========================================================================================================
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from jitcdde import jitcdde, y, t
import symengine

import function_library as sf
import plotting_functions as ptf
import simulation_parameters as p

plt.style.use('dark_background')

'''
The user plugs in a node and a component of that node; the pair (node, variable), and the function returns the position
of the requested state in the network wide state vector.
'''
def idx(node, var):
    return (node)*p.N + var + 3

'''
Building out the network dynamical system formula f, in the format desired by Jitcdde. 
'''
zoh_states = [symengine.Symbol(f"u_{k}") for k in range(p.n * p.N + 3)]

# Target state at y0,y1,y2 in the state vector
f = [
    p.sigma * (y(1,t) - y(0,t)),
    y(0,t) * (p.r - y(2,t))- y(1,t),
    y(0,t) * y(1,t) - p.b * y(2,t)
]

for i in range(p.n):
    sums = [0.0] *p.N
    
    for v in range(p.N): 
        coupling_term = 0.0
        for j in range(p.n):
            inner_prod = (p.Gamma[v][0] * y(idx(j, 0), t - p.tau) +
                          p.Gamma[v][1] * y(idx(j, 1), t - p.tau) +
                          p.Gamma[v][2] * y(idx(j, 2), t - p.tau))
            coupling_term += p.A[i][j] * inner_prod
        
        sums[v] = p.c * coupling_term
        
        if sf.check_pinned(i,p.P): 
            sums[v] -= p.D[i][i] * (zoh_states[idx(i, v)] - y(v,t)) # [-d*(sample_state - target_state)]
            
    # Local dynamics
    xi = y(idx(i, 0), t)
    yi = y(idx(i, 1), t)
    zi = y(idx(i, 2), t)
    
    dx = p.sigma * (yi - xi) + sums[0]
    dy = xi * (p.r - zi) - yi + sums[1]
    dz = xi * yi - p.b * zi + sums[2]
    
    f.extend([dx, dy, dz])


'''
Jitcdde Configuration
'''
DDE = jitcdde(f, control_pars=zoh_states)
DDE.set_integration_parameters(rtol=p.relative_tolerance, atol=p.absolute_tolerance)
np.random.seed(42) 
sample_state = p.initial_conditions
DDE.constant_past(p.initial_conditions)
DDE.set_parameters(*sample_state)
DDE.step_on_discontinuities()


'''
Getting Results
'''
data = []
etc_error = []
etc_envelope = []
event_times = []

times = np.arange(DDE.t, DDE.t + p.sim_length, p.time_step)
for i, time in enumerate(times):

    state = DDE.integrate(time)
    data.append(state)

    updated_sample_state = sf.ETC_control_update_check(state, sample_state, time, p.P, p.D, p.alpha, p.beta)
    


    etc_error.append(sf.compute_etc_error_measure(state, sample_state, p.D, p.P))
    etc_envelope.append(p.alpha*np.exp(-p.beta*time))
    
    if not np.array_equal(updated_sample_state, sample_state):
        sample_state = updated_sample_state
        DDE.set_parameters(*sample_state)
        event_times.append(time)


    if i % 100 == 0:  
        print(f"Simulated up to t = {time:.2f}")

data = np.vstack(data)

plot_data = data
plot_times = times

sync_error = sf.compute_sync_error(plot_data[:, 3:], p.n, p.N)
ptf.frames_plot(plot_times, plot_data, p.n)
ptf.quad_plot(sync_error, etc_error, plot_times, etc_envelope, event_times, plot_data, p.n)


