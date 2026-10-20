# ====================================================
# PLOTTING FUNCTIONS (ptf)
# ====================================================
import matplotlib.pyplot as plt
import numpy as np
import function_library as sf
'''
Synchronization error and ETC error over time - plot
'''
def plot_errors(sync_error, etc_error, plot_times, etc_envelope):
    # Create 2 subplots sharing the same X (Time) axis
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Top Plot: Synchronization Error
    ax1.plot(plot_times, sync_error, color='red', lw=1.0)
    ax1.set_ylabel("Total Sync Error")
    ax1.set_title("Network Synchronization Over Time")
    ax1.grid(True, alpha=0.2)
    
    # Bottom Plot: ETC Error (|sample_state - state|)
    ax2.plot(plot_times, etc_error, color='cyan', lw=1.0)
    ax2.plot(plot_times, etc_envelope, color='yellow', lw=1.0)
    ax2.set_xlabel("Time")
    ax2.set_ylabel("ETC Error Plot")
    ax2.set_title("Sampling Error")
    ax2.grid(True, alpha=0.2)
    
    plt.tight_layout()
    plt.show()


import matplotlib.pyplot as plt
import numpy as np

def quad_plot(sync_error, etc_error, plot_times, etc_envelope, event_times, plot_data, num_nodes):
    fig = plt.figure(figsize=(16, 10))
    
    # Top Left: Synchronization Error
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.plot(plot_times, sync_error, color='red', lw=1.0)
    ax1.set_ylabel("Total Sync Error")
    ax1.set_title("Network Synchronization Over Time")
    ax1.grid(True, alpha=0.2)
    
    # Bottom Left: Sampling Error
    ax2 = fig.add_subplot(2, 2, 3, sharex=ax1)
    ax2.plot(plot_times, etc_error, color='cyan', lw=1.0)
    ax2.plot(plot_times, etc_envelope, color='yellow', lw=1.0)
    ax2.set_xlabel("Time")
    ax2.set_ylabel("ETC Error Plot")
    ax2.set_title("Sampling Error")
    ax2.grid(True, alpha=0.2)
    
    # Top Right: Inter-event Times
    ax3 = fig.add_subplot(2, 2, 2)
    if len(event_times) > 1:
        event_times_arr = np.array(event_times)
        inter_event_times = np.diff(event_times_arr)
        ax3.scatter(event_times_arr[1:], inter_event_times, color='magenta', s=8, alpha=0.7)
    ax3.set_xlabel("Time")
    ax3.set_ylabel("Inter-Event Time")
    ax3.set_title("Triggering Intervals Over Time")
    ax3.grid(True, alpha=0.2)
    
    # Bottom Right: 2D Phase Portrait (X-Z plane)
    ax4 = fig.add_subplot(2, 2, 4)
    
    # Filter data for the first 4 seconds
    mask = plot_times <= (plot_times[0] + 0.5)
    early_data = plot_data[mask]
    
    # Plot leader FIRST with a thicker, slightly dashed/transparent line
    ax4.plot(early_data[:, 0], early_data[:, 2], 
             color='#70c1b3', lw=1.2, alpha=0.8, linestyle='--', label='Leader (Target)')
    
    # Plot followers with higher opacity and fixed labels
    colors = ['#ff6666', '#66ff66', '#6666ff', '#ffff66', '#ff66ff']
    for i in range(num_nodes):
        c = colors[i % len(colors)]
        idx_x = 3 * i + 3
        idx_z = 3 * i + 5
        ax4.plot(early_data[:, idx_x], early_data[:, idx_z], 
                 color=c, lw=0.6, alpha=0.7, label=f'Node {i}')
        
    ax4.set_title("Phase Portrait (0-0.5s)")
    ax4.set_xlabel("x")
    ax4.set_ylabel("z")
    ax4.grid(False) 
    
    # Rebuilt legend logic
    handles, labels = ax4.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax4.legend(by_label.values(), by_label.keys(), loc='upper left', fontsize='small')
    
    plt.tight_layout()
    plt.show()


'''
Pass fail boundary plot
'''
def plot_theoretical_boundary(A, L, rho, Gamma, l_scalar):
    c_vals = np.linspace(0.1, 5.0, 50)
    tau_bounds = []
    
    for c_test in c_vals:
        t_max = sf.compute_tau_max(A, L, rho, Gamma, c_test, l_scalar)
        tau_bounds.append(t_max)
        
    tau_bounds = np.array(tau_bounds)
    valid = tau_bounds > 0
    
    # --- NEW: Safety check to prevent the crash ---
    if not np.any(valid):
        print("\n--- Plotting Warning ---")
        print("No positive theoretical max delays were found in this range.")
        print("The theorem's conditions are too conservative for these network parameters (tau_max <= 0).")
        print(f"Max calculated tau: {np.max(tau_bounds):.4f}")
        return
    # ----------------------------------------------
    
    plt.figure(figsize=(8, 6))
    plt.fill_between(c_vals[valid], 0, tau_bounds[valid], 
                     color='green', alpha=0.3, label='Theoretical Guarantee (Pass)')
    plt.plot(c_vals[valid], tau_bounds[valid], color='lime', lw=2)
    
    plt.xlabel('Coupling Strength (c)')
    plt.ylabel('Maximum Permissible Delay ($\\tau_{max}$)')
    plt.title('Theoretical Feasibility Boundary')
    plt.xlim(min(c_vals), max(c_vals))
    plt.ylim(0, max(tau_bounds[valid]) * 1.2)
    plt.legend()
    plt.grid(True, alpha=0.2)
    plt.show()
import numpy as np
import matplotlib.pyplot as plt


def _gaussian_smooth(a, sigma):
    """Gaussian smoothing along axis 0 with edge padding (no lag)."""
    if sigma <= 0 or len(a) < 2:
        return a.copy()
    r = int(np.ceil(3 * sigma))
    x = np.arange(-r, r + 1)
    w = np.exp(-0.5 * (x / sigma) ** 2)
    w /= w.sum()
    padded = np.pad(a, [(r, r), (0, 0)], mode='edge')
    return np.stack([np.convolve(padded[:, d], w, mode='valid')
                     for d in range(a.shape[1])], axis=1)


def _rolling_max(a, r):
    """Max over a window of +/- r frames along axis 0."""
    if r <= 0:
        return a.copy()
    out = np.empty_like(a)
    for k in range(len(a)):
        out[k] = a[max(0, k - r):k + r + 1].max(axis=0)
    return out


def frames_plot(plot_times, plot_data, num_nodes,
                frames=240, overlap=0.8, ghost=True,
                follow=True, zoom_pad=0.15, min_span_frac=0.02,
                smooth=0.04, lock_aspect=False, fps=30):
    """
    Scrubbable phase portrait (x vs z) with overlapping time windows, a
    smooth follow camera, and autoplay.

    Keys
    ----
    space      : play / pause
    + / -      : faster / slower playback
    right / n  : next frame           left / p  : previous frame
    pagedown   : jump forward         pageup    : jump back
    home / end : first / last frame
    f          : toggle follow camera <-> full view

    stride = span/frames and window = stride/(1-overlap), so with the
    defaults each window covers 5 frames' worth of time (2% of the run).

    plot_times    : (N,) array of time stamps, ascending
    plot_data     : (N, M) array; cols 0,1,2 = leader x,y,z
                    cols 3i+3, 3i+4, 3i+5 = node i x,y,z
    frames        : number of windows
    overlap       : 0.0 = windows just touch, 0.9 = heavy smear
    ghost         : draw the full trajectory faintly underneath
    follow        : start with the follow camera on
    zoom_pad      : margin around the window's points, as a fraction of
                    the window's own extent
    min_span_frac : smallest view size, as a fraction of the full range
    smooth        : camera smoothing time as a fraction of the whole run
                    (0 = snap to every window, 0.04 = 4% of the run).
                    Independent of `frames`, so adding frames doesn't
                    make the camera twitchier.
    lock_aspect   : keep x and z at the same data scale while following
    fps           : initial playback speed
    """
    plot_times = np.asarray(plot_times, dtype=float)
    plot_data = np.asarray(plot_data, dtype=float)

    if plot_times.shape[0] != plot_data.shape[0]:
        raise ValueError("plot_times and plot_data must have the same length")
    if not 0.0 <= overlap < 1.0:
        raise ValueError("overlap must be in [0, 1)")
    if smooth < 0:
        raise ValueError("smooth must be >= 0")
    if np.any(np.diff(plot_times) < 0):
        raise ValueError("plot_times must be ascending")

    t0, t1 = plot_times[0], plot_times[-1]
    span = t1 - t0
    stride = span / frames
    width = stride / (1.0 - overlap)
    starts = t0 + stride * np.arange(frames)
    his = np.minimum(starts + width, t1)

    # Index slices instead of boolean masks: cheap even with many frames
    i0s = np.searchsorted(plot_times, starts, side='left')
    i1s = np.searchsorted(plot_times, his, side='right')

    x_cols = [0] + [3 * i + 3 for i in range(num_nodes)]
    z_cols = [2] + [3 * i + 5 for i in range(num_nodes)]

    colors = ['#ff6666', '#66ff66', '#6666ff', '#ffff66', '#ff66ff']

    # ---- Full-view limits -------------------------------------------------
    xmin, xmax = plot_data[:, x_cols].min(), plot_data[:, x_cols].max()
    zmin, zmax = plot_data[:, z_cols].min(), plot_data[:, z_cols].max()
    xpad = 0.05 * (xmax - xmin) or 1.0
    zpad = 0.05 * (zmax - zmin) or 1.0
    full_xlim = (xmin - xpad, xmax + xpad)
    full_zlim = (zmin - zpad, zmax + zpad)

    min_half_x = 0.5 * max(min_span_frac * (xmax - xmin), 1e-9)
    min_half_z = 0.5 * max(min_span_frac * (zmax - zmin), 1e-9)

    # ---- Camera path ------------------------------------------------------
    sigma = smooth * frames          # smoothing expressed in frames
    reach = int(np.ceil(2 * sigma))  # lookahead for the zoom envelope

    centers = np.full((frames, 2), np.nan)
    for k in range(frames):
        if i1s[k] > i0s[k]:
            seg = plot_data[i0s[k]:i1s[k]]
            xs, zs = seg[:, x_cols], seg[:, z_cols]
            centers[k] = [0.5 * (xs.min() + xs.max()),
                          0.5 * (zs.min() + zs.max())]

    valid = ~np.isnan(centers[:, 0])
    if not valid.any():
        centers[:] = [0.5 * sum(full_xlim), 0.5 * sum(full_zlim)]
    else:
        idx = np.arange(frames)
        for d in range(2):
            centers[:, d] = np.interp(idx, idx[valid], centers[valid, d])

    centers = _gaussian_smooth(centers, sigma)

    halves = np.empty((frames, 2))
    for k in range(frames):
        if i1s[k] > i0s[k]:
            seg = plot_data[i0s[k]:i1s[k]]
            hx = np.abs(seg[:, x_cols] - centers[k, 0]).max()
            hz = np.abs(seg[:, z_cols] - centers[k, 1]).max()
        else:
            hx, hz = min_half_x, min_half_z
        halves[k] = [max(hx * (1 + zoom_pad), min_half_x),
                     max(hz * (1 + zoom_pad), min_half_z)]

    # The camera starts widening *before* a big window arrives, then eases
    # back in afterwards, instead of popping out on the frame it's needed.
    envelope = _gaussian_smooth(_rolling_max(halves, reach), sigma)
    halves = np.maximum(envelope, halves)

    # ---- Figure -----------------------------------------------------------
    fig, ax = plt.subplots(figsize=(16, 10))

    if ghost:
        ax.plot(plot_data[:, 0], plot_data[:, 2],
                color='#70c1b3', lw=0.5, alpha=0.15, zorder=0)
        for i in range(num_nodes):
            ax.plot(plot_data[:, 3 * i + 3], plot_data[:, 3 * i + 5],
                    color=colors[i % len(colors)], lw=0.4, alpha=0.12, zorder=0)

    leader_line, = ax.plot([], [], color='#70c1b3', lw=1.2, alpha=0.9,
                           linestyle='--', label='Leader (Target)')
    node_lines = [
        ax.plot([], [], color=colors[i % len(colors)], lw=0.8, alpha=0.9,
                label=f'Node {i}')[0]
        for i in range(num_nodes)
    ]
    leader_head, = ax.plot([], [], 'o', color='#70c1b3', ms=6)
    node_heads = [
        ax.plot([], [], 'o', color=colors[i % len(colors)], ms=5)[0]
        for i in range(num_nodes)
    ]

    ax.set_xlim(*full_xlim)
    ax.set_ylim(*full_zlim)
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    ax.grid(False)
    ax.legend(loc='upper left', fontsize='small')

    jump = max(10, frames // 20)
    state = {'k': 0, 'follow': follow, 'playing': False, 'fps': fps}

    def apply_camera(k):
        if not state['follow']:
            ax.set_xlim(*full_xlim)
            ax.set_ylim(*full_zlim)
            return
        cx, cz = centers[k]
        hx, hz = halves[k]
        if lock_aspect:
            bbox = ax.get_window_extent()
            ratio = bbox.height / bbox.width if bbox.width else 1.0
            if hz / hx < ratio:
                hz = hx * ratio
            else:
                hx = hz / ratio
        ax.set_xlim(cx - hx, cx + hx)
        ax.set_ylim(cz - hz, cz + hz)

    def draw_frame(k):
        seg = plot_data[i0s[k]:i1s[k]]

        leader_line.set_data(seg[:, 0], seg[:, 2])
        for i, line in enumerate(node_lines):
            line.set_data(seg[:, 3 * i + 3], seg[:, 3 * i + 5])

        if len(seg):
            leader_head.set_data([seg[-1, 0]], [seg[-1, 2]])
            for i, h in enumerate(node_heads):
                h.set_data([seg[-1, 3 * i + 3]], [seg[-1, 3 * i + 5]])
        else:
            leader_head.set_data([], [])
            for h in node_heads:
                h.set_data([], [])

        apply_camera(k)

        mode = "follow" if state['follow'] else "full"
        play = f"playing {state['fps']:.0f} fps" if state['playing'] else "paused"
        ax.set_title(f"Phase Portrait   frame {k + 1}/{frames}"
                     f"   t = {starts[k]:.3f}-{his[k]:.3f} s"
                     f"   [{mode}, {play}]")
        fig.canvas.draw_idle()

    timer = fig.canvas.new_timer(interval=int(1000 / fps))

    def tick():
        state['k'] = (state['k'] + 1) % frames
        draw_frame(state['k'])

    timer.add_callback(tick)

    def on_key(event):
        k = state['k']
        if event.key == ' ':
            state['playing'] = not state['playing']
            (timer.start if state['playing'] else timer.stop)()
        elif event.key in ('+', '='):
            state['fps'] = min(state['fps'] * 1.5, 120)
            timer.interval = int(1000 / state['fps'])
        elif event.key == '-':
            state['fps'] = max(state['fps'] / 1.5, 2)
            timer.interval = int(1000 / state['fps'])
        elif event.key in ('right', 'n'):
            state['k'] = (k + 1) % frames
        elif event.key in ('left', 'p'):
            state['k'] = (k - 1) % frames
        elif event.key == 'pageup':
            state['k'] = (k - jump) % frames
        elif event.key == 'pagedown':
            state['k'] = (k + jump) % frames
        elif event.key == 'home':
            state['k'] = 0
        elif event.key == 'end':
            state['k'] = frames - 1
        elif event.key == 'f':
            state['follow'] = not state['follow']
        else:
            return
        draw_frame(state['k'])

    our_keys = {'left', 'right', 'n', 'p', 'pageup', 'pagedown',
                'home', 'end', 'f', ' ', '+', '=', '-'}
    for key in [k for k in plt.rcParams if k.startswith('keymap.')]:
        plt.rcParams[key] = [b for b in plt.rcParams[key] if b not in our_keys]

    fig.canvas.mpl_connect('key_press_event', on_key)
    fig._frames_timer = timer  # keep a reference so it isn't garbage-collected
    fig.tight_layout()
    draw_frame(0)
    plt.show()

    return fig, ax