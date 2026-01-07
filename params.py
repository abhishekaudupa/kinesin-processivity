import matplotlib as mpl

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['mathtext.fontset'] = 'cm'

dx              = 8             # nm
D               = 0.1           # nm^2 us^-1
T               = 300           # K
kb              = 1.380649e-8   # nm^2 ug us^-2 K^-1
V0              = kb*T          # nm^2 ug us^-2
s               = 0.31          # slope of potential landscape

# Rates
k_atp_bind      = 4e-1          # M^-1 us^-1
k_adp_release   = 200e-6        # us^-1
k_atp_hyd       = 600e-6        # us^-1

# Simulation parameters.
dt              = 1e-1          # us
N_iter          = int(3e6)      # 300 ms
