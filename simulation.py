import params
import kinesin
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def compute_speed_distribution(n_reps, make_plot=True, atp_conc=1e-3):
    speeds = []
    plot_choice = np.random.choice(range(n_reps))

    tarr = 0
    poslarr = []
    posrarr = []

    for i in range(n_reps):
        print(f'Iteration {i}')
        kin = kinesin.Kinesin(params.N_iter, 0, params.dx, 'D', 'E', params.D, params.dt, 
                              params.kb, params.T, params.s, atp_conc, params.k_atp_bind, 
                              params.k_atp_hyd, params.k_adp_release)
        for _ in range(params.N_iter):
            kin.update()
    
        [tarr_i, poslarr_i, posrarr_i, histl_i, histr_i] = kin.data()
        total_time = (params.N_iter * params.dt * 1e-6)
        speed = ((poslarr_i[-1] + posrarr_i[-1]) - (poslarr_i[0] + posrarr_i[0])) / (2 * total_time) # nm s^-1
        speeds.append(speed)

        if i == plot_choice and make_plot:
            print(f"Plot choice is {i}")
            tarr = tarr_i
            poslarr = poslarr_i
            posrarr = posrarr_i

    # Calculate average speed.
    speeds = np.array(speeds)
    l = np.arange(len(speeds))
    mean_speed = np.mean(speeds)
    speed_sd = np.std(speeds)

    if make_plot:
        # Plot trajectory
        plt.figure()
        plt.plot(tarr, poslarr, lw=1, color='red', label='Left')
        plt.plot(tarr, posrarr, lw=1, color='black', label='Right')
        plt.xlabel('Time, ms')
        plt.ylabel('Position, nm')
        speed = (poslarr[-1] + posrarr[-1]) - (poslarr[0] + posrarr[0]) / (2 * total_time) # nm s^-1
        plt.title(f'Speed = {speed:.2f} nm/s')
        plt.tight_layout()
        plt.savefig(f'traj_t{int(atp_conc*1000)}_s{int(params.s*100)}.png', dpi=300, bbox_inches='tight')

        # Plot separation
        plt.figure()
        plt.plot(tarr, np.abs(poslarr-posrarr), lw=1)
        plt.xlabel('Time, ms')
        plt.ylabel('Separation, nm')
        plt.ylim(-1, 15)
        plt.tight_layout()
        plt.savefig(f'sep_t{int(atp_conc*1000)}_s{int(params.s*100)}.png', dpi=300, bbox_inches='tight')

        # Plot speed distribution
        plt.figure()
        plt.hist(speeds, bins='auto')
        plt.xlabel('Speeds')
        plt.ylabel('Frequency')
        plt.title(f'Speed Distribution \nMean speed = {mean_speed:.2f} nm/s, Standard Deviation = {speed_sd:.2f} nm/s')
        plt.tight_layout()
        plt.savefig(f'dist_t{int(atp_conc*1000)}_s{int(params.s*1000)}.png', dpi=300, bbox_inches='tight')

    return (mean_speed, speed_sd)

def plot_speed_vs_conc():
    atp_concs = np.array([1e-3, 1e-4, 1e-5, 1e-6])

    mean_speeds = []
    std_devs = []
    for atp_conc in atp_concs:
        (mean_spd, sd) = compute_speed_distribution(3, False)
        mean_speeds.append(mean_spd)
        std_devs.append(sd)

    plt.figure()
    plt.errorbar(atp_concs*1000, mean_speeds, yerr=std_devs, fmt='o', capsize=4)
    plt.xlabel('ATP Concentration, mM')
    plt.ylabel('Kinesin motion speed, nm/s')
    plt.tight_layout()
    plt.savefig(f'velerb.png', dpi=300, bbox_inches='tight')

compute_speed_distribution(100)
plot_speed_vs_conc()
