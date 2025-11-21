import numpy as np
import matplotlib.pyplot as plt
from initialize import *
import gas
import plottingClass as plots

M_ball = 10
M_gas = 0.1
N_particles = int(1e3)
init_ball_pos = np.array([0.5,0.5])
init_ball_vel = np.array([0,-0.1])

def stepper(r_arr, v_arr, ball_pos, ball_vel, R_ball, step_length):
    new_r_arr = np.zeros((r_arr.shape[0], 2))
    new_v_arr = np.zeros((r_arr.shape[0], 2))
    dv_arr = np.zeros((r_arr.shape[0], 2)) 
    for i in range(r_arr.shape[0]):
       r_new, v_new, dv = gas.evolve_position(r_arr[i], v_arr[i], ball_pos, ball_vel, R_ball, 1, step_length, M_gas, M_ball)
       new_r_arr[i,0] += r_new[0]
       new_r_arr[i,1] += r_new[1]
       new_v_arr[i,0] += v_new[0]
       new_v_arr[i,1] += v_new[1]
       dv_arr[i,0] += dv[0]
       dv_arr[i,1] += dv[1]
    return new_r_arr, new_v_arr, dv_arr

r_arr_init = initialize.make_particles_pos(N_particles)
v_arr_init = initialize.make_particles_vel(N_particles)

def run_simulation(r_arr_init, v_arr_init, ball_pos_init, ball_vel_init, R_ball, step_length, sim_time):
    """
    Runs the n-body simulation
    Returns:
        average_particle_energy_hist : (N_steps+1,) array
        ball_energy_hist             : (N_steps+1,) array
        ball_pos_hist                : (N_steps+1, 2) array
        ball_vel_hist                : (N_steps+1, 2) array
        force_hist                   : (N_steps+1, 2) array  (drag force on the ball)
    """
    # number of steps
    N_steps = int(sim_time / step_length)

    # history arrays
    average_particle_energy_hist = np.zeros(N_steps + 1)
    ball_energy_hist = np.zeros(N_steps + 1)
    ball_pos_hist = np.zeros((N_steps + 1, 2))
    ball_vel_hist = np.zeros((N_steps + 1, 2))
    force_hist = np.zeros((N_steps + 1, 2))

    # local copies 
    r_arr = r_arr_init.copy()
    v_arr = v_arr_init.copy()
    ball_pos = ball_pos_init.copy()
    ball_vel = ball_vel_init.copy()

    # --- initial state at step 0 ---
    ball_pos_hist[0] = ball_pos
    ball_vel_hist[0] = ball_vel

    speeds2_init = v_arr[:, 0]**2 + v_arr[:, 1]**2
    average_particle_energy_hist[0] = 0.5 * M_gas * np.mean(speeds2_init)
    ball_energy_hist[0] = 0.5 * M_ball * np.dot(ball_vel, ball_vel)
    force_hist[0] = np.array([0.0, 0.0])

    # --- time stepping ---
    for step in range(1, N_steps + 1):
        # evolves all gas particles by one step and get total dv per particle
        r_new, v_new, dv_arr = stepper(r_arr, v_arr, ball_pos, ball_vel, R_ball, step_length)

        # net change in gas momentum
        delta_p_gas = M_gas * np.sum(dv_arr, axis=0)

        # by momentum conservation, ball gets opposite impulse
        net_impulse_ball = -delta_p_gas

        # updates projectile (box size = 1.0)
        ball_pos_new, ball_vel_new, ball_accel = gas.update_projectile(
            ball_pos, ball_vel, net_impulse_ball, M_ball, 1.0, step_length
        )

        # records histories
        ball_pos_hist[step] = ball_pos_new
        ball_vel_hist[step] = ball_vel_new
        force_hist[step] = M_ball * ball_accel  # F = m a

        # energies
        speeds2 = v_new[:, 0]**2 + v_new[:, 1]**2
        average_particle_energy_hist[step] = 0.5 * M_gas * np.mean(speeds2)
        ball_energy_hist[step] = 0.5 * M_ball * np.dot(ball_vel_new, ball_vel_new)

        # advances state for next step
        r_arr = r_new
        v_arr = v_new
        ball_pos = ball_pos_new
        ball_vel = ball_vel_new

    return average_particle_energy_hist, ball_energy_hist, ball_pos_hist, ball_vel_hist, force_hist


if __name__ == "__main__":
    # --- simulation parameters ---
    R_ball = 0.1
    dt = 0.01          # time step
    sim_time = 20.0    # total simulation time
    
    # --- initial gas particle state ---
    r_arr_init = initialize.make_particles_pos(N_particles)
    v_arr_init = initialize.make_particles_vel(N_particles)
    
    # --- run the simulation ---
    avg_E_hist, ball_E_hist, ball_pos_hist, ball_vel_hist, force_hist = run_simulation(
        r_arr_init,
        v_arr_init,
        init_ball_pos,
        init_ball_vel,
        R_ball,
        dt,
        sim_time,
    )
    
    # --- builds time axis for plotting ---
    time = np.linspace(0, sim_time, len(avg_E_hist))
    
    # --- energies vs time ---
    total_E = avg_E_hist + ball_E_hist

    plt.figure()
    plt.plot(time, avg_E_hist, label="avg particle KE")
    plt.plot(time, ball_E_hist, label="ball KE")
    plt.plot(time, total_E, label="total KE")
    plt.xlabel("time")
    plt.ylabel("energy")
    plt.title("Kinetic energies vs time")
    plt.legend()
    plt.grid(True)

   # ball speed over time
    ball_speed = np.linalg.norm(ball_vel_hist, axis=1)

    plt.figure()
    plt.plot(time, ball_speed)
    plt.xlabel("time")
    plt.ylabel("ball speed")
    plt.title("Ball speed vs time")
    plt.grid(True)

    # drag force magnitude vs time
    drag_mag = np.linalg.norm(force_hist, axis=1)

    plt.figure()
    plt.plot(time, drag_mag)
    plt.xlabel("time")
    plt.ylabel("drag force magnitude")
    plt.title("Drag force vs time")
    plt.grid(True)

    plt.show()
