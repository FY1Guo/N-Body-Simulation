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

def run_simulation(r_arr_init, v_arr_init, ball_pos_init, ball_vel_init,  R_ball, step_length, sim_time):
    N_steps = int(sim_time/step_length) #N steps
    i = 0
    ball_accel_hist = np.zeros((N_steps + 1, 2))
    ball_pos_hist = np.zeros((N_steps + 1, 2))
    ball_vel_hist = np.zeros((N_steps + 1, 2))
    average_particle_energy_hist = np.zeros(N_steps + 1)
    ball_energy_hist = np.zeros(N_steps + 1)
    force_hist = np.zeros((N_steps + 1, 2))

    while i <= N_steps:
        r_new, v_new, dv_list = stepper(r_arr_init, v_arr_init, ball_pos_init, ball_vel_init, R_ball, step_length) #update gas particles 
        net_impulse_x = np.sum(M_gas * dv_list[0:,0])
        net_impulse_y = np.sum(M_gas * dv_list[0:,1])
        
        #we can estimate the force on the ball using the impulse: delta p / delta t 
        force_hist[i,0] += net_impulse_x / step_length
        force_hist[i,1] += net_impulse_y / step_length

        net_impulse = np.array([net_impulse_x, net_impulse_y])
        ball_r_new, ball_v_new, ball_accel = gas.update_projectile(ball_pos_init, ball_vel_init, net_impulse, M_ball, 1, step_length)
        # print(ball_v_new) 
        energy = []
        for k in range(N_particles):
            v_mag = v_new[k][0] ** 2 + v_new[k][1] ** 2
            e_mag = 1/2 * M_gas * v_mag
            energy.append(e_mag)
        average_particle_energy_hist[i] += np.mean(energy)
        ball_energy_hist[i] += 1/2 * M_ball * (ball_v_new[0] ** 2 + ball_v_new[1] ** 2)
        
        ball_accel_hist[i,0] += ball_accel[0]
        ball_accel_hist[i,1] += ball_accel[1]

        r_arr_init = r_new #update frame
        v_arr_init = v_new
        ball_pos_init = ball_r_new
        ball_vel_init = ball_v_new
        i += 1

    return avg_E_hist, ball_E_hist, ball_pos_hist, ball_vel_hist, force_hist

if __name__ == "__main__":
    # --- simulation parameters ---
    R_ball = 0.1
    dt = 0.01          # time step
    sim_time = 20.0    # total simulation time
    
    # --- initial gas particle state ---
    r_arr_init = initialize.make_particles_pos(N_particles)
    v_arr_init = initialize.make_particles_vel(N_particles)
    
    # --- runs the simulation ---
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
    # run_simulation fills N_steps+1 entries where N_steps = sim_time/dt
    time = np.linspace(0, sim_time, len(avg_E_hist))
    
    # --- plots ---
    plotter = plots(time)
    
    # energies vs time
    total_E = avg_E_hist + ball_E_hist
    
    plotter.total_energy_vs_time(total_E)
    plotter.kinetic_energy_vs_time(ball_E_hist, scope="ball")
    plotter.kinetic_energy_vs_time(avg_E_hist, scope="particles")
    
    # force vs velocity and force vs time
    plotter.force_vs_velocity(force_hist, ball_vel_hist)
    plotter.force_vs_time(force_hist)
    
    # velocity vs time
    plotter.velocity_vs_time(ball_vel_hist)
    
    # trajectories
    plotter.trajectory_ball(ball_pos_hist)
