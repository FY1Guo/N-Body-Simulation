import numpy as np
import matplotlib.pyplot as plt
from initialize import *
import gas

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
       new_v_arr[i,0] += v_new[1]
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

    while i <= N_steps:
        r_new, v_new, dv_list = stepper(r_arr_init, v_arr_init, ball_pos_init, ball_vel_init, R_ball, step_length) #update gas particles 
        net_impulse_x = np.sum(M_gas * dv_list[0:,0])
        net_impulse_y = np.sum(M_gas * dv_list[0:,1])
        net_impulse = np.array([net_impulse_x, net_impulse_y])
        ball_r_new, ball_v_new, ball_accel = gas.update_projectile(ball_pos_init, ball_vel_init, net_impulse, M_ball, 1, step_length)
        print(ball_v_new) 
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
        ball_vel_init = ball_vel_init
        i += 1

    return average_particle_energy_hist, ball_energy_hist 

print(run_simulation(r_arr_init, v_arr_init, init_ball_pos, init_ball_vel, 0.1, 1, 20))


