import numpy as np
import matplotlib.pyplot as plt
import initialize

def stepper(r_arr, v_arr, ball_pos, R_ball, step_length):
    new_r_arr, new_v_arr = np.zeros((r_arr.shape[0], 2))
    for i in range(r_arr.shape[0]):
       r_new, v_new, dv = evolve_position(r_arr[i], v_arr[i], ball_pos, R_ball, 1, step_length)
       new_r_arr[i,0], new_r_arr[i.1] += r_new[0], r_new[1]
       new_v_arr[i,0]. new_v_arr[i,1] += v_new[0], v_new[1]
    return new_r_arr, new_v_arr, dv

def run_simulation(r_arr_init, v_arr_init, ball_pos_init, R_ball, step_length, sim_time):
    N_steps = sim_time/step_length #N steps
    dv_list = []
    i=0
    while i <= N_steps:
        r_new, v_new, dv = stepper(r_arr_init, v_arr_init, ball_pos, R_ball, step_length)
        

