import numpy as np

epsilon = 1e-10

def get_ball_collisions(r0, v, ball_pos, R_ball, t_max):
    #Check for collisions with the ball which occur in time less than t_max
    return #position, velocity, time when collision occurs. Return [] if no collision occurs    

def get_wall_collisions(r0, v, box_size, t_max):
    x0, y0 = r0
    vx, vy = v
    L = boxsize
    
    # Set x0 and y0 in bounds
    x0 = min(max(x0, epsilon), boxsize-epsilon)
    y0 = min(max(y0, epsilon), boxsize-epsilon)
    
    if vx==0:
        t_right = np.inf
        t_left = np.inf
    else:
        t_right = (L-x0)/vx
        t_left = -x0/vx
    if vy==0:
        t_top = np.inf
        t_bot = np.inf
    else:
        t_top = (L-y0)/vy
        t_bot = -y0/vy
    
    new_coordinates = np.array([(L-epsilon, y0+vy*t_right), (epsilon, y0+vy*t_left), (x0+vx*t_top, L-epsilon), (x0+vx*t_bot,epsilon)])
    new_velocity = np.array([(-vx, vy), (-vx,vy), (vx, -vy), (vx,-vy)])
    
    times = np.array([t_right, t_left, t_top, t_bot])
    valid_time_indices = np.where((times>0) & (times<=t_max))
    valid_times =  times[valid_time_indices]
    if len(valid_times)==0:
        return []
    time_idx = np.argmin(valid_times)
    time = valid_times[time_idx]
    collision_idx = valid_time_indices[0][time_idx]
    
    return new_coordinates[collision_idx,:], new_velocity[collision_idx,:], time

def evolve_position(r0, v, ball_pos, R_ball, box_size, step_length):
    t_remaining = step_length
    while time_remaining > 0:
        ball_intersections = get_ball_intersections(r0, v, ball_pos, R_ball, t_remaining)
        wall_collisions = get_wall_collisions(r0, v, box_size, t_remaining)
        if len(ball_intersections)>0:
            v, r0, delta_t = ball_intersections
        elif len(wall_collisions)>0:
            v, r0, delta_t = wall_collisions
        else:
            r0 = r0 + v*t
            break
        time_remaining -= delta_t
