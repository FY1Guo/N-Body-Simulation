import numpy as np

def get_ball_v_new(r_col, v, ball_pos, R_ball):
    #Compute new velocity after collision with ball

    n = r_col - ball_pos
    n_hat = n / np.linalg.norm(n)
    v_normal = np.dot(v, n_hat) * n_hat
    v_new = v - 2 * v_normal
    return v_new


def get_ball_collisions(r0, v, ball_pos, R_ball, t_max):
    #Check for collisions with the ball which occur in time less than t_max

    # Relative coordinates to the ball center
    dr = r0 - ball_pos

    # Solve |dr + v*t|^2 = R_ball^2
    a = np.dot(v, v)
    b = 2 * np.dot(dr, v)
    c = np.dot(dr, dr) - R_ball**2

    if a == 0:
        return []
    
    discriminant = b**2 - 4*a*c
    if discriminant < 0:
        return []  # No real roots, no collision
    sqrt_disc = np.sqrt(discriminant)
    t1 = (-b - sqrt_disc) / (2*a)
    t2 = (-b + sqrt_disc) / (2*a)

    collision_times = [t for t in (t1, t2) if 0 < t <= t_max]
    if not collision_times:
        return []
    t_collision = min(collision_times)
    r_col = r0 + v * t_collision
    v_col = get_ball_v_new(r_col, v, ball_pos, R_ball)
    
    #Compute change in velocity so we can compute total change in the ball's momentum due to collisions.
    delta_v = v_col-v
    return [r_col, v_col, t_collision, delta_v]
    #position, velocity, time when collision occurs. Return [] if no collision occurs
    
    
def get_ball_delta_v(r0, v, ball_pos, R_ball):
    return
    

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
    dv_total = np.array([0,0])
    while time_remaining > 0:
        ball_intersections = get_ball_collisions(r0, v, ball_pos, R_ball, t_remaining)
        wall_collisions = get_wall_collisions(r0, v, box_size, t_remaining)
        if len(ball_intersections)>0:
            v, r0, delta_t, delta_v = ball_intersections
            dv_total += delta_v
        elif len(wall_collisions)>0:
            v, r0, delta_t = wall_collisions
        else:
            r0 = r0 + v*t
            break
        time_remaining -= delta_t
   return r0, v, dv_total
