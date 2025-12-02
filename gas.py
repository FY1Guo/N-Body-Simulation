import numpy as np


def get_v_ball_new(r_col, v, ball_pos, v_ball, m, M):
    # Compute new velocity after collision with ball

    n = r_col - ball_pos
    n_hat = n / np.linalg.norm(n)

    rel = v - v_ball
    u_n = np.dot(rel, n_hat) * n_hat

    coeff_p = 2 * M / (m + M)
    #coeff_b = 2 * m / (m + M)

    v_p_new = v - coeff_p * u_n
    #v_b_new = v_ball + coeff_b * u_n

    delta_v = v_p_new - v

    return v_p_new, delta_v #v_b_new


def get_ball_collisions(r0, v, ball_pos, R_ball, t_max, m_gas, M_ball):
    # Check for collisions with the ball which occur in time less than t_max

    # Relative coordinates to the ball center
    dr = r0 - ball_pos

    # Solve |dr + v*t|^2 = R_ball^2
    a = np.dot(v, v)
    b = 2 * np.dot(dr, v)
    c = np.dot(dr, dr) - R_ball**2

    if a == 0:
        return []

    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        return []  # No real roots, no collision
    sqrt_disc = np.sqrt(discriminant)
    t1 = (-b - sqrt_disc) / (2 * a)
    t2 = (-b + sqrt_disc) / (2 * a)

    collision_times = [t for t in (t1, t2) if 0 < t <= t_max]
    if not collision_times:
        return []
    t_collision = min(collision_times)
    r_col = r0 + v * t_collision
    v_col, delta_v = get_v_ball_new(r_col, v, ball_pos, R_ball, m_gas, M_ball)
    
    # Set new position to the edge of the ball.
    relative_pos = r_col - ball_pos
    r_new = ball_pos + relative_pos*R_ball/np.linalg.norm(relative_pos)

    # Compute change in velocity so we can compute total change in the ball's momentum due to collisions.
    return [r_new, v_col, t_collision, delta_v]
    # position, velocity, time when collision occurs. Return [] if no collision occurs


def get_wall_collisions(r0, v, box_size, t_max):
    x0, y0 = r0
    vx, vy = v
    L = box_size

    # Set x0 and y0 in bounds
    x0 = min(max(x0, 0), L)
    y0 = min(max(y0, 0), L)

    if vx == 0:
        t_right = np.inf
        t_left = np.inf
    else:
        t_right = (L - x0) / vx
        t_left = -x0 / vx
    if vy == 0:
        t_top = np.inf
        t_bot = np.inf
    else:
        t_top = (L - y0) / vy
        t_bot = -y0 / vy

    new_coordinates = np.array(
        [
            (L, y0 + vy * t_right),
            (0, y0 + vy * t_left),
            (x0 + vx * t_top, L),
            (x0 + vx * t_bot, 0),
        ]
    )
    new_velocity = np.array([(-vx, vy), (-vx, vy), (vx, -vy), (vx, -vy)])

    times = np.array([t_right, t_left, t_top, t_bot])
    valid_time_indices = np.where((times > 0) & (times <= t_max))
    valid_times = times[valid_time_indices]
    if len(valid_times) == 0:
        return []
    time_idx = np.argmin(valid_times)
    time = valid_times[time_idx]
    collision_idx = valid_time_indices[0][time_idx]

    return new_coordinates[collision_idx, :], new_velocity[collision_idx, :], time


def evolve_position(r0, v, ball_pos, v_ball, R_ball, box_size, step_length, m_gas, M_ball):
    """
    Evolves the position of a single gas particle through a time step.

    Returns:
    r0 : (2,) array
        Position of the particle at the end of the time step
    v : (2,) array
        Velocity of the particle at the end of the time step
    dv_total : (2,) array
        Total change in velocity of the gas particle due to ball collisions
    """
    t_remaining = step_length
    dv_total = np.array([0, 0], dtype = float)
    while t_remaining > 0:
        ball_intersections = get_ball_collisions(r0, v, ball_pos, R_ball, t_remaining, m_gas, M_ball)
        wall_collisions = get_wall_collisions(r0, v, box_size, t_remaining)
        if len(ball_intersections) > 0:
            v, r0, delta_t, delta_v = ball_intersections
            dv_total += delta_v
        elif len(wall_collisions) > 0:
            v, r0, delta_t = wall_collisions
        else:
            r0 = r0 + v * t_remaining
            break
        t_remaining -= delta_t
    return r0, v, dv_total

def update_projectile(ball_pos, ball_vel, impulse_from_particles, M_ball, box_size, dt):
    """
    Update the ball velocity and position for one time step.

    Returns:
    new_pos : (2,) array
    new_vel : (2,) array
    acceleration : (2,) array
    """
    
    delta_v = impulse_from_particles / M_ball # change in v
    new_vel = ball_vel + delta_v # updated velocity
    new_pos = ball_pos + new_vel * dt # updated position
    for i in range(2):

        # left wall
        if new_pos[i] < 0:
            new_pos[i] = -new_pos[i]
            new_vel[i] = -new_vel[i]

        # right wall
        elif new_pos[i] > box_size:
            new_pos[i] = 2*box_size - new_pos[i]
            new_vel[i] = -new_vel[i]

    #acceleration = (new_vel - ball_vel) / dt
    acceleration = delta_v / dt
    return new_pos, new_vel, acceleration

