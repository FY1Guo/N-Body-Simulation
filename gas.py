import numpy as np

def get_ball_collisions(r0, v, ball_pos, R_ball, t_max):
    #Check for collisions with the ball which occur in time less than t_max
    return #position, velocity, time when collision occurs. Return [] if no collision occurs
    
def get_ball_delta_v(r0, v, ball_pos R_ball):
    

def get_wall_collisions(r0, v, box_size, t_max):
    #Check for collisions with the wall which occur in time less than t_max
    return #Positions, velocity, time when collision occurs
    
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
