import numpy as np
import pytest
import gas
import initialize
from main import stepper, M_ball, M_gas


""" Wall collisions"""

def test_wall_collision_left_and_right():
    r0 = np.array([0.1, 0.5])
    v  = np.array([-1.0, 0.0])
    box = 1.0
    t_max = 1.0

    r_new, v_new, t_hit = gas.get_wall_collisions(r0, v, box, t_max)

    assert t_hit > 0
    assert np.isclose(r_new[0], 0.0) or np.isclose(r_new[0], box) # Check that we hit a vertical wall: x should be 0 or box
    assert np.allclose(v_new, np.array([-v[0], v[1]])) # Velocity in x should flip sign, y unchanged

""" Update projectile zero impulse """

def test_update_projectile_zero_impulse_no_change():
    ball_pos = np.array([0.5, 0.5])
    ball_vel = np.array([1.0, 0.0])
    impulse = np.array([0.0, 0.0])
    M_ball = 10.0
    dt = 0.01
    box = 1.0

    new_pos, new_vel, acc = gas.update_projectile(ball_pos, ball_vel, impulse, M_ball, box, dt)

    assert np.allclose(new_vel, ball_vel)
    # new_pos should just be old_pos + v * dt
    assert np.allclose(new_pos, ball_pos + ball_vel * dt)


""" Evolve position free motion """

def test_evolve_position_free_motion():
    r0 = np.array([0.5, 0.5])
    v0 = np.array([0.2, -0.1])
    ball_pos = np.array([10.0, 10.0])   # Far away, no collisions
    ball_vel = np.array([0.0, 0.0])
    R_ball   = 0.1
    box      = 1.0
    dt       = 0.05

    r_new, v_new, dv = gas.evolve_position(
        r0, v0, ball_pos, ball_vel, R_ball,
        box, dt, M_gas, M_ball
    )

    assert np.allclose(r_new, r0 + v0 * dt)
    assert np.allclose(v_new, v0)
    assert np.allclose(dv, np.zeros(2))
