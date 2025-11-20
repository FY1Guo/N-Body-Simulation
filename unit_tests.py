import numpy as np
from gas import get_wall_collisions, update_projectile, evolve_position

def test_wall_collision_left_and_right():
    r0 = np.array([0.1, 0.5])
    v  = np.array([-1.0, 0.0])
    box = 1.0
    t_max = 1.0

    t_hit, wall = get_wall_collisions(r0, v, box, t_max)

    assert t_hit > 0
    assert wall in ("left", "right")


def test_update_projectile_zero_impulse_no_change():
    ball_pos = np.array([0.5, 0.5])
    ball_vel = np.array([1.0, 0.0])
    impulse = np.array([0.0, 0.0])
    M_ball = 10.0
    dt = 0.01
    box = 1.0

    new_pos, new_vel = update_projectile(ball_pos, ball_vel, impulse, M_ball, box, dt)

    assert np.allclose(new_vel, ball_vel)
    # new_pos should just be old_pos + v * dt
    assert np.allclose(new_pos, ball_pos + ball_vel * dt)


def test_evolve_position_no_collision_free_motion():
    r0 = np.array([0.5, 0.5])
    v = np.array([0.2, -0.1])

    ball_pos = np.array([2.0, 2.0])  # far away
    ball_vel = np.array([0.0, 0.0])
    R_ball = 0.1
    box_size = 10.0
    dt = 0.05
    m_gas = 1.0
    M_ball = 1000.0

    r_new, v_new, dv = evolve_position(r0, v, ball_pos, ball_vel, R_ball, box_size, dt, m_gas, M_ball)

    assert np.allclose(r_new, r0 + v * dt)
    assert np.allclose(v_new, v)        # no collision
    assert np.allclose(dv, np.zeros_like(v))  # no impulse on ball
