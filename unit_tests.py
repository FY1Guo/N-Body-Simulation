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
