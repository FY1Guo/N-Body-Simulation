import numpy as np

class initialize:
    def make_particles_vel(N, mu_vel=464, sigma_vel=50):
        """
        Input:
        N - number of particles

        Output:
        2xN Array of ball velocity in x-y basis
        """

        initial_angles = np.random.uniform(0,2*np.pi, size = N) 
        initial_magnitudes = np.random.normal(loc = mu_vel,scale = sigma_vel, size = N)
        out_arr = np.zeros((N,2))
        out_arr[:,0] = initial_magnitudes*np.cos(initial_angles)
        out_arr[:,1] = initial_magnitudes*np.sin(initial_angles)
        
        return out_arr
    
    def make_particles_pos(N, box_size, ball_pos, ball_radius):
        """
        Input:
        N - number of particles

        Output:
        Array of particle positions randomly distributed within the box
        """
        
        initial_x = np.random.uniform(0,1, size = N)*box_size
        initial_y = np.random.uniform(0,1, size = N)*box_size
        pos = np.array([initial_x, initial_y]).T
        for i in range(N):
            while np.sum((ball_pos-pos[i])**2)<ball_radius:
                pos[i] = np.random.uniform(0,1, size = 2)*box_size
        #out_arr = np.zeros((N,2))
        #for i in range(N):
        #    out_arr[i,0] += initial_x[i]
        #    out_arr[i,1] += initial_y[i]

        return pos
