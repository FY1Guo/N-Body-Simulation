import numpy as np

class initialize:
    def make_particles_vel(N):
        """
        Input:
        N - number of particles

        Output:
        2xN Array of ball velocity in x-y basis
        """

        initial_angles = np.asarray([np.random.uniform(0,2*np.pi) for i in range(N)]) 
        initial_magnitudes = np.asarray([np.random.normal(loc = 464,scale = 50) for i in range(N)])
        out_arr = np.zeros((N,2))
        for i in range(N):
            out_arr[i,0] += initial_magnitudes[i] * np.cos(initial_angles[i]) 
            out_arr[i,1] += initial_magnitudes[i] * np.sin(initial_angles[i])

        return out_arr
    
    def make_particles_pos(N, box_size, ball_pos, ball_radius):
        """
        Input:
        N - number of particles

        Output:
        Array of particle positions randomly distributed within the box
        """
        
        initial_x = np.asarray([np.random.uniform(0,1) for i in range(N)])*box_size
        initial_y = np.asarray([np.random.uniform(0,1) for i in range(N)])*box_size
        for i in range(N):
            pos_i = np.array([initial_x[i], initial_y[i]])
            while np.sum((ball_pos-pos_i)**2)<ball_radius:
                initial_x[i] = np.random.uniform(0,1)*box_size
                initial_y[i] = np.random.uniform(0,1)*box_size
                pos_i = np.array([initial_x[i], initial_y[i]])
        out_arr = np.zeros((N,2))
        for i in range(N):
            out_arr[i,0] += initial_x[i]
            out_arr[i,1] += initial_y[i]

        return out_arr
