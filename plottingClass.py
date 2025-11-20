import matplotlib.pyplot as plt

class plotting:
    def __init__(self, time):
        """
        'time' is a 1D array from initial time to end time with step size dt
        """
        
        self.time = time 
        return
    
    def trajectory_ball(self, traj):
        """
        This shows the trajectory of the ball over the entire simulation
        
        'traj' is the ball trajectory in 2D array form
        """
        
        plt.plot(traj[:,0], traj[:,1])
        plt.title("Ball Trajectory for Full Simulation")
        plt.xlabel("x Position")
        plt.ylabel("y Position")
        plt.grid()
        plt.show()
        return
    
    def trajectory_fluid_particle(self, traj, timespan="full"):
        """
        This shows the trajectory of a single fluid particle over either a
        single time step or the full simulation.
        
        timespan = "full" for the full simulation
        timespan = "step" for a single timestep
        
        'traj' is the particle trajectory in 2D array form
        """
        
        plt.plot(traj[:,0], traj[:,1])
        if timespan == "step":
            plt.title("Single Fluid Particle Trajectory for 1 Time Step")
        elif timespan == "full":
            plt.title("Single Fluid Particle Trajectory for Full Simulation")
        plt.xlabel("x Position")
        plt.ylabel("y Position")
        plt.grid()
        plt.show()
        return
    
    def velocity_vs_time(self, vel):
        """
        This shows the velocity of the ball over time
        
        'vel' is the velocity in 2D array form
        """
        
        vel_mag = (vel[:,0]**2 + vel[:,1]**2)**0.5
        plt.plot(self.time, vel_mag, label="full magnitude")
        plt.plot(self.time, vel[:,0], label="x magnitude")
        plt.plot(self.time, vel[:,1], label="y magnitude")
        plt.title("Ball Velocity vs Time")
        plt.xlabel("Time")
        plt.ylabel("Velocity")
        plt.legend()
        plt.grid()
        plt.show()
        return
        
    def force_vs_time(self, force):
        """
        This shows the force of drag on the ball over time
        
        'force' is the drag force in 1D array form
        """
        
        plt.plot(self.time, force)
        plt.title("Drag Force on Ball vs Time")
        plt.xlabel("Time")
        plt.ylabel("Drag Force")
        plt.grid()
        plt.show()
        return
    
    def force_vs_velocity(self, force, vel):
        """
        This shows the force of drag on the ball vs velocity of the ball
        
        'force' is the drag force in 1D array form
        'vel' is the velocity in 2D array form
        """
        
        vel_mag = (vel[:,0]**2 + vel[:,1]**2)**0.5
        plt.plot(vel_mag, force)
        plt.title("Drag Force vs Velocity Magnitude")
        plt.xlabel("Velocity Magnitude")
        plt.ylabel("Drag Force")
        plt.grid()
        plt.show()
        return
    
    def total_energy_vs_time(self, energy):
        """
        This shows the total system energy over time
        
        'energy' is the total system energy in 1D array form
        """
        
        plt.plot(self.time, energy)
        plt.title("System Energy vs Time")
        plt.xlabel("Time")
        plt.ylabel("Energy")
        plt.grid()
        plt.show()
        return
    
    def kinetic_energy_vs_time(self, k_energy, scope="ball"):
        """
        This shows kinetic energy of the ball over time
        
        scope = "ball" for ball kinetic energy
        scope = "particles" for the average kinetic energy of the particles
        
        'k_energy' is the kinetic energy in 1D array form
        """
        
        plt.plot(self.time, k_energy)
        if scope == "ball":
            plt.title("Ball Kinetic Energy vs Time")
        elif scope == "particles":
            plt.title("Average Kinetic Energy of Particles vs Time")
        plt.xlabel("Time")
        plt.ylabel("Energy")
        plt.grid()
        plt.show()
