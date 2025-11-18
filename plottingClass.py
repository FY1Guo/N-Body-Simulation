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
        plt.show()
        return
    
    def trajectory_fluid_particle(self, traj, timespan="step", step_select=0):
        """
        This shows the trajectory of a single fluid particle over either a
        single time step or the full simulation.
        
        timespan = "step" for a single timestep
        timespan = "full" for the full simulation
        
        'traj' is the particle trajectory in 2D array form
        'step_select' is which step to use if plotting a single step
        """
        if timespan == "step":
            plt.plot(traj[step_select,0], traj[step_select,1])
            plt.title("Single Fluid Particle Trajectory for 1 Time Step")
            plt.xlabel("x Position")
            plt.ylabel("y Position")
            plt.show()
        elif timespan == "full":
            plt.plot(traj[:,0], traj[:,1])
            plt.title("Single Fluid Particle Trajectory for Full Simulation")
            plt.xlabel("x Position")
            plt.ylabel("y Position")
            plt.show()
        return
    
    def velocity_vs_time(self, vel):
        """
        'vel' is the velocity in 2D array form
        """
        vel_mag = (vel[:,0]**2 + vel[:,1]**2)**0.5
        plt.plot(self.time, vel_mag)
        plt.title("Ball Velocity Magnitude vs Time")
        plt.xlabel("Time")
        plt.ylabel("Velocity Magnitude")
        plt.show()
        return
        
    def force_vs_time(self, force):
        """
        'force' is the drag force in 1D array form
        """
        plt.plot(self.time, force)
        plt.title("Drag Force on Ball vs Time")
        plt.xlabel("Time")
        plt.ylabel("Drag Force")
        plt.show()
        return
    
    def force_vs_velocity(self, force, vel):
        """
        'force' is the drag force in 1D array form
        'vel' is the velocity in 2D array form
        """
        vel_mag = (vel[:,0]**2 + vel[:,1]**2)**0.5
        plt.plot(vel_mag, force)
        plt.title("Drag Force vs Velocity Magnitude")
        plt.xlabel("Velocity Magnitude")
        plt.ylabel("Drag Force")
        plt.show()
        return
    
    def energy_vs_time(self, energy):
        """
        'energy' is the total system energy in 1D array form
        """
        plt.plot(self.time, energy)
        plt.title("System Energy vs Time")
        plt.xlabel("Time")
        plt.ylabel("Energy")
        plt.show()
        return
