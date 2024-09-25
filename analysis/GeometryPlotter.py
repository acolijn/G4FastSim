import json
import numpy as np
import matplotlib.pyplot as plt

class GeometryPlotter:
    def __init__(self, geometry_file):
        """Constructor that loads the geometry from a JSON file."""
        self.geometry_file = geometry_file
        self.geometry_data = self.load_geometry(geometry_file)

    def load_geometry(self, geometry_file):
        """Loads the geometry from a JSON file."""
        with open(geometry_file, 'r') as file:
            return json.load(file)

    def draw_geometry(self):
        """Main function to draw the entire geometry."""
        fig, ax = plt.subplots()
        world = self.geometry_data['world']
        #ax.set_xlim(-world['size']*1000, world['size']*1000)
        #ax.set_ylim(-world['size']*1000, world['size']*1000)
        ax.set_xlim(-200,200)
        ax.set_ylim(-200,200)
        for volume in self.geometry_data['volumes']:
            print(f"Drawing volume: {volume['name']}")
            self.draw_volume(volume, ax)
        plt.show()

    def draw_volume(self, volume, ax):
        """Draws an individual volume based on its shape."""
        shape = volume['shape']
        if shape == 'tubs':
            self.draw_tubs(volume, ax)
        elif shape == 'box':
            self.draw_box(volume, ax)
        elif shape == 'sphere':
            self.draw_sphere(volume, ax)
        elif shape == 'union':
            self.draw_union(volume, ax)
        elif shape == 'subtraction':
            self.draw_subtraction(volume, ax)
        # Add more shapes if needed

    def draw_tubs(self, volume, ax):
        """Draws a 'tubs' shape."""
        dims = volume['dimensions']
        rMin = dims['rMin']
        rMax = dims['rMax']
        z = dims['z']
        placement = volume['placement']
        self.plot_cylinder(rMax, rMin, z, placement, ax)

    def draw_box(self, volume, ax):
        """Draws a 'box' shape."""
        dims = volume['dimensions']
        x = dims['x']
        y = dims['y']
        z = dims['z']
        placement = volume['placement']
        self.plot_box(x, y, z, placement, ax)

    def draw_sphere(self, volume, ax):
        """Draws a 'sphere' shape."""
        dims = volume['dimensions']
        rMin = dims['rMin']
        rMax = dims['rMax']
        placement = volume['placement']
        self.plot_circle(rMax, placement, ax)

    def draw_union(self, volume, ax):
        """Draws a 'union' of components."""
        for component in volume['components']:
            self.draw_volume(component, ax)
        self.apply_rotation(volume, ax)

    def draw_subtraction(self, volume, ax):
        """Handles 'subtraction' shape drawing."""
        # Subtract components visually (not detailed here)
        for component in volume['components']:
            self.draw_volume(component, ax)
        self.apply_rotation(volume, ax)

    def plot_cylinder(self, rMax, rMin, height, placement, ax):
        """Plots a cylinder (tubs) shape."""
        x, y = self.apply_translation(placement['x'], placement['y'])
        circle_outer = plt.Circle((x, y), rMax, color='blue', fill=False)
        circle_inner = plt.Circle((x, y), rMin, color='red', fill=False)
        ax.add_artist(circle_outer)
        if rMin > 0:
            ax.add_artist(circle_inner)

    def plot_box(self, x, y, z, placement, ax):
        """Plots a box shape."""
        x0, y0 = self.apply_translation(placement['x'], placement['y'])
        ax.add_patch(plt.Rectangle((x0 - x/2, y0 - y/2), x, y, fill=None))

    def plot_circle(self, rMax, placement, ax):
        """Plots a circle shape (for spheres)."""
        x, y = self.apply_translation(placement['x'], placement['y'])
        circle = plt.Circle((x, y), rMax, color='blue', fill=False)
        ax.add_artist(circle)

    def apply_translation(self, x, y):
        """Applies translation to the given coordinates."""
        return x, y

    def apply_rotation(self, volume, ax):
        """Applies rotation to a volume if it has rotation."""
        if 'rotation' in volume['placement']:
            rotation = volume['placement']['rotation']
            angle_z = rotation.get('z', 0.0)
            angle_y = rotation.get('y', 0.0)
            angle_x = rotation.get('x', 0.0)
            # Apply rotation matrix logic

    def rotation_matrix(self, angle_z, angle_y, angle_x):
        """Returns a combined rotation matrix for X, Y, Z rotations."""
        rz = np.array([[np.cos(angle_z), -np.sin(angle_z)], [np.sin(angle_z), np.cos(angle_z)]])
        ry = np.array([[np.cos(angle_y), 0, np.sin(angle_y)], [0, 1, 0], [-np.sin(angle_y), 0, np.cos(angle_y)]])
        rx = np.array([[1, 0, 0], [0, np.cos(angle_x), -np.sin(angle_x)], [0, np.sin(angle_x), np.cos(angle_x)]])
        return rz @ ry @ rx  # Combine all rotations

    def handle_union_rotation(self, volume, ax):
        """Handles the union of rotations in multiple components."""
        if 'rotation' in volume['placement']:
            rotation = volume['placement']['rotation']
            matrix = self.rotation_matrix(rotation['x'], rotation['y'], rotation['z'])
            # Apply the rotation to the components within the union
