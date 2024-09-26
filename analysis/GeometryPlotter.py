import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R

class GeometryPlotter:
    def __init__(self, geometry_file):
        """Initialize the GeometryPlotter with a geometry JSON file."""
        with open(geometry_file, 'r') as file:
            self.geometry_data = json.load(file)

    def plot_geometry(self, ax=None, view="xy"):
        """Plot the geometry on the provided axis."""
        if ax is None:
            fig, ax = plt.subplots()

        ax.set_ylim(-150, 150)
        ax.set_xlim(-150, 150)
        # Plot each volume recursively, starting with the world volume
        world_volume = [vol for vol in self.geometry_data['volumes'] if vol.get('parent') == 'World'][0]
        self.plot_volume(world_volume, ax, view, parent_transform=None)

        # Set axis labels based on the view
        if view == "xy":
            plt.xlabel("X (mm)")
            plt.ylabel("Y (mm)")
        elif view == "rz":
            plt.xlabel("Radius (mm)")
            plt.ylabel("Z (mm)")

        plt.show()

    def plot_volume(self, volume, ax, view, parent_transform=None):
        """Plot an individual volume, applying transformations from parent volumes."""
        shape_type = volume['shape']

        # Apply parent transformation (position and rotation)
        volume_transform = self.get_volume_transform(volume)
        if parent_transform:
            volume_transform = self.combine_transforms(parent_transform, volume_transform)

        # Plot based on volume type
        if shape_type in ('tubs', 'box', 'sphere'):
            self.plot_simple_volume(volume, ax, view, volume_transform)
        elif shape_type in ('union', 'subtraction'):
            self.plot_composite_volume(volume, ax, view, volume_transform)

        # Recursively plot child volumes
        for child_volume in self.get_child_volumes(volume['name']):
            self.plot_volume(child_volume, ax, view, parent_transform=volume_transform)

    def plot_composite_volume(self, volume, ax, view, volume_transform):
        """Handle composite volumes with union or subtraction."""
        components = volume['components']

        # Collect component coordinates, applying only relative transformations
        all_coords = []
        for component in components:
            # Apply the local placement (relative to the composite object) but not the final placement
            coords = self.get_component_coords(component, view)
            all_coords.append((component, coords))

        # Apply final transformation (rotation + translation) to the entire composite
        transformed_coords = []
        for component, coords in all_coords:
            transformed_coords.append(self.apply_placement(component, coords, volume_transform))

        # Plot the final combined object with the overall transformation
        # for coords in transformed_coords:
        #     for coord in coords:
        #         if len(coord) == 2:
        #             x, y = coord
        #             z = 0
        #         elif len(coord) == 3:
        #             x, y, z = coord
        #         else:
        #             raise ValueError("Coordinate format not supported.")

        #         if view == "xy":
        #             self.draw_circle(ax, x, y, 50, {})  # Example, replace with actual shape drawing
        #         elif view == "rz":
        #             r = np.sqrt(x**2 + y**2)
        #             self.draw_rectangle(ax, r, z, 50, 100, {}, view)

    def get_component_coords(self, component, view):
        """Get coordinates of a component relative to its parent composite volume."""
        placement = component.get('placement', {})
        x = placement.get('x', 0)
        y = placement.get('y', 0)
        z = placement.get('z', 0)
        
        # Return the appropriate coordinate based on the view (xy or rz)
        if view == "xy":
            return [(x, y, z)]
        elif view == "rz":
            r = np.sqrt(x**2 + y**2)
            return [(r, z)]

    def get_volume_transform(self, volume):
        """Return the transformation (translation + rotation) of the volume."""
        placement = volume.get('placement', {})
        x = placement.get('x', 0)
        y = placement.get('y', 0)
        z = placement.get('z', 0)
        rotation = placement.get('rotation', {'x': 0, 'y': 0, 'z': 0})

        return {'position': np.array([x, y, z]), 'rotation': rotation}

    def combine_transforms(self, parent_transform, child_transform):
        """Combine parent and child transformations."""
        # Apply parent's rotation to child's position
        child_position = self.apply_rotation(child_transform['position'], parent_transform['rotation'])
        combined_position = parent_transform['position'] + child_position

        # Combine rotations (Euler angles)
        combined_rotation = {
            'x': parent_transform['rotation']['x'] + child_transform['rotation']['x'],
            'y': parent_transform['rotation']['y'] + child_transform['rotation']['y'],
            'z': parent_transform['rotation']['z'] + child_transform['rotation']['z']
        }

        return {'position': combined_position, 'rotation': combined_rotation}

    def apply_rotation(self, coords, rotation_angles):
        """Apply 3D rotation using Euler angles."""
        rotation = R.from_euler('xyz', [rotation_angles['x'], rotation_angles['y'], rotation_angles['z']], degrees=True)
        return rotation.apply(coords)

    def apply_placement(self, volume, coords, transform):
        """Apply the global transformation (rotation + translation) to the entire object."""
        x_offset, y_offset, z_offset = transform['position']

        # Apply rotation if needed and handle 2D or 3D coordinates
        transformed_coords = []
        for coord in coords:
            if len(coord) == 2:  # If only (x, y) is provided
                x, y = coord
                z = 0  # Default z to 0 for 2D cases
            elif len(coord) == 3:  # If (x, y, z) is provided
                x, y, z = coord
            else:
                raise ValueError(f"Invalid coordinate format: {coord}")
            
            # Apply rotation
            rotated_coord = self.apply_rotation(np.array([x, y, z]), transform['rotation'])
            transformed_coords.append(rotated_coord)

        # Apply translation after rotation
        final_coords = [(x + x_offset, y + y_offset, z + z_offset) for x, y, z in transformed_coords]
        return final_coords


    def get_child_volumes(self, parent_name):
        """Retrieve child volumes of a given parent."""
        return [vol for vol in self.geometry_data['volumes'] if vol.get('parent') == parent_name]

    def plot_simple_volume(self, volume, ax, view, transform):
        """Plot a simple volume like a tub, box, or sphere."""
        shape = volume['shape']
        coords = transform['position']

        if shape == "tubs":
            radius = volume['dimensions']['rMax']
            height = volume['dimensions']['z']
            x, y, z = coords
            if view == "xy":
                self.draw_circle(ax, x, y, radius, self.get_color_properties(volume))
            elif view == "rz":
                r = np.sqrt(x**2 + y**2)
                self.draw_rectangle(ax, r, z, radius*2, height, self.get_color_properties(volume), view)

        elif shape == "box":
            x_length = volume['dimensions']['x']
            y_length = volume['dimensions']['y']
            z_length = volume['dimensions']['z']
            x, y, z = coords
            if view == "xy":
                self.draw_rectangle(ax, x, y, x_length, y_length, self.get_color_properties(volume))
            elif view == "rz":
                radius = np.sqrt(x**2 + y**2)
                self.draw_rectangle(ax, radius, z, x_length, z_length, self.get_color_properties(volume), view)

    def draw_circle(self, ax, x, y, radius, properties):
        """Draw a circle (for tubs) on the given axis."""
        circle = plt.Circle((x, y), radius, **properties)
        ax.add_patch(circle)

    def draw_rectangle(self, ax, x, y, width, height, properties, view):
        """Draw a rectangle (for box) on the given axis."""
        rect = plt.Rectangle((x - width / 2, y - height / 2), width, height, **properties)
        ax.add_patch(rect)

    def get_color_properties(self, volume):
        """Convert color from list to dictionary suitable for matplotlib."""
        color = volume.get('color', [0.5, 0.5, 0.5, 0.5])
        color = [0.5, 0.5, 0.5, 0.1]
        return {'color': (color[0], color[1], color[2]), 'alpha': color[3]}
