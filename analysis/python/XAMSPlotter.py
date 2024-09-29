import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Rectangle, Circle
from shapely.geometry import Polygon
from shapely.ops import unary_union

class XAMSPlotter:
    def __init__(self, geometry):
        """
        Initializes the XAMSPlotter with the provided geometry data.

        Args:
            geometry (dict): A dictionary containing the geometry data. 
                             It should include a 'volumes' key with a list of volume dictionaries.

        Attributes:
            geometry_data (dict): Stores the provided geometry data.
            volumes (dict): A dictionary mapping volume names to their respective volume data.
        """
        #with open(geometry_file, 'r') as file:
        self.geometry_data = geometry
        self.volumes = {volume['name']: volume for volume in self.geometry_data['volumes']}


    def plot_geometry(self, ax=None, view="xy"):
        """
        Plots the geometry of the simulation setup.

        Parameters:
        ax (matplotlib.axes.Axes, optional): The axes on which to plot. If None, a new figure and axes are created.
        view (str, optional): The view for the plot, either "xy" or "rz". Defaults to "xy".

        The function plots various components of the simulation setup including:
        - Outer cryostat
        - Inner cryostat
        - Gaseous xenon
        - Liquid xenon
        - PTFE bucket
        - PTFE reflection cylinder
        - NaI
        - Collimator
        - Source

        The axis labels are set based on the view parameter:
        - "xy": x (mm) and y (mm)
        - "rz": r (mm) and z (mm)

        The plot is displayed using plt.show().
        """
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_ylim(-150, 150)
            ax.set_xlim(-150, 150)

        # if a lead shield is present, plot it
        if 'LeadShield' in self.volumes:
            self.plot_cylinder(ax, view, "LeadShield")


        # plot outer cryostat
        self.plot_cylinder(ax, view, "OuterCryostat")
        # plot inner cryostat
        self.plot_cylinder(ax, view, "InnerCryostat")
        # plot gaseous xenon
        self.plot_cylinder(ax, view, "GaseousXenon")
        # plot liquid xenon
        self.plot_cylinder(ax, view, "LiquidXenon")
        # plot PTFE bucket
        self.plot_cylinder(ax, view, "PTFEBucket")
        if 'FieldShape' in self.volumes:
            self.plot_cylinder(ax, view, "FieldShape")
        # plot PTFE reflection cylinder
        self.plot_cylinder(ax, view, "PTFEReflectionCylinder")
        # plot NaI
        if 'NaI' in self.volumes:
            self.plot_cylinder(ax, view, "NaI")
        # plot collimator
        if 'Collimator' in self.volumes:
            self.plot_collimator(ax, view)
        # plot source
        if 'SourceSphere' in self.volumes:
            self.plot_source(ax, view)

        # Set axis labels based on the view
        if view == "xy":
            ax.set_xlabel("$x$ (mm)")
            ax.set_ylabel("$y$ (mm)")
        elif view == "rz":
            ax.set_xlabel("r (mm)")
            ax.set_ylabel("z (mm)")

        #plt.show()

    def plot_source(self, ax, view):
        """
        Plots the source sphere on the given matplotlib axis.

        Parameters:
        ax (matplotlib.axes.Axes): The axis on which to plot the source sphere.
        view (str): The view type, currently supports "xy" for plotting in the xy-plane.

        Returns:
        None
        """
        volume     = self.volumes['SourceSphere']

        if view == "xy":
            x = volume['placement']['x']
            y = volume['placement']['y']
            r = volume['dimensions']['rMax']
            circle = Circle((x, y), r, edgecolor='black', facecolor='green', alpha=0.7)
            ax.add_patch(circle)
        elif view == "rz":
            x = volume['placement']['x']
            y = volume['placement']['y']
            z = volume['placement']['z']
            r = np.sqrt(x**2 + y**2)
            rSource = volume['dimensions']['rMax']
            circle = Circle((r, z), rSource, edgecolor='black', facecolor='green', alpha=0.7)
            ax.add_patch(circle)


    def plot_collimator(self, ax, view):
        """
        Plots the collimator on the given axes in the specified view.
        Parameters:
        ax (matplotlib.axes.Axes): The axes on which to plot the collimator.
        view (str): The view in which to plot the collimator. Currently supports "xy".
        Returns:
        None
        """
        volume     = self.volumes['Collimator']
        if 'components' in volume:
            components = {component['name']: component for component in volume['components']}
     
        if view == "xy":
            #print("Plotting xy view : Collimator")
            # left side of collimator block
            x = volume['placement']['x'] - components['CollimatorBlock']['dimensions']['x']/2.
            y = volume['placement']['y'] - components['CollimatorBlock']['dimensions']['z']/2.
            width = components['CollimatorBlock']['dimensions']['x']/2.-components['CollimatorHole']['dimensions']['rMax']
            height = components['CollimatorBlock']['dimensions']['z']
            rectangle = Rectangle((x, y), width, height, edgecolor='black', facecolor='grey', alpha=0.7)
            ax.add_patch(rectangle)
            # right side of collimator block
            x = volume['placement']['x'] + components['CollimatorHole']['dimensions']['rMax']
            y = volume['placement']['y'] - components['CollimatorBlock']['dimensions']['z']/2.

            rectangle = Rectangle((x, y), width, height, edgecolor='black', facecolor='grey', alpha=0.7)
            ax.add_patch(rectangle)
        elif view == "rz":
            #print("Plotting rz view : Collimator")
            x = volume['placement']['x']
            y = volume['placement']['y']
            z = volume['placement']['z']
            r = np.sqrt(x**2 + y**2)
            height = components['CollimatorBlock']['dimensions']['x']/2
            rHole = components['CollimatorHole']['dimensions']['rMax']
            rectangle = Rectangle((r-components['CollimatorBlock']['dimensions']['z']/2., z-height/2.), components['CollimatorBlock']['dimensions']['z'], height/2-rHole, edgecolor='black', facecolor='grey', alpha=0.7)
            ax.add_patch(rectangle)
            rectangle = Rectangle((r-components['CollimatorBlock']['dimensions']['z']/2., z+rHole), components['CollimatorBlock']['dimensions']['z'], height/2-rHole, edgecolor='black', facecolor='grey', alpha=0.7)
            ax.add_patch(rectangle)

    def plot_cylinder(self, ax, view, name):
        """
        Plots a cylindrical volume on the given axis in the specified view.
        Parameters:
        ax (matplotlib.axes.Axes): The axis on which to plot the cylinder.
        view (str): The view in which to plot the cylinder. Currently supports "xy".
        name (str): The name of the volume to plot. Should be one of "OuterCryostat" or "InnerCryostat".
        Raises:
        SystemExit: If the name is invalid.
        Notes:
        - For "OuterCryostat", the inner radius (rMin) is set to the outer radius of the 'Vacuum' volume.
        - For "InnerCryostat", the inner radius (rMin) is set to the outer radius of the 'GaseousXenon' volume.
        - The function uses the 'draw_ring' method to plot an annulus representing the cylinder.
        """
        volume     = self.volumes[name]
        if 'components' in volume:
            components = {component['name']: component for component in volume['components']}
     
        rMin = 0
        rMax = 0
        alpha = 0.5
        edgecolor = None

        if view == "xy":
            #print("Plotting xy view :", name)
            if name == "OuterCryostat":
                rMin = self.volumes['Vacuum']['dimensions']['rMax']
                rMax = components['OuterCryostatCylinder']['dimensions']['rMax']
                facecolor = 'grey'
                edgecolor = 'black'
                alpha = 0.8
            elif name == "InnerCryostat":
                rMin = self.volumes['GaseousXenon']['dimensions']['rMax'] 
                rMax = components['InnerCryostatCylinder']['dimensions']['rMax']
                facecolor = 'grey'
                edgecolor = 'black'
                alpha = 0.8
            elif name == "LeadShield":
                rMin = self.volumes[name]['dimensions']['rMin']
                rMax = self.volumes[name]['dimensions']['rMax']
                facecolor = 'grey'
                edgecolor = 'black'
                alpha = 0.8
            elif name == "GaseousXenon":
                vtmp = self.volumes['PTFEBucket']
                ptfeComponents = {component['name']: component for component in vtmp['components']}
                rMin = ptfeComponents['PTFEWall']['dimensions']['rMax']
                rMax = self.volumes[name]['dimensions']['rMax']
                facecolor = 'lightblue'
                alpha = 0.3
            elif name == "LiquidXenon":
                rMin = 0.0
                rMax = self.volumes[name]['dimensions']['rMax']
                facecolor = 'darkblue'
                alpha = 0.3
            elif name == "PTFEBucket":
                rMin = components['PTFEWall']['dimensions']['rMin']
                rMax = components['PTFEWall']['dimensions']['rMax']
                facecolor = 'white'
                edgecolor = 'black'
                alpha = 0.8
            elif name == "PTFEReflectionCylinder":
                rMin = self.volumes[name]['dimensions']['rMin']
                rMax = self.volumes[name]['dimensions']['rMax']
                facecolor = 'white'
                edgecolor = 'black'
                alpha = 0.8
            elif name == "NaI":
                rMin = 0.
                rMax = self.volumes[name]['dimensions']['rMax']
                facecolor = 'green'
                edgecolor = 'black'
                alpha = 0.2
            elif name == "FieldShape":
                return
            else:
                print("Invalid name")
                exit(-1)

            xc = volume['placement']['x']
            yc = volume['placement']['y']
            # plot the x-y view of the cylinder. draw an annulus 360deg with rMin and rMax and centr (xc,yc)
            self.draw_ring(ax, (xc, yc), rMin, rMax, angle_range=(0, 360), edgecolor=edgecolor, facecolor=facecolor, alpha=alpha)
        
        elif view == 'rz':
            #print("Plotting rz view :", name)
            x = []
            y = []
            if name == "OuterCryostat" or name == "InnerCryostat":
                rMin = 0.0
                rMax = components[name+'Cylinder']['dimensions']['rMax']
                dz = components[name+'Cylinder']['dimensions']['z']
                mainCylinder = self.create_rectangle(0, -dz/2, rMax, dz)
                
                rMin = 0.0
                rMax = components[name+'BottomFlange']['dimensions']['rMax']
                dz = components[name+'BottomFlange']['dimensions']['z']
                z = components[name+'BottomFlange']['placement']['z']
                bottomFlange = self.create_rectangle(0, z-dz/2, rMax, dz)

                rMin = 0.0
                rMax = components[name+'TopFlange']['dimensions']['rMax']
                dz = components[name+'TopFlange']['dimensions']['z']
                z = components[name+'TopFlange']['placement']['z']
                topFlange = self.create_rectangle(0, z-dz/2, rMax, dz)
                
                union_shape =unary_union([mainCylinder, bottomFlange, topFlange])
                z_offset = 0
                if name == "OuterCryostat":
                    
                    # subtract the vacuum
                    rMin = 0.0
                    rMax = self.volumes['Vacuum']['dimensions']['rMax']
                    dz = self.volumes['Vacuum']['dimensions']['z']
                    z = self.volumes['Vacuum']['placement']['z']
                    vacuum = self.create_rectangle(0, z-dz/2, rMax, dz)
                    
                    union_shape = union_shape.difference(vacuum)
                else:
                    #subtract the gaseous xenon
                    rMin = 0.0
                    rMax = self.volumes['GaseousXenon']['dimensions']['rMax']
                    dz = self.volumes['GaseousXenon']['dimensions']['z']
                    z = self.volumes['GaseousXenon']['placement']['z']
                    gaseousXenon = self.create_rectangle(0, z-dz/2, rMax, dz)
                    union_shape = union_shape.difference(gaseousXenon)

                x,y=union_shape.exterior.xy

                z_offset = self.volumes[name]['placement']['z']
                x = np.array(x)
                y = np.array(y)+z_offset
                ax.fill(x, y, color='lightgrey', alpha=0.6)  # Fill the shape with color
                ax.plot(x, y, linewidth=1., color='black', alpha=0.6)	
            elif name == "FieldShape":
                rMin = self.volumes[name]['dimensions']['rMin']
                rMax = self.volumes[name]['dimensions']['rMax']
                dz_ring = self.volumes[name]['dimensions']['z']  # Thickness of a ring
                num_rings = self.volumes[name]['repetitions']['count']  # Number of rings
                dz_spacing = self.volumes[name]['repetitions']['dz']  # Spacing between rings
        
                z_offset = self.volumes['InnerCryostat']['placement']['z']+self.volumes['LiquidXenon']['placement']['z']+self.volumes[name]['placement']['z']

                coppercolor = '#B87333'
                for i in range(num_rings):
                    z = z_offset + i * dz_spacing
                    rectangle = Rectangle((rMin, z - dz_ring / 2), rMax - rMin, dz_ring,
                                  edgecolor=None, facecolor=coppercolor, alpha=0.4)
                    ax.add_patch(rectangle)
            elif name == "LeadShield":
                rMin = self.volumes[name]['dimensions']['rMin']
                rMax = self.volumes[name]['dimensions']['rMax']
                dz = self.volumes[name]['dimensions']['z']
                z_offset = self.volumes[name]['placement']['z']
                rectangle = Rectangle((rMin, -dz/2.+z_offset), rMax-rMin, dz, edgecolor='black', facecolor='grey', alpha=0.8)
                ax.add_patch(rectangle)
            elif name == "PTFEBucket":
                rMin = components['PTFEWall']['dimensions']['rMin']
                rMax = components['PTFEWall']['dimensions']['rMax']
                dz = components['PTFEWall']['dimensions']['z']
                mainCylinder = self.create_rectangle(rMin, -dz/2, rMax-rMin, dz)
                
                rMin = 0.0
                rMax = components['PTFEBottom']['dimensions']['rMax']
                dz = components['PTFEBottom']['dimensions']['z']
                z = components['PTFEBottom']['placement']['z']
                bottomFlange = self.create_rectangle(0, z-dz/2, rMax, dz)

                union_shape =unary_union([mainCylinder, bottomFlange])

                z_offset = self.volumes[name]['placement']['z']
                z_offset += self.volumes['InnerCryostat']['placement']['z']

                x,y=union_shape.exterior.xy
                x = np.array(x)
                y = np.array(y)+z_offset
                ax.fill(x, y, color='white', alpha=0.9)
                ax.plot(x, y, linewidth=1., color='black', alpha=0.9)

            elif name == "GaseousXenon":
                rMin = 0.0
                rMax = self.volumes[name]['dimensions']['rMax'] 
                height = self.volumes[name]['dimensions']['z']
                z_offset = self.volumes['InnerCryostat']['placement']['z']
                rectangle = Rectangle((0, -height/2.+z_offset), rMax, height, edgecolor=None, facecolor='lightblue', alpha=0.3)
                ax.add_patch(rectangle)

            elif name == "LiquidXenon":
                rMin = self.volumes[name]['dimensions']['rMin']
                rMax = self.volumes[name]['dimensions']['rMax'] 
                height = self.volumes[name]['dimensions']['z']
                z_offset = self.volumes['InnerCryostat']['placement']['z']+self.volumes[name]['placement']['z']
                rectangle = Rectangle((0, -height/2.+z_offset), rMax, height, edgecolor=None, facecolor='darkblue', alpha=0.3)
                ax.add_patch(rectangle)
            elif name == "PTFEReflectionCylinder":
                rMin = self.volumes[name]['dimensions']['rMin']
                rMax = self.volumes[name]['dimensions']['rMax'] 
                height = self.volumes[name]['dimensions']['z']
                z_offset = self.volumes['InnerCryostat']['placement']['z']+self.volumes['LiquidXenon']['placement']['z']+self.volumes[name]['placement']['z']
                rectangle = Rectangle((rMin, -height/2.+z_offset), rMax-rMin, height, edgecolor='grey', facecolor='white', alpha=0.9)
                ax.add_patch(rectangle)
            elif name == "NaI":
                r = np.sqrt(self.volumes[name]['placement']['x']**2+self.volumes[name]['placement']['y']**2)
                rMin = self.volumes[name]['dimensions']['rMin']
                rMax = self.volumes[name]['dimensions']['rMax'] 
                height = self.volumes[name]['dimensions']['z']
                z_offset = self.volumes[name]['placement']['z']
                rectangle = Rectangle((r-rMax, -height/2.+z_offset), 2*rMax, height, edgecolor='black', facecolor='green', alpha=0.2)
                ax.add_patch(rectangle)
            else:
                print("Invalid name")
                exit(-1)

    # Function to create a Shapely polygon for a rectangle
    def create_rectangle(self, x, y, width, height):
        return Polygon([(x, y), (x + width, y), (x + width, y + height), (x, y + height)])


    def draw_ring(self, ax, center, inner_radius, outer_radius, angle_range=(0, 360), **kwargs):
        """
        Draws a ring (annulus) on a plot using Wedge.

        Parameters:
        - ax: The matplotlib axis to draw on.
        - center: (x, y) tuple for the center of the ring.
        - inner_radius: The inner radius of the ring.
        - outer_radius: The outer radius of the ring.
        - angle_range: (start_angle, end_angle) of the ring in degrees. Default is full 360.
        - kwargs: Additional keyword arguments to pass to the Wedge (e.g., color, alpha).
        """
        wedge = Wedge(center, outer_radius, angle_range[0], angle_range[1], width=outer_radius - inner_radius, **kwargs)
        ax.add_patch(wedge)
