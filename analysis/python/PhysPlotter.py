import uproot
import numpy as np
import awkward as ak
from matplotlib.patches import Circle, Rectangle
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from XAMSPlotter import XAMSPlotter

from RunManager import RunManager
from mendeleev import element

def is_jagged(array):
    """Check if the given array is a jagged array."""
    return isinstance(array, ak.highlevel.Array) and isinstance(array.layout, ak.contents.ListOffsetArray)


class PhysPlotter:
    def __init__(self, run_id, first_only=True):
        """
        Initializes the analyzer with the given file path.

        Args:
            file_path (str): The path to the ROOT file.
            label (str, optional): The label to use for the plot.
        """

        manager = RunManager("../../run/rundb.json")

        self.file_paths = manager.get_output_root_files(run_id, first_only=first_only)
        self.settings = manager.get_run_settings(run_id, convert_units=True)
        self.geometry = manager.get_geometry(run_id)   
        self.materials = manager.get_materials(run_id)

        self.raw = None
        self.data = {}

        print(f"Initialized PhysPlotter with run_id={run_id}")	
        self.load_data()

    def load_data(self):
        """
        Load only the first file from the list of file paths: they all contain the physics data.

        Args:
            file_paths (str or list): Path or list of paths to the ROOT files.

        Raises:
            FileNotFoundError: If any file is not found.
        """

        if isinstance(self.file_paths, str):
            self.file_paths = [self.file_paths]
    
        if len(self.file_paths) == 0:
            raise FileNotFoundError("No files found")
            exit(-1)
        else:
            self.data = uproot.open(self.file_paths[0])["gam"].arrays(library="ak")

        print(f"Data loaded from {self.file_paths}")

    def plot_sigma(self, material=None, ax=None):
        """
        Plots the absorption cross-section ($\sigma$) for a given material.
        Parameters:
        -----------
        material : str, optional
            The material for which the absorption cross-section is to be plotted.
            If None, the function will print an error message and return.
        ax : matplotlib.axes.Axes, optional
            The axes on which to plot. If None, a new figure and axes will be created.
        Returns:
        --------
        ax : matplotlib.axes.Axes
            The axes with the plotted data.
        Notes:
        ------
        The function plots the following processes:
        - Compton scattering (dashed blue line)
        - Photoelectric effect (dash-dot red line)
        - Total absorption (solid black line)
        The x-axis represents the energy in MeV (log scale), and the y-axis represents
        the cross-section in barns (log scale).
        """
        if ax is None:
            fig, ax = plt.subplots()
            fig.set_size_inches(5, 5)

        if material is None:
            # exit with error: 'mat' is not defined
            print("Error: 'mat' is not defined")
            return
        
        d = self.data[self.data['mat'] == material]
        cut = (d['proc'] == "compton")
        ax.plot(d['e'][cut], d['att'][cut], linestyle='--', color="blue", label="Compton")	
        cut = (d['proc'] == "phot")
        ax.plot(d['e'][cut], d['att'][cut], linestyle='-.', color="red", label="Photoelectric")
        cut = (d['proc'] == "tot")
        ax.plot(d['e'][cut], d['att'][cut], color="black", label="Total")

        ax.legend(frameon=False)

        ax.set_xlabel("Energy (MeV)")
        ax.set_ylabel("$\\sigma$ (barn)")
        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.set_title(f"$\gamma$ absorption for {material}")

        return ax
    
    def plot_attenuation(self,material=None, ax=None, **kwargs):
        """
        Plots the attenuation length for a given material.
        Parameters:
        -----------
        material : str, optional
            The material for which the attenuation length is to be plotted. If not provided, the function will exit with an error.
        ax : matplotlib.axes.Axes, optional
            The axes on which to plot. If not provided, a new figure and axes will be created.
        **kwargs : dict, optional
            Additional keyword arguments to customize the plot. Supported keys:
            - 'linecolor': str, default 'black'
                The color of the plot line.
            - 'linestyle': str, default '-'
                The style of the plot line.
        Returns:
        --------
        ax : matplotlib.axes.Axes
            The axes with the plotted attenuation length.
        Notes:
        ------
        The function expects the data to be available in `self.data` with columns 'mat', 'proc', 'e', and 'att'.
        The 'proc' column should have the value 'att' for attenuation data.
        """


        if ax is None:
            fig, ax = plt.subplots()
            fig.set_size_inches(5, 5)

        if material is None:
            # exit with error: 'mat' is not defined
            print("Error: 'mat' is not defined")
            return
        
        if material not in self.data['mat']:
            print(f"Error: Material {material} not found in the data")
            return
        
        linecolor = kwargs.get("linecolor", "black")
        linestyle = kwargs.get("linestyle", "-")
        
        d = self.data[self.data['mat'] == material]
        cut = (d['proc'] == "att")
        ax.plot(d['e'][cut], d['att'][cut], linestyle=linestyle, color=linecolor, label=material)	

        ax.set_xlabel("Energy (MeV)")
        ax.set_ylabel("Attenuation length (cm)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend(frameon=False)
        ax.set_title(f"Linear attenuation for {material}")

        return ax


        