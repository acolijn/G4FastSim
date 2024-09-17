import uproot
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from RunManager import RunManager
from mendeleev import element

def is_jagged(array):
    """Check if the given array is a jagged array."""
    return isinstance(array, ak.highlevel.Array) and isinstance(array.layout, ak.contents.ListOffsetArray)

class Geant4Analyzer:
    def __init__(self, run_id, label="", first_only=True):
        """
        Initializes the analyzer with the given file path.

        Args:
            file_path (str): The path to the ROOT file.
            label (str, optional): The label to use for the plot.
        """

        manager = RunManager("../run/rundb.json")

        self.file_paths = manager.get_output_root_files(run_id, first_only=first_only)
        self.settings = manager.get_run_settings(run_id, convert_units=True)
        self.label = ""

        if label == "":
            if 'ion' in self.settings['gps_settings']:
                txt = self.settings['gps_settings']['ion']
                Z = txt.split(" ")[0]
                A = txt.split(" ")[1]
                element_symbol = element(int(Z)).symbol
                label_text = 'ion: $^{{{:2d}}}${:s}'.format(int(A), element_symbol)
                self.label = label_text
            else:
                self.label = ""
        else:
            self.label = label

        self.raw = None
        self.data = {}

        print(f"Initialized Geant4Analyzer with run_id={run_id}, label={self.label}")	
        self.load_data()

    def load_data(self):
        """
        Loads the raw data from the file or list of files.

        Args:
            file_paths (str or list): Path or list of paths to the ROOT files.

        Raises:
            FileNotFoundError: If any file is not found.
        """

        print(f"Loading data from {self.file_paths}")
        if isinstance(self.file_paths, str):
            self.file_paths = [self.file_paths]

        data_list = []

        for file_path in self.file_paths:
            try:
                print(f"Loading {file_path}")
                root = uproot.open(file_path)
                data = root["ev"].arrays(library="ak")  # Use awkward array
                data_list.append(data)
            except FileNotFoundError:
                raise FileNotFoundError(f"File not found: {file_path}")

        if data_list:
            # Concatenate awkward arrays
            self.raw = ak.concatenate(data_list, axis=0)

            # Add derived variables
            self.raw['r'] = np.sqrt(self.raw['xh']**2 + self.raw['yh']**2)

            print(f"Data loaded from {len(self.file_paths)} files")
        else:
            print("No data loaded")

        print(f"Data loaded from {self.file_paths}")


    def load_data_obsolete(self):
        """
        Loads the raw data from the file.

        Raises:
            FileNotFoundError: If the file is not found.
        """
        if isinstance(self.file_paths, str):
            self.file_paths = [self.file_paths]

        for file in self.file_paths:

            root = uproot.open(file)
            self.raw = root["ev"].arrays()
            # add derived variables
            # radius
            self.raw['r'] = np.sqrt(self.raw['xh']**2 + self.raw['yh']**2)

        print(f"Data loaded from {self.file_paths}")

    def preprocess_data(self, cut=None, cut_hit=None):
        """
        Preprocesses the raw data by applying filters and converting it to a format suitable for analysis.

        Args:
            cut (array-like, optional): The primary cut to apply to the data.
            cut_hit (array-like, optional): Additional cut for jagged arrays.
        
        Raises:
            ValueError: If the data is not loaded. Call load_data() first.
        """
        if self.raw is None:
            raise ValueError("Data not loaded. Call load_data() first.")
        
        if cut is None:
            cut = cut
        else:
            cut = cut(self.raw)
            
        if cut_hit is None:
            cut_hit = cut  & (self.raw['eh'] > 0. )
        else:
            cut_hit = cut_hit(self.raw) & cut & (self.raw['eh'] > 1.)

        for field in self.raw.fields:
            data_field = self.raw[field]
            if is_jagged(data_field):
                # Flatten jagged arrays
                
                # make sure that you do not apply the cuts on the hits on the other fields   
                if ( (field == 'edet') or (field == 'ndet') or (field == 'ncomp') or (field == 'nphot') ):
                    data_field = data_field[cut]
                else:
                    data_field = ak.flatten(data_field[cut_hit])
            else:
                data_field = data_field[cut]

            self.data[field] = ak.to_numpy(data_field)


    def plot_histogram(self, variable, ax=None, bins=50, range=None, show=True):
        """
        Plots a histogram of the given variable.

        Args:

            variable (str): The variable to plot.
            ax (matplotlib.axes.Axes, optional): The axis to plot on. If None, a new figure is created.
            bins (int or array-like, optional): The number of bins or bin edges.
            range (tuple, optional): The range of the histogram.
            show (bool, optional): Whether to display the plot.

        Returns:
            matplotlib.axes.Axes: The axis object.

        Raises:
            ValueError: If the variable is not found in the preprocessed data.
        """ 
        if variable not in self.data:
            raise ValueError(f"Variable '{variable}' not found in preprocessed data.")

        if ax is None:
            fig, ax = plt.subplots()

        # use the event weights for the event variables, otherwise use the hit weights
        weights = self.data['w'] if len(self.data['w']) == len(self.data[variable]) else self.data['wh']

        hist, _ = np.histogram(self.data[variable], weights=np.exp(weights), bins=bins, range=range)
        print("integral =",np.sum(hist))

        ax.hist(self.data[variable], weights=np.exp(weights), bins=bins, range=range, histtype='step', label=self.label)

        if variable == 'r':
            ax.set_xlabel('radius (mm)')
        elif (variable == 'xp') or (variable == 'xh'):
            ax.set_xlabel('x (mm)')
        elif (variable == 'yp') or (variable == 'yh'):
            ax.set_xlabel('y (mm)')
        elif (variable == 'zp') or (variable == 'zh'):
            ax.set_xlabel('z (mm)')
        elif (variable == 'eh') or (variable == 'e'):
            ax.set_xlabel('Energy (keV)')
        else:
            ax.set_xlabel(variable)
    
        ax.set_ylabel('Counts')

        if show:
            ax.legend(frameon=False)
            plt.show()

        return ax

