#
# make an analyzer class. It should read the root files that were create
#
import uproot
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt

class Geant4Analyzer:
    def __init__(self, file_path, label=""):
        self.file_path = file_path
        self.label = label
        self.data = None
        self.preprocessed_data = {}

        self.load_data()

    def load_data(self):
        file = uproot.open(self.file_path)
        self.data = file["ev"].arrays()
        print(f"Data loaded from {self.file_path}")

    def preprocess_data(self):
        
        if self.data is None:
            raise ValueError("Data not loaded. Call load_data() first.")
        
        cut = (self.data['ncomp'] + self.data['nphot'] == 1) & (self.data['type'] == 0)
        cute = (self.data['eh'] > 1.)

        self.preprocessed_data['x'] = ak.to_numpy(ak.flatten(self.data['xh'][cut & cute]))
        self.preprocessed_data['y'] = ak.to_numpy(ak.flatten(self.data['yh'][cut & cute]))
        self.preprocessed_data['z'] = ak.to_numpy(ak.flatten(self.data['zh'][cut & cute]))
        self.preprocessed_data['e'] = ak.to_numpy(ak.flatten(self.data['eh'][cut & cute]))

    def plot_histogram(self, variable, bins=50, alpha=0.5, show=True):
        if variable not in self.preprocessed_data:
            raise ValueError(f"Variable '{variable}' not found in preprocessed data.")

        fig, ax = plt.subplots()
        ax.hist(self.preprocessed_data[variable], bins=bins, alpha=alpha, label=self.label)
        ax.set_xlabel(variable)
        ax.set_ylabel('Counts')
        ax.legend()
        ax.set_title(f'Histogram of {variable}')

        if show:
            plt.show()

        return ax