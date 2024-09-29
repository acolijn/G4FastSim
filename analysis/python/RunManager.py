import json
import pandas as pd
import os
import shutil
from mendeleev import element

class RunManager:
    def __init__(self, config_file="config.json"):
        """
        Initializes a RunManager object.

        Args:
            config_file (str): The path to the configuration file. Default is "config.json".
        """
        self.config_file = config_file
        self.load_config()

    def load_config(self):
        """
        Loads the configuration from the specified config file.
        If the 'status' field is missing for any run, it is initialized to 'active'.
        The modified config is saved and a DataFrame is created from the 'runs' field.
        """
        with open(self.config_file, 'r') as file:
            self.config = json.load(file)
        # Initialize status field if not present
        for run in self.config.get("runs", []):
            if "status" not in run:
                run["status"] = "active"
        self.save_config()  # Save the config if it was modified
        self.runs_df = pd.DataFrame(self.config.get("runs", []))

    def save_config(self):
        """
        Saves the configuration data to a JSON file.
        
        This method writes the configuration data stored in the `config` attribute
        to a JSON file specified by the `config_file` attribute. The data is written
        with an indentation level of 4 spaces.
        """
        with open(self.config_file, 'w') as file:
            json.dump(self.config, file, indent=4)

    def display_all_runs(self, include_deleted=False, detector_type=None):
        """
        Display all runs in the RunManager, with optional filtering by detector type.

        Args:
            include_deleted (bool): If True, include deleted runs in the output. 
                                    If False (default), exclude deleted runs.
            detector_type (str): If provided, filter runs that have this detector type
                                as specified in the geometry JSON file.

        Returns:
            pandas.DataFrame: A DataFrame containing the filtered runs to be displayed.
        """
        # Filter the runs by 'deleted' status
        if include_deleted:
            filtered_runs = self.runs_df
        else:
            filtered_runs = self.runs_df[self.runs_df['status'] != 'deleted']

        # If detector_type is specified, filter by that
        if detector_type:
            filtered_indices = []
            for index, run in filtered_runs.iterrows():
                run_id = run['id']
                geometry = self.get_geometry(run_id)
                if geometry:
                    # Check if the 'detector' key exists in the geometry
                    detector = geometry.get('detector')
                    if detector and detector == detector_type:
                        filtered_indices.append(index)
            
            filtered_runs = filtered_runs.loc[filtered_indices]
        
        # if particle is 'ion', get the ion name
        for index, run in filtered_runs.iterrows():
            # if particle is 'ion', get the ion name
            if run['particle'] == 'ion':
                settings = self.get_run_settings(run['id'])
                txt = settings['gps_settings']['ion']
                Z = txt.split(" ")[0]
                A = txt.split(" ")[1]
                element_symbol = element(int(Z)).symbol
                filtered_runs['ion'] = '{:s}{:2d}'.format(element_symbol, int(A))

            # get the geometry description
            geometry = self.get_geometry(run['id'])
            filtered_runs.loc[index,'geo description'] = geometry['description']


        # if detector is 'xams' drop a few columns that are not relevant
        if detector_type == 'xams':
            filtered_runs = filtered_runs.drop(columns=['fastSimulation','maxEnergy', 'maxScatters','outputFile', 'numJobs', 'sourceVolume'])
        return filtered_runs

    def get_output_root_files(self, run_id, first_only=False):
        """
        Retrieves the paths of the output ROOT files associated with the given run ID.

        Args:
            run_id (int): The ID of the run.
            first_only (bool): If True, return only the first ROOT file. If False, return all ROOT files.

        Returns:
            list or str or None: A list of paths to the output ROOT files if they exist, otherwise an empty list.
                                If first_only is True, returns the path to the first ROOT file as a string.
        """
        run = self.get_run_by_id(run_id)
        if run and run.get("status") != "deleted":
            output_dir = run.get("outputDir")
            if output_dir:
                root_files = [os.path.join(output_dir, file) for file in os.listdir(output_dir) if file.endswith(".root")]
                if first_only:
                    return root_files[0] if root_files else None
                return root_files
        return [] if not first_only else None


    def get_run_settings(self, run_id, convert_units=False):
        """
        Retrieves the settings for a specific run.

        Args:
            run_id (int): The ID of the run.
            convert_to_mm (bool): If True, convert all distance-related settings to millimeters.

        Returns:
            dict or None: The settings for the run if found, None otherwise.
        """
        run = self.get_run_by_id(run_id)
        if run and run.get("status") != "deleted":
            settings_file = os.path.join(run["outputDir"], "settings.json")
            if os.path.exists(settings_file):
                with open(settings_file, 'r') as file:
                    settings = json.load(file)
                if convert_units:
                    settings = self.convert_units(settings)
                #print(settings)
                
                return settings
        return None
    
    def get_geometry(self, run_id):
        """
        Retrieves the geometry for a specific run.

        Args:
            run_id (int): The ID of the run.

        Returns:
            dict or None: The geometry for the run if found, None otherwise.
        """

        settings = self.get_run_settings(run_id)

        if settings:
            # get the geometryFileName from the settings
            geometryFileName = settings["detector_configuration"]["geometryFileName"]
            # check if the file geometryFileName exists
            if os.path.exists(geometryFileName):
                with open(geometryFileName, 'r') as file:
                    geometry = json.load(file)
                return geometry
            else:
                print(f"Geometry file {geometryFileName} not found.")
                return None
            
        return None

    def convert_units(self, settings):
        """
        Convert all units in the given settings dictionary

        Args:
            settings (dict): The settings dictionary.

        Returns:
            dict: The settings dictionary with all distances converted to mm and keV.
        """
        conversion_factors = {
            "mm": 1,
            "cm": 10,
            "m": 1000,
            "keV": 1,
            "MeV": 1000
        }
        for key, value in settings.items():
            if isinstance(value, str) and any(unit in value for unit in conversion_factors):
                for unit, factor in conversion_factors.items():
                    #print(key, unit, factor)
                    if value.endswith(unit):
                        number = float(value.replace(unit, "").strip())
                        settings[key] = number * factor
                        break
        return settings

    def get_run_by_id(self, run_id):
        """
        Retrieves a run from the configuration by its ID.

        Args:
            run_id (int): The ID of the run to retrieve.

        Returns:
            dict or None: The run with the specified ID if found, otherwise None.
        """
        for run in self.config.get("runs", []):
            if run["id"] == run_id:
                return run
        return None

    def delete_run(self, run_id):
        """
        Deletes a run with the given run_id.

        Args:
            run_id (int): The ID of the run to be deleted.

        Returns:
            None

        Raises:
            None
        """
        run = self.get_run_by_id(run_id)
        if run:
            run["status"] = "deleted"
            self.save_config()
            output_dir = run.get("outputDir")
            if output_dir and os.path.exists(output_dir):
                shutil.rmtree(output_dir)
            print(f"Run {run_id} marked as deleted and its output has been removed.")
        else:
            print(f"Run {run_id} not found.")

