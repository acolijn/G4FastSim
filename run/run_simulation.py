import json
import os
import shutil
import subprocess
import datetime
import argparse
import random

class PathManager:
    """
    A class that manages the paths and directories for a simulation.

    Args:
        project_base_dir (str): The base directory of the project.
        output_dir (str, optional): The output directory for the simulation. If not provided, a default directory will be used.

    Attributes:
        project_base_dir (str): The base directory of the project.
        timestamp (str): The current timestamp in the format "%Y%m%d_%H%M%S".
        output_base_dir (str): The base directory for the simulation output.
        output_dir (str): The output directory for the simulation.
        jobs_dir (str): The directory for storing job files.
        logs_dir (str): The directory for storing log files.
    """

    def __init__(self, project_base_dir, output_dir=None):
        self.project_base_dir = project_base_dir
        self.timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        self.output_base_dir = output_dir if output_dir else "/data/xenon/acolijn/G4FastSim/"
        self.output_dir = self.create_output_dir()
        self.jobs_dir = os.path.join(self.output_dir, "jobs")
        self.logs_dir = os.path.join(self.output_dir, "logs")
        self.create_directory_structure()

    def create_output_dir(self):
        output_dir = os.path.join(self.output_base_dir, f"{self.timestamp}")
        os.makedirs(output_dir, exist_ok=True)
        return output_dir

    def create_directory_structure(self):
        os.makedirs(self.jobs_dir, exist_ok=True)
        os.makedirs(self.logs_dir, exist_ok=True)

##############################################################################################

def load_settings(json_file):
    """
    Load settings from a JSON file.
    
    Parameters:
        json_file (str): The path to the JSON file.
        
    Returns:
        dict: The loaded settings as a dictionary.
    """
    if os.path.exists(json_file):
        with open(json_file, 'r') as file:
            return json.load(file)
    else:
        return {"runs": []}

def save_settings(settings, json_file):
    """
    Save the given settings dictionary to a JSON file.

    Parameters:
    - settings (dict): The settings dictionary to be saved.
    - json_file (str): The path to the JSON file where the settings will be saved.
    """
    with open(json_file, 'w') as file:
        json.dump(settings, file, indent=4)

def generate_gps_settings(gps_settings):
    """
    Generate GPS commands based on the given gpsSettings.
    """
    commands = [
        f"/gps/particle {gps_settings['particle']}",
        f"/gps/energy {gps_settings['energy']}",
        f"/gps/pos/type {gps_settings['posType']}",
        f"/gps/pos/centre {gps_settings['posCentre']}"
    ]
    
    if gps_settings['posType'] == "Volume":
        if 'posRadius' in gps_settings:
            commands.append(f"/gps/pos/radius {gps_settings['posRadius']}")
        if 'posHalfz' in gps_settings:
            commands.append(f"/gps/pos/halfz {gps_settings['posHalfz']}")
        if 'posConfine' in gps_settings:
            commands.append(f"/gps/pos/confine {gps_settings['posConfine']}")
    elif gps_settings['posType'] == "Point":
        if 'direction' in gps_settings:
            commands.append(f"/gps/direction {gps_settings['direction']}")
    
    if 'angType' in gps_settings:
        commands.append(f"/gps/ang/type {gps_settings['angType']}")
    
    return "\n".join(commands)

def generate_geometry_settings(geometry_settings):
    """
    Generate the geometry commands based on the provided geometry_settings.
    """
    commands = [
        f"/detector/setOuterCryostatRadius {geometry_settings['outerCryostatRadius']}",
        f"/detector/setOuterCryostatHeight {geometry_settings['outerCryostatHeight']}",
        f"/detector/setOuterCryostatWallThickness {geometry_settings['outerCryostatWallThickness']}",
        f"/detector/setInnerCryostatRadius {geometry_settings['innerCryostatRadius']}",
        f"/detector/setInnerCryostatHeight {geometry_settings['innerCryostatHeight']}",
        f"/detector/setInnerCryostatWallThickness {geometry_settings['innerCryostatWallThickness']}",
        f"/detector/setFiducialRadius {geometry_settings['fiducialRadius']}",
        f"/detector/setFiducialHeight {geometry_settings['fiducialHeight']}",
    ]
    return "\n".join(commands)

def generate_run_settings(run_settings, path_manager, job_id):
    """
    Generate the run-related commands based on the provided run_settings.
    """
    output_file_name = os.path.join(path_manager.output_dir, f"{run_settings['outputFileName']}_{job_id}")
    
    commands = [
        f"/run/setOutputFileName {output_file_name}",
        f"/run/setFastSimulation {run_settings['fastSimulation']}",
        f"/run/setNumberOfScatters {run_settings['numberOfScatters']}",
        f"/run/setMaxEnergy {run_settings.get('maxEnergy', '2 MeV')}",
        f"/run/printProgress {run_settings['printProgress']}",
    ]
    
    return "\n".join(commands)

def generate_run_control(beam_on, random_seed1, random_seed2):
    """
    Generate the run section commands for the macro file.
    
    Args:
        beam_on (int): Number of particles to simulate.
        random_seed1 (int): Seed for randomization.
        random_seed2 (int): Second seed for randomization.
    
    Returns:
        str: A string containing the run section commands.
    """
    commands = [
        f"/random/setSeeds {random_seed1} {random_seed2}",
        f"/run/beamOn {beam_on}"
    ]
    return "\n".join(commands)

def generate_mac_file(settings, path_manager, beam_on, random_seed1, job_id):
    """
    Generates a macro file for simulation based on the provided settings.

    Args:
        settings (dict): A dictionary containing simulation settings.
        path_manager (PathManager): Path management instance.
        beam_on (int): Number of beam particles to simulate.
        random_seed1 (int): Random seed for the simulation.
        job_id (int): Job identifier.

    Returns:
        str: The path to the generated macro file.
    """
    # Generate different sections of the macro file
    gps_commands = generate_gps_settings(settings["gpsSettings"])
    geometry_commands = generate_geometry_settings(settings["geometry_settings"])
    run_commands = generate_run_settings(settings["run_settings"], path_manager, job_id)
    run_section = generate_run_control(beam_on, random_seed1, random_seed1 + 1)
    
    # Combine all the sections into the final macro content
    mac_content = "\n".join([
        f"/control/verbose {settings['verbose']}",
        f"/run/verbose {settings['verbose']}",
        f"/tracking/verbose {settings['verbose']}",
        "# define the detector geometry",
        geometry_commands,
        "/run/initialize",
        gps_commands,
        run_commands,
        run_section
    ])

    # Write the generated content to a file
    mac_file = os.path.join(path_manager.jobs_dir, f"run_macro_{job_id}.mac")
    with open(mac_file, 'w') as file:
        file.write(mac_content)
    
    return mac_file

def submit_job(mac_file, path_manager, job_name="G4Job"):
    """
    Submits a Geant4 job to the batch queue using a job submission system (e.g., SLURM, PBS).
    Modify this function according to the specifics of your batch system.
    """
    # Move the macro file to the jobs directory
    # shutil.move(mac_file, path_manager.jobs_dir)

    # Create a condor submit file
    submit_file = os.path.join(path_manager.jobs_dir, f"submit_{os.path.basename(mac_file).replace('.mac', '.submit')}")
    log_file = os.path.join(path_manager.logs_dir, f"{job_name}.log")
    script_file = os.path.join(path_manager.jobs_dir, f"submit_{os.path.basename(mac_file).replace('.mac', '.sh')}")

    submit_content = f"""
executable = {script_file}
output = {log_file}
error = {log_file}
log = {log_file}

## Can use "el7", "el8", or "el9" for UseOS or you can specify your own
## SingularityImage but an OS must be specified and in string quotations.
+UseOS                  = "el9"
## This job can run up to 4 hours. Can choose "express", "short", "medium", or "long".
+JobCategory            = "short"
queue
    """
    with open(submit_file, 'w') as file:
        file.write(submit_content)

    # Create job script
    script_content = f"""#!/bin/bash
source /user/z37/.bashrc
conda activate g4
cd {path_manager.jobs_dir}
/user/z37/g4/G4FastSim/build/G4FastSim {mac_file}
"""
    with open(script_file, 'w') as file:
        file.write(script_content)

    # chmod +x the script file
    os.chmod(script_file, 0o755)
    
    # Submit the job
    subprocess.run(["condor_submit", submit_file])

def run_simulation(mac_file, path_manager):
    """
    Run the simulation using the specified macro file and path manager.

    Args:
        mac_file (str): The path to the macro file.
        path_manager (PathManager): An instance of the PathManager class.

    Returns:
        None
    """
    executable = os.path.join(path_manager.project_base_dir, "build", "G4FastSim")
    print(executable, mac_file)
    subprocess.run([executable, mac_file])

def parse_arguments():
    """
    Parse command line arguments for running Geant4 simulation.

    Returns:
        argparse.Namespace: Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description="Run Geant4 simulation with specified settings.")
    parser.add_argument("-json", "--json_file", required=True, help="Path to the JSON settings file.")
    parser.add_argument("-n", "--beam_on", type=int, required=True, help="Number of events to simulate.")
    parser.add_argument("-o", "--output_dir", default=None, help="Optional output directory. Defaults to a parameter-based directory in '../output'.")
    parser.add_argument("-rundb", "--rundb_file", default="rundb.json", help="Path to the run database file.")
    parser.add_argument("-jobs", "--num_jobs", type=int, default=1, help="Number of jobs to submit.")
    parser.add_argument("--batch", action="store_true", help="Submit jobs to batch queue.")
    parser.add_argument("--base_dir", default="/user/z37/g4/G4FastSim", help="Base directory of the project.")
    return parser.parse_args()

def initialize_paths(args):
    """
    Initializes the paths for the simulation.

    Args:
        args (Namespace): The command line arguments.

    Returns:
        PathManager: An instance of the PathManager class.
    """
    return PathManager(args.base_dir, args.output_dir)

def prepare_settings(args, path_manager):
    """
    Prepare the simulation settings.

    Args:
        args (object): The command line arguments.
        path_manager (object): The path manager object.

    Returns:
        dict: The prepared simulation settings.
    """
    settings = load_settings(args.json_file)
    random_seed = random.randint(0, 1000000)
    settings['randomSeed'] = random_seed
    settings['beamOn'] = args.beam_on
    settings_file = os.path.join(path_manager.output_dir, "settings.json")

    with open(settings_file, 'w') as json_file:
        json.dump(settings, json_file, indent=4)
    return settings

def execute_jobs(args, settings, path_manager):
    """
    Executes multiple simulation jobs.

    Args:
        args (object): Command-line arguments.
        settings (dict): Simulation settings.
        path_manager (object): Path manager object.

    Returns:
        None
    """
    for job_id in range(args.num_jobs):
        mac_file = generate_mac_file(settings, path_manager, args.beam_on // args.num_jobs, settings["randomSeed"] + job_id * 10, job_id)
        if args.batch:
            submit_job(mac_file, path_manager, f"job_{job_id}")
        else:
            run_simulation(mac_file, path_manager)

def update_master_rundb(rundb, settings, path_manager, args):
    """
    Update the master runs database with a new run based on the provided settings.

    Args:
        rundb (dict): The master run database dictionary.
        settings (dict): The settings dictionary containing GPS settings and other simulation parameters.
        path_manager (PathManager): An instance of the PathManager class.
        args (argparse.Namespace): The command-line arguments.

    Returns:
        None

    """
    gps_settings = settings["gpsSettings"]
    run_settings = settings["run_settings"]
    new_run = {
        "id": f"run_{len(rundb['runs']) + 1:02}",
        "particle": gps_settings["particle"],
        "energy": gps_settings["energy"],
        "fastSimulation": run_settings["fastSimulation"],
        "maxScatters": run_settings.get("numberOfScatters", 1),
        "maxEnergy": run_settings.get("maxEnergy", "2 MeV"),
        "sourceVolume": gps_settings["posConfine"] if gps_settings["posType"] == "Volume" else "",
        "outputDir": path_manager.output_dir,
        "outputFile": run_settings['outputFileName'],
        "numEvents": args.beam_on,
        "numJobs": args.num_jobs,
        "randomSeed": settings["randomSeed"],
        "settingsFile": "settings.json",
        "status": "active"
    }
    rundb['runs'].append(new_run)
    save_settings(rundb, args.rundb_file)

def main():
    """
    Main function that executes the simulation.

    This function performs the following steps:
    1. Parses the command line arguments.
    2. Loads the settings from the specified configuration file.
    3. Initializes the path manager.
    4. Prepares the simulation settings.
    5. Executes the simulation jobs.
    6. Updates the master configuration file.

    Note: This function assumes that the necessary helper functions are defined elsewhere.

    Args:
        None

    Returns:
        None
    """
    # 1. Parse command line arguments
    args = parse_arguments()
    # 2. Load the master configuration file
    rundb = load_settings(args.rundb_file)
    # 3. Initialize the path manager
    path_manager = initialize_paths(args)
    # 4. Prepare the simulation settings
    settings = prepare_settings(args, path_manager)
    # 5. Execute the simulation jobs
    execute_jobs(args, settings, path_manager)
    # 6. Update the master configuration file
    update_master_rundb(rundb, settings, path_manager, args)

if __name__ == "__main__":
    main()
