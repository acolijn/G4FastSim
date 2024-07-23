import json
import os
import subprocess
import datetime
import shutil
import argparse
import random


def load_settings(json_file):
    if os.path.exists(json_file):
        with open(json_file, 'r') as file:
            return json.load(file)
    else:
        return {"runs": []}

def save_settings(settings, json_file):
    with open(json_file, 'w') as file:
        json.dump(settings, file, indent=4)

def generate_gps_commands(gps_settings):
    commands = []
    commands.append(f"/gps/particle {gps_settings['particle']}")
    commands.append(f"/gps/energy {gps_settings['energy']}")
    commands.append(f"/gps/pos/type {gps_settings['posType']}")
    if 'posShape' in gps_settings:
        commands.append(f"/gps/pos/shape {gps_settings['posShape']}")
    commands.append(f"/gps/pos/centre {gps_settings['posCentre']}")
    if gps_settings['posType'] == "Volume":
        if 'posRadius' in gps_settings:
            commands.append(f"/gps/pos/radius {gps_settings['posRadius']}")
        if 'posHalfz' in gps_settings:
            commands.append(f"/gps/pos/halfz {gps_settings['posHalfz']}")
        if 'posConfine' in gps_settings:
            commands.append(f"/gps/pos/confine {gps_settings['posConfine']}")
        commands.append(f"/gps/ang/type {gps_settings['angType']}")
    elif gps_settings['posType'] == "Point":
        if 'direction' in gps_settings:
            commands.append(f"/gps/direction {gps_settings['direction']}")
    return "\n".join(commands)

def generate_mac_file(settings, mac_template, output_dir, output_filename, beam_on, random_seed1=12345):
    with open(mac_template, 'r') as template_file:
        template = template_file.read()

    gps_commands = generate_gps_commands(settings["gpsSettings"])


    # Generate random seeds
    random_seed2 = random_seed1 + 1
    # Substitute placeholders in the template with actual settings
    mac_content = template.format(
        verbose=settings["verbose"],
        outerCryostatRadius=settings["outerCryostatRadius"],
        outerCryostatHeight=settings["outerCryostatHeight"],
        outerCryostatWallThickness=settings["outerCryostatWallThickness"],
        innerCryostatRadius=settings["innerCryostatRadius"],
        innerCryostatHeight=settings["innerCryostatHeight"],
        innerCryostatWallThickness=settings["innerCryostatWallThickness"],
        fiducialRadius=settings["fiducialRadius"],
        fiducialHeight=settings["fiducialHeight"],
        gpsCommands=gps_commands,
        outputFileName=os.path.join(output_dir, output_filename),
        fastSimulation=settings["fastSimulation"],
        maxEnergy=settings.get("maxEnergy", "2 MeV"),  # Use default if not provided
        numEvents=beam_on,
        numberOfScatters=settings.get("maxScatters", 1),
        printProgress=settings["printProgress"],
        beamOn=beam_on,
        randomSeed1=random_seed1,
        randomSeed2=random_seed2
    )

    mac_file = os.path.join(output_dir, "run_macro.mac")
    with open(mac_file, 'w') as file:
        file.write(mac_content)
    
    # Create README.txt with re-run instructions
    readme_content = f"To re-run this simulation, use the following command:\n\n./G4FastSim {mac_file}\n"
    readme_file = os.path.join(output_dir, "README.txt")
    with open(readme_file, 'w') as file:
        file.write(readme_content)
    
    return mac_file

def run_simulation(mac_file):
    executable = os.path.join(os.path.dirname(__file__), "../build", "G4FastSim")
    print(executable, mac_file)
    subprocess.run([executable, mac_file])

def get_run_by_id(config, run_id):
    for run in config['runs']:
        if run['id'] == run_id:
            return run
    return None

def regenerate_mac_from_settings(output_dir):
    settings_file = os.path.join(output_dir, "settings.json")
    if os.path.exists(settings_file):
        settings = load_settings(settings_file)
        mac_template = os.path.join(os.path.dirname(__file__), "../macros/run_macro_template.mac")
        beam_on = settings['beamOn']
        output_filename = settings['outputFileName']
        random_seed = settings['randomSeed']
        mac_file = generate_mac_file(settings, mac_template, output_dir, output_filename, beam_on, random_seed)
        return mac_file
    else:
        raise FileNotFoundError(f"Settings file not found: {settings_file}")

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run Geant4 simulation with specified settings.")
    parser.add_argument("-json", "--json_file", help="Path to the JSON settings file.")
    parser.add_argument("-n", "--beam_on", type=int, help="Number of events to simulate.")
    parser.add_argument("-o", "--output_dir", default=None, help="Optional output directory. Defaults to a parameter-based directory in '../output'.")
    parser.add_argument("-config", "--config_file", default="config.json", help="Path to the master configuration file.")
    parser.add_argument("-id", "--run_id", help="ID of the run to re-execute.")
    args = parser.parse_args()

    if args.run_id is None and (args.json_file is None or args.beam_on is None):
        parser.error("Either --run_id or both --json_file and --beam_on must be specified.")

    # Load master configuration file
    config_file = args.config_file
    config = load_settings(config_file)

    # Paths
    base_dir = os.path.dirname(os.path.abspath(__file__))
    mac_template = os.path.join(base_dir, "../macros/run_macro_template.mac")
    output_base_dir = "/data/xenon/acolijn/G4FastSim/"

    if args.run_id:
        # Re-run based on the run ID
        run = get_run_by_id(config, args.run_id)
        if run is None:
            print(f"No run found with ID {args.run_id}")
            return

        output_dir = run["outputDir"]
        mac_file = regenerate_mac_from_settings(output_dir)
        run_simulation(mac_file)
    else:
        # Load settings
        settings = load_settings(args.json_file)

        # Add the beamOn field to the settings
        settings['beamOn'] = args.beam_on

        # Generate a dynamic output directory based on key parameters
        gps_settings = settings["gpsSettings"]
        particle = gps_settings["particle"]
        energy = gps_settings["energy"].replace(" ", "")  # Remove spaces
        simulation_type = settings["fastSimulation"]

        # Generate a timestamp
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        # Append timestamp to output directory
        if args.output_dir:
            output_dir = os.path.join(args.output_dir, timestamp)
        else:
            #output_dir_name = f"{timestamp}_{particle}_{energy}_{simulation_type}"
            output_dir = os.path.join(output_base_dir, timestamp)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Generate a dynamic output filename
        output_filename = f"{settings['outputFileName']}"

        # Save the updated settings to a new JSON file in the output directory
        random_seed = random.randint(1, 1000000)
        settings['randomSeed'] = random_seed
        with open(os.path.join(output_dir, "settings.json"), 'w') as json_file:
            json.dump(settings, json_file, indent=4)

        # Generate .mac file
        mac_file = generate_mac_file(settings, mac_template, output_dir, output_filename, args.beam_on, random_seed)

        # Run simulation
        run_simulation(mac_file)

        # Debugging: Print new run information before adding to config
        new_run = {
            "id": f"run_{len(config['runs']) + 1:02}",
            "particle": particle,
            "energy": energy,
            "fastSimulation": simulation_type,
            "maxScatters": settings.get("numberOfScatters", 1),
            "maxEnergy": settings.get("maxEnergy", "2 MeV"),
            "sourceVolume": gps_settings["posConfine"] if gps_settings["posType"] == "Volume" else "",
            "outputDir": output_dir,
            "outputFile": output_filename,
            "settingsFile": "settings.json",	
            "randomSeed": random_seed,
            "numEvents": args.beam_on,
        }
        print("Adding new run to config:", new_run)

        # Update the master configuration file
        config['runs'].append(new_run)
        save_settings(config, config_file)

if __name__ == "__main__":
    main()
