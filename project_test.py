import subprocess
import argparse
import multiprocessing
import platform
import os

if platform.system() == "Windows":
    PYTHON_CMD = "python"
else:
    PYTHON_CMD = "python3"

def run_command(command, cwd=None, requires_linux=False):
    if platform.system() == "Windows" and requires_linux==True:
        return subprocess.run(["wsl"] + command, cwd=cwd, check=True)
    else:
        return subprocess.run(command, cwd=cwd, check=True)

def generate_execution_command(app, params):
    return ["mpirun", "-n", "1", "--map-by", f"node:PE={multiprocessing.cpu_count() // 2}", "--bind-to", "core", f"./build/{app}", params]

###################################################################################################
def run_tests_PhdUtility():
    run_command(["make", "test", "-j"], cwd="PhdUtility", requires_linux=True)

def install_PhdUtility():
    run_command(["make", "install"], cwd="PhdUtility", requires_linux=True)

###################################################################################################
def build_FermionCommute():
    run_command(["make"], cwd="cpp/FermionCommute", requires_linux=True)

###################################################################################################
def build_continuum():
    run_command(["make", "-j"], cwd="cpp/ContinuumSystem", requires_linux=True)

def run_continuum(i):
    if i == 0:
        param_name = "no_coulomb.config"
    elif i == 1:
        param_name = "normal_screening.config"
    elif i == 2:
        param_name = "weak_screening.config"
    elif i == 3:
        param_name = "strong_attraction.config"
    else:
        raise ValueError("Continuum execution: Invalid index")
    
    run_command(generate_execution_command("ContinuumSystem", f"test_params/{param_name}"), cwd="cpp/ContinuumSystem", requires_linux=True)

###################################################################################################
def build_hubbard():
    run_command(["make", "-j"], cwd="cpp/Hubbard", requires_linux=True)

def run_hubbard(i):
    if i == 0:
        param_name = "sc_cdw.config"
    elif i == 1:
        param_name = "sc.config"
    elif i == 2:
        param_name = "cdw.config"
    elif i == 3:
        param_name = "afm.config"
    else:
        raise ValueError("Hubbard execution: Invalid index")
    
    run_command(generate_execution_command("HubbardMeanField", f"test_params/{param_name}"), cwd="cpp/Hubbard", requires_linux=True)

###################################################################################################
def manual_data_check(app, i=None):
    if i is not None:
        run_command([PYTHON_CMD, os.path.join("phd_plot_scripts", "tests", f"{app}.py"), f"{i}"])
    else:
        run_command([PYTHON_CMD, os.path.join("phd_plot_scripts", "tests", f"{app}.py")])
    is_good = input("Does the data look good? (yes/no): ")
    
    if is_good != "yes" and is_good != "y" and is_good != "Yes" and is_good != "Y":
        print("Data result is bad.")
        raise SystemExit

###################################################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Project Test Orchestration")
    parser.add_argument("--step", 
                        choices=['all', 'PhdUtility', 'Hubbard', 'Continuum'], 
                        default='all',
                        help="Select which step to run")
    args = parser.parse_args()
    
    if args.step == 'all':
        run_tests_PhdUtility()
        install_PhdUtility()
        build_FermionCommute()
        
        build_hubbard()
        for i in range(4):
            run_hubbard(i)    
            manual_data_check("hubbard", i)
        
        build_continuum()
        for i in range(4):
            run_continuum(i)
            manual_data_check("continuum", i)
    elif args.step == 'PhdUtility':
        run_tests_PhdUtility()
    elif args.step == 'Hubbard':
        build_hubbard()
        for i in range(4):
            run_hubbard(i)    
            manual_data_check("hubbard", i)
    elif args.step == 'Continuum':
        build_continuum()
        for i in range(4):
            run_continuum(i)
            manual_data_check("continuum", i)

    print("Testing successful!")