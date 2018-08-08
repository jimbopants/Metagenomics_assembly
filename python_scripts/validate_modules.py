#!/software/anaconda3.6/bin/python
"""
When executed, this script tries to load all of the modules specififed in
in my assembly scripts and prints a console message for each module.
I have no idea if Quest frequently changes module availability,
but this could be useful to verify that everything works before sending jobsself.

It also checks the locations of shared binaries used by our lab (quast, spades etc.)
to make sure they all still exist at the anticipated locations.
"""

# Imports
import subprocess
import os
import shared_utilities as ll

def validation_procedure():
    params = ll.read_config()
    failures = 0
    failed_modules = []
    failures, failed_modules = check_modules(params, failures, failed_modules)
    failures, failed_modules = check_binaries(params, failures, failed_modules)
    recap_failures(failures, failed_modules)

def check_modules(params, failures, failed_modules):
    """ Validate module paths
    """
    print ("\n...Checking Available Modules...\n")
    for module in params["modules"]:
        statement = "module load {}".format(module)
        print("Loading {}......".format(module))
        err = subprocess.check_output([statement], stderr=subprocess.STDOUT, shell=True)
        if err:
            print(err)
            failures += 1
            failed_modules.append(module)
            print("Can't find module {}".format(module))
    return failures, failed_modules

def check_binaries(params, failures, failed_modules):
    print ("\n\n...Checking Binary paths...")
    for binary in binary_paths:
        print("Checking {}......".format(binary))
        if os.path.exists(binary):
            pass
        else:
            failures += 1
            failed_modules.append(binary)
            print("Can't find binary {}, did it move?\n".format(binary))
    return failures, failed_modules

def recap_failures(failures, failed_modules):
    """ Print results
    """
    if failures == 0:
        print("\n\nAll Modules/binaries Available\n")
    else:
        print("\n\nSome Modules/binaries Missing")
        print("Total failed modules: {}".format(failures))
        for i in failed_modules:
            print("Missing: {}".format(i))

if __name__ == __main__:
    validation_procedure()
