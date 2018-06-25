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

# Modules:
module_paths = ["idba/2016_12_longread",
                "bowtie2/2.2.6",
                "samtools/1.2",
                "megahit/1.0.6.1",
                "python/anaconda3.6",
                "boost/1.56.0"
                ]

# Binaries:
binary_paths = [
"/projects/b1052/shared/quast-4.6.3/quast.py",
"/projects/b1052/shared/SPAdes-3.12.0-Linux/bin"
]

# Check Modules
failures = 0
failed_modules = []
print ("\n...Checking Available Modules...\n")
for module in module_paths:
    statement = "module load {}".format(module)
    print("Loading {}......".format(module))

    err = subprocess.check_output([statement], stderr=subprocess.STDOUT, shell=True)
    if err:
        print(err)
        failures += 1
        failed_modules.append(module)
        print("Can't find module {}".format(module))


# Check binary paths:
print ("\n\n...Checking Binary paths...")
for binary in binary_paths:
    print("Checking {}......".format(binary))
    if os.path.exists(binary):
        pass
    else:
        failures += 1
        failed_modules.append(binary)
        print("Can't find binary {}, did it move?\n".format(binary))

# Recap failures:
if failures == 0:
    print("\n\nAll Modules/binaries Available\n")
else:
    print("\n\nSome Modules/binaries Missing")
    print("Total failed modules: {}".format(failures))
    for i in failed_modules:
        print("Missing: {}".format(i))
