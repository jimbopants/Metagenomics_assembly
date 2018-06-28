# Written by Jim Griffin
# 6/20/18

# Imports:
import subprocess
import os
import argparse
import sys

# Cleanup intermediate files
# assume directories look like:
# assembled /
#   megahit /
#       S1...
#   idba /
#       S1 /
#           {crap}, contigs, scaffolds, log
