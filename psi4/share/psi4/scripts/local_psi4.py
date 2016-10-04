#!/usr/bin/env python
import sys
import os
import argparse

parser = argparse.ArgumentParser(description="A hybrid C++/Python quantum chemistry module.")
parser.add_argument("-i", "--input", default="input.dat", help="Input file name. Default input.dat.")
parser.add_argument("-o", "--output", help="Redirect output elsewhere.\n"
                                           "Default filename.out if input is filename"
                                           "filename.out if input is filename.in\n"
                                           "output.dat if input is input.dat\n")

parser.add_argument("-v", "--verbose", action='store_true', help="Print a lot of information.")
parser.add_argument("-V", "--version", action='store_true', help="Print version information.")

parser.add_argument("-d", "--debug", action='store_true', help="Flush the outfile at every print statement.")
parser.add_argument("-k", "--skip-preprocessor", action='store_true', help="Skips input preprocessing. Expert mode.")
parser.add_argument("-m", "--messy", action='store_true', help="Leave temporary files after the run is completed.")
parser.add_argument("-r", "--restart", action='store_true', help="Number to be used instead of process id.")

parser.add_argument("-s", "--scratch", help="Psi4 scratch directory to use.")
parser.add_argument("-a", "--append", help="Append results to output file. Default Truncate first")
parser.add_argument("-l", "--psidatadir", help="Specify where to look for the Psi data directory. Overrides PSIDATADIR.")
parser.add_argument("-n", "--nthread", default=1, help="Number of threads to use (overrides OMP_NUM_THREADS)")
parser.add_argument("-p", "--prefix", help="Prefix name for psi files. Default psi")

# For plugins
parser.add_argument("--new-plugin", help="Creates a new directory with files for writing a "
                                         "new plugin. You can specify an additional argument "
                                         "that specifies a template to use, for example "
                                         "--new-plugin name +mointegrals")
parser.add_argument("--new-plugin-makefile", help="Creates Makefile that can be used to compile"
                                                  "plugins. The Makefile is placed in the current"
                                                  "directory.")

# print("Environment Variables\n");
# print("     PSI_SCRATCH           Directory where scratch files are written.")
# print("                           Default: $TMPDIR (or /tmp/ when not set)")
# print("                           This should be a local, not network, disk")

#parser.print_help()
args, unknown = parser.parse_known_args()
args = args.__dict__ # Namespace object seems silly

# Insert the python path

cmake_install_prefix = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + '..'
#cmake_libdir =  "@CMAKE_INSTALL_LIBDIR@" # CMAKE_INSTALL_LIBDIR
#cmake_datadir =  "@CMAKE_INSTALL_DATADIR@" # CMAKE_INSTALL_DATADIR
#if "CMAKE_INSTALL_LIBDIR" in cmake_libdir:
#    raise ImportError("Psi4 was not installed correctly!")

# Bake in default PSIDATADIR
#data_dir = cmake_install_prefix + os.path.sep + cmake_datadir + os.path.sep + "psi4"
data_dir = "/Users/daniel/Gits/dgas_psi4/psi4/share/psi4"


# Replace input/output if unknown kwargs
if len(unknown) > 0:
    args["input"] = unknown[0]
elif len(unknown) > 1:
    args["output"] = unknown[1]
elif len(unknown) > 2:
    raise KeyError("Too many unknown arguments: %s" % str(unknown))

# Figure out output arg
if args["output"] is None:
    if args["input"] == "input.dat":
        args["output"] = "output.dat"
    elif args["input"].endswith(".in"):
        args["output"] = args["input"].replace(".in", ".out")
    else:
        args["output"] = args["input"] + ".dat"

if not os.path.isfile(args["input"]):
    raise KeyError("The file %s does not exist." % args["input"])


# Figure out psidata dir
if "PSIDATADIR" in os.environ.keys():
    data_dir = os.environ["PSIDATADIR"]

if args["psidatadir"] is not None:
    data_dir = os.path.abspath(args["psidatadir"])

if not os.path.isdir(data_dir):
    raise Exception("PSIDATADIR (%s) is not a directory" % data_dir)

os.environ["PSIDATADIR"] = data_dir


### Actually import psi4 and apply setup ###

# Import installed psi4
#lib_path = cmake_install_prefix + os.path.sep + cmake_libdir
lib_path = "/Users/daniel/Gits/dgas_psi4"
sys.path.insert(1, lib_path)
import psi4

psi4.psi4core.set_environment("PSIDATADIR", data_dir)

# Read input
with open(args["input"]) as f:
    content = f.read()

# Preprocess
if not args["skip_preprocessor"]:
    content = psi4.process_input(content, psi4_imported=False)

# Handle Verbose
if args["verbose"]:
    print('-' * 74)
    print('Parsed Psithon:')
    print(content)
    print('-' * 74)

# Handle Messy
if args["messy"]:
    import atexit
    for handler in atexit._exithandlers:
        if handler == psi4.psi4core.clean:
            atexit._exithandlers.remove(handler)


# Run the program!
exec(content)


