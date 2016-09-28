import sys
# Modify this to the required path
sys.path.insert(1, '/Users/daniel/Gits/dgas_psi4/objdir/stage/usr/local/lib/')
import psi4core

# Init psi4core
psi4core.initialize()
psi4core.set_memory(int(512e6)) # Set to 512 MB

# Cleanup psi4core at exit
import atexit
atexit.register(psi4core.set_legacy_molecule, None)
atexit.register(psi4core.finalize)
