# get the runtime arguments
import sys
inputfile = sys.argv[1]
outputfile = sys.argv[2]

# load pymatgen
from pymatgen.io.cif import CifParser

# read structure from input and write the primitive cell
CifParser(inputfile).get_structures()[0].get_primitive_structure().to(filename=outputfile)
