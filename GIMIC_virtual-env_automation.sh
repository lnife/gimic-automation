#!/bin/bash
# GIMIC Calculation Setup Script
# This script automates the setup and execution of Gauge-Including Magnetically Induced Current calculations

# Define path to the existing environment's activate script
# CUSTOMIZE THIS: Update with your GIMIC environment location
ACTIVATE_SCRIPT="/path/to/gimic_env/bin/activate"

echo

# Step 1: Check if the environment exists and activate it
if [ -f "$ACTIVATE_SCRIPT" ]; then
    # Display activation message with green color
    echo -e "\e[38;2;1;144;83mActivating gimic virtual environment...\e[0m"
    # Source the environment to activate it
    source "$ACTIVATE_SCRIPT"
else
    # Display error message if environment not found
    echo -e "\e[3;38;2;104;219;187mActivation script not found at $ACTIVATE_SCRIPT\e[0m"
    # Exit the script with error status
    exit 1
fi

# Step 2: Initiate GIMIC Calculations
# Display section header with custom color
echo -e "\e[3;38;2;139;96;87mInitiating GIMIC Calculations\n---------------------------------------------\n\e[0m"

# Create Gaussian2gimic.py in current directory
# This Python script converts Gaussian output to GIMIC input format
cat << 'EOFX' > Gaussian2gimic.py  #Don't alter this part unless you want use turbo2gimic.py or lsdalton2gimic or etc. then paste code of those python files from next line to EOFX (in this code it's line 530) 
#! /usr/bin/env python




import BasisSet
import os.path
from optparse import OptionParser
import numpy
import math
import re

Table=['Xx', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn']

class FCHK (object) :

    def __init__ (self, fchkname='gaussian.fchk'):
        self.filename = fchkname

    def get_property(self, property):
        """
        Read a property name "property" in the fchk file
        return a numpy vector
        """
        f = open(self.filename, 'r')
        lines = f.readlines()
        f.close()
        # find the block "property"
        first = 0
        last = 0
        cols = {"R":5, "I":6}
        numpytype = {"R":"d", "I":"int"}
        stype = "R"
        size = 0
        for i, l in enumerate(lines):
            if property in l:
                l = l[len(property):]
                if "N=" in l:
                    # read array of data
                    pattern = r"\d+"
                    size = int(re.findall(pattern, l)[0])
                    typepattern = r"   [IR]   "
                    try:
                        stype = re.findall(typepattern, l)[0].strip()
                    except IndexError:
                        raise Exception("The type of data in fchk is not recognized")
                    first = i+1
                    nlines = (size-1)//cols[stype]+1
                    last = first+nlines
                    break
                else:
                    # single data
                    typepattern = r"   [IR]   "
                    try:
                        stype = re.findall(typepattern, l)[0].strip()
                    except IndexError:
                        raise Exception("The type of data in fchk is not recognized")
                    if stype == 'I':
                        pattern = r"\d+"
                        return int(re.findall(pattern, l)[0])
                    elif stype == 'R':
                        pattern = r"-?\d+\.\d*[ED]?[+\-]?\d*"
                        return float(re.findall(pattern, l)[0])
        lines = lines[first:last]
        if len(lines) == 0:
            return None
        # read the data
        data = []
        for line in lines:
            data.extend([float (fl) for fl in line.split()])
        data = numpy.array(data, numpytype[stype])
        if data.size != size:
            raise Exception("The number of data recovered [%i] is not the same as the size written in the file [%i]"%(data.size, size))
        return data

    def get_density(self):
        """
        Get the non-perturbed density
        """
        data = self.get_property("Total SCF Density")
        density = None
        if isinstance(data, numpy.ndarray):
            nbasis = self.get_property("Number of basis functions")
            index = 0
            density = numpy.zeros((nbasis, nbasis))
            for i in range(nbasis):
                for j in range(i+1):
                    density[i, j] = data[index]
                    index += 1

            # fill the other half of the density matrix
            density = density+density.transpose()-numpy.diag(density.diagonal())
        return density

    def get_spindensity(self):
        """
        Get the spin density
        """
        data = self.get_property("Spin SCF Density")
        density = None
        if isinstance(data, numpy.ndarray):
            nbasis = self.get_property("Number of basis functions")
            index = 0
            density = numpy.zeros((nbasis, nbasis))
            for i in range(nbasis):
                for j in range(i+1):
                    density[i, j] = data[index]
                    index += 1

            # fill the other half of the density matrix
            density = density+density.transpose()-numpy.diag(density.diagonal())
        return density

    def get_density_magnetic(self):
        """
        Get the perturbed density by magnetic field
        """
        data = self.get_property("Magnetic Field P1 (GIAO)")
        density = None
        if isinstance(data, numpy.ndarray):
            nbasis =  self.get_property("Number of basis functions")
            index = 0
            density = numpy.zeros((3, nbasis, nbasis))
            for idir in range(3):
                for i in range(nbasis):
                    for j in range(i+1):
                        density[idir, i, j] = data[index]
                        index += 1

                # fill the other half of the density matrix
                density[idir] = density[idir]-density[idir].transpose()
        return density

    def get_basisset(self):
        """
        read the basis set
        return a list of AO objects
        """
        nshell = self.get_property("Number of contracted shells")
        shelltype = self.get_property("Shell types")
        primitive_pershell = self.get_property("Number of primitives per shell")
        shell2atom = self.get_property("Shell to atom map")
        primitive_exp = self.get_property("Primitive exponents")
        primitive_coeff = self.get_property("Contraction coefficients")
        primitive_coeff_p = self.get_property("P(S=P) Contraction coefficients")
        atomicnumbers = self.get_property("Atomic numbers")
        # basis set list
        basisset = []
        start = 0
        for ishell in range(nshell):
            nprim = primitive_pershell[ishell]
            # get the exponents and coeffs for all the primitive of the shell
            exponents = primitive_exp[start:start+nprim]
            coefficients = primitive_coeff[start:start+nprim]
            coefficients_p = None
            if primitive_coeff_p is not None:
                if numpy.nonzero(primitive_coeff_p[start:start+nprim])[0].size !=0 :
                    coefficients_p = primitive_coeff_p[start:start+nprim]
            iatom = shell2atom[ishell]
            an = int(atomicnumbers[iatom-1])
            type = int(shelltype[ishell])
            shell = BasisSet.SHELL(iatom, an, type, exponents=exponents, coefficients=coefficients, coefficients_p=coefficients_p)
            basisset.append(shell)
            start = start + nprim
        if len(basisset) == 0:
            basisset =  None
        else:
            basisset =  BasisSet.BasisSet(basisset)
#            basisset.check_contractedcoeff()
        return basisset

class Gaussian (object) :

    def __init__ (self, logname='gaussian.log',w=None) :
        self.filename = logname

    def get_density(self):
        """
        Get the density matrix for alpha and beta electrons
        Need pop=regular or full keyword
        """
        f = open(self.filename, 'r')
        lines = f.readlines()
        f.close()
        pattern = r"\w* Density Matrix"
        size = 0
        densities = []
        for iline, line in enumerate(lines):
            if "basis functions," in line:
                size = int(re.findall(r"\d+", line)[0])
            elif size > 0 and re.search(pattern, line) is not None:
                ntab = (size-1) // 5 + 1
                nlines = (size+1) * ntab - 5*ntab*(ntab-1)//2
                densities.append(self.read_integrals(string=lines[iline+1:iline+1+nlines], size=size, symmetric=True))
        if len(densities) == 0:
            return None
        densities = numpy.array(densities)
        if densities.shape[0] == 1:
            densities = densities[0]
        return densities

    def get_density_magnetic(self):
        """
        Get the perturbed density by magnetic field for alpha and beta electrons
        Need IOP(10/33=2) keyword
        """
        f = open(self.filename, 'r')
        lines = f.readlines()
        f.close()
        pattern = r"P1 \w+ \(asymm, AO basis\)"
        size = 0
        densities = []
        for iline, line in enumerate(lines):
            if "basis functions," in line:
                size=int(re.findall(r"\d+", line)[0])
            elif size > 0 and re.search(pattern, line) is not None:
                ntab = (size-1) // 5 + 1
                nlines = (size+1) * ntab - 5*ntab*(ntab-1)//2
                densities.append(self.read_integrals(string=lines[iline+1:iline+1+nlines], size=size, symmetric=False))
        if len(densities) == 0:
            return None
        return numpy.array(densities)


    def read_integrals(self, string, size, symmetric=None):
        """
        Get integrals of the form:
                1             2             3             4             5
      1  0.100000D+01
      2  0.219059D+00  0.100000D+01
      3  0.000000D+00  0.000000D+00  0.100000D+01
      4  0.000000D+00  0.000000D+00  0.000000D+00  0.100000D+01
      5  0.000000D+00  0.000000D+00  0.000000D+00  0.000000D+00  0.100000D+01
      6  0.184261D+00  0.812273D+00  0.000000D+00  0.000000D+00  0.000000D+00
        """
        # read the data (triangular matrix divided into different block)
        integral = numpy.zeros((size,size))
        iblock = 0
        while len(string) > 0:
            end = size+1-iblock*5 #skip the first line that are the label of the  columns
            block = string[1:end]
            string = string[end:]
            for irow, blockline in enumerate(block):
                pattern = r"-?\d+\.\d*[ED]?[+\-]?\d*"
                alist = re.findall(pattern, blockline)
                line=[float(fl.replace('D','E')) for fl in alist]
                integral[irow+iblock*5, iblock*5:iblock*5+len(line)] = line
            iblock = iblock+1

        if symmetric is None:
            return integral
        elif symmetric:
            # fill the other half of the matrix if symmetric = True
            integral = integral+integral.transpose()-numpy.diag(integral.diagonal())
        else:
            # fill the other half of the matrix if antisymmetric (symmetric=False)
            integral = integral-integral.transpose()

        return integral

    def get_basisset(self):
        """
        read the basis set given by gfprint in log file
        return a dictionary with the basis set information of each atom
        The key is the label
        The value is a list made of tuple for each shell
        Shell[0] is the type of orbital
        Shell[1] is a numpy array with coeff and exp for each primitive
                 One line per primitive with an exponent and the coefficients
        """
        f = open(self.filename, 'r')
        lines = f.readlines()
        f.close()
        start = 0
        end = 0
        # get the block with the basis set informations
        for i, l in enumerate(lines):
            if "AO basis set" in l:
                start=i+1
            elif (start != 0) and ((l.strip().split()[0] != "Atom") and (l.strip().split()[0][:2] != "0.")) : # End of the basis set
                end = i
                break
        lines = lines[start:end]
        if len(lines) == 0:
            return None
        basisset = []
        while len(lines) > 0 :
            # line should look like:
            # Atom C1       Shell     1 S   6     bf    1 -     1          1.000000000000          2.000000000000          8.000000000000
            info = lines[0].split()
            label = info[1]
            # get iatom and an from label
            iatom = int(re.findall(r"\d+",label)[0])
            symbol = re.findall(r"[a-zA-Z]+",label)[0]
            an = Table.index(symbol)
            nprim = int(info[5])
            mult = int(info[9]) - int(info[7])+1
            typeoforbital = BasisSet.mult2type[mult]
            coord = [float(x) for x in re.findall(r"[+-]?\d+\.\d*", lines[0])]
            coord = numpy.array(coord)
            # Define a block that correspond to a shell
            block = lines[1:1+nprim]
            lines = lines[1+nprim:]
            # read the block of coeff and exp
            data = []
            for blockline in block:
                aline = [float(fl.replace('D', 'E')) for fl in blockline.split()]
                data.append(aline)
            data = numpy.array(data)
            exponents = data[:, 0]
            coefficients = data[:, 1]
            coefficients_p = None
            if data.shape[-1] == 3:
                coefficients_p = data[:, 2]
            shell = BasisSet.SHELL(iatom, an, typeoforbital, exponents=exponents, coefficients=coefficients, coefficients_p=coefficients_p, coord=coord)
            basisset.append(shell)

        if len(basisset) == 0:
            basisset =  None
        else:
            basisset = BasisSet.BasisSet(basisset)
        return basisset

usage = "usage: %prog [options] "
parser = OptionParser(usage=usage)

parser.add_option("-i", "--inputfile", dest="inputfile",
                  help="Input log or fchk file from an NMR calculation. For log file, you need to have specified the following keywords: gfprint pop=regular iop(10/33=2). gfprint: to print the database. pop=regular or pop=full: to print the density matrix. iop(10/33=2) to print the perturbed density matrices. NOTE: for OPENSHELL systems, only the log file can be used (with the three keywords above).",
                  action="store", type="string",default="none")

parser.add_option("-t", "--turbomole", dest="turbomole",
                  help="Arrange the XDENS like with Turbomole instead of like with CFOUR. Default: False (like CFOUR)",
                  action="store_true", default=False)

parser.add_option("", "--XDENS", dest="XDENS",
                  help="XDENS file to compare with the generated one. Default: none",
                  action="store", type="string", default="none")

(options, args) = parser.parse_args()

#check that the user introduce  one input filename
if options.inputfile == "none":
    parser.error("No inputfile given")

if not os.path.isfile(options.inputfile):
    parser.error("The file ['%s'] does not exist" % options.inputfile)



# open the fchk or log file
if os.path.splitext(options.inputfile)[1] == ".fchk":
    file=FCHK(options.inputfile)
    # Read the density matrices
    # Only for closeshell. For opensell, use log file
    density = file.get_density()
    if density is None:
        raise Exception("Error while reading the density matrix.")
    spindensity = file.get_spindensity()
    if spindensity is not None:
        parser.error("Read the log file for openshell systems.\nUse pop=regular gfprint iop(10/33=2) keywords.")
    density_m = file.get_density_magnetic()
    if density_m is None:
        raise Exception("Error while reading the perturbed density matrices.\nCheck that NMR keyword has been used.")
    # create the 'densities' numpy array to collect all of them
    densities = numpy.zeros((1, 4) +density.shape)
    densities[0, 0] = density
    densities[0, 1:] = density_m
elif os.path.splitext(options.inputfile)[1] == ".log":
    file = Gaussian(options.inputfile)
    density = file.get_density()
    if density is None:
        raise Exception("Error while reading the density matrix.\nCheck that pop=regular or pop=full keyword has been used.")
    density_m = file.get_density_magnetic()
    if density_m is None:
        raise Exception("Error while reading the perturbed density matrices.\nCheck that iop(10/33=2) keyword has been used.")
    # create the 'densities' numpy array to collect all of them
    if density.ndim == 2:# close shell
        densities = numpy.zeros((1, 4) +density.shape)
        densities[0, 0] = density
        densities[0, 1:] = density_m
    else:#open shell
        k = density.shape[-1]
        densities = numpy.zeros((2, 4, k, k))
        densities[0, 0] = density[0]
        densities[1, 0] = density[1]
        density_m = density_m.reshape((2, 3, k, k))
        densities[0, 1:] = density_m[0]
        densities[1, 1:] = density_m[1]
else:
    parser.error("The inputfile is not a fchk or log file")


#read the basis set
basisset = file.get_basisset()
if basisset is None:
    parser.error("No basis set information found.\nIf you read a log file, check that gfprint keyword has been used.")

basisset = basisset.split_SP()

# prepare a transformation matrix from Spherical-harmonic to Cartesian orbitals
# create the block diagonal matrix of transformation
if options.turbomole:
    CartOrdering = "turbomole"
else:
    CartOrdering = "cfour"

if basisset.spherical:
    tmat = basisset.Cart2Spher(square=False, real=True, order=CartOrdering, normalized=False)
else:
    tmat = basisset.Cart2Cart(order1=CartOrdering, order2="gaussian", normalized1=False, normalized2=True)
if tmat is None:
    raise Exception("Error when creating the transformation matrix")

nbasisset = basisset.to_Cartesian()

# arrange the order of the basis function like CartOrdering
# S for all atoms, P for all atoms, D for all atoms, ...
aos = nbasisset.to_AOs(CartOrdering)
if options.turbomole:
    index = aos.index_sort_by_type()
else:
    index = aos.index_sort_by_atom()
aos = aos.basis_from_index(index)

if densities.shape[2] != tmat.shape[1]:
    raise Exception("The shape of the density matrix [%i] is not the same as the shape of the transformation matrix [%i]"%(densities.shape[2], tmat.shape[1]))
if len(index) != tmat.shape[0]:
    raise Exception("The length of the list of index [%i] is not the same as the shape of the transformation matrix [%i]"%(len(index), tmat.shape[0]))
tmat = tmat[index, :]#new, old
if options.turbomole:
    nbasisset = nbasisset.basis_from_index(nbasisset.index_sort_by_type())
else:
    nbasisset = nbasisset.basis_from_index(nbasisset.index_sort_by_atom())

# apply the transformation to 'densities' array
# the transformation reorganized the basis functions
# with the shells organized by types
# with the angular momenta organized like CartOrdering
nspin = densities.shape[0]
ndensities = numpy.zeros((nspin, 4, tmat.shape[0], tmat.shape[0]))
for ispin in range(nspin):
    for idir in range(4):# 0, Bx, By, Bz
        ndensities[ispin, idir] = numpy.dot(numpy.dot(tmat, densities[ispin, idir]), tmat.transpose())

# ATTENTION values in XDENS for Bx, By, Bz are 2 times bigger for closeshell
# and 4 times biggen for openshell
if nspin == 1:
    ndensities[:, 1:, :, :] *= 2
else:
    ndensities[:, 1:, :, :] *= 4

# open outputfile
outfile=open("XDENS", "w")

# write densities matrices on file
# write in a fortran way
for ispin in range(nspin):
    for idir in range(4):# 0, Bx, By, Bz
        for j in range(ndensities.shape[3]):
            for i in range(ndensities.shape[2]):
                outfile.write("%16.8e\n"%(ndensities[ispin, idir, i, j]))
        outfile.write("\n")
outfile.close()

# write the MOL file with coordinates and basis set
coordinates = None
if hasattr(file, "get_property"):
    coordinates = file.get_property("Current cartesian coordinates").reshape((-1, 3))
nbasisset.write_MOL(filename='MOL', coords=coordinates, turbomole=options.turbomole)


if options.XDENS != 'none':
    if not  os.path.isfile(options.XDENS):
        parser.error("The file ['%s'] does not exist"%(options.XDENS))
    f = open(options.XDENS, 'r')
    lines = f.readlines()
    f.close()
    data = []
    for line in lines:
        if len(line.strip()) != 0:
            data.append(float(line.strip()))
    data = numpy.array(data)
    k = nbasisset.nbasisfunctions
    data = data.reshape((-1, 4, k, k))
    ndensities_o = data.copy()
    # fortran to python arrangement
    for ispin in range(nspin):
        for idir in range(4):
            ndensities_o[ispin, idir] = data[ispin, idir].transpose()
    # check the densities matrices
    for ispin in range(nspin):
        for idir in range(4):
            print("Spin %s, %s"%("alpha + beta" if nspin == 1 else["alpha", "beta"][ispin]\
, "Unperturbed" if idir == 0 else "Magnetic field %i"%(idir)))
            print("        AO orbitals            :   GAUSSIAN     OTHER     DIFF    ")
            for i, ao1 in enumerate(aos.basis):
                for j, ao2 in enumerate(aos.basis):
                    print("%13s -- %13s : %12.8f %12.8f %12.8f"%(ao1.label, ao2.label, ndensities[ispin, idir, j, i], ndensities_o[ispin, idir, j, i], ndensities[ispin, idir, j, i] -ndensities_o[ispin, idir, j, i]))
            print("Maximum Error: %12.8f"%(numpy.amax(numpy.abs(ndensities[ispin, idir] -ndensities_o[ispin, idir]))))
EOFX

# Create BasisSet.py in current directory
# This module handles basis set operations needed for GIMIC calculations
cat << 'EOFX' > BasisSet.py




import numpy
import math
import re

Table=['Xx', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn']

multiplicity = {0:1, 1:3, -1:4, -2:5, 2:6, -3:7, 3:10, -4:9, 4:15, -5:11, 5:21, -6:13, 6:28}
mult2type = {1:0, 3:1, 4:-1, 5:-2, 6:2, 7:-3, 10:3, 9:-4, 15:4, 11:-5, 21:5, 13:-6, 28:6}
shelltypes = {0:"S", 1:"P", 2:"D", 3:"F", 4:"G", 5:"H", 6:"I", -2:"D", -3:"F", -4:"G", -5:"H", -6:"I", -1:"SP"}
orbitaltypes = ["S", "SP", "P", "D", "F", "G", "H", "I"]

angularmomenta_turbomole = {}
angularmomenta_turbomole["S"] = [(0, 0, 0)]
angularmomenta_turbomole["P"] = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
#angularmomenta_turbomole["SP"]=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
angularmomenta_turbomole["D"] = [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)]
angularmomenta_turbomole["F"] = [(3, 0, 0), (0, 3, 0), (0, 0, 3), (2, 1, 0), (2, 0, 1), (1, 2, 0), (0, 2, 1), (1, 0, 2), (0, 1, 2), (1, 1, 1)]
angularmomenta_turbomole["G"] = [(4, 0, 0), (0, 4, 0), (0, 0, 4), (3, 1, 0), (3, 0, 1), (1, 3, 0), (0, 3, 1), (1, 0, 3), (0, 1, 3), (2, 2, 0), (2, 0, 2), (0, 2, 2), (2, 1, 1), (1, 2, 1), (1, 1, 2)]
angularmomenta_turbomole["H"] = [(5, 0, 0), (0, 5, 0), (0, 0, 5), (4, 1, 0), (4, 0, 1), (1, 4, 0), (0, 4, 1), (1, 0, 4), (0, 1, 4), (3, 2, 0), (3, 0, 2), (2, 3, 0), (0, 3, 2), (2, 0, 3), (0, 2, 3), (3, 1, 1), (1, 3, 1), (1, 1, 3), (2, 2, 1), (2, 1, 2), (1, 2, 2)]
angularmomenta_turbomole["I"] = [(6, 0, 0), (0, 6, 0), (0, 0, 6), (5, 1, 0), (5, 0, 1), (1, 5, 0), (0, 5, 1), (1, 0, 5), (0, 1, 5), (4, 2, 0), (4, 0, 2), (2, 4, 0), (0, 4, 2), (2, 0, 4), (0, 2, 4), (4, 1, 1), (1, 4, 1), (1, 1, 4), (3, 3, 0), (3, 0, 3), (0, 3, 3), (3, 2, 1), (3, 1, 2), (2, 3, 1), (1, 3, 2), (2, 1, 3), (1, 2, 3), (2, 2, 2)]

angularmomenta_gaussian = {}
angularmomenta_gaussian["S"] = [(0, 0, 0)]
angularmomenta_gaussian["P"] = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
#angularmomenta_gaussian["SP"]=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
angularmomenta_gaussian["D"] = [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)]
angularmomenta_gaussian["F"] = [(3, 0, 0), (0, 3, 0), (0, 0, 3), (1, 2, 0), (2, 1, 0), (2, 0, 1), (1, 0, 2), (0, 1, 2), (0, 2, 1), (1, 1, 1)]
angularmomenta_gaussian["G"] = [(0, 0, 4), (0, 1, 3), (0, 2, 2), (0, 3, 1), (0, 4, 0), (1, 0, 3), (1, 1, 2), (1, 2, 1), (1, 3, 0), (2, 0, 2), (2, 1, 1), (2, 2, 0), (3, 0, 1), (3, 1, 0), (4, 0, 0)]
angularmomenta_gaussian["H"] = [(0, 0, 5), (0, 1, 4), (0, 2, 3), (0, 3, 2), (0, 4, 1), (0, 5, 0), (1, 0, 4), (1, 1, 3), (1, 2, 2), (1, 3, 1), (1, 4, 0), (2, 0, 3), (2, 1, 2), (2, 2, 1), (2, 3, 0), (3, 0, 2), (3, 1, 1), (3, 2, 0), (4, 0, 1), (4, 1, 0), (5, 0, 0)]
angularmomenta_gaussian["I"] = [(0, 0, 6), (0, 1, 5), (0, 2, 4), (0, 3, 3), (0, 4, 2), (0, 5, 1), (0, 6, 0), (1, 0, 5), (1, 1, 4), (1, 2, 3), (1, 3, 2), (1, 4, 1), (1, 5, 0), (2, 0, 4), (2, 1, 3), (2, 2, 2), (2, 3, 1), (2, 4, 0), (3, 0, 3), (3, 1, 2), (3, 2, 1), (3, 3, 0), (4, 0, 2), (4, 1, 1), (4, 2, 0), (5, 0, 1), (5, 1, 0), (6, 0, 0)]

angularmomenta_cfour = {}
angularmomenta_cfour["S"] = [(0, 0, 0)]
angularmomenta_cfour["P"] = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
#angularmomenta_cfour["SP"]=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
angularmomenta_cfour["D"] = [(2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 2, 0), (0, 1, 1), (0, 0, 2)]
angularmomenta_cfour["F"] = [(3, 0, 0), (2, 1, 0), (2, 0, 1), (1, 2, 0), (1, 1, 1), (1, 0, 2), (0, 3, 0), (0, 2, 1), (0, 1, 2), (0, 0, 3)]
angularmomenta_cfour["G"] = [(4, 0, 0), (3, 1, 0), (3, 0, 1), (2, 2, 0), (2, 1, 1), (2, 0, 2), (1, 3, 0), (1, 2, 1), (1, 1, 2), (1, 0, 3), (0, 4, 0), (0, 3, 1), (0, 2, 2), (0, 1, 3), (0, 0, 4)]
angularmomenta_cfour["H"] = [(5, 0, 0), (4, 1, 0), (4, 0, 1), (3, 2, 0), (3, 1, 1), (3, 0, 2), (2, 3, 0), (2, 2, 1), (2, 1, 2), (2, 0, 3), (1, 4, 0), (1, 3, 1), (1, 2, 2), (1, 1, 3), (1, 0, 4), (0, 5, 0), (0, 4, 1), (0, 3, 2), (0, 2, 3), (0, 1, 4), (0, 0, 5)]
angularmomenta_cfour["I"] = [(6, 0, 0), (5, 1, 0), (5, 0, 1), (4, 2, 0), (4, 1, 1), (4, 0, 2), (3, 3, 0), (3, 2, 1), (3, 1, 2), (3, 0, 3), (2, 4, 0), (2, 3, 1), (2, 2, 2), (2, 1, 3), (2, 0, 4), (1, 5, 0), (1, 4, 1), (1, 3, 2), (1, 2, 3), (1, 1, 4), (1, 0, 5), (0, 6, 0), (0, 5, 1), (0, 4, 2), (0, 3, 3), (0, 2, 4), (0, 1, 5), (0, 0, 6)]

Fields = ["x", "y", "z"]

mqn = {}
mqn["S"] = [0]
#mqn["SP"]=[0, 1, -1, 0]
mqn["P"] = [1, -1, 0]
mqn["D"] = [0, 1, -1, 2, -2]
mqn["F"] = [0, 1, -1, 2, -2, 3, -3]
mqn["G"] = [0, 1, -1, 2, -2, 3, -3, 4, -4]
mqn["H"] = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5]
mqn["I"] = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6]


def coeff_Cart2Spher(m, l, real=False, normalized=True):
    r"""
    Coefficient of transformation between a normalized Cartesian orbitals and a normalized spherical-harmonic ones
    \chi^sphe_{l, m} = \sum_{lx, ly, lz} c_{l, m, lx, ly, lz) \chi^cart_{lx, ly, lz}
    Param m: m value for normalized spherical-harmonic orbital
    Param l: l vector with lx, ly, lz values for normalized Cartesian orbital
    Param real: transform into real spherical-harmonic instead of complex (default)
    Param normalized: whether or not to use normalized Cartesian orbitals
    """
    lx, ly, lz = l
    l = sum(l)
    j = (lx + ly - abs(m))
    if j%2 == 1:
        return 0.0
    j = j // 2
    N = math.sqrt(float(math.factorial(2*lx) * math.factorial(2*ly) * math.factorial(2*lz) * math.factorial(l) * math.factorial(l-abs(m))) /\
                 float(math.factorial(2*l) * math.factorial(lx) * math.factorial(ly) * math.factorial(lz) * math.factorial(l+abs(m))))
    N = N / (2**l * math.factorial(l))
    fact1 = 0.0
    for i in range((l-abs(m))//2+1):
        bin1 = math.factorial(l) / (math.factorial(i) * math.factorial(l - i))
        try:
            bin2 = math.factorial(i) / (math.factorial(j) * math.factorial(i - j))
        except ValueError:
            bin2 = 0.0
        N2 = (-1)**i * math.factorial(2*l-2*i) / math.factorial(l - abs(m) - 2*i)
        fact1 = fact1 + bin1*bin2*N2
    fact2 = 0.0
    for k in range(j+1):
        bin3 = math.factorial(j) / (math.factorial(k) * math.factorial(j - k))
        try:
            bin4 = math.factorial(abs(m)) / (math.factorial(lx-2*k) * math.factorial(abs(m) - lx + 2*k))
        except ValueError:
            bin4 = 0.0
        N3 = (abs(m)-lx+2*k)
#        if N3 % 2 == 0:
#            N3 = (-1)**(N3//2)
#        else:
#            N3 = numpy.sign(m)*1j *( (-1)**((N3-1)//2))
        if not real:
            # N3 = (-1)**(\pm N3/2)
            fp = [1, 1j, -1, -1j]
            fm = [1, -1j, -1, 1j]
        else:
            # N3 = [ (-1)**(N3/2) + (-1)**(-N3/2)] / \sqrt{2}  for m>0
            # N3 = [ (-1)**(N3/2) - (-1)**(-N3/2)] / \sqrt{2}i for m<0
            fp = [math.sqrt(2.0), 0, -math.sqrt(2.0), 0]
            fm = [0, math.sqrt(2.0), 0, -math.sqrt(2.0)]
        if m > 0:
            N3 = fp[N3 % 4]
        elif m < 0:
            N3 = fm[N3 % 4]
        else:
#            fp = [1, 1j, -1, -1j]
#            N3 = fp[N3 % 4]
            N3 = 1.0 # the formula can be simplified into 1.0
        fact2 = fact2 + bin3 * bin4 * N3
    coeff = N * fact1 * fact2
    if not normalized:
        coeff = coeff * math.sqrt(float(math.factorial(lx)*math.factorial(ly)*math.factorial(lz) * 2**l) / float(math.factorial(2*lx) * math.factorial(2*ly) * math.factorial(2*lz)))
    return coeff


class SHELL(object):
    """
    Class to describe one shell
    A shell is described by a contraction of GTO
    A shell is quantified by l quantum number
    """

    def __init__(self, iatom, atnum, l, exponents, coefficients, coefficients_p=None, coord=None):
        """
        Param iatom: atom number (one-based)
        Param atnum: atomic number
        Param l: int to represent the type of atomic orbitals
                    the integer correspond to the 'l' quantum number with a sign
                    + sign means Cartesian function
                    - sign means Spherical-harmonic function
                    -1 means SP orbital
        Param exponents: numpy array with exponents of the contraction of GTO
        Param coefficients: numpy array with coefficients of the contraction of GTO
        Param coefficients_p: numpy array with coefficients of p orbitals for SP type
        Param coord: numpy array with the coordinates of the atoms in Angstrom
        """
        self.iatom = iatom
        self.atnum = atnum
        self.type = l
        self.exponents = exponents
        self.coefficients = coefficients
        self.coefficients_p = coefficients_p
        self.coord = coord
        self.nprimitive = self.exponents.size
        self.multiplicity = multiplicity[self.type]
        self.shelltype = shelltypes[self.type]
        self.label = Table[self.atnum]+"%i@"%(self.iatom)+self.shelltype

    def copy(self):
        newSHELL = SHELL(self.iatom, self.atnum, self.type, self.exponents, self.coefficients, self.coefficients_p, self.coord)
        return newSHELL

    def to_Cartesian(self):
        """
        Return a new SHELL where the orbitals are Cartesian
        """
        type = self.type
        if type < -1:
            type = -type
        newSHELL = SHELL(self.iatom, self.atnum, type, self.exponents, self.coefficients, self.coefficients_p, self.coord)
        return newSHELL

    def to_Spherical(self):
        """
        Return a new SHELL where the orbitals are Spherical-harmonic
        """
        type = self.type
        if type > 1:
            type = -type
        newSHELL = SHELL(self.iatom, self.atnum, type, self.exponents, self.coefficients, self.coefficients_p, self.coord)
        return newSHELL


    @staticmethod
    def mat_Cart2Spher(shelltype, square=False, real=False, order="gaussian", normalized=True):
        r"""
        Transformation matrices from normalized Cartesian orbitals to normalized spherical-harmonic ones
        for one type of orbital, one shell
        \chi^sphe_ {l, m}= \sum_ {lx, ly, lz} c_ {l, m, lx, ly, lz)\chi^cart_ {lx, ly, lz}
        Lines: normalized cartesian orbitals. Columns:normalized spherical-harmonic orbitals
        Param shelltype: one of orbitaltypes
        Param square: return square matrices with blank lines instead of rectangular matrices(default)
        Param real: transform into real spherical-harmonic instead of complex(default)
        Param order: order of the angularmomenta. Either "gaussian", "turbomole" or "cfour"
        Param normalized:whether or not to use normalized Cartesian orbitals
        Return: c_cart, spher
        """
        mqns = mqn[shelltype]
        angmoms = eval("angularmomenta_"+order)[shelltype]
        if real:
            dtype = "float"
        else:
            dtype = "complex"
        if square:
            c = numpy.zeros((len(angmoms), len(angmoms)), dtype)
        else:
            c = numpy.zeros((len(angmoms), len(mqns)), dtype)
        for il, l in enumerate(angmoms):
            for im, m in enumerate(mqns):
                c[il, im] = coeff_Cart2Spher(m, l, real=real, normalized=normalized)
        return c

    @staticmethod
    def Overlap_Cart(shelltype, order="gaussian", normalized=True):
        """
        Overlap matrix between Cartesian atomic orbitals of same l
        Param shelltype: one of orbitaltypes
        Param order: order of the angularmomenta. Either "gaussian", "turbomole" or "cfour"
        Param normalized: whether or not to use normalized Cartesian orbitals
        """
        angmoms = eval("angularmomenta_"+order)[shelltype]
        ov = numpy.zeros((len(angmoms), len(angmoms)))
        for il, l1 in enumerate(angmoms):
            for jl, l2 in enumerate(angmoms):
                if len(numpy.flatnonzero((numpy.array(l1)+numpy.array(l2))%2)) > 0:
                    ov[il, jl] = 0.0
                else:
                    num1 = math.factorial(l1[0]+l2[0]) * math.factorial(l1[1]+l2[1]) * math.factorial(l1[2]+l2[2])
                    den1 = math.factorial((l1[0]+l2[0])//2) * math.factorial((l1[1]+l2[1])//2) * math.factorial((l1[2]+l2[2])//2)
                    num2 = math.factorial(l1[0]) * math.factorial(l1[1]) * math.factorial(l1[2]) *\
                        math.factorial(l2[0]) * math.factorial(l2[1]) * math.factorial(l2[2])
                    den2 = math.factorial(2*l1[0]) * math.factorial(2*l1[1]) * math.factorial(2*l1[2]) *\
                        math.factorial(2*l2[0]) * math.factorial(2*l2[1]) * math.factorial(2*l2[2])
                    ov[il, jl] = float(num1) / float(den1) * math.sqrt(float(num2)/float(den2))
                    if not normalized:
                        f1 = math.sqrt(float(math.factorial(2*l1[0]) * math.factorial(2*l1[1]) * math.factorial(2*l1[2])) / float((math.factorial(l1[0])*math.factorial(l1[1])*math.factorial(l1[2])) * 2**(sum(l1))))
                        f2 = math.sqrt(float(math.factorial(2*l2[0]) * math.factorial(2*l2[1]) * math.factorial(2*l2[2])) / float((math.factorial(l2[0])*math.factorial(l2[1])*math.factorial(l2[2])) * 2**(sum(l2))))
                        ov[il, jl] = ov[il, jl] * f1 * f2
        return ov

    @staticmethod
    def mat_Spher2Cart(shelltype, square=False, real=False, order="gaussian", normalized=True):
        r"""
        Transformation matrices from normalized spherical-harmonic orbitals to normalized Cartesian ones
        for one type of orbital
        \chi^cart_{lx, ly, lz} = \sum_{l <= lx+ly+lz} c-1_{l, m, ly, lz, lz)\chi^sphe_{l, m}
        c-1 = c^dagger S since c^dagger S c = 1 = c-1 c
        Lines: normalized spherical-harmonic orbitals. Columns:normalized cartesian orbitals.
        Param shelltype: one of orbitaltypes
        Param square: return square matrices with blank lines instead of rectangular matrices(default)
        Param real: transform into real spherical-harmonic instead of complex(default)
        Param order: order of the angularmomenta. Either "gaussian", "turbomole" or "cfour"
        Param normalized: whether or not to use normalized Cartesian orbitals
        Return: c-1_spher, cart
        """
        C = SHELL.mat_Cart2Spher(shelltype, square=square, real=real, order=order, normalized=normalized)
        if numpy.dtype("complex") == C.dtype:
            C = C.conjugate()
        S = SHELL.Overlap_Cart(shelltype, order=order, normalized=normalized)
        Cm1 = numpy.dot(C.transpose(), S)
        return Cm1

    @staticmethod
    def mat_Cart2Cart(shelltype, order1="gaussian", order2="gaussian", normalized1=True, normalized2=True):
        """
        Transformation matrices from normalized Cartesian orbitals to normalized Cartesian ones
        for one type of orbital
        Lines: normalized cartesian orbitals 1. Columns:normalized cartesian orbitals 2.
        The arrangement of the angular momenta of the orbitals can be different between the two sets
        Param shelltype: one of orbitaltypes
        Param order1, order2: order of the angularmomenta. Either "gaussian", "turbomole" or "cfour"
        Param normalized1, normalized2: whether or not to use normalized Cartesian orbitals
        Return: c_cartold, cartnew
        """
        angmoms1 = eval("angularmomenta_"+order1)[shelltype]
        angmoms2 = eval("angularmomenta_"+order2)[shelltype]
        c = numpy.zeros((len(angmoms1), len(angmoms2)))
        for il, l1 in enumerate(angmoms1):
            fct = 1.0
            if (not normalized1) and normalized2:
                fct = math.sqrt(float((math.factorial(l1[0])*math.factorial(l1[1])*math.factorial(l1[2])) * 2**(sum(l1))) / float(math.factorial(2*l1[0]) * math.factorial(2*l1[1]) * math.factorial(2*l1[2])))
            elif normalized1 and (not normalized2):
                fct = math.sqrt(float(math.factorial(2*l1[0]) * math.factorial(2*l1[1]) * math.factorial(2*l1[2])) / float((math.factorial(l1[0])*math.factorial(l1[1])*math.factorial(l1[2])) * 2**(sum(l1))))
            jl = angmoms2.index(l1)
            c[il, jl] = fct
        return c

    @staticmethod
    def mat_Spher2Spher(shelltype, real1=False, real2=False):
        """
        Transformation matrices from normalized spherical-harmonic orbitals to normalized spherical-harmonic ones
        for one type of orbital
        Lines: normalized spherical-harmonic orbitals 1. Columns:normalized spherical-harmonic orbitals 2.
        Both sets can be real or complex
        Param shelltype:one of orbitaltypes
        Param real1, real1:transform into real spherical-harmonic instead of complex(default)
        Return: c_spherold, sphernew
        """
        mqns = mqn[shelltype]
        if real1 == real2:
            c = numpy.eye(len(mqns), len(mqns))
        else:
            # transformation matrix from imaginary to real
            ci2r = numpy.zeros((len(mqns), len(mqns)), "complex")
            for im, m1 in enumerate(mqns):
                if m1 == 0:
                    ci2r[im, im] = 1.0
                elif m1 > 0:
                    ci2r[im, im] = 1.0/math.sqrt(2.0)
                    jm = mqns.index(-m1)
                    ci2r[jm, im] = 1.0/math.sqrt(2.0)
                else:
                    ci2r[im, im] = 1.0j/math.sqrt(2.0)
                    jm = mqns.index(-m1)
                    ci2r[jm, im] = -1.0j/math.sqrt(2.0)
            if real2:
                c = ci2r
            else:
                c = ci2r.conjugate().transpose()
        return c

    def split_SP(self):
        """
        Split a SP shell in a S and a P shells
        return a list of two shells out of the SP shell
               or a list with one shell as a copy
        """
        shells = []
        type = self.type
        if type == -1:
            shells.append(SHELL(self.iatom, self.atnum, 0, exponents=self.exponents, coefficients=self.coefficients, coefficients_p=None, coord=self.coord))
            shells.append(SHELL(self.iatom, self.atnum, 1, exponents=self.exponents, coefficients=self.coefficients_p, coefficients_p=None, coord=self.coord))
        else:
            shells.append(self.copy())
        return shells

    def to_AO(self, order):
        """
        Convert a shell to a list of AOs
        Param order:order of the angularmomenta. Either "gaussian" or "turbomole"
        """
        AOs = []
        type = self.type
        if type == -1:#special for SP
            AOs.append(AO(self.iatom, self.atnum, -1, (0, 0, 0), exponents=self.exponents, coefficients=self.coefficients, coord=self.coord))
            AOs.append(AO(self.iatom, self.atnum, -1, (1, 0, 0), exponents=self.exponents, coefficients=self.coefficients_p, coord=self.coord))
            AOs.append(AO(self.iatom, self.atnum, -1, (0, 1, 0), exponents=self.exponents, coefficients=self.coefficients_p, coord=self.coord))
            AOs.append(AO(self.iatom, self.atnum, -1, (0, 0, 1), exponents=self.exponents, coefficients=self.coefficients_p, coord=self.coord))
        else:
            for im in range(self.multiplicity):
                if type < -1:
                    m = mqn[self.shelltype][im]
                else:
                    m = eval("angularmomenta_"+order)[self.shelltype][im]
                AOs.append(AO(self.iatom, self.atnum, type, m, exponents=self.exponents, coefficients=self.coefficients, coord=self.coord))
        return AOs


class AO(object):
    """
    Class to describe a AO function
    An AO function is described by a contraction of GTO
    An AO is quantified by l and m quantum numbers
    """

    def __init__(self, iatom, atnum, l, m, exponents, coefficients, coord=None):
        """
        Param iatom:atom number(one-based)
        Param atnum:atomic number
        Param l:int to represent the type of atomic orbitals
                    the integer correspond to the 'l' quantum number with a sign
                    + sign means Cartesian function
                    - sign means Spherical-harmonic function
                    -1 means SP orbital
        Param m:m quantum number either a int(spherical-harmonic) or tuple(cartesian)
        Param exponents:numpy array with exponents of the contraction of GTO
        Param coefficients:numpy array with coefficients of the contraction of GTO
        Param coord:numpy array with the coordinates of the atoms in Angstrom
        """
        self.iatom = iatom
        self.atnum = atnum
        self.type = l
        self.m = m
        self.exponents = exponents
        self.coefficients = coefficients
        self.coefficients_p = None
        self.coord = coord
        self.nprimitive = self.exponents.size
        self.multiplicity = 1
        self.shelltype = shelltypes[self.type]
        self.label = Table[self.atnum]+"%i@"%(self.iatom)+self.shelltype
        if isinstance(self.m, int):
            self.label = self.label + "%i"%(self.m)
        elif len(self.m) == 3:
            power = ""
            for idir in range(3):
                for _ in range(self.m[idir]):
                    power += Fields[idir]
            self.label = self.label + power
        else:
            raise Exception("m quantum number for AO is not correct")


    def copy(self):
        newAO = AO(self.iatom, self.atnum, self.type, self.m, self.exponents, self.coefficients, self.coord)
        return newAO

class BasisSet(object):

    def __init__(self, SHELLs):
        """
        Param SHELLs:list of SHELL objects
        """

        self.basis = SHELLs
        self.nshells = len(self.basis)
        self.maptype = []
        self.mapatom = []
        self.mapatnum = []
        self.nbasisfunctions = 0
        for shell in self.basis:
            self.maptype.append(shell.type)
            self.mapatom.append(shell.iatom)
            self.mapatnum.append(shell.atnum)
            self.nbasisfunctions += shell.multiplicity
        self.maptype = numpy.array(self.maptype)
        self.mapatom = numpy.array(self.mapatom)
        self.mapatnum = numpy.array(self.mapatnum)
        self.atoms = sorted(set(self.mapatom))# list of atoms
        self.atnums = sorted(set(self.mapatnum))#list of atnums
        self.natoms = len(self.atoms)
        # check the type of orbital
        self.spherical = None
        if numpy.all(self.maptype >= -1):# are all Cartesian orbitals
            self.spherical = False
        elif numpy.all(self.maptype <= 1):# are all Spherical-harmonic orbitals
            self.spherical = True
        else:
            raise Exception("Basis set contains both Cartesian and Spherical-harmonic functions")

    def to_Cartesian(self):
        """
        Return a new basis set where all the orbitals are Cartesian
        """
        newbasis = [shell.to_Cartesian() for shell in self.basis]
        newbasis = BasisSet(newbasis)
        return newbasis

    def to_Spherical(self):
        """
        Return a new basis set where all the orbitals are Spherical-harmonic
        """
        newbasis = [shell.to_Spherical() for shell in self.basis]
        newbasis = BasisSet(newbasis)
        return newbasis

    def split_SP(self):
        """
        Return a new basis set where the SP functions have been converted into S and P shells
        """
        if not isinstance(self.basis[0], SHELL):
            raise Exception("This method is only available for basis set made of shells")
        newbasis = []
        for shell in self.basis:
            newbasis.extend(shell.split_SP())
        newbasis = BasisSet(newbasis)
        return newbasis

    def to_AOs(self, order):
        """
        Return a new basis set where the list is made of AO functions instead of shells
        Param order:order of the angularmomenta. Either "gaussian" or "turbomole"
        """
        if isinstance(self.basis[0], SHELL):
            newbasis = []
            for shell in self.basis:
                newbasis.extend(shell.to_AO(order=order))
        else:
            newbasis = [ao.copy() for ao in self.basis]
        newbasis = BasisSet(newbasis)
        return newbasis

    def applyTransform(self, command):
        """
        Param command:python command to apply to every shell
        """
        index1 = 0
        index2 = 0
        m = numpy.zeros((10000, 10000))
        # Build the block-diagonal transformation matrix m
        for shell in self.basis:
            shelltype = shell.shelltype
            data = eval(command.replace("@T", shelltype))
            end1 = index1 + data.shape[0]
            end2 = index2 + data.shape[1]
            m[index1:end1, index2:end2] = data
            index1 = end1
            index2 = end2
        return m[:index1, :index2]


    def Cart2Spher(self, square=False, real=False, order="gaussian", normalized=True):
        """
        Transformation matrices from normalized Cartesian orbitals to normalized spherical-harmonic ones
        Lines:normalized cartesian orbitals. Columns:normalized spherical-harmonic orbitals
        Param square:return square matrices with blank lines instead of rectangular matrices(default)
        Param real:transform into real spherical-harmonic instead of complex(default)
        Param order:order of the angularmomenta. Either "gaussian" or "turbomole"
        Param normalized:whether or not to use normalized Cartesian orbitals
        Return:c_cart, spher
        """
        command = "SHELL.mat_Cart2Spher('@T', square="+str(square)+", real="+str(real)+",order='"+order+"', normalized="+str(normalized) +")"
        return self.applyTransform(command)

    def Spher2Cart(self, square=False, real=False, order="gaussian", normalized=True):
        """
        Transformation matrices from normalized spherical-harmonic orbitals to normalized Cartesian ones
        c-1 = c^dagger S since c^dagger S c = 1 = c-1 c
        Lines:normalized spherical-harmonic orbitals. Columns:normalized cartesian orbitals.
        Param square:return square matrices with blank lines instead of rectangular matrices(default)
        Param real:transform into real spherical-harmonic instead of complex(default)
        Param order:order of the angularmomenta. Either "gaussian" or "turbomole"
        Param normalized:whether or not to use normalized Cartesian orbitals
        Return:c-1_spher, cart
        """
        command = "SHELL.mat_Spher2Cart('@T', square="+str(square)+", real="+str(real)+",order='"+order+"', normalized="+str(normalized) +")"
        return self.applyTransform(command)

    def Cart2Cart(self, order1="gaussian", order2="gaussian", normalized1=True, normalized2=True):
        """
        Transformation matrices from normalized Cartesian orbitals to normalized Cartesian ones
        Lines:normalized cartesian orbitals 1. Columns:normalized cartesian orbitals 2.
        The arrangement of the angular momenta of the orbitals can be different between the two sets
        Param order1, order2:order of the angularmomenta. Either "gaussian" or "turbomole"
        Param normalized1, normalized2:whether or not to use normalized Cartesian orbitals
        Return:c_cartold, cartnew
        """
        command = "SHELL.mat_Cart2Cart('@T',order1='"+order1+"',order2='"+order2+"', normalized1="+str(normalized1)+", normalized2="+str(normalized2)+")"
        return self.applyTransform(command)

    def Spher2Spher(self, real1=False, real2=False):
        """
        Transformation matrices from normalized spherical-harmonic orbitals to normalized spherical-harmonic ones
        Lines:normalized spherical-harmonic orbitals 1. Columns:normalized spherical-harmonic orbitals 2.
        Both sets can be real or complex
        Param shelltype:one of orbitaltypes
        Param real1, real1:transform into real spherical-harmonic instead of complex(default)
        Return:c_spherold, sphernew
        """
        command = "SHELL.mat_Spher2Spher('@T', real1="+str(real1)+", real2="+str(real2)+")"
        return self.applyTransform(command)

    def get_index_for_atom(self, iatom):
        """
        return the index of the SHELLs associated to an atom 'iatom'
        Param iatom:atom number (one-based)
        """
        index = numpy.where(self.mapatom == iatom)[0]
        return index

    def get_index_for_type(self, type):
        """
        return the index of the SHELLs associated to a type
        Param type:int to represent the type of atomic orbitals
                    the integer correspond to the 'l' quantum number with a sign
                    + sign means Cartesian function
                    - sign means Spherical-harmonic function
                    -1 means SP orbital
        """
        index = numpy.where(self.maptype == type)[0]
        return index

    def index_sort_by_type(self):
        """
        Sort the SHELLs by type (slow axis) and by atom (fast axis)
        """
        mapatom = self.mapatom
        # convert self.maptype into index array using 'orbitaltypes' order
        maptype = numpy.array([orbitaltypes.index(shelltypes[t]) for t in self.maptype])
        # first sort along mapatom
        index1 = numpy.argsort(mapatom, kind="mergesort")
        # sort maptype with these index
        maptype = maptype[index1]
        # then sort along maptype
        index2 = numpy.argsort(maptype, kind="mergesort")
        index = index1[index2]
        return index

    def index_sort_by_atom(self):
        """
        Sort the SHELLs by atom (slow axis) and by type (fast axis)
        """
        mapatom = self.mapatom
        # convert self.maptype into index array using 'orbitaltypes' order
        maptype = numpy.array([orbitaltypes.index(shelltypes[t]) for t in self.maptype])
        # first sort along maptype
        index1 = numpy.argsort(maptype, kind="mergesort")
        # sort mapatom with these index
        mapatom = mapatom[index1]
        # then sort along mapatom
        index2 = numpy.argsort(mapatom, kind="mergesort")
        index = index1[index2]
        return index

    def basis_from_index(self, index):
        """
        Return a new basis set with the SHELLs in index
        Param index:list of index
        """
        basis = BasisSet([self.basis[i].copy() for i in index])
        return basis


    def get_atombasis(self, iatom):
        """
        Return a new basis set with the basis functions associated to an atom iatom
        """
        return self.basis_from_index(self.get_index_for_atom(iatom))

    def get_nborbitals(self):
        """
        Return a list with the number orbitals for each 'orbitaltypes'
        """
        typelist = [[k for k, v in list(shelltypes.items()) if v == t] for t in orbitaltypes]
        types = []
        for ts in typelist:
            index = []
            for t in ts:
                index.extend(self.get_index_for_type(t))
            types.append(len(index))
        return types

    def write(self, file):
        """
        write the basis set in file
        """
        for atnum in self.atnums:
            iatom = self.mapatom[numpy.where(self.mapatnum == atnum)][0]
            atombasis = self.get_atombasis(iatom)
            file.write("****\n")
            file.write("%s\n"%(Table[int(atnum)]))
            file.write("Nb_shells %i\n"%(atombasis.nshells))
            for shell in atombasis.basis:
                file.write("%s  %i \n"%(shell.shelltype, shell.nprimitive))
                coeffs = shell.coefficients
                coeffs_p = shell.coefficients_p
                exp = shell.exponents
                if coeffs_p is not None:
                    for a, b, c in zip(exp, coeffs, coeffs_p):
                        file.write("  %18.10E  %18.10E  %18.10E\n"%(a, b, c))
                else:
                    for a, b in zip(exp, coeffs):
                        file.write("  %18.10E  %18.10E\n"%(a, b))

    def write_MOL(self, filename, coords=None, turbomole=False):
        """
        write the basis set in file
        Param filename:filename to write into
        Param coords:coordinates of the atoms
        """
        if turbomole:
            HEADER = """INTGRL        1    0    1    0    0    0    0    0    0
TURBOMOLE
              Generated by Gaussian2gimic
"""
        else:
            HEADER = """INTGRL        1    0    1    0    0    0    0    0    0
CFOUR
              Generated by Gaussian2gimic
"""
        file = open(filename, "w")
        file.write(HEADER)
        file.write("%i    0            0.10E-08              0    0\n"%(self.natoms))
        file.write("9999.00      3.00\n")
        for iatom  in self.atoms:
            atombasis = self.get_atombasis(iatom)
            nborbitals = [x for x in atombasis.get_nborbitals() if x != 0]# keep the nb of orbitals that are non-zero
            atnum = atombasis.atnums[0]
            if coords is not None:
                coord = coords[iatom-1]
            else:
                coord = atombasis.basis[0].coord
                if coord is None:
                    raise Exception("No atomic coordinate found for atom %i"%(iatom))
            file.write("%2.1f  %3i %3i"%(float(atnum), 1, len(nborbitals)))
            for norb in nborbitals:
                file.write(" %3i"%(norb))
            file.write("\n")
            file.write("%s %i %19.12f %19.12f %19.12f\n"%(Table[int(atnum)], 1, coord[0], coord[1], coord[2]))
            for shell in atombasis.basis:
                file.write("   %3i %3i\n"%(shell.nprimitive, 1))
                coeffs = shell.coefficients
                coeffs_p = shell.coefficients_p
                exp = shell.exponents
                if coeffs_p is not None:
                    for a, b, c in zip(exp, coeffs, coeffs_p):
                        file.write("  %16.10f  %16.10f  %16.10f\n"%(a, b, c))
                else:
                    for a, b in zip(exp, coeffs):
                        file.write("  %16.10f  %16.10f\n"%(a, b))
        file.close()

    @staticmethod
    def read_MOL(filename):
        """
        read the basis set from MOL file
        Param filename:filename of the MOL file
        """
        file = open(filename, "r")
        lines = file.readlines()
        file.close()
#        natoms = int(lines[3].split()[0])
        basisset = []
        lines = lines[5:]
        iatom = 1
        # remove empty lines
        lines = [l for l in lines if len(l.strip()) != 0]
        while len(lines) > 0:
            data = lines[0].strip().split()
            atnum = int(float(data[0]))
            nborbitals = [int(x) for x in data[3:]]
            coord = [float(x) for x in lines[1].split()[2:]]
            lines = lines[2:]
            for itype, nborbital in enumerate(nborbitals):
                for _ in range(nborbital):
                    nprim = int(lines[0].split()[0])
                    block = lines[1:1+nprim]
                    lines = lines[1+nprim:]
                    data = []
                    for line in block:
                        data.append([float(x) for x in line.split()])
                    data = numpy.array(data)
                    exponents = data[:, 0]
                    coefficients = data[:, 1]
                    try:
                        coefficients_p = data[:, 2]
                    except IndexError:
                        coefficients_p = None
                    # FIXME: itype does not say if Cart or Spherical and no SP
                    shell = SHELL(iatom, atnum, itype, exponents, coefficients, coefficients_p, coord=coord)
                    basisset.append(shell)
            iatom += 1
        if len(basisset) == 0:
            basisset = None
        else:
            basisset = BasisSet(basisset)
        return basisset

    def check_contractedcoeff(self):
        r"""
        Check that=
        \sum_p, q d_p\mu d_q\mu  2^(l+3/2)\sqrt(a_p^l+3/2 a_q^l+3/2) / (a_p+a_q)^l+3/2
        """
        for shell in self.basis:
            type = shell.type
            if type == -1:
                type = [0, 1]
            else:
                type = [abs(type)]
            coeffs = shell.coefficients
            coeffs_p = shell.coefficients_p
            exponents = shell.exponents
            if coeffs_p is None:
                coefficients = coeffs[numpy.newaxis, :]
            else:
                coefficients = numpy.vstack((coeffs, coeffs_p))
            # loop over all set of coefficients for each l
            # (only important for sp shell)
            value = []
            for il, l in enumerate(type):
                aa = numpy.sqrt(numpy.outer(exponents, exponents) **(l+1.5))
                dd = numpy.outer(coefficients[il], coefficients[il])
                apa = numpy.add.outer(exponents, exponents) **(l+1.5)
                value.append((2**(l+1.5) * dd * aa / apa).sum())
            print(shell.label, ["%.8f"%(v) for v in value])

EOFX

# Confirm Python files have been created successfully
echo -e "\e[38;2;142;180;139mPython files created successfully in the current directory.\e[0m"
echo

# Step 3: Creating MOL and XDENS Files
# Search for input files in the current directory

# Look for .fchk files first (preferred format)
fchk_file=$(ls *.fchk 2>/dev/null | head -n 1)

if [ -n "$fchk_file" ]; then
    # If .fchk file is found, process it
    echo -e "\e[38;2;142;180;139mFound \".fchk\" file: $fchk_file\nCreating MOL and XDENS Files...\n\e[0m"
    # Run Python converter on the .fchk file
    python Gaussian2gimic.py -i "$fchk_file"
else
    # If no .fchk file, look for .log files
    log_file=$(ls *.log 2>/dev/null | head -n 1)
    if [ -n "$log_file" ]; then
        # If .log file is found, process it instead
        echo -e "\e[38;2;142;180;139mNo \".fchk\" found. Found \".log\" file: $log_file\nCreating MOL and XDENS Files...\n\e[0m"
        # Run Python converter on the .log file with appropriate flag
        python Gaussian2gimic.py --input="$log_file"
    else
        # If neither file type is found, display error and exit
        echo -e "\e[3;38;2;250;103;117mError: No \".fchk\" or \".log\" file found in the current directory\e[0m"
        exit 1
    fi
fi

# Step 4: Perform a dry run to check configuration
# This validates the input without running the full calculation
echo -e "\e[38;2;142;180;139mPerforming GIMIC dry run...\e[0m"
gimic --dryrun

# Check if gimic.inp was created successfully
if [ ! -f gimic.inp ]; then
    # Exit if configuration file wasn't generated
    echo -e "\e[3;38;2;250;103;117mError: Failed to generate gimic.inp file\e[0m"
    exit 1
fi

echo

# Step 5: Creating CML coordinate file for visualization software
echo -e "\e[38;2;142;180;139m---------------------------------------------\nConverting mol.xyz file to mol.cml/mol-bhor.cml file\n\e[0m"
function molecule() {
    obabel -ixyz mol.xyz -O mol.cml
    awk '{FS="\""; OFS="\"";
         if ($1 ~ "<atom id") {
             if ($5 ~ "spinMultiplicity")
                 { print $1, $2, $3, $4, $5, $6, $7, $8/0.526, $9, $10/0.526, $11, $12/0.526, $13 }
             else  { print $1, $2, $3, $4, $5, $6/0.526, $7, $8/0.526, $9, $10/0.526, $11 }
             }
         else print $0; }' mol.cml > mol-bohr.cml
}


# It takes an XYZ file as an argument:

molecule mol.xyz

echo -e "\e[38;2;142;180;139mCoordinates converted and saved to mol-bohr.cml\e[0m"

# Step 6: Run the main GIMIC calculation
echo -e "\e[38;2;142;180;139mRunning main GIMIC calculation...\e[0m"
# Execute GIMIC with the input file and redirect output to gimic.out
gimic gimic.inp > gimic.out

# Check if the calculation completed successfully
if [ $? -eq 0 ]; then
    echo -e "\e[38;2;142;180;139mGIMIC run complete. Output written to gimic.out\e[0m"
else
    echo -e "\e[3;38;2;250;103;117mError: GIMIC calculation failed. Check gimic.out for details.\e[0m"
fi

# Step 7: Clean up by deactivating the virtual environment
echo -e "\e[3;38;2;139;96;87m---------------------------------------------\nDeactivating environment...\n\e[0m"
# Exit the virtual environment
deactivate

# Final completion message
echo -e "\e[38;2;1;144;83mDone. Environment closed.\e[0m"

