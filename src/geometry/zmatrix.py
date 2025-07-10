#  Planck
#  Copyright (C) 2024 Hemanth Haridas, University of Utah
#  Contact: hemanthhari23@gmail.com
# 
#  This program is free software: you can redistribute it and/or modify it under
#  the terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 3 of the License, or a later version.
# 
#  This program is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#  PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License along with
#  this program.  If not, see <http://www.gnu.org/licenses/>.

from planck.src.exceptions.base import IllDefinedGeometryError
from planck.src.geometry.base import BaseMolecule
from planck.src.symmetry import symm
from planck.src.helpers import maths
from planck.src.helpers import tables
import numpy

class Molecule(BaseMolecule):
    """
    A class representing a molecule, inheriting from the BaseMolecule class. This class processes molecular
    structures provided as a string, extracting details such as molecular charge, multiplicity, number of 
    atoms, atomic symbols, and coordinates.

    Attributes:
    -----------
    atomicnumbers (list): A list to store the atomic numbers.
    atoms (list)        : A list to store the atomic symbols of the molecule.
    charge (int)        : The molecular charge.
    coords (list)       : A list to store the atomic coordinates as lists of floats.
    multi (int)         : The multiplicity of the molecule.
    natoms (int)        : The number of atoms in the molecule.

    Methods:
    --------
    geometry(structure: str) -> None:
        Parses the molecular structure string to populate the molecular attributes such as charge, multiplicity, 
        number of atoms, atomic symbols, and coordinates.
    """
    
    def geometry(self, structure: str, use_symm: bool = True) -> None:  
        """
        Builds the molecular structure from a Z-matrix input.

        Args:
        -----
            structure (str): A Z-matrix formatted string defining the molecular structure.
                - The first line contains the molecular charge and multiplicity, separated by a space.
                - Subsequent lines define the atoms and their relative positions.

        Attributes Initialized:
        -----------------------
            - atoms (list of str): List of atomic symbols in the molecule.
            - coords (lists of float): List of Cartesian coordinates for each atom.
              The first atom is always placed at the origin [0, 0, 0].
            - charge (int): Net charge of the molecule, extracted from the first line.
            - multi (int): Spin multiplicity of the molecule, extracted from the first line.
        """
        self.atoms         = []
        self.coords        = []
        self.atomicnumbers = []
        
        # Split the input structure into lines, ignoring empty ones
        _structure  = [line for line in structure.splitlines() if line.strip()]
        self.natoms = len(_structure) - 1
        
        # Extract charge and multiplicity from the first line
        self.charge = int(_structure[0].split()[0])
        self.multi  = int(_structure[0].split()[1])
        
        if self.natoms < 1:
            raise IllDefinedGeometryError(message="No atoms have been defined. Please define atleast one atom in the input block.")
            return 
        
        # Append the first atom to the list of atoms
        self.atoms.append(_structure[1].split()[0])
        self.atomicnumbers.append(tables.atomic_numbers[self.atoms[0]]) # Append atomic numbers
        self.coords.append([0, 0, 0]) # Because first atom is always placed at origin
        
        if self.natoms == 1:
            self.coords        = numpy.array(self.coords).flatten()
            self.atomicnumbers = numpy.array(self.atomicnumbers).flatten()
            return None
        
        # Now read the second atom
        self.atoms.append(_structure[2].split()[0])
        self.atomicnumbers.append(tables.atomic_numbers[self.atoms[1]]) # Append atomic numbers
        _distance = float(_structure[2].split()[-1]) # Second atom is always placed on x-axis and will be connected to first
        self.coords.append([_distance, 0, 0])
        
        if self.natoms == 2:
            self.coords        = numpy.array(self.coords).flatten()
            self.atomicnumbers = numpy.array(self.atomicnumbers).flatten()
            return None
                    
        # Now read the third line
        self.atoms.append(_structure[3].split()[0])
        self.atomicnumbers.append(tables.atomic_numbers[self.atoms[2]]) # Append atomic numbers
        _distance = float(_structure[3].split()[2])     # Get the distance to the atom bonded
        _angle    = float(_structure[3].split()[4])     # Get the angle to the next nearest neighbour
        _bonded_a = int(_structure[3].split()[1]) - 1   # Get the atom to which it is bonded
        _angle_a  = int(_structure[3].split()[3]) - 1   # Get the next nearest neighbour
        
        # Calculate the bond vector along x-axis
        _bond_vec = numpy.array(self.coords[_angle_a]) - numpy.array(self.coords[_bonded_a])
        _bond_vec = _bond_vec / numpy.linalg.norm(_bond_vec)

        # Rotate this along z-axis to required angle and place the atom
        _x_coord = self.coords[_bonded_a][0] + _distance * _bond_vec[0] * numpy.cos(_angle * numpy.pi/180) - _distance * _bond_vec[1] * numpy.sin(_angle * numpy.pi/180)# Calculate the position in XY plane
        _y_coord = self.coords[_bonded_a][1] + _distance * _bond_vec[0] * numpy.sin(_angle * numpy.pi/180) + _distance * _bond_vec[1] * numpy.cos(_angle * numpy.pi/180)# Calculate the position in XY plane
        self.coords.append([_x_coord, _y_coord, 0])
        
        if self.natoms == 3:
            self.coords        = numpy.array(self.coords).flatten()
            self.atomicnumbers = numpy.array(self.atomicnumbers).flatten()
            return None
        
        # Now process remaining lines
        for _line in _structure[4:]:
            self.atoms.append(_line.split()[0])
            self.atomicnumbers.append(tables.atomic_numbers[_line.split()[0]]) # Append atomic numbers
            _distance = float(_line.split()[2])     # Get the distance to the atom bonded
            _angle    = float(_line.split()[4])     # Get the angle to the next nearest neighbour
            _dihedral = float(_line.split()[6])     # Get the dihedral to the first atom
            _bonded_a = int(_line.split()[1]) - 1   # Get the atom to which it is bonded
            _angle_a  = int(_line.split()[3]) - 1   # Get the next nearest neighbour
            _dihed_a  = int(_line.split()[5]) - 1   # Get the first atom
        
            # Calculate the bond vector
            _bond_vec = numpy.array(self.coords[_angle_a]) - numpy.array(self.coords[_bonded_a])
            _bond_vec = _bond_vec / numpy.linalg.norm(_bond_vec)  
            
            # Calculate the angle vector
            _angle_vec = numpy.array(self.coords[_angle_a]) - numpy.array(self.coords[_dihed_a])
            _angle_vec = _angle_vec / numpy.linalg.norm(_angle_vec)   
            
            # Calculate the normal plane
            _norm_vec = numpy.cross(_bond_vec, _angle_vec)
            
            # Rotate the  bond_vec about normal plane by angle
            _ini_vec = _distance * _bond_vec
            _rot_mat  = maths.rotation_matrix(axis = _norm_vec, angle = _angle)
            _rot_bond = numpy.dot(_rot_mat, _ini_vec)
            
            # Rotate the rotated bond about bond_vec by dihedral
            _rot_mat  = maths.rotation_matrix(axis = _ini_vec, angle = _dihedral)
            _rot_bond = numpy.dot(_rot_mat, _rot_bond)
            
            # Calculate the new position
            _x_coord = self.coords[_bonded_a][0] + _rot_bond[0]
            _y_coord = self.coords[_bonded_a][1] + _rot_bond[1]
            _z_coord = self.coords[_bonded_a][2] + _rot_bond[2]
            
            self.coords.append([_x_coord, _y_coord, _z_coord])
            
        # Flatten the list for easier processing
        self.coords        = numpy.array(self.coords)
        self.atomicnumbers = numpy.array(self.atomicnumbers).flatten()