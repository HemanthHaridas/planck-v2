import numpy

from planck.geometry.base import BaseMolecule
from planck.base.errors import IllDefinedGeometryError
from planck.helpers import tables


class Cartesian(BaseMolecule):

    # A class representing a molecule, inheriting from the BaseMolecule class .
    # This class processes molecular structures provided as a string,
    # extracting details such as molecular charge, multiplicity, number
    # of atoms, atomic symbols, and coordinates.

    # Attributes:
    # -----------
    # atomicnumbers(list): A list to store the atomic numbers.
    # atoms(list): A list to store the atomic symbols of the molecule.
    # charge(int): The molecular charge.
    # coords(list): A list to store the atomic coordinates as lists of floats.
    # multi(int): The multiplicity of the molecule.
    # natoms(int): The number of atoms in the molecule.

    # Methods:
    # --------
    # geometry(structure: str) -> None:
    # Parses the molecular structure string to populate the molecular
    # attributes such as charge, multiplicity, number of atoms, atomic symbols,
    # and coordinates.

    def geometry(self, structure: str):
        # Split the input structure into lines, ignoring empty ones
        _structure = [line for line in structure.splitlines() if line.strip()]

        # Extract charge and multiplicity from the first line
        self.charge = int(_structure[0].split()[0])
        self.multiplicity = int(_structure[0].split()[1])

        # Calculate the number of atoms
        # Subtract 1 to exclude the first line
        self.natoms = len(_structure) - 1

        # Leave this exception unhandled
        # Should stop the execution
        if self.natoms < 1:
            raise IllDefinedGeometryError(
                message="No atoms have been defined in the coordinate section."
            )

        # Initialize lists to store atomic symbols and coordinates
        self.atoms = []
        self.atomicnumbers = []
        self.coords = []

        # Process each atom line to extract atomic symbols and coordinates
        for _line in _structure[1:]:
            _atom = _line.split()[0]  # Extract the atomic symbol
            _coords = [float(_coord) for _coord in _line.split()[
                1:]]  # Extract and convert coordinates
            self.atoms.append(_atom)  # Append atomic symbol to atoms list
            self.coords.append(_coords)  # Append coordinates to coords list
            # Append atomic numbers
            self.atomicnumbers.append(tables.atomic_numbers[_atom])

        # Flatten the list for easier processing
        self.coords = numpy.array(self.coords).flatten()
        self.atomicnumbers = numpy.array(self.atomicnumbers).flatten()

    def build(self, atoms: list[str], coords: list[list[float]],
              charge: int, multiplicity: int):
        self.atoms = atoms
        self.coords = coords.to_list()
        self.charge = charge
        self.multiplicity = multiplicity
        self.atomicnumbers = [tables.atomic_numbers[_atom]
                              for _atom in self.atoms]

        # Flatten the list for easier processing
        self.coords = numpy.array(self.coords).flatten()
        self.atomicnumbers = numpy.array(self.atomicnumbers).flatten()
