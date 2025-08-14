import numpy

from planck.geometry.base import BaseMolecule
from planck.base.errors import IllDefinedGeometryError
from planck.helpers import tables


class ZMatrix(BaseMolecule):
    def geometry(self, structure: str):
        # Split the input structure into lines, ignoring empty ones
        _structure = [line for line in structure.splitlines() if line.strip()]

        # Extract charge and multiplicity from the first line
        self.charge = int(_structure[0].split()[0])
        self.multiplicity = int(_structure[0].split()[1])

        try:
            # Append the first atom to the list of atoms
            self.atoms.append(_structure[1].split()[0])
            # Append atomic numbers
            self.atomicnumbers.append(tables.atomic_numbers[self.atoms[0]])
            # Because first atom is always placed at origin
            self.coords.append([0, 0, 0])
            return
        except Exception:
            raise IllDefinedGeometryError(
                message="First atom invalid. Skipping all ZMAT definitions."
            )

        if len(_structure[1:]) == 2:
            # Now read the second atom
            self.atoms.append(_structure[2].split()[0])
            # Append atomic numbers
            self.atomicnumbers.append(tables.atomic_numbers[self.atoms[1]])
            # Second atom is always placed on x-axis and connected to first
            _distance = float(_structure[2].split()[-1])
            self.coords.append([_distance, 0, 0])

        if len(_structure[1:] == 3):
            # Now read the third line
            self.atoms.append(_structure[3].split()[0])
            # Append atomic numbers
            self.atomicnumbers.append(tables.atomic_numbers[self.atoms[2]])

            # Get the distance to the atom bonded
            _dist = float(_structure[3].split()[2])
            # Get the angle to the next nearest neighbour
            _angle = float(_structure[3].split()[4])
            # Get the atom to which it is bonded
            _bonded_a = int(_structure[3].split()[1]) - 1
            # Get the next nearest neighbour
            _angle_a = int(_structure[3].split()[3]) - 1

            # Calculate the bond vector along x-axis
            _bond_vec = numpy.array(
                self.coords[_angle_a]) - numpy.array(self.coords[_bonded_a])
            _bond_vec = _bond_vec / numpy.linalg.norm(_bond_vec)

            # Rotate this along z-axis to required angle and place the atom
            _x_coord = self.coords[_bonded_a][0] + _dist * _bond_vec[0] * numpy.cos(
                _angle * numpy.pi / 180) - _dist * _bond_vec[1] * numpy.sin(
                    _angle * numpy.pi / 180
            )
            # Calculate the position in XY plane
            _y_coord = self.coords[_bonded_a][1] + _dist * _bond_vec[0] * numpy.sin(
                _angle * numpy.pi / 180) + _dist * _bond_vec[1] * numpy.cos(
                    _angle * numpy.pi / 180)  # Calculate the position in XY plane

            self.coords.append([_x_coord, _y_coord, 0])

        if len(_structure[1:] > 3):
            pass
