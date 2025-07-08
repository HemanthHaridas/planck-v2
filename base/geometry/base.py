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

from abc import ABC, abstractmethod
import numpy

class BaseMolecule(ABC):
    """
    Abstract Base Class (ABC) for representing a molecule.

    This class serves as a blueprint for creating specific molecule implementations. 
    It requires subclasses to define the `geometry` method, which should describe 
    the geometric properties of the molecule, such as its spatial structure, coordinates, 
    or connectivity.

    Methods:
        geometry():
            Abstract method that must be implemented by subclasses to define 
            the geometric properties of the molecule.
    """

    @abstractmethod
    def geometry(self) -> None:
        """
        Define the geometric properties of the molecule.

        This method should be implemented in subclasses to provide specific details
        about the molecule's geometry. The implementation can include attributes 
        such as atomic positions, bond lengths, angles, and molecular shape.

        Raises:
            NotImplementedError: If the method is not implemented in a subclass.
        """
        pass