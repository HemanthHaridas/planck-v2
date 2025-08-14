from abc import ABC, abstractmethod


class BaseMolecule(ABC):

    # Abstract Base Class (ABC) for representing a molecule.

    # This class is a blueprint for creating specific molecule implementations.
    # It requires subclasses to define the `geometry`
    # method, which should describe the geometric properties of the molecule,
    # such as its spatial structure, coordinates, or connectivity.

    # Methods:
    #     geometry():
    #         Abstract method that must be implemented by subclasses
    #         to define the geometric properties of the molecule.

    @abstractmethod
    def geometry(self) -> None:

        # Define the geometric properties of the molecule.

        # This method should be implemented in subclasses to provide
        # specific details about the molecule's geometry. The implementation
        # can include attributes such as atomic positions, bond lengths,
        # angles, and molecular shape.

        # Raises:
        #     NotImplementedError: If the method is not implemented
        #     in a subclass.
        pass
