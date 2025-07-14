from abc import ABC, abstractmethod


class BaseMolecule(ABC):
    """
    Abstract base class representing a generic molecular object.

    This class serves as a blueprint for concrete molecular representations
    by enforcing the implementation of geometry-related methods.

    Methods
    -------
    geometry()
        Abstract method that must be implemented to define or return
        the molecular geometry in derived classes.
    """

    @abstractmethod
    def geometry(self):
        """
        Define or return the molecular geometry.

        This method must be implemented by subclasses. It typically provides
        information such as atomic coordinates, bonding structure, or
        spatial configuration relevant to molecular modeling tasks.
        """
        pass
