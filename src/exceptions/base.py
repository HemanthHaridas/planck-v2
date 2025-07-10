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

class ChargeMultiplicityError(Exception):
    """
    Exception raised for errors related to invalid charge and multiplicity combinations.

    This custom exception is used to handle cases where the combination of charge 
    and multiplicity is physically or logically invalid in quantum chemistry calculations.

    Attributes
    ----------
    message : str
        A description of the error.

    Methods
    -------
    __init__(message: str)
        Initializes the exception with a specific error message.
    """
    def __init__(self, message: str):
        """
        Initializes the ChargeMultiplicityError with a custom error message.

        Parameters
        ----------
        message : str
            A description of the error, providing details about the invalid combination.
        """
        self.message = message
        super().__init__(self.message)

class IllDefinedGeometryError(Exception):
    """
    Exception raised for errors related to invalid geometry.

    This custom exception is used to handle cases where the geometry is ill defined 
    in quantum chemistry calculations.

    Attributes
    ----------
    message : str
        A description of the error.

    Methods
    -------
    __init__(message: str)
        Initializes the exception with a specific error message.
    """
    def __init__(self, message: str):
        """
        Initializes the IllDefinedGeometryError with a custom error message.

        Parameters
        ----------
        message : str
            A description of the error, providing details about the invalid combination.
        """
        self.message = message
        super().__init__(self.message)