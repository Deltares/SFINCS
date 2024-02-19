from typing import List
from src.config.check_parameters import CheckParameters


class ResultChecks(object):
    """Checks that should be done for comparison run."""

    def __init__(self):
        self.__file: str = ""
        self.__type: str = ""
        self.__parameters: List[CheckParameters] = []

    @property
    def file(self) -> str:
        """file to compare results"""
        return self.__file

    @file.setter
    def file(self, value: str):
        self.__file = value

    @property
    def type(self) -> str:
        """type of result file for comparison"""
        return self.__type

    @type.setter
    def type(self, value: str):
        self.__type = value

    @property
    def parameters(self) -> List[CheckParameters]:
        """return parameters to check with this file."""
        return self.__parameters
