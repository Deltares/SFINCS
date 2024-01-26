class CheckParameters(object):
    """Parameters to use when comparing values."""
    def __init__(self):
        self.__name: str = ""
        self.__location: str = ""
        self.__tolerance_absolute: float = 0

    @property
    def name(self) -> str:
        """Name of the parameter"""
        return self.__name

    @name.setter
    def name(self, value: str):
        self.__name = value

    @property
    def location(self) -> str:
        """check location"""
        return self.__location

    @location.setter
    def location(self, value: str):
        self.__location = value

    @property
    def tolerance_absolute(self) -> float:
        """Absolute difference tolerated for comparison"""
        return self.__tolerance_absolute

    @tolerance_absolute.setter
    def tolerance_absolute(self, value: float):
        self.__tolerance_absolute = value
