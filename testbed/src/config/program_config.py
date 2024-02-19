class ProgramConfig():
    """Configurations for programs(engines) to execute for testcases."""

    def __init__(self):
        self.__name: str = ""
        self.__path: str = ""

    @property
    def name(self) -> str:
        """Program name"""
        return self.__name

    @name.setter
    def name(self, value: str):
        self.__name = value

    @property
    def path(self) -> str:
        """Path for program"""
        return self.__path

    @path.setter
    def path(self, value: str):
        self.__path = value
