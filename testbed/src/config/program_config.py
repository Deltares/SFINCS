class ProgramConfig(object):
    def __init__(self):
        self.__name: str = ""
        self.__path: str = ""
        self.__local_dir: str = ""

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

    @property
    def local_dir(self) -> str:
        """WorkDirectory for program"""
        return self.__local_dir

    @local_dir.setter
    def local_dir(self, value: str):
        self.__local_dir = value