# store settings from parameters
from src.config.credentials import Credentials


class TestSettings(object):
    def __init__(self):
        self.__config_file: str
        self.__credentials: Credentials
        self.__server_base_url: str
        self.__filter: str

    @property
    def config_file(self) -> str:
        """config file to run testcases from"""
        return self.__config_file

    @config_file.setter
    def config_file(self, value: str):
        self.__config_file = value

    @property
    def credentials(self) -> Credentials:
        """credentials for cloud"""
        return self.__credentials

    @credentials.setter
    def credentials(self, value: Credentials):
        self.__credentials = value

    @property
    def server_base_url(self) -> str:
        """Server url for cloud storage of cases/dependencies"""
        return self.__server_base_url

    @server_base_url.setter
    def server_base_url(self, value: str):
        self.__server_base_url = value

    @property
    def filter(self) -> str:
        """Filter testcases to not run all within the config."""
        return self.__filter

    @filter.setter
    def filter(self, value: str):
        self.__filter = value
