# store settings from parameters
from src.config.credentials import Credentials


class TestSettings(object):
    """General settings for the tests to run."""
    __test__ = False

    def __init__(self):
        self.__credentials: Credentials
        self.__server_base_url: str = ""
        self.__run_type: str = ""
        self.__case_path: str
        self.__reference_path: str
        self.__engine_path: str
        self.__result_path: str

    @property
    def credentials(self) -> Credentials:
        """credentials for cloud"""
        return self.__credentials

    @credentials.setter
    def credentials(self, value: Credentials):
        self.__credentials = value

    @property
    def server_base_url(self) -> str:
        """Server url to cloud storage of cases/dependencies"""
        return self.__server_base_url

    @server_base_url.setter
    def server_base_url(self, value: str):
        self.__server_base_url = value

    @property
    def case_path(self) -> str:
        """Path to the local testcase folder to run the program with."""
        return self.__case_path

    @case_path.setter
    def case_path(self, value: str):
        self.__case_path = value

    @property
    def reference_path(self) -> str:
        """Path to local references to compare against results."""
        return self.__reference_path

    @reference_path.setter
    def reference_path(self, value: str):
        self.__reference_path = value

    @property
    def engine_path(self) -> str:
        """Path to engine directory for the testbench."""
        return self.__engine_path

    @engine_path.setter
    def engine_path(self, value: str):
        self.__engine_path = value

    @property
    def result_path(self) -> str:
        """Path to engine directory for the testbench."""
        return self.__result_path

    @result_path.setter
    def result_path(self, value: str):
        self.__result_path = value

    @property
    def run_type(self) -> str:
        """Path to engine directory for the testbench."""
        return self.__run_type

    @run_type.setter
    def run_type(self, value: str):
        self.__run_type = value
