from typing import List
from src.config.program_config import ProgramConfig
from src.config.result_checks import ResultChecks
from src.config.test_settings import TestSettings
from src.config.graph_parameters import GraphParameters


class TestCaseConfig(object):
    """Testcase configuration to run and check the results"""

    __test__ = False

    def __init__(self, test_settings: TestSettings):
        self.__path: str
        self.__max_run_time: float = 300
        self.__program_config: ProgramConfig
        self.__graph_parameters: GraphParameters
        self.__result_checks: List[ResultChecks] = []
        self.__name__: str = ""
        self.__test_settings = test_settings
        self.__success: bool = True

    @property
    def name(self) -> str:
        """Testcase name"""
        return self.__name__

    @name.setter
    def name(self, value: str):
        self.__name__ = value

    @property
    def program_config(self) -> ProgramConfig:
        """Testcase program configuration"""
        return self.__program_config

    @program_config.setter
    def program_config(self, value: ProgramConfig):
        self.__program_config = value

    @property
    def graph_parameters(self) -> GraphParameters:
        """Configuration for plotting results in a graph."""
        return self.__graph_parameters

    @graph_parameters.setter
    def graph_parameters(self, value: GraphParameters):
        self.__graph_parameters = value

    @property
    def max_run_time(self) -> float:
        """Testcase maximum run time for execution"""
        return self.__max_run_time

    @max_run_time.setter
    def max_run_time(self, value: float):
        self.__max_run_time = value

    @property
    def path(self) -> str:
        """Path to testcase folder to run program in"""
        return self.__path

    @path.setter
    def path(self, value: str):
        self.__path = value

    @property
    def result_checks(self) -> List[ResultChecks]:
        """files to check"""
        return self.__result_checks

    @property
    def test_settings(self) -> TestSettings:
        """get general settings for testcase"""
        return self.__test_settings

    @property
    def success(self) -> bool:
        """Check if issues occured within the testcase"""
        return self.__success

    @success.setter
    def success(self, value: bool):
        self.__success = value
