from src.config.program_config import ProgramConfig
from src.config.result_checks import ResultChecks


class TestCaseConfig(object):
    def __init__(self):
        self.__name: str = ""
        self.__program_config: ProgramConfig = ""
        self.__result_checks: ResultChecks() = []
        self.__max_run_time: float = 300

    @property
    def name(self) -> str:
        """Testcase name"""
        return self.__name

    @name.setter
    def name(self, value: str):
        self.__name = value

    @property
    def program_config(self) -> ProgramConfig:
        """Testcase name"""
        return self.__program_config

    @program_config.setter
    def program_config(self, value: ProgramConfig):
        self.__program_config = value

    @property
    def max_run_time(self) -> float:
        """Testcase name"""
        return self.__max_run_time

    @max_run_time.setter
    def max_run_time(self, value: float):
        self.__max_run_time = value
