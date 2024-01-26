from src.utils.config_parser import ConfigParser
from src.suite.test_settings import TestSettings

# interpret settings and run testcases defined by xml config


class TestRunner():
    def __init__(self, settings: TestSettings):
        self.__settings = settings

    def run(self):
        config_parser = ConfigParser(self.__settings)
        config_parser.validate()
        #config_parser.parse_config()

    def run_test_case(self):
        pass