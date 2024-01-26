from typing import Dict, List, Any
from lxml import etree
from src.suite.test_settings import TestSettings
from src.config.program_config import ProgramConfig
from src.config.test_case_config import TestCaseConfig
# parse test case configurations


class ConfigParser(object):
    "Parse xml configurations for the runner"

    def loop(self, dictionary: Dict[str, Any], key: str) -> List:
        "Create a dictionary list for looping."
        if key in dictionary:
            if isinstance(dictionary[key], list):
                return dictionary[key]
            elif isinstance(dictionary[key], dict):
                return list(dictionary[key].values())
            else:
                return [dictionary[key]]
        else:
            return []

    def __init__(self, settings: TestSettings):
        self.__xml_path = settings.config_file
        self.__xml_file = etree.parse(self.__xml_path)
        self.__program_configs: List[ProgramConfig] = []
        self.__test_case_config: List[TestCaseConfig] = []

    def validate(self):
        """Validate the xml to be run against the schema"""
        xmlschema_path = etree.parse("configs/xsd/testbed.xsd")
        xmlschema = etree.XMLSchema(xmlschema_path)

        result = xmlschema.validate(self.__xml_file)

        if (result):
            print("XML configuration validation correct.")
        else:
            print("XML configuration validation error.")
            assert result, xmlschema.error_log.filter_from_errors()[0]

    def parse_config(self):
        """Parse the specified configuration from paramater"""
        for programs in self.loop(self.__xml_file, "programs"):
            for program in programs["program"]:
                program_config = self.parse_program(program)
                self.__program_configs.append(program_config)

    def parse_program(self, program):
        """Parse program list for testcases within the configuration"""
        program_config = ProgramConfig()
        if "name" in program:
            program_config.name = str(program["name"][0])
        if "location" in program:
            if "from" in program["location"]:
                program_config.local_dir = str(program["location"]["from"]["txt"])
        if "path" in program:
            program_config.path = str(program["path"]["txt"])
        return program_config
