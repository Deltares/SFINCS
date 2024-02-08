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

    def branch(self, xml_tree: etree._ElementTree, prefix: str) -> Dict[str, Any]:
        new_tree = {"txt": xml_tree.text}

        for ch in xml_tree.getchildren():
            if type(ch.tag) == str:
                branch_name = ch.tag.replace(prefix, "")

                if branch_name not in new_tree:
                    new_tree[branch_name] = []

                new_tree[branch_name].append(self.branch(ch, prefix))

        for key, val in xml_tree.attrib.items():
            new_tree[key.replace(prefix, "")] = [val]

        return new_tree

    def __maketree__(self, path):
        parser = etree.XMLParser(remove_blank_text=True, attribute_defaults=False)
        xml_file: etree._ElementTree = etree.parse(path, parser)

        xml_file.xinclude()
        root_node: etree._Element = xml_file.getroot()

        prefix = "{%s}" % root_node.nsmap

        parsed_tree = self.branch(root_node, prefix)
        return [xml_file, parsed_tree]

    def __init__(self, config_file: str):
        tree_values = self.__maketree__(config_file)
        self.__xml_file: etree._ElementTree = tree_values[0]
        self.__parsed_tree: Dict[str, Any] = tree_values[1]
        self.__program_configs: List[ProgramConfig] = []
        self.__test_case_config: List[TestCaseConfig] = []

    def validate(self):
        """Validate the xml by using the testbed schema"""
        xmlschema_path = etree.parse("configs/xsd/testbed.xsd")
        xmlschema = etree.XMLSchema(xmlschema_path)

        result = xmlschema.validate(self.__xml_file)

        if (result):
            print("XML configuration validation correct.")
        else:
            print("XML configuration validation error.")
            assert result, xmlschema.error_log.filter_from_errors()[0]

    def generate_testcases(self):
        testcase_names = []
        for testcases in self.loop(self.__parsed_tree, "testCases"):
            for testcase in testcases["testCase"]:
                if "name" in testcase:
                    testcase_names.append(str(testcase["name"][0]))
        print(testcase_names)
        return testcase_names

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
