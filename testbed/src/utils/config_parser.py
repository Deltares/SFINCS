from typing import Dict, List, Any
from lxml import etree
from src.config.program_config import ProgramConfig
from src.config.test_case_config import TestCaseConfig
from src.config.result_checks import ResultChecks
from src.config.check_parameters import CheckParameters
from src.config.test_settings import TestSettings


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
            if isinstance(ch.tag, str):
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

    @property
    def test_case_config(self):
        """Retrieve list of parse test case configurations"""
        return self.__test_case_config

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

    def generate_testcases(self, config_file):
        """Get all the testcase names within the xml_file"""
        testcase_names = []
        for testcases in self.loop(self.__parsed_tree, "testCases"):
            for testcase in testcases["testCase"]:
                if "name" in testcase:
                    testcase_names.append(config_file + " - " + str(testcase["name"][0]))
        return testcase_names

    def parse_config(self, test_settings: TestSettings, case_filter: str = ""):
        """Parse the specified configuration from paramater"""
        for config in self.loop(self.__parsed_tree, "config"):
            self.parse_settings(config, test_settings)
        for programs in self.loop(self.__parsed_tree, "programs"):
            for program in programs["program"]:
                program_config = self.parse_program(program)
                self.__program_configs.append(program_config)
        for testcases in self.loop(self.__parsed_tree, "testCases"):
            for testcase in testcases["testCase"]:
                if "name" in testcase:
                    if case_filter not in str(testcase["name"][0]):
                        continue
                test_case = self.parse_testcase(testcase, test_settings)
                self.__test_case_config.append(test_case)

    def parse_settings(self, config, test_settings: TestSettings):
        if "serverUrl" in config:
            if test_settings.server_base_url == "":
                test_settings.server_base_url = str(config["serverUrl"][0]["txt"])
        if "localPaths" in config:
            local_path = config["localPaths"][0]
            if "casesDir" in local_path:
                test_settings.case_path = str(local_path["casesDir"][0]["txt"])
            if "referenceDir" in local_path:
                test_settings.reference_path = str(local_path["referenceDir"][0]["txt"])
            if "enginesDir" in local_path:
                test_settings.engine_path = str(local_path["enginesDir"][0]["txt"])
            if "resultsDir" in local_path:
                test_settings.result_path = str(local_path["resultsDir"][0]["txt"])

    def parse_program(self, program):
        """Parse program for testcases within the configuration"""
        program_config = ProgramConfig()
        if "name" in program:
            program_config.name = str(program["name"][0])
        if "path" in program:
            program_config.path = str(program["path"][0]["txt"])
        return program_config

    def parse_testcase(self, testcase, test_settings: TestSettings):
        """"Parse testcase and link to program config"""
        test_case = TestCaseConfig(test_settings)
        if "name" in testcase:
            test_case.name = str(testcase["name"][0])
        if "ref" in testcase:
            ref_name = str(testcase["ref"][0])
            test_case.program_config = self.find_program_reference(ref_name)
        if "path" in testcase:
            test_case.path = testcase["path"][0]["txt"]
        if "maxRunTime" in testcase:
            test_case.max_run_time = float(testcase["maxRunTime"][0]["txt"])
        for checks in self.loop(testcase, "checks"):
            for file in checks["file"]:
                parsed_checks = self.__parse_checks__(file)
                test_case.result_checks.append(parsed_checks)
        return test_case

    def __parse_checks__(self, checks):
        """parse file checks for processing results of testcase"""
        result_checks = ResultChecks()
        if "name" in checks:
            result_checks.file = checks["name"][0]
        if "type" in checks:
            result_checks.type = checks["type"][0]
        for parameters in self.loop(checks, "parameters"):
            for parameter in parameters["parameter"]:
                parsed_parameters = self.__parse_parameters__(parameter)
                result_checks.parameters.append(parsed_parameters)
        return result_checks

    def __parse_parameters__(self, parameter):
        """parse parameters associated with result file"""
        check_parameters = CheckParameters()
        if "name" in parameter:
            check_parameters.name = parameter["name"][0]
        if "location" in parameter:
            check_parameters.location = parameter["location"][0]
        if "toleranceAbsolute" in parameter:
            check_parameters.tolerance_absolute = float(parameter["toleranceAbsolute"][0])
        return check_parameters

    def find_program_reference(self, name):
        """find program for testcase by referenced name in program list"""
        for program in self.__program_configs:
            if program.name == name:
                return program
        assert False, "Failed to find program referenced by testcase"
