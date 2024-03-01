from typing import Dict, List, Any
from lxml import etree
from src.config.program_config import ProgramConfig
from src.config.test_case_config import TestCaseConfig
from src.config.result_checks import ResultChecks
from src.config.check_parameters import CheckParameters
from src.config.test_settings import TestSettings
from src.config.graph_parameters import GraphParameters


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
        """Create a dictionary from parsed elements in the xml for iteration"""
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
        """Create a parser and return parsed tree for xml file with xinclude"""
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
        self.__testcase_config: List[TestCaseConfig] = []

    @property
    def testcase_config(self):
        """Retrieve list of parse test case configurations"""
        return self.__testcase_config

    def validate(self, xml_file):
        """Validate the xml by using the testbed schema"""
        xmlschema_path = etree.parse("configs/xsd/testbed.xsd")
        xmlschema = etree.XMLSchema(xmlschema_path)

        result = xmlschema.validate(self.__xml_file)

        if result:
            print(f"{xml_file} validated succesfully.")
        else:
            print(f"{xml_file} gave a validation error.")
            assert result, xmlschema.error_log.filter_from_errors()[0]

    def parse_config(self, test_settings: TestSettings, case_filter: str = ""):
        """Parse the specified configuration from paramater"""
        for config in self.loop(self.__parsed_tree, "config"):
            self.parse_settings(config, test_settings)
        for programs in self.loop(self.__parsed_tree, "programs"):
            for program in programs["program"]:
                program_config = self.parse_program(program)
                self.__program_configs.append(program_config)
        for testcases in self.loop(self.__parsed_tree, "testCases"):
            for testcase_xml in testcases["testCase"]:
                if "name" in testcase_xml:
                    if case_filter not in str(testcase_xml["name"][0]):
                        continue
                testcase = self.parse_testcase(testcase_xml, test_settings)
                self.__testcase_config.append(testcase)
        for plot_settings in self.loop(self.__parsed_tree, "plotSettings"):
            for plot_setting in plot_settings["plotSetting"]:
                self.parse_plot_settings(plot_setting)

    def parse_settings(self, config, test_settings: TestSettings):
        """Parse general test settings for testcases."""
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

    def parse_testcase(self, testcase_xml, test_settings: TestSettings):
        """"Parse testcase and link to program config"""
        testcase = TestCaseConfig(test_settings)
        if "name" in testcase_xml:
            testcase.name = str(testcase_xml["name"][0])
        if "ref" in testcase_xml:
            ref_name = str(testcase_xml["ref"][0])
            testcase.program_config = self.find_program_reference(ref_name)
        if "path" in testcase_xml:
            testcase.path = testcase_xml["path"][0]["txt"]
        if "maxRunTime" in testcase_xml:
            testcase.max_run_time = float(testcase_xml["maxRunTime"][0]["txt"])
        for checks in self.loop(testcase_xml, "checks"):
            for file in checks["file"]:
                parsed_checks = self.parse_checks(file)
                testcase.result_checks.append(parsed_checks)
        return testcase

    def parse_checks(self, checks):
        """parse file checks for processing results of testcase"""
        result_checks = ResultChecks()
        if "name" in checks:
            result_checks.file = checks["name"][0]
        if "type" in checks:
            result_checks.type = checks["type"][0]
        for parameters in self.loop(checks, "parameters"):
            for parameter in parameters["parameter"]:
                parsed_parameters = self.parse_parameters(parameter)
                result_checks.parameters.append(parsed_parameters)
        return result_checks

    def parse_parameters(self, parameter):
        """parse parameters associated with result file"""
        check_parameters = CheckParameters()
        if "name" in parameter:
            check_parameters.name = parameter["name"][0]
        if "location" in parameter:
            check_parameters.location = parameter["location"][0]
        if "toleranceAbsolute" in parameter:
            check_parameters.tolerance_absolute = float(parameter["toleranceAbsolute"][0])
        return check_parameters

    def parse_plot_settings(self, plot_settings):
        graph_parameters = GraphParameters()
        if "observations" in plot_settings:
            if plot_settings["observations"][0].lower() == "true":
                graph_parameters.observations = True
        if "ref" in plot_settings:
            testcase = self.find_testcase_reference(str(plot_settings["ref"][0]))
            testcase.graph_parameters = graph_parameters
        if "his" in plot_settings:
            his = plot_settings["his"][0]
            graph_parameters.his = True
            if "location" in his:
                graph_parameters.his_loc = self.get_float_array(str(his["location"][0]))
            if "variable" in his:
                graph_parameters.his_var = str(his["variable"][0])
        if "map" in plot_settings:
            plot_map = plot_settings["map"][0]
            if "mapType" in plot_map:
                graph_parameters.map = str(plot_map["mapType"][0])
            if "variable" in plot_map:
                graph_parameters.map1D_var = str(plot_map["variable"][0])
            if "yloc" in plot_map:
                graph_parameters.map1D_yloc = str(plot_map["yloc"][0])
            if "time" in plot_map:
                graph_parameters.map1D_t = self.get_float_array(str(plot_map["time"][0]))
            if "variable2D" in plot_map:
                graph_parameters.map2D_var = str(plot_map["variable2D"][0])
            if "clim" in plot_map:
                graph_parameters.clim_2D = self.get_float_array(plot_map["clim"][0])
        if "axes" in plot_settings:
            axes = plot_settings["axes"][0]
            if "xlimit" in axes:
                graph_parameters.xlim = self.get_float_array(str(axes["xlimit"][0]))
            if "ylimit" in axes:
                graph_parameters.ylim = self.get_float_array(str(axes["ylimit"][0]))
            if "xlabel" in axes:
                graph_parameters.xlabel = str(axes["xlabel"][0])
            if "ylabel" in axes:
                graph_parameters.ylabel = str(axes["ylabel"][0])

    def find_program_reference(self, name):
        """find program for testcase by referenced name in program list"""
        for program in self.__program_configs:
            if program.name == name:
                return program
        assert False, "Failed to find program referenced by testcase"

    def find_testcase_reference(self, name):
        """find plot settings for testcase by referenced name in testcase list"""
        for testcase in self.testcase_config:
            if testcase.name == name:
                return testcase
        assert False, "Failed to find testcase referenced by latex plotter"

    def get_float_array(self, value):
        stripped_value = value.strip("[]")
        array = stripped_value.split(",")
        float_array = []
        for number in array:
            float_array.append(float(number))
        return float_array
