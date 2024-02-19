import pytest_check as check
from src.utils.config_parser import ConfigParser
from src.config.test_settings import TestSettings
from src.suite.test_runner import TestRunner


def test_run_testxml(load_xmls,
                     test_settings: TestSettings):
    config_parser = ConfigParser(load_xmls)
    config_parser.validate()
    config_parser.parse_config(test_settings)
    run_test_cases(config_parser, test_settings)
    assert True


def test_run_testcase(load_xml_testcases: ConfigParser,
                      test_settings: TestSettings,
                      request):
    config_parser = load_xml_testcases
    config_parser.validate()
    case_filter = request.node.callspec.id.split(" - ")[-1]
    config_parser.parse_config(test_settings, case_filter)
    run_test_cases(config_parser, test_settings)
    assert True


def run_test_cases(config_parser: ConfigParser,
                   test_settings: TestSettings):
    for test_case in config_parser.test_case_config:
        test_runner = TestRunner(test_case, test_settings)
        test_runner.build_program()
        result = test_runner.run_programs()
        check.is_not_instance(result, Exception, msg="Error occured while executing testcase.")
        print(f"finished running {test_case.name}")
