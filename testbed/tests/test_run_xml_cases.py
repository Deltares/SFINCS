from typing import List
import pytest_check as check
from src.config.test_case_config import TestCaseConfig
from src.suite.test_runner import TestRunner


def test_run_testxml(load_xmls: List[TestCaseConfig]):
    """Run an entire xml test configuration
    This test is for a local execution."""
    for test_case in load_xmls:
        run_test_cases(test_case)
    assert True


def test_run_testcase(load_xml_testcases: TestCaseConfig):
    """Run a single test case defined in a xml configuration."""
    run_test_cases(load_xml_testcases)
    assert True


def run_test_cases(test_case: TestCaseConfig):
    """Execute test case by their test configuration."""
    test_runner = TestRunner(test_case)
    test_runner.build_program()
    result = test_runner.run_programs()
    check.is_not_instance(result, Exception, msg=f"{test_case.name} had an error occur while executing.")
    test_runner.move_results()
