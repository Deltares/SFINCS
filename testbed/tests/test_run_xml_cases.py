from typing import List
import pytest_check as check
from src.config.test_case_config import TestCaseConfig
from src.suite.test_runner import TestRunner
from src.suite.comparison_runner import ComparisonRunner
from src.suite.result_plotter import ResultPlotter


def test_run_testxml(load_xmls: List[TestCaseConfig]):
    """Run an entire xml test configuration
    This test is for a local execution."""
    # TODO: add download for testcase
    for testcase in load_xmls:
        run_testcase(testcase)
        post_process(testcase)
    assert True


def test_run_testcase(load_xml_testcases: TestCaseConfig):
    """Run a single test case defined in a xml configuration."""
    # TODO: add download for testcase
    run_testcase(load_xml_testcases)
    post_process(load_xml_testcases)
    assert True


def run_testcase(testcase: TestCaseConfig):
    """Execute test case by their test configuration."""
    test_runner = TestRunner(testcase)
    test_runner.build_program()
    result = test_runner.run_program()
    check.is_not_instance(result, Exception, msg=f"{testcase.name} had an error occur while executing.")
    test_runner.move_results()


def post_process(testcase: TestCaseConfig):
    """Do all post processing actions for the results of a test run."""
    if not testcase.success:
        print("Skipped post processing due to unsuccesfull run.")
        return
    comparison_runner = ComparisonRunner(testcase)
    comparison_runner.build_comparison_paths()
    if not comparison_runner.check_result_files():
        check.fail("One or more result files specified in the xml do not exist.")
        return

    comparison_runner.compare_results()
    result_plotter = ResultPlotter(testcase)
    result_plotter.sfincs_postprocess()
