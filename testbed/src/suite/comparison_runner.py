import os
from typing import List
from src.config.test_case_config import TestCaseConfig
from src.config.test_settings import TestSettings
from src.config.result_checks import ResultChecks

class ComparisonRunner():
    """Compare results against reference after a run is completed."""

    def __init__(self, testcase: TestCaseConfig):
        self.__testcase_config = testcase
        self.__test_settings: TestSettings = testcase.test_settings
        self.__result_checks: List[ResultChecks] = testcase.result_checks
        self.__reference_path: str
        self.__result_path: str

    def build_comparison_paths(self):
        reference_folder = self.__test_settings.reference_path
        testcase_path = self.__testcase_config.path
        self.__reference_path = os.path.join(reference_folder, testcase_path)

        result_folder = self.__test_settings.result_path
        self.__result_path = os.path.join(result_folder, testcase_path)

    def check_result_files(self) -> bool:
        """Check if all the files specified to check in the xml exist."""
        for result_check in self.__result_checks:
            check_file = os.path.join(self.__result_path, result_check.file)
            if not os.path.exists(check_file):
                return False
        return True

    def compare_results(self):
        if not self.__test_settings.compare:
            print("comparison was skipped due to argument.")
        
