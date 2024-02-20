import subprocess
import os
from src.config.test_case_config import TestCaseConfig
from src.config.test_settings import TestSettings


class TestRunner():
    """Run the tests parsed by the xml configuration parser"""
    __test__ = False

    def __init__(self, test_case_config: TestCaseConfig):
        self.__test_case_config: TestCaseConfig = test_case_config
        self.__test_settings: TestSettings = test_case_config.test_settings
        self.__execute_cmd: str
        self.__working_directory: str
        self.__timeout: float
        self.__output: str

    def build_program(self):
        """Setup the runner environment for a program configuration."""
        path_to_engine_folder = self.__test_settings.engine_path
        path_to_program = self.__test_case_config.program_config.path
        self.__execute_cmd = os.path.join(path_to_engine_folder, path_to_program)

        path_to_testcase_folder = self.__test_settings.case_path
        path_to_testcase = self.__test_case_config.path
        self.__working_directory = os.path.join(path_to_testcase_folder, path_to_testcase)

        self.__timeout = self.__test_case_config.max_run_time

        self.__output = os.path.join("./logs", f"{self.__test_case_config.name}.txt")
        # test_settings.local_path
        # find location of testcases
        # find location of engine from testcase location
        # setup environment
        # execute engine within testcase

    def run_programs(self):
        """Run program with environment and arguments"""
        try:
            with open(self.__output, 'w', encoding="utf-8") as open_output_file:
                completed_process = subprocess.run(
                    self.__execute_cmd,
                    stdout=open_output_file,
                    check=True,
                    # env=program_env,
                    cwd=self.__working_directory,
                    timeout=self.__timeout,
                )
        except Exception as exception:
            self.__test_case_config.success = False
            return exception
        return completed_process

    def move_results(self):
        """Move results from test case to another folder for processing"""
        if not self.__test_case_config.success:
            return
        try:
            destination = self.__working_directory.replace(
                self.__test_settings.case_path, self.__test_settings.result_path)
            self.remove_result_directory(destination)
            os.makedirs(destination)
            for checks in self.__test_case_config.result_checks:
                file = os.path.join(self.__working_directory, checks.file)
                destinationfile = os.path.join(destination, checks.file)
                os.replace(file, destinationfile)
        except OSError as exception:
            self.__test_case_config.success = False
            return exception

    def remove_result_directory(self, folder):
        """Empty and remove an existing folder."""
        if os.path.isdir(folder):
            for filename in os.listdir(folder):
                file = os.path.join(folder, filename)
                os.remove(file)
            os.removedirs(folder)
