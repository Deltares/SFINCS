import os
from typing import List
from src.config.credentials import Credentials
from src.config.test_settings import TestSettings
from src.utils.config_parser import ConfigParser


def pytest_addoption(parser):
    """Parse commandline parameters"""
    parser.addoption(
        "--config", action="store", default="noconfig",
        help="Specifies an xml config to use for running."
    )
    parser.addoption(
        "--reference", action="store_true",
        help="Does a reference run and skips results processing."
    )
    parser.addoption(
        "--username", action="store",
        help="Your username to use for downloading testcases and references."
    )
    parser.addoption(
        "--password", action="store",
        help="Your password to use for downloading testcases and references."
    )


def setup_credentials():
    """Setup a set of credentials from commandline or enviornment"""
    credentials = Credentials()
    credentials.name = "python"
    # ToDo: parse these parameters either from environment or arguments
    credentials.username = "test"
    credentials.password = "test"
    return credentials


def init_settings():
    """Use standardized settings for tests."""
    settings = TestSettings()
    settings.credentials = setup_credentials()
    settings.server_base_url = "https://s3.deltares.nl"
    return settings


def pytest_generate_tests(metafunc):
    """Generate dynamic testcase from configuration for test pytest tests"""
    xml_file = metafunc.config.getoption('--config')
    test_settings = init_settings()
    xml_files = find_xml_files(xml_file)
    if 'load_xmls' in metafunc.fixturenames:
        grouped_testcases = []
        for xml_file_path in xml_files:
            testcases = parse_testcases_from_xml([xml_file_path], test_settings)
            grouped_testcases.append(testcases)
        metafunc.parametrize('load_xmls', grouped_testcases, ids=xml_files)
    if 'load_xml_testcases' in metafunc.fixturenames:
        testcases = parse_testcases_from_xml(xml_files, test_settings)
        metafunc.parametrize('load_xml_testcases', testcases)


def parse_testcases_from_xml(xml_files: List[str], test_settings: TestSettings):
    """Returns a list of parsed test case configurations from xml."""
    testcases = []
    for xml_file_path in xml_files:
        config_parser = ConfigParser(xml_file_path)
        config_parser.validate(xml_file_path)
        config_parser.parse_config(test_settings)
        testcases.extend(config_parser.test_case_config)
    return testcases


def find_xml_files(xml_config_filter):
    """return configurations based on configs folder"""
    xml_files = []
    if xml_config_filter != "noconfig":
        # if a configuration is specified return it as a list.
        return [xml_config_filter]
    config_folder = "configs"
    files = os.listdir(config_folder)
    for filename in files:
        if filename.endswith('.xml'):
            xml_file = os.path.join(config_folder, filename)
            xml_files.append(xml_file)
    return xml_files
