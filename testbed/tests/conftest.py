import pytest
import os
from src.config.credentials import Credentials
from src.suite.test_settings import TestSettings
from src.utils.config_parser import ConfigParser

def pytest_addoption(parser):
    """Parse commandline parameters"""
    parser.addoption(
        "--config", action="store", default="noconfig",
        help="Specifies an xml config to use for running?"
    )


@pytest.fixture(scope="module")
def xml_config(request):
    """Get specific parameter for xml config filter"""
    return request.config.getoption("--config")


@pytest.fixture(scope="module")
def setup_credentials():
    """Setup a set of credentials from commandline or enviornment"""
    credentials = Credentials()
    credentials.name = "python"
    # ToDo: parse these parameters either from environment or arguments
    credentials.username = "test"
    credentials.password = "test"
    return credentials


@pytest.fixture(scope="module")
def init_settings(setup_credentials):
    """Use standardized settings for tests."""
    settings = TestSettings()
    settings.credentials = setup_credentials
    settings.server_base_url = "https://s3.deltares.nl"
    # ToDo: add filter to parameter from commandline or environment
    settings.filter = ""
    return settings


def pytest_generate_tests(metafunc):
    xml_file = metafunc.config.getoption('--config')
    xml_files = []
    if xml_file == "noconfig":
        xml_files = load_xmls()
    else:
        xml_files.append(xml_file)
    if 'load_xmls' in metafunc.fixturenames:
        metafunc.parametrize('load_xmls', xml_files)
    if 'load_xml_testcases' in metafunc.fixturenames:
        testcases = parse_testcases_from_xml(xml_files)
        metafunc.parametrize('load_xml_testcases', testcases)


def parse_testcases_from_xml(xml_files):
    testcases = []
    for xml_file_path in xml_files:
        config_parser = ConfigParser(xml_file_path)
        testcases = config_parser.generate_testcases()
    return testcases
    #return ["parent 1 - test 1", "parent 2 - test 2"]


def load_xmls():
    xml_files = []
    config_folder = "configs"
    for filename in os.listdir(config_folder):
        if filename.endswith('.xml'):
            xml_files.append(os.path.join(config_folder, filename))
    return xml_files
