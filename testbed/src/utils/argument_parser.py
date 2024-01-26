import getpass
from typing import Any, Optional
from argparse import ArgumentParser, Namespace
from src.suite.test_settings import TestSettings
from src.config.credentials import Credentials


class CmdArgumentParser():
    """Parse arguments from the commandline."""

    def parse_arguments_to_settings(cls) -> TestSettings:
        """Parse arguments and set them in settings for runner settings."""
        settings = TestSettings()
        parser = cls.__create_argument_parser()
        args: Namespace = parser.parse_args()
        credentials = cls.__get_credentials(args)

        # Additionally, extra TeamCity messages will be produced.
        server_base_url = cls.__get_argument_value("server_base_url", args) or ""

        config_file = cls.__get_argument_value("config", args)
        settings.credentials = credentials
        settings.server_base_url = server_base_url
        settings.config_file = config_file
        return settings

    @classmethod
    def __get_argument_value(
        cls,
        name: str,
        args: Namespace,
        is_interactive: bool = False,
        secret_value: bool = False,
    ) -> Optional[Any]:
        return_value = None

        if hasattr(args, name):
            return_value = getattr(args, name)

        if not return_value and is_interactive:
            if secret_value:
                return getpass.getpass(f"{name} : ")

            return input(f"{name}")

        return return_value

    @classmethod
    def __get_credentials(cls, args: Namespace) -> Credentials:
        credentials = Credentials()
        credentials.name = "commandline"

        is_interactive = cls.__get_argument_value("interactive", args) or False
        credentials.username = (
            cls.__get_argument_value("username", args, is_interactive) or ""
        )
        if credentials.username == "":
            print(
                'No username on commandline. add "-i" to enable interactive input'
            )

        credentials.password = (
            cls.__get_argument_value("password", args, is_interactive, True) or ""
        )
        if credentials.password == "":
            print(
                'No password on commandline. add "-i" to enable interactive input'
            )

        return credentials

    @classmethod
    def __create_argument_parser(cls) -> ArgumentParser:
        parser = ArgumentParser(
            description="TestBed, test runner for testcases.",
            add_help=True,
        )

        # Compulsory arguments
        parser.add_argument(
            "--config",
            default="",
            help="Path to config file",
            dest="config",
            required=True,
        )
        
        # Optional arguments
        parser.add_argument(
            "--filter",
            default="",
            help="Specify what tests to run based on filter",
            dest="filter",
        )
        parser.add_argument(
            "--log-level",
            default="",
            help="CRITICAL, ERROR, WARNING, INFO or DEBUG",
            dest="loglevel",
        )
        parser.add_argument(
            "--teamcity",
            action="store_true",
            help="Turns on specific TeamCity logging.",
            dest="teamcity",
        )
        parser.add_argument(
            "--server-base-url",
            help="e.g. SVN, S3, Git LFS",
            default="https://repos.deltares.nl/repos/DSCTestbench/trunk",
            required=False,
            dest="server_base_url",
        )
        parser.add_argument(
            "--username",
            help="Server username (MinIO).",
            default=None,
            # required=True,
            dest="username",
        )
        parser.add_argument(
            "--password",
            help="Server password (MinIO).",
            default=None,
            # required=True,
            dest="password",
        )
        parser.add_argument(
            "-i",
            "--interactive",
            action="store_true",
            help="Must be True to enable username/password via keyboard.",
            dest="interactive",
        )

        return parser
