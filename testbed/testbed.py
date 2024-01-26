from src.suite.test_runner import TestRunner
from src.utils.argument_parser import CmdArgumentParser


if __name__ == "__main__":
    settings = CmdArgumentParser().parse_arguments_to_settings()
    test_runner = TestRunner(settings)
    test_runner.run()
