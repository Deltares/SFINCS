from typing import Dict, List
from src.config.check_parameters import CheckParameters

class ResultChecks(object):
    """Checks that should be done for comparison run."""

    def __init__(self):
        self.__file: str = ""
        self.__type: str = ""
        self.__parameters: Dict[str, List[CheckParameters]] = {}