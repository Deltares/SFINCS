from typing import List


class GraphParameters():
    """Values to plot on graphics for latex documentation."""

    def __init__(self):
        self.__his: bool = False
        self.__his_var: str = "point_zs"
        self.__his_loc: List[int] = None

        self.__map: str = ""
        self.__map1D_t: List[int] = None
        self.__map1D_var: str = "zs"
        self.__map1D_yloc: int = 0

        self.__observations: bool = False
        self.__datalabel: str = "Analytical"

        self.__title: str
        self.__xlim: List[int] = None
        self.__ylim: List[int] = None
        self.__xlabel: str
        self.__ylabel: str

        self.__title_2D: str
        self.__xlim_2D: List[int]
        self.__ylim_2D: List[int]
        self.__xlabel_2D: str
        self.__ylabel_2D: str
        self.__clim_2D: List[int] = [-5, 5]
        self.__map2D_var: str = "zsmax"

    @property
    def his(self) -> bool:
        """Plot his file"""
        return self.__his

    @his.setter
    def his(self, value: bool):
        self.__his = value

    @property
    def his_var(self) -> str:
        """Variable to plot for his"""
        return self.__his_var

    @his_var.setter
    def his_var(self, value: str):
        self.__his_var = value

    @property
    def his_loc(self) -> bool:
        """Locations to plot from his"""
        return self.__his_loc

    @his_loc.setter
    def his_loc(self, value: bool):
        self.__his_loc = value

    @property
    def map(self) -> str:
        """Plot map type 1D or 2D"""
        return self.__map

    @map.setter
    def map(self, value: str):
        self.__map = value

    @property
    def map1D_t(self) -> List[int]:
        """Time series for plot"""
        return self.__map1D_t

    @map1D_t.setter
    def map1D_t(self, value: List[int]):
        self.__map1D_t = value

    @property
    def map1D_var(self) -> str:
        """1D variable to use for plotting"""
        return self.__map1D_var

    @map1D_var.setter
    def map1D_var(self, value: str):
        self.__map1D_var = value

    @property
    def map1D_yloc(self) -> int:
        return self.__map1D_yloc

    @map1D_yloc.setter
    def map1D_yloc(self, value: int):
        self.__map1D_yloc = value

    @property
    def observations(self) -> bool:
        return self.__observations

    @observations.setter
    def observations(self, value: bool):
        self.__observations = value

    @property
    def datalabel(self) -> str:
        return self.__datalabel

    @datalabel.setter
    def datalabel(self, value: str):
        self.__datalabel = value

    @property
    def title(self) -> str:
        return self.__title

    @title.setter
    def title(self, value: str):
        self.__title = value

    @property
    def xlim(self) -> List[int]:
        return self.__xlim

    @xlim.setter
    def xlim(self, value: List[int]):
        self.__xlim = value

    @property
    def ylim(self) -> List[int]:
        return self.__ylim

    @ylim.setter
    def ylim(self, value: List[int]):
        self.__ylim = value

    @property
    def xlabel(self) -> str:
        return self.__xlabel

    @xlabel.setter
    def xlabel(self, value: str):
        self.__xlabel = value

    @property
    def ylabel(self) -> str:
        return self.__ylabel

    @ylabel.setter
    def ylabel(self, value: str):
        self.__ylabel = value

    @property
    def title_2D(self) -> str:
        return self.__title_2D

    @title_2D.setter
    def title_2D(self, value: str):
        self.__title_2D = value

    @property
    def xlim_2D(self) -> List[int]:
        return self.__xlim_2D

    @xlim_2D.setter
    def xlim_2D(self, value: List[int]):
        self.__xlim_2D = value

    @property
    def ylim_2D(self) -> List[int]:
        return self.__ylim_2D

    @ylim_2D.setter
    def ylim_2D(self, value: List[int]):
        self.__ylim_2D = value

    @property
    def xlabel_2D(self) -> str:
        return self.__xlabel_2D

    @xlabel_2D.setter
    def xlabel_2D(self, value: str):
        self.__xlabel_2D = value

    @property
    def ylabel_2D(self) -> str:
        return self.__ylabel_2D

    @ylabel_2D.setter
    def ylabel_2D(self, value: str):
        self.__ylabel_2D = value

    @property
    def clim_2D(self) -> List[int]:
        return self.__clim_2D

    @clim_2D.setter
    def clim_2D(self, value: List[int]):
        self.__clim_2D = value

    @property
    def map2D_var(self) -> str:
        return self.__map2D_var

    @map2D_var.setter
    def map2D_var(self, value: str):
        self.__map2D_var = value
