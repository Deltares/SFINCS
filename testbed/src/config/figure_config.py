from typing import List
import netCDF4
import matplotlib.pyplot as plt
class FigureConfig:
    def __init__(self):
        self.__observations: netCDF4.Dataset
        self.__reference: netCDF4.Dataset
        self.__output: netCDF4.Dataset
        self.__figure_plot: plt.Figure
        self.__total_points: int
        self.__columns: int
        self.__rows: int
        self.__locations: int
        self.__var: str

    @property
    def observations(self) -> netCDF4.Dataset:
        """Dataset that has observation data."""
        return self.__observations

    @observations.setter
    def observations(self, value: netCDF4.Dataset):
        self.__observations = value

    @property
    def reference(self) -> netCDF4.Dataset:
        """Dataset of a stable reference."""
        return self.__reference

    @reference.setter
    def reference(self, value: netCDF4.Dataset):
        self.__reference = value

    @property
    def output(self) -> netCDF4.Dataset:
        """Dataset of the executed testcase."""
        return self.__output

    @output.setter
    def output(self, value: netCDF4.Dataset):
        self.__output = value

    @property
    def figure_plot(self) -> plt.Figure:
        """Figure to plot data"""
        return self.__figure_plot

    @figure_plot.setter
    def figure_plot(self, value: plt.Figure):
        self.__figure_plot = value

    @property
    def total_points(self) -> int:
        """Total amount of locations to get data from"""
        return self.__total_points

    @total_points.setter
    def total_points(self, value: int):
        self.__total_points = value

    @property
    def columns(self) -> int:
        """Amount of columns for the figure."""
        return self.__columns

    @columns.setter
    def columns(self, value: int):
        self.__columns = value

    @property
    def rows(self) -> int:
        """Amount of rows for the figure."""
        return self.__rows

    @rows.setter
    def rows(self, value: int):
        self.__rows = value

    @property
    def locations(self) -> int:
        """Locations to check in current plot"""
        return self.__locations

    @locations.setter
    def locations(self, value: int):
        self.__locations = value

    @property
    def var(self) -> str:
        """The variable to plot from netCDF"""
        return self.__var

    @var.setter
    def var(self, value: str):
        self.__var = value
