import os
import shutil
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import xarray
import pandas as pd
from scipy import interpolate
from src.config.test_case_config import TestCaseConfig
from src.config.figure_config import FigureConfig


class ResultPlotter:
    """Plot a comparison graph into a latex file for specified data."""

    def __init__(self, testcase_config: TestCaseConfig):
        self.__testcase_config = testcase_config
        self.__test_settings = testcase_config.test_settings
        self.__graph_parameters = testcase_config.graph_parameters
        self.figure_folder: str
        self.factor = 3600

    def sfincs_postprocess(self):
        """Postprocess data from run with reference and/or observations."""
        self.figure_folder = os.path.join(self.__test_settings.figure_path,
                                          self.__testcase_config.path)
        self.prepare_directory(self.figure_folder)

        if self.__graph_parameters.his:
            figure = self.prepare_dataset("sfincs_his.nc")
            self.prepare_plot(figure, "his")
            self.plot_dataset(figure)

        if self.__graph_parameters.map == "1D":
            figure = self.prepare_dataset("sfincs_map.nc")
            self.prepare_plot(figure, self.__graph_parameters.map)
            self.plot_dataset(figure)

        if self.__graph_parameters.map == "2D":
            figure = self.prepare_dataset("sfincs_map.nc")
            self.plot_2D_map(figure)

    def prepare_dataset(self, file: str):
        """Open datasets and create a figure config for plotting."""

        figure = FigureConfig()
        reference_path = os.path.join(self.__test_settings.reference_path, self.__testcase_config.path)
        if self.__graph_parameters.observations:
            observation_file = os.path.join(reference_path, "observations.nc")
            figure.observations = nc.Dataset(observation_file)

        reference_file = os.path.join(reference_path, file)
        figure.reference = nc.Dataset(reference_file)

        case_path = os.path.join(self.__test_settings.result_path, self.__testcase_config.path)
        outputfile = os.path.join(case_path, file)
        figure.output = nc.Dataset(outputfile)
        return figure

    def prepare_plot(self, figure: FigureConfig, type: str):
        """Prepare dimensions of figure and create it."""
        if type == "his":
            total_points = self.identify_locations_his(figure)
            figure.var = self.__graph_parameters.his_var
        elif type == "1D":
            total_points = len(self.__graph_parameters.map1D_t)
            figure.locations = self.__graph_parameters.map1D_t
            figure.var = self.__graph_parameters.map1D_var
        figure.total_points = total_points

        if total_points > 3:
            columns = 2
        elif total_points > 10:
            columns = 3
        else:
            columns = 1
        # Compute Rows required
        rows = total_points // columns
        # Add a row if necessary
        if total_points % columns != 0:
            rows += 1
        figure.columns = columns
        figure.rows = rows

        size = self.set_size(width="a4", fraction=0.9, subplots=(rows, columns))
        figure.figure_plot = plt.figure(figsize=size)

    def identify_locations_his(self, figure: FigureConfig):
        """Get locations to plot for his file by provided or netCDF input"""
        if self.__graph_parameters.his_loc:
            total_points = len(self.__graph_parameters.his_loc)
            figure.locations = self.__graph_parameters.his_loc
        else:
            total_points = figure.output["station_id"].size
            figure.locations = figure.output["station_id"][:]
        return total_points

    def plot_dataset(self, figure: FigureConfig):
        """Plot dataset by input from FigureConfig object"""
        for i in range(figure.total_points):
            current_loc = figure.locations[i]

            axes = figure.figure_plot.add_subplot(figure.rows, figure.columns, i + 1)
            axes.grid()
            self.set_limits(figure, axes)

            if self.__graph_parameters.his:
                self.plot_axes_his(axes, figure, current_loc - 1)
            else:
                self.plot_axes_1D(axes, figure, current_loc)

            # Set labels only for outer subplots
            if (i + 1) / figure.columns > figure.rows - 1:
                axes.set_xlabel(self.__graph_parameters.xlabel)
            else:
                axes.xaxis.set_ticklabels([])
            if (i + 1) % figure.columns != 0 or figure.columns == 1:
                axes.set_ylabel(self.__graph_parameters.ylabel)
            else:
                axes.yaxis.set_ticklabels([])

        if i == figure.total_points - 1:
            axes.legend(loc='upper right')
        plt.suptitle(self.__testcase_config.name)

        figure_path = os.path.join(self.figure_folder,
                                   self.__testcase_config.name + '.png')

        plt.savefig(figure_path, dpi=300)

    def plot_axes_his(self, axes: plt.Axes, figure: FigureConfig, current_loc: int):
        """Plot his data on axes."""
        axes.plot(figure.reference["time"][:] / self.factor,
                  figure.reference[figure.var][:, current_loc],
                  color='C2',
                  label=figure.reference.getncattr("Build-Revision"))

        axes.plot(figure.output["time"][:] / self.factor,
                  figure.output[figure.var][:, current_loc],
                  color='C1',
                  label=figure.output.getncattr("Build-Revision"))

        if self.__graph_parameters.observations:
            axes.plot(figure.observations["time"][:] / self.factor,
                      figure.observations[figure.var][:, current_loc],
                      color='k', linestyle="--",
                      label=self.__graph_parameters.datalabel)

        axes_title = str(nc.chartostring(figure.output["station_name"][current_loc, :]))
        axes.set_title(axes_title.strip())

    def plot_axes_1D(self, axes: plt.Axes, figure: FigureConfig, current_loc: int):
        """Plot 1D data on axes"""
        ref_index = np.argmin(np.abs(np.array(figure.reference["time"]) - current_loc))  
        axes.plot(np.transpose(figure.reference["x"][:]),
                  np.squeeze(figure.reference[figure.var][ref_index, self.__graph_parameters.map1D_yloc, :]),
                  color='C2',
                  label=figure.reference.getncattr("Build-Revision"))  

        index = np.argmin(np.abs(np.array(figure.output["time"]) - current_loc))
        axes.plot(np.transpose(figure.output["x"][:]),
                  np.squeeze(figure.output[figure.var][index, self.__graph_parameters.map1D_yloc, :]),
                  color='C1',
                  label=figure.output.getncattr("Build-Revision"))

        if self.__graph_parameters.observations:
            obs_index = np.argmin(np.abs(np.array(figure.observations["time"]) - current_loc))
            axes.plot(figure.observations["x"][:],
                      np.squeeze(figure.observations[figure.var][obs_index, :, self.__graph_parameters.map1D_yloc]),
                      color='k', linestyle="--",
                      label=self.__graph_parameters.datalabel)
        # Set labels
        axes.set_title(f"t = {float(figure.output['time'][index]):.2f} s")

    def set_limits(self, figure: FigureConfig, axes: plt.Axes):
        """Set x and y limit on axes for plot"""
        if self.__graph_parameters.xlim:
            xlim = self.__graph_parameters.xlim
        else:
            time = figure.output["time"][:]
            xlim = [np.min(time) / self.factor, np.max(time / self.factor)]

        if self.__graph_parameters.ylim:
            ylim = self.__graph_parameters.ylim
        else:
            plot_netcdf_var = figure.output[figure.var][:, :]
            ylim = [np.min(plot_netcdf_var), np.max(plot_netcdf_var)]

        plt.setp(axes, xlim=xlim, ylim=ylim)

    def create_table(self, figure):
        # TODO: create table for his, 1D and 2D!
        pass

    def prepare_xarray(self, figure: FigureConfig, file: str):
        """Open 2D datasets as xarray"""
        reference_path = os.path.join(self.__test_settings.reference_path, self.__testcase_config.path)
        reference_file = os.path.join(reference_path, file)
        figure.reference = xarray.open_dataset(reference_file, decode_times=False)

        case_path = os.path.join(self.__test_settings.result_path, self.__testcase_config.path)
        outputfile = os.path.join(case_path, file)
        figure.output = xarray.open_dataset(outputfile, decode_times=False)

    def set_2D_limits(self, figure: FigureConfig, axes: plt.Axes):
        """Set x and y limit on axes for 2D type plot"""
        if self.__graph_parameters.xlim_2D:
            xlim = self.__graph_parameters.xlim_2D
        else:
            x = figure.output["x"]
            xlim = [np.nanmin(x) / self.factor, np.nanmax(x) / self.factor]

        if self.__graph_parameters.ylim_2D:
            ylim = self.__graph_parameters.ylim_2D
        else:
            y = figure.output["y"]
            ylim = [np.nanmin(y) / self.factor, np.nanmax(y) / self.factor]

        plt.setp(axes, xlim=xlim, ylim=ylim)

    def plot_2D_map(self, figure: FigureConfig):
        '''Plot a 2d map'''

        figure.rows = rows = 1
        figure.columns = columns = 2
        figure.var = self.__graph_parameters.map2D_var
        self.factor = 1000

        self.prepare_xarray(figure, "sfincs_map.nc")

        size = self.set_size(width="a4", fraction=0.9, subplots=(rows, columns))
        fig = plt.figure(figsize=size)

        for i in range(figure.columns):
            axes = fig.add_subplot(rows, columns, i + 1)

            self.set_2D_limits(figure, axes)

            if i == 0:
                dataset = figure.reference
                title = figure.reference.attrs["Build-Revision"].replace('$', '').strip()
            else:
                dataset = figure.output
                title = figure.output.attrs["Build-Revision"].replace('$', '').strip()
            pcm = self.plot_2D_dataset(dataset, axes, figure.var)

            # Set labels
            axes.set_aspect('equal')
            axes.set_title(title)
            axes.set_xlabel(self.__graph_parameters.xlabel_2D)
            # Set labels only for left subplots
            if i == 0:
                axes.set_ylabel(self.__graph_parameters.ylabel_2D)

        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.04, 0.7])
        cbar = fig.colorbar(pcm, cax=cbar_ax)
        cbar.set_label(f"{figure.output[figure.var].long_name} in {figure.output[figure.var].units}", rotation=270)
        plt.suptitle(self.__testcase_config.name)

        figure_path = os.path.join(self.figure_folder,
                                   self.__testcase_config.name + '_2d.png')
        plt.savefig(figure_path, dpi=300)

    def plot_2D_dataset(self, dataset: nc.Dataset, axes: plt.Axes, variable: str):
        """Plot dataset for a 2d map"""
        var = dataset[variable]
        if var.ndim == 3:
            var = np.nanmax(var, axis=0)
        x = dataset["x"]
        y = dataset["y"]
        pcm = axes.pcolormesh(x / self.factor, y / self.factor, var,
                              vmin=self.__graph_parameters.clim_2D[0],
                              vmax=self.__graph_parameters.clim_2D[1])
        return pcm

    def set_size(self, width, fraction=1, subplots=(1, 1)):
        """Set figure dimensions to avoid scaling in LaTeX.

        Parameters
        ----------
        width: float or string
                Document width in points, or string of predined document type
        fraction: float, optional
                Fraction of the width which you wish the figure to occupy
        subplots: array-like, optional
                The number of rows and columns of subplots.
        Returns
        -------
        fig_dim: tuple
                Dimensions of figure in inches
        """
        if width == 'a4':
            width_pt = 594
        else:
            width_pt = width

        # Width of figure (in pts)
        fig_width_pt = width_pt * fraction
        # Convert from pt to inches
        inches_per_pt = 1 / 72.27

        # Golden ratio to set aesthetic figure height
        # https://disq.us/p/2940ij3
        golden_ratio = (5**.5 - 1) / 2

        # Figure width in inches
        fig_width_in = fig_width_pt * inches_per_pt
        # Figure height in inches
        fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])
        if subplots[0] == 1:
            fig_height_in += fig_height_in

        return (fig_width_in, fig_height_in)

    def rmse(self, predictions, targets):
        return np.sqrt(np.mean((predictions-targets)**2))

    def prepare_directory(self, folder: str):
        """Delete directory and contents and create a new empty directory."""
        if os.path.exists(folder):
            shutil.rmtree(folder, ignore_errors=True)
        os.makedirs(folder)
