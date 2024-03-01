import os
import shutil
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import xarray
import pandas as pd
from scipy import interpolate
#from typing import List
from src.config.test_case_config import TestCaseConfig
from src.config.figure_config import FigureConfig

# create an object for latex and plot configuration.
# overwrite standard with parsed xml values
# try to store as little data in xml as possible.
# connect the latex data to the testcase
# rewrite resultplotter to use the object for it's configuration


class ResultPlotter:
    """Plot a comparison graph into a latex file for specified data."""

    def __init__(self, testcase_config: TestCaseConfig):
        self.__testcase_config = testcase_config
        self.__test_settings = testcase_config.test_settings
        self.__graph_parameters = testcase_config.graph_parameters
        self.factor = 3600

    def sfincs_postprocess(self):
        """Postprocess data from run with reference and/or observations."""
        if self.__graph_parameters.his:
            figure = self.prepare_dataset("sfincs_his.nc")
            self.prepare_plot(figure, "his")
            self.plot_dataset(figure)

        if self.__graph_parameters.map != "":
            figure = self.prepare_dataset("sfincs_map.nc")
            self.prepare_plot(figure, self.__graph_parameters.map)

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
            current_loc = figure.locations[i] - 1

            axes = figure.figure_plot.add_subplot(figure.rows, figure.columns, i + 1)
            axes.grid()
            self.set_limits(figure, axes)

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

        figure_folder = os.path.join(self.__test_settings.figure_path,
                                     self.__testcase_config.path)
        self.prepare_directory(figure_folder)
        figure_path = os.path.join(figure_folder,
                                   self.__testcase_config.name + '.png')

        plt.savefig(figure_path, dpi=300)

    def set_limits(self, figure: FigureConfig, axes: plt.Axes):
        """Set x and y limit on axes for plot"""
        # TODO: find out if necessary specific xlim and ylim for 1d plots
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

    def create_table(self):
        # Create table with version, total_runtime and RMSE between stable release and latest
        if self.config["stable_version"]:
            if self.config["revision_nr"]:
                initial_data = {"Release version":["Rev : " + str(self.config["revision_nr"]),output.getncattr("Build-Revision")], "Total runtime [s]": [float(stable["total_runtime"][:]),float(output["total_runtime"][:])]}
            else:
                initial_data = {"Release version":[stable.getncattr("Build-Revision"),output.getncattr("Build-Revision")], "Total runtime [s]": [float(stable["total_runtime"][:]),float(output["total_runtime"][:])]}
            df = pd.DataFrame(initial_data, columns=['Release version', 'Total runtime [s]'])

            for i in range(0,Tot):
                if self.config["his_loc"]:
                    #python starts numbering at 0, therefore -1
                    iloc = self.config["his_loc"][i] -1
                else:
                    iloc= i

                #interpolate stable release onte latest timeseries ...
                f = interpolate.interp1d(stable["time"][:],stable[var][:,iloc],fill_value="extrapolate")
                stable_var_new = f(np.asarray(output["time"][:]))

                RMSE = self.rmse(output[var][:,iloc],stable_var_new)
                df_new=pd.DataFrame({'RMSE loc:' + str(nc.chartostring(output["station_name"][iloc,:])).strip() : [np.nan, RMSE]})
                df = pd.concat([df, df_new], axis=1)

            texname = os.path.join(self.testbedloc,"models",model.testid, model.runid, model.runid + '.tex')
            if Tot > 3:
                df.set_index('Release version',inplace=True)
                df = df.transpose()
                df.to_latex(buf=texname, index=True, float_format="%.3f",na_rep="-")
            else:
                df.to_latex(buf=texname, index=False, float_format="%.3f",na_rep="-")

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
