# Running SFINCS

You can run sfincs using the `sfincs` command. It assumes the input file are in the current working directory. This assumes the sfincs executable is in your path and the corresponding libraries are in the library path (the directory containing the libraries should be in the `LD_LIBRARY_PATH` variable in linux, in `DYLD_LIBRARY_PATH` under OSX).

``` shell
    $ sfincs
```

You can also run sfincs using a docker container. For docker we have two versions, the cpu version (compiled with gfortran) and the gpu version (compiled with the nvidia hpc sdk).

``` shell
    $ docker pull deltares/sfincs-cpu
    $ docker run -v$(pwd):/data deltares/sfincs-cpu
```

If you are on a cluster that supports singularity, you can start sfincs using the singularity run command.

``` shell
    $ singularity run -B$(pwd):/data --nv docker://deltares/sfincs-gpu
```

# Building the linux version

Make sure you have the following tools installed:
- compiler (gnu fortran or intel for CPU only, nvidia HPC SDK for GPU version)
- autotools (autoconf, automake, libtool, m4, make, also collectively available as build-essentials in many linux distributions)

This is a short summary of the installation and building procedure. See the general INSTALL file for detailed instructions.

To prepare the source code directory for building you need to run autoreconf once like this:

```
    $ autoreconf -vif
```

At this point you should be able to build the software. Please consult `./configure --help` for extra options, such as compiler selection (`FC=gfortran`) and the install location (`--prefix`). If you want to use openmp and openacc the options are checked against the corresponding C compiler. So if you change the compiler make sure you also spcify the corresponding C compiler (for example CC=gcc FC=gfortran or CC=nvc FC=nvfortran). You might need to pass `--with-pic` or `--disable-shared` on certain platforms.

```
    $ ./configure
    $ make
```

This should give you an sfincs executable in the src directory and a libsfincs.la that you can link to.
To install it into your system you can use make install.

```
    $ make install
```

See the [autobook](https://www.sourceware.org/autobook/autobook/autobook_toc.html) for documentation on the build system.
