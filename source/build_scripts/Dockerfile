FROM ubuntu:jammy
ENV DEBIAN_FRONTEND=noninteractive
RUN apt clean && apt autoclean && apt update --fix-missing && apt upgrade -y && apt install -y libnetcdf-dev build-essential autoconf automake libtool gfortran gdb tzdata m4
COPY . /usr/src/sfincs
WORKDIR /usr/src/sfincs
# Install software into /usr/local by default
# No OpenACC for GNU/CPU build
# RUN autoreconf -ivf && ./autogen.sh && FFLAGS="-w -fallow-argument-mismatch -O2" && FCFLAGS="-w -fallow-argument-mismatch -O2 -fopenmp -O3" ./configure --disable-openacc && make && make install
RUN chmod -R 777 autogen.sh
# -fallow-argument-mismatch needed for https://github.com/Unidata/netcdf-fortran/issues/212
ENV FCFLAGS="-fopenmp -O3 -fallow-argument-mismatch -w"
ENV FFLAGS="-fopenmp -O3 -fallow-argument-mismatch -w"
RUN autoreconf -ivf && ./autogen.sh && ./configure --disable-openacc && make && make install

# Make a new container and copy /usr/local so we only have the source
FROM ubuntu:jammy
ENV DEBIAN_FRONTEND=noninteractive
# install netcdf
RUN apt update && apt install -y libnetcdf-dev build-essential
COPY --from=0 /usr/local /usr/local
VOLUME /data
WORKDIR /data
CMD ["sfincs"]
