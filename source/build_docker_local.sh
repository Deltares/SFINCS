cd /mnt/c/work/checkouts/git/sfincs/source

find . -name \*.m4|xargs dos2unix && find . -name \*.ac|xargs dos2unix && find . -name \*.am|xargs dos2unix
find . -name \*.f90|xargs dos2unix
find . -name \*.F90|xargs dos2unix
find . -name \*.am|xargs dos2unix
find . -name \*.sh|xargs dos2unix

docker build -f Dockerfile.cpu . -t mvanormondt/sfincs_cpu
docker push mvanormondt/sfincs_cpu:latest
