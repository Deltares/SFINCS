cd /mnt/c/users/leijnse/repos/SFINCS/source

find . -name \*.m4|xargs dos2unix && find . -name \*.ac|xargs dos2unix && find . -name \*.am|xargs dos2unix
find . -name \*.f90|xargs dos2unix
find . -name \*.F90|xargs dos2unix
find . -name \*.am|xargs dos2unix
find . -name \*.sh|xargs dos2unix

docker build -f Dockerfile_cpu_local . -t leynse/sfincs-cpu
docker push leynse/sfincs-cpu
