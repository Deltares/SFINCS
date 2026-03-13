cd /mnt/c/work/checkouts/git/sfincs/source

find . -name \*.m4|xargs dos2unix && find . -name \*.ac|xargs dos2unix && find . -name \*.am|xargs dos2unix
find . -name \*.f90|xargs dos2unix
find . -name \*.F90|xargs dos2unix
find . -name \*.am|xargs dos2unix
find . -name \*.sh|xargs dos2unix

docker build -f Dockerfile.gpu.25.5.ccall . -t leynse/sfincs-gpu:latest_11032026 > build.log 2>&1