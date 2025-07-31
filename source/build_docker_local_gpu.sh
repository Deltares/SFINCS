cd /mnt/c/work/checkouts/git/sfincs/source

find . -name \*.m4|xargs dos2unix
find . -name \*.ac|xargs dos2unix
find . -name \*.am|xargs dos2unix
find . -name \*.f90|xargs dos2unix
find . -name \*.F90|xargs dos2unix
find . -name \*.sh|xargs dos2unix

docker build -f Dockerfile.gpu.25.5.ccall . -t mvanormondt/sfincs-gpu:coldeze_combo_ccall_255 > build.log 2>&1
docker push mvanormondt/sfincs-gpu:coldeze_combo_ccall_255

#docker build -f Dockerfile.gpu.cc75 . -t mvanormondt/sfincs-gpu:coldeze_combo_cc75
#docker push mvanormondt/sfincs-gpu:coldeze_combo_cc75

#docker build -f Dockerfile.gpu.cc86 . -t mvanormondt/sfincs-gpu:coldeze_combo_cc86
#docker push mvanormondt/sfincs-gpu:coldeze_combo_cc86

#docker build -f Dockerfile.gpu.cc90 . -t mvanormondt/sfincs-gpu:coldeze_combo_cc90
#docker push mvanormondt/sfincs-gpu:coldeze_combo_cc90

