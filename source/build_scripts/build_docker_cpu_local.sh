cd /mnt/c/users/leijnse/repos/SFINCS/source

docker build -f ./build_scripts/Dockerfile_cpu_local . -t leijnse/sfincs-cpu > build.log 2>&1
docker push leijnse/sfincs-cpu
