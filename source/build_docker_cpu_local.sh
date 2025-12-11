cd /mnt/c/users/leijnse/repos/SFINCS/source

docker build -f Dockerfile_cpu_local . -t leynse/sfincs-cpu
docker push leynse/sfincs-cpu
