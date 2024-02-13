cd /mnt/c/work/checkouts/git/sfincs/source
docker build -f Dockerfile.gpu . -t mvanormondt/sfincs-gpu
docker push mvanormondt/sfincs-gpu:latest
