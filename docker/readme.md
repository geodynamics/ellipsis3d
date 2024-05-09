This container hosts a built version of ellipsis3d.

docker run -it --rm -v $HOME/ellipsis3d:/home/ellipsis3d_user/work geodynamics/ellipsis3d

This command will start the ellipsis3d docker image and give you terminal access. Any changes made in the /home/ellipsis3d_user/work directory will be reflected on the host machine at home/ellipsis3d.

OpenDX is no longer maintained, so the docker container has not been verified against the test data.
