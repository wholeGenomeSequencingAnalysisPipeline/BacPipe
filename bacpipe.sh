#!/bin/bash 
#############################################################################################################
xhost +
echo "Paste path to input folder" && read input_folder
echo 
echo 
echo docker run -it --volume $input_folder/:/mahmed/mnt/mydata --env="DISPLAY" --env="QT_X11_NO_MITSHM=1" --volume="/tmp/.X11-unix:/tmp/.X11-unix:rw" mahmed/bacpipe python ./Pipeline.py unix

docker run -it --volume $input_folder/:/mahmed/mnt/mydata --env="DISPLAY" --env="QT_X11_NO_MITSHM=1" --volume="/tmp/.X11-unix:/tmp/.X11-unix:rw" mahmed/bacpipe python ./Pipeline.py unix
