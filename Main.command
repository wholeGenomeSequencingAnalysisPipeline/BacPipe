#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd "${DIR}"

if [ $1 == "install" ];then
	#if [ $2 == "" ];then
	if [ -z "$2" ];then
		echo Please add a PATH for prokka, e.g. ./BacPipe.v1.run install /PATH/TO/Prokka
	else
		mkdir -p $2
		echo Please insert \sudo password
		sudo ./install.bash unix $2
	fi
else
	if [ $1 == "run" ];then
		python ./Pipeline.py unix
	else
		echo To install type: ./Main.command install or ./BacPipe.v?.?_unix.run install
		echo To run type: ./Main.command run or ./BacPipe.v?.?_unix.run run
	fi
fi
