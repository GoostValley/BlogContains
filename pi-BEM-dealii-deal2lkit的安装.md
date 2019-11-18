---
title: 'pi_BEM,dealii,deal2lkit的安装'
date: 2019-11-02 19:00:04
tags: 
 - BEM 
 - dealii
categories:
 - BEM
 - dealii
---

[pi-BEM](https://github.com/mathLab/pi-BEM) 是一个高阶边界元求解函数库，集成了快速多极子算法、网格细化、并行求解等功能，实现优雅。其基于dealii以及deal2lkit进行开发，并且在安装dealii时还需要安装一系列第三方库，安装过程较为繁琐。

## Ubuntu 18.04安装

安装的过程分为三个步骤：

* 安装dealii及相应的第三方库
* 安装deal2lkit
* 安装pi-BEM

### dealii及第三方库的安装

dealii以及相应的第三方库的安装是最难的。在这里我卡壳了好几天。在项目的github主页中提供了三种安装方式：

1. 通过candi脚本安装

2. 通过spark安装

3. 手动从头安装

经过几天尝试，第二种和第三种我还是没有成功，通过candi脚本安装一开始我遇到了报错，以为失败了，可能是由于网络原因第三方库下载失败，后来再次尝试运行脚本，顺利地完成了安装。



[candi](https://github.com/koecher/candi)安装极为简单，只需要按照github主页上所说，下载文件，进入candi目录，然后运行脚本即可。

```bash
 git clone https://github.com/dealii/candi
 cd candi
 ./candi.sh -j8
```



安装前首先要确保自己已经安装好相应的环境，我的环境中已经安装好了Intel Parrallel Xe，在开始安装前，candi会提示安装依赖的库，如果不确定是否已经安装，可以运行一下，对于ubuntu18.04，提示的安装命令为：

```bash
sudo apt-get install build-essential lsb-release wget \
   automake autoconf gfortran \
   openmpi-bin openmpi-common libopenmpi-dev cmake subversion git \
   libblas-dev liblapack-dev libblas3 liblapack3 \
   libsuitesparse-dev libtool libboost-all-dev \
   splint tcl tcl-dev environment-modules qt4-dev-tools
```



pi-BEM要求额外安装一些第三方库，如petsc、Opencascad等，可以通过修改candi.cfg文件来更改安装的第三方库，注释掉不需要安装的库即可。如果在安装完之后发现哪个库忘了安装，取消注释并重新运行candi.sh即可。

```bash
# These packages determine the active components of deal.II:
#PACKAGES="${PACKAGES} once:adolc"
#PACKAGES="${PACKAGES} once:arpack-ng"
#PACKAGES="${PACKAGES} once:assimp"
#PACKAGES="${PACKAGES} once:nanoflann"
PACKAGES="${PACKAGES} once:opencascade"
PACKAGES="${PACKAGES} once:parmetis"
#PACKAGES="${PACKAGES} once:sundials"
#PACKAGES="${PACKAGES} once:superlu_dist"
PACKAGES="${PACKAGES} once:hdf5"
PACKAGES="${PACKAGES} once:p4est"
PACKAGES="${PACKAGES} once:trilinos"
PACKAGES="${PACKAGES} once:petsc"
#PACKAGES="${PACKAGES} once:slepc"
PACKAGES="${PACKAGES} dealii"
```



安装的过程旷日持久，所以为了加快编译速度，务必记得加上-j n。安装极为缓慢，反反复复安装了好几个小时……



最后dealii便安装完成啦，dealii以及相应的第三方库安装在/home/host_name/deal.ii-candi/中。最后根据命令行中的提示设置环境变量：

```bash
source /home/host_name/deal.ii-candi/configuration/deal.II-v9.1.1
source /home/host_name/deal.ii-candi/configuration/enable.sh
```

## deal2lkit安装

从[deal2lkit](https://github.com/mathLab/deal2lkit)主页上下载项目，并进行编译，由于之前已经在环境变量中设置好了dealii的安装位置，所以此处编译不需要声明dealii的安装路径：

```bash
git clone https://github.com/mathlab/deal2lkit.git
mkdir build
cd build
cmake ..
make
```

然后在环境变量中添加deal2lkit：

```bash
export D2K_DIR=/home/ghost/ykzhengWorkspace/learn_pi_bem/deal2lkit9.1.1/
export DEAL2LKIT_DIR=/home/ghost/ykzhengWorkspace/learn_pi_bem/deal2lkit9.1.1/
```



## pi-BEM 安装

dealii以及deallkit安装成功之后，pi-BEM安装便是小菜一碟：

```bash
git clone https://github.com/mathLab/pi-BEM.git
cd pi-BEM
mkdir build
cd build
cmake ../
make -j4
```

在安装完成后可以尝试着运行一下

```
cd build
mpirun -np 6 ./bem_fma_2d
mpirun -np 6 ./bem_fma_3d
```

至此pi-BEM便安装成功了。



## Docker 安装

在进行Docker安装之前，先了解一下Docker的一些基本命令：

```bash
docker search foo             查找镜像foo 
docker images                 列出本地已有的镜像 
docker ps -a                  列出本地所有的容器
docker run -it foo /bin/bash  运行镜像foo
docker start foo              启动容器foo
docker attach foo             链接到容器foo
```

在pi-BEM主页上的dealii Docker镜像并不能使用，这是因为目前pi-BEM已经基于dealii9.1，而提供的dealii版本为9.0。因此，在项目的Docker文件夹下，我发现作者在deal2lkit Docker镜像的基础上写了一个Docker file，镜像中没有编辑器很不方便，我在镜像中添加了vim编辑器。

```
FROM mathlab/deal2lkit:v9.1.1-debugrelease

MAINTAINER luca.heltai@gmail.com

# pi-BEM master image
RUN git clone https://github.com/mathLab/pi-BEM/ &&\
    mkdir pi-BEM/build && cd pi-BEM/build &&\
    cmake -DCMAKE_BUILD_TYPE=DebugRelease \
	  -GNinja \
          ../ && \
    ninja -j4 

# Need this for Travis!
USER root
RUN apt-get update && apt-get -y install vim
```

创建Dockerfile后，编译Docker镜像，编译pi-BEM镜像。

```bash
docker build -t pi_bem .
```

然后运行该镜像

```bash
docker run -it pi-BEM /bin/bash
```

进入镜像后，由于mpich没有添加到环境变量，无法执行mpirun命令，因此，修改.bashrc文件，讲mpich添加到系统路径。

```bash
export PATH=/usr/local/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/mpich-3.3-j5u4l3i4w5xjawupwn4gsrb43tg6wntz/bin:$PATH
export MANPATH=/usr/local/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/mpich-3.3-j5u4l3i4w5xjawupwn4gsrb43tg6wntz/man:$MANPATH
```

之后要运行pi-BEM只需启动容器即可。

```bash
docker ps -a #获取pi-BEM容器的编号 foo
docker start foo
docker attach foo
```

