---
title: deal.II学习--Step-3
date: 2019-12-15 16:23:24
tags:
- dealii
categories:
- dealii
mathjax: true
type: "picture"
layout: post
---

# 简介

## 有限元法简介

之前的例子涉及网格与稀疏模式，在这个例子中完整地运用有限元法进行计算。Step-3中求解右端项非0的泊松方程，且其边界条件为0：

\begin{aligned}
-\Delta u & = f \quad\quad & \text{in} \quad\Omega,\\\\
u & = 0                   & \text{on} \quad\partial\Omega.
\end{aligned}

计算域为单位正方形域：$\Omega=[0,1]^2$,在Step-1中已经学习了相应网格的生成。在这个示例程序中，考虑较为简单的情况：$f(\mathbf{x})=1$,在Step-4中会考虑更加复杂的情形。



在有限元法中，对解$u$进行近似，首先需要对方程进行处理，在方程的两边左乘测试函数$\varphi$，并在整个域上积分：

$$-\int_{\Omega}\varphi\Delta u=\int_{\Omega}\varphi f.$$

由格林第一公式可将上式化为：

$$\int_{\Omega}\nabla\varphi \cdot \nabla u-\int_{\partial\Omega}\varphi \mathrm{n} \cdot \nabla u=\int_{\Omega}\varphi f.$$



> 在此记录一下格林第一公式的证明:
>
> 假设$u(x,y,z)$，$v(x,y,z) \in C^2(\bar\Omega)$,设$F=u\nabla v$根据高斯公式有
>
> $$\iiint_{\Omega}\nabla F=\iiint_{\Omega} \nabla \cdot (u\nabla v)=\iint_{\partial\Omega }u \nabla v \cdot \vec n.$$
>
> 积分中的项$\nabla \cdot (u\nabla v)=\nabla u \nabla v+u\Delta v$，代入上式中，即可得到格林第一公式：
>
> $$\iiint_{\Omega}u \Delta v=\iint_{\partial \Omega}u \frac{\partial v}{\partial n} -\iiint_{\Omega} \nabla u\cdot \nabla v.$$



测试函数$\varphi$必须满足同样的边界条件，在此为$\varphi=0$，从而方程可以转化为

$$(\nabla \varphi,\nabla u)=(\varphi,f),$$

其中记号$(a,b)=\int_{\Omega}ab$。那么问题转化为求一个函数$u$，使得对于任意的测试函数，都要满足以上的方程。



在大部分情况下，通过计算机找到这样的函数是无法实现的，但是我们可以作近似$u_h(\mathbf{x})=\sum_j U_j \varphi_j(\mathbf{x})$,其中$U_j$是展开系数，为未知量，即Step-2中所说的自由度DoF；$\varphi_i(\mathbf{x})$为有限元法中的型函数。为了定义型函数，需要做以下的步骤：

* 必须在网格上定义型函数。在Step-1以及Step-2中已经讨论了相应的操作。
* 在有限元法中需要定义参考坐标系。在deal.II中对于一维，定义在$[0,1]$，对于二维，定义在$[0,1]^2$，对于三维，定义在$[0,1]^3$，即均为单位长度的线、面、体。在Step-2中，我们已经定义过FE_Q<2>，表示二维拉格朗日单元。其中，最为简单的便是线性单元FE_Q<2>(1)，表示多项式的阶数为1，即线性的。
* 类DoFHandler枚举了网格上所有的自由度，并将其与相应的单元进行对应。
* 指定在单元上要采用怎样的型函数，在deal.II中默认采用线性函数。



通过这些步骤，便可以得到型函数$\varphi_i$，并且型函数代入，便可以得到方程的离散形式：

$$(\nabla \varphi_i,\nabla u_h)=(\varphi_i,f), \quad \quad i=0\dots N-1.$$



在C或者C++中，下标都是从0开始的。这个方程可以写成一个线性方程组



\begin{aligned}
(\nabla \varphi_i,\nabla u_h) & = (\nabla \varphi_i,\nabla[\sum_j U_j \varphi_j])\\\\ 
& = \sum_j(\nabla \varphi_i,\nabla[U_j \varphi_j])\\\\
& = \sum_j(\nabla \varphi_i,\nabla \varphi_j) U_j.
\end{aligned}



线性方程组可以写成：

$$AU=F$$

其中矩阵$A$以及右端向量$F$定义为

$$A_{ij}=(\nabla\varphi_i,\nabla\varphi_j),$$

$$F_i=(\varphi_i,f).$$

## 测试函数的左乘与右乘

假如将在方程两边右乘测试函数，那么最终得到的线性方程组为

$$U^TA=F^T$$

其中$F^T$为列向量。对上述矩阵方程进行转置：

$$A^TU=F$$

其中，$A=A^T$是对称矩阵。不过在大多数情况下$A$并不是对称矩阵，在此只是为了说明左乘右乘的问题对此进行简化。大多数情况下在方程左乘测试函数只是习惯，无论左乘右乘最终的结果都是正确的。



## 计算矩阵以及右端向量

现在我们已经知道了线性方程组是如何得到的，那么接下来要探究其中的细节。

* 线性方程组中的矩阵$A$为稀疏矩阵，$U$和$F$为向量，在示例程序中会展示如何求解线性方程组。

* 我们需要求解积分。在此采用数值积分的方式，即，最终积分的值是单元中各个点的函数值的加权和。于是，将计算域$\Omega$分解为各个单元：

  $$A_{ij}=(\nabla \varphi_i,\nabla \varphi_j)=\sum_{K \in T} \int_K \nabla \varphi_i \cdot \nabla \varphi_j ,$$

  $$F_i=(\varphi_i,f)=\sum_{K \in T}\int_K \varphi_i f.$$

  在每个单元中数值积分的求取为：

  $$A_{ij}^K=\int_K \nabla \varphi_i \cdot \nabla \varphi_j \approx\sum_q \nabla \varphi_i(\mathbf{x}_q^K) \cdot \nabla \varphi_j (\mathbf{x}_q^K) \omega_q^k,$$

  $$F_i^K=\int_K \varphi_i f \approx \sum_q \varphi_i (\mathbf{x}_q^K) f(\mathbf{x}_q^K) \omega_q^K.$$

  ​



# 示例程序

