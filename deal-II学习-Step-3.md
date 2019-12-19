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

  其中$\mathbf{x}_q^K$是单元$K$中的第$q$个积分点，$\omega_q^K$是积分权重。通常参考单元与形函数之间的对应方式是一直地，默认为MappingQ1，不过也可以显式地指定其它的形式。在单元中的积分点以及积分权重的设定在Quadrature类中实现。通常数值积分的值可以等于理论积分的值，因为所有需要积分的项均为多项式，在此一般用高斯积分，在deal.II中由QGauss类实现。

* 随后需要定义一个能在单元$k$中计算$\varphi_i(\mathbf{x}_q^K)$的类，在此由类FEValues实现。FEValues需要输入的参数有：一个有限元实例，实例中描述了参考单元中的形函数$\varphi$；一个数值积分实例，实例中指明了积分点、积分权重；以及一个匹配实例（默认采用MappingQ1)，该实例将提供形函数的微分积分，数值积分点等积分所需要的信息。


FEValues是合成线性方程组的核心过程。可以这样看待这一过程：类FiniteElement以及衍生的类描述了形函数，即保证了每个点上都有相应的函数值。需要形函数是因为我们需要对函数进行积分。但是，计算机只能处理离散信息，所以需要使用数值积分。由Mapping类匹配实际单元与参考单元，由Quadrature类提供数值积分点，通过对函数值加权求和得到数值积分的值。更为简单地说，我们需要在有限个点上拥有形函数值以及导数、高斯积分权重以及法向量等，便可以求解这个问题。FEValues类便是整合这些信息的类。



这三个类的工作也可以通过程序实现，并整合相应的信息。但是deal.II中的FEValues提供了大量优化，实现简洁、快速。



最后采用线性求解器求解线性方程组，并通过DataOut类输出结果。



## 关于程序实现

虽然这是用有限元法实现的最简单的问题，但是在示例程序中展示了有限元法程序的基本结构，这也可以作为求解其它问题的一个基石模板，程序中类的结构如下：

```c++
class Step3
{
  public:
    Step3 ();
    void run ();
  private:
    void make_grid ();
    void setup_system ();
    void assemble_system ();
    void solve ();
    void output_results () const;
    Triangulation<2>     triangulation;
    FE_Q<2>              fe;
    DoFHandler<2>        dof_handler;
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double>       solution;
    Vector<double>       system_rhs;
};
```



在代码中努力实现封装的目的，即尽可能地隐藏类内部的细节，使其无法从类外部访问。



首先介绍类中的成员变量：需要网格Triangulation类以及自由度分配DoFHandler类的实例对象、以及一个有限元对象来描述形函数。接着还需要与线性代数有关的对象：包括矩阵、右端向量、解向量以及一个描述稀疏模式的矩阵。这便是所有在程序运行始终都需要的类了。相比之下，FEValues类只需要在组成线性方程组的时候用到，便只在相关的函数中声明，在使用结束后便销毁。



接着，介绍一下成员函数。这些函数组成了一个基本的结构，在之后的教程中都会用到：

* make_grid(): 在这个函数中保存了网格Triangulation。在之后的教程中，还会处理边界条件、几何等相关信息。
* setup_system(): 这个函数处理求解这个问题有关的线性方程组相关的所有数据结构。首先，它要初始化DoFHandler类，并且处理与线性方程组有关的各个类。将这个函数与处理网格的函数分开，是因为在时域问题中，网格每隔几步才更新一次；而建立线性方程组的数据结构在程序一开始就完成，因此将其分离开来。
* assemble_system(): 计算线性方程组中矩阵与右端向量的系数。
* solve(): 计算线性方程组的解，由于目前的矩阵较为简单，所以solve()函数的实现也较为简单。但是随着问题逐渐复杂，求解器的设计也会相应变得复杂（在之后的Step-20、Step-22以及Step-31会涉及到）。
* output_results(): 最后，输出计算的结果。比如，需要将计算的结果以一定的格式组织，从而输入一些可视化软件作图。或者只输出关心的物理量，忽略其他数据：比如热交换问题只关心热量流动、机翼问题关心空气摩擦系数、桥梁的最大载荷或者最简单的——某一个点处解的值。所以在程序的最后加上这个模块来方便后处理。



上述的所有模块组合成为一个函数——run()。在run()函数中，它由其所属的类调用，并调用所有相关的函数来求解问题。将函数封装进run()函数，相比于在main()函数中逐个调用有关模块更好，这是因为可以在类中做出改动，而不影响类之外的代码。



## 变量类型提醒

deal.II在命名空间内定义了一系列数据类型。特别地，在这个程序中可以看到type::global_dof_index的变量类型，这是一个整数类型，用来表示自由度的数目，即DoFHandler在整个网格上得到的未知量个数。对于目前的程序，可能得到DoF的总数为几千到几百万个。因此，需要一个足够大的数来存储DoF的数目。

unsigned int 的范围为0到40亿，在deal.II 7.3版本之前，type::global_dof_index都是这个类型。不过，现在deal.II已经支持非常大规模的计算，超过了40亿的范围。因此对于type::global_dof_index做了设定，默认类型为unsigned int，如果超过表示范围，那么类型更改为unsigned long long int。



知道变量的类型有助于我们了解变量的数量级，此外，也可以避免可能出现的错误：比如将types::global_dof_index赋值给types::subdomain_id，因为两者都是unsigned int类型，所以编译器不会报错，但是可能出现错误。



在实际的代码编写过程中，这个类型在组成总矩阵要用到，比如创建一个二维Q1单元，那么在局部会有一个$4\times 4$的矩阵，随后需要将这个矩阵添加到总矩阵中。那么便需要知道局部坐标下的点在总体坐标下的下标是多少，实现的方式如下：

```c++
cell->get_dof_indices (local_dof_indices);
```

其中，下标local_dof_indices声明方式如下:

```c++
std::vector<types::global_dof_index> local_dof_indices (fe.dofs_per_cell);
```

  

# 示例程序

首先需要引入相关的库，有些库在之前已经有过介绍：

```c++
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
```

提供局部坐标下单元中DoF的信息：

```c++
#include <deal.II/dofs/dof_accessor.h>
```

在每个单元上进行数值积分：

```c++
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
```

处理边界条件的值：

```c++
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
```

接下来引入的库与求解线性方程组有关，需要在每个单元用到向量以及满矩阵来合成局部的矩阵，然后在合成全局的稀疏矩阵。在示例程序中用CG法求解线性方程组，所以需要对矩阵进行预处理。

```c++
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
```

最后是输出文件以及命名空间

```c++
using namespace dealii;
```

## Step-3中的类

与前面面向过程的程序不同，这个程序采用面向对象的风格编写，所有的细节都封装在类中。类中的函数各自执行相应的有限元法中的功能。在这里首先需要声明的是一个调用所有函数的类似于"main"函数的函数，然后需要定义各种成员变量。



类的公有部分非常简单：包括一个构造函数以及一个run()函数。run()函数在此的功能类似于"main"，即按顺序调用其它函数，而其它的函数都封装在私有部分。

```c++
class Step3
{
public:
  Step3();
  void run();
```



随后是一系列的私有函数，在之前已经讨论过其相应的功能。这些函数无需从外部调用，因此设为类的私有成员。

```c++
private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results() const;
```



最后是相应的成员变量，其中有关于网格的Triangulation，单元的类型以及自由度的处理：

```c++
Triangulation<2> triangulation;
FE_Q<2>          fe;
DoFHandler<2>    dof_handler;
```

总矩阵以及稀疏模式：

```c++
SparsityPattern      sparsity_pattern;
SparseMatrix<double> system_matrix;
```

解向量以及右端向量：

```c++
  Vector<double> solution;
  Vector<double> system_rhs;
```



### Step3::Step3

构造函数需要完成对象的初始化。在此具体只有两个内容：

* 设定单元的类型。在此设为线性单元。
* 将自由度的dof_handler与网格Triangulation关联起来，这里并没有涉及到具体的网格，只有当调用distribute_dofs()时，才会真正地根据网格生成自由度信息。

### Step3::make_grid

首先需要生成网格，在Step-1中已经介绍过。在这里，首先生成$[-1,1]\times [-1,1]$的正方形域，然后全局细化五次，那么，总共有$2^5\times 2^5=1024$个单元。



需要注意的是这里用的n_active_cells()，指的是最终不再需要细化的网格，如果用的是n_cell()，那么还会包括一些粗的网格。

```c++
void Step3::make_grid()
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(5);
  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
}
```



### Step3::setup_system

随后需要枚举所有的自由度，并建立矩阵以及向量对象来存储数据。与Step-2中相同，利用DoFHanler::distribute_dofs()，可以枚举出所有的自由度。



由于在此采用线性单元，那么每一个节点有一个自由度，一共有$33\times 33$个节点，所以自由度的数目为1089。

```c++
void Step3::setup_system()
{
  dof_handler.distribute_dofs(fe);
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;
```



与Step-2中相同，首先建立动态的稀疏模式，随后复制到SparsityPattern中。值得注意的值SparityPattern本身只存储非零项的位置，不存储矩阵的值，矩阵值存储再SparseMatrix中。

```c++
DynamicSparsityPattern dsp(dof_handler.n_dofs());
DoFTools::make_sparsity_pattern(dof_handler, dsp);
sparsity_pattern.copy_from(dsp);
```



涉及稀疏模式这个类的目的在于让不同的矩阵可以公用一个稀疏模式。在大规模问题中，建立稀疏模式是极为耗时的，因此这样的处理方式在时间空间上都提高了效率。

```c++
system_matrix.reinit(sparsity_pattern);
```

最后设定解向量以及右端向量的大小：

```c++
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
```



### Step3::assemble_system

在这个步骤中需要计算矩阵以及右端向量中的系数，这是有限元法的核心步骤，在简介中已经有所介绍。



循环遍历所有单元，通过数值积分计算系数，便可以得到矩阵及右端向量里的系数。问题的关键在于得到实际坐标系下单元中积分点的位置以及形函数的值。不过实际上提供的是参考坐标系下的结果。实际上，在程序中并未提供相应的接口，并没有办法直接查询某个单元形函数的值以及积分点的位置。解决这个问题的办法是定义一个Mapping类，用来将参考单元与实际的单元做匹配。



所以现在便讲解完了所需的三种类：单元、数值积分以及mapping。类FEValues将三个类的信息综合起来，如果给出了三个类的对象，那么FEValues便可以返回实际网格中在积分点上形函数的值与梯度。最终的函数如下：



```c++
void Step3::assemble_system()
{
```

首先需要指定在每个单元中进行积分的函数。在此采用高斯积分，对于二维情况，每个方向2个积分点，一共有四个积分点，对于一维情况来说，积分多项式的阶数可以达到3阶，对于目前的问题已经足够。

```c++
QGauss<2> quadrature_formula(fe.degree + 1);
```



随后初始化fe_values类，为每个单元计算积分点以及权重，随后，隐式地使用Q1 mapping的方式匹配参考单元与实际单元。最后，需要指出需要的结果，如形函数在积分点处的值、梯度以及高斯积分的权重等信息。



在程序的实现中出现了一系列的判断信号，用来指明需要更新的物理量，比如update_values指出需要更新形函数的值；update_gradients为需要更新梯度；update_JxW_values为需要更新高斯积分的权重以及Jacob矩阵的乘积。



```c++
FEValues<2> fe_values(fe,
                      quadrature_formula,
                      update_values | update_gradients | update_JxW_values);
```

这种实现方式的好处是可以指明在每个单元中到底需要什么样的输出信息。相比于计算所有信息，有选择地计算需要地信息可以大大减少计算量。在此update_values | update_gradients | update_JxW_values是一种典型地C语言中的写法，在计算机内数字用二进制表示，假设update_values=0b00001=1, update_gradients=0b00010=2, update_JxW_values=0b10000=16,在计算机内相应的位上为0表示关闭相应位上的功能，为1表示开启该功能，通过按位取或之后得到最终的结果为update_values | update_gradients | update_JxW_values = 0b10011 = 19，最终通过查询结果相关的位上是0还是1即可知道需要输出什么信息。



在此定义两个常量，一个是每个单元中的自由度个数，另一个是每个单元中的积分点数目。

```c++
const unsigned int dofs_per_cell = fe.dofs_per_cell;
const unsigned int n_q_points    = quadrature_formula.size();
```



随后初始化局部坐标系下的矩阵以及右端向量，为其设定空间大小。

```c++
FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
Vector<double>     cell_rhs(dofs_per_cell);
```



从局部坐标系到整体坐标系，需要有相应的对应关系，因此需要一个向量来存储局部坐标系下的点在整体坐标系下的下标，注意到这个向量的类型是types::global_dof_index，在之前已经讲过这个类型。

```c++
std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
```



现在循环所有单元，DoFHandler中的迭代器与之前Triangulation的迭代器使用方法类似。这里用的是const auto &，这是因为在此没有对变量做修改，因此可以常量引用；在Step-1中修改了cell是否细化的标记，所以不加const。

```c++
for (const auto &cell : dof_handler.active_cell_iterators())
  {
```



现在在一个单元中操作，需要计算每个单元的形函数、梯度以及Jacob矩阵、积分点等。由于每个单元中这些值都依赖于单元的几何条件，因此对每个单元都要独立计算。

```c++
fe_values.reinit(cell);
```

将局部的矩阵与右端向量置零：

```c++
cell_matrix = 0;
cell_rhs    = 0;
```

循环遍历单元中的所有点：

```c++
for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
  {
```

对于拉普拉斯问题，每个单元中的矩阵系数是形函数i和j的梯度的积分。在数值积分中，直接对值进行加权求和，因此只需要将点i点j的梯度相乘，并乘上权重与Jacob矩阵：

```c++
for (unsigned int i = 0; i < dofs_per_cell; ++i)
  for (unsigned int j = 0; j < dofs_per_cell; ++j)
    cell_matrix(i, j) +=
      (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
       fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
       fe_values.JxW(q_index));           // dx
```

随后对右端向量做类似的操作：

```c++
  for (unsigned int i = 0; i < dofs_per_cell; ++i)
    cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                    1 *                                 // f(x_q)
                    fe_values.JxW(q_index));            // dx
}
```

完成局部矩阵的生成，随后查找单元中的点在整体坐标系中的下标：

```c++
cell->get_dof_indices(local_dof_indices);
```

循环所有单元中的点，将局部矩阵以及右端向量的值转移到整体矩阵与右端向量中：

```c++
for (unsigned int i = 0; i < dofs_per_cell; ++i)
  for (unsigned int j = 0; j < dofs_per_cell; ++j)
    system_matrix.add(local_dof_indices[i],
                      local_dof_indices[j],
                      cell_matrix(i, j));
for (unsigned int i = 0; i < dofs_per_cell; ++i)
  system_rhs(local_dof_indices[i]) += cell_rhs(i);
```

现在需要添加边界值，事实上，如果没有Dirichlet边界条件，Laplace方程的解不是唯一的，因为将求得的解加上任意常数，仍然满足Laplace方程。



接下来需要获取边界上的自由度，一起相应位置的形函数。为了简化代码，这里在插值处边界值而非在边界处设置形函数。在deal.II中，VectorTools::interpolate_boundary_values()可以完成这个任务，所需的参数有：DoFHandler；需要插值处理的边界；边界条件的函数；以及输出对象。



在许多情况下，我们只想在某一段特定的边界设置边界条件的值。比如，流体力学中常常在进流段与出流段设定边界条件；或者在计算物体变形时固定端或者自由端。那么对于不同的边界便设定不同的数字，这样对于不用的边界便采用不同的边界条件函数。所有边界的类型默认设置为0，如果没有更改，那么就会保持默认值；因此如果设定类型为0的边界进行处理，那么便会对整段边界进行处理。如果要对不同的边界区别处理，便要对其命名不同的数字。



描述边界值的函数是类Function的对象或者Function的衍生类的对象。Function::ZeroFunction便是衍生类之一，设置边界值处处为0。自由度与边界值之间的匹配由类std::map完成。

```c++
std::map<types::global_dof_index, double> boundary_values;
VectorTools::interpolate_boundary_values(dof_handler,
                                         0,
                                         Functions::ZeroFunction<2>(),
                                         boundary_values);
```



设置好了边界值以及边界函数，那么便可以相应地修改总矩阵：

```c++
  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);
```



### Step3::solve

这个函数求解线性方程组。由于方程组规模较大，不适合用直接姐发如高斯消去法或者LU分解法求解，在此用CG算法求解。这个例子中未知量的个数为1089，对于有限元法来说，这是一个非常小的数字，一般来说，未知量的个数在百万量级，在这种量级上，直接法无法开销巨大，无法求解，因此只能选择如CG这样的迭代法。



首先需要确定求解器的终止条件，在此通过SolverControl类实现，这里参数设置的意思是，最多进行1000次迭代计算，或者在残差小于$10^{12}$时停止迭代：
```c++
SolverControl solver_control(1000, 1e-12);
```



随后便要设置求解器，SolverCG默认的数据类型是Vector，这里没有进行设置，则取默认设置：

```c++
SolverCG<> solver(solver_control);
```



最后求解方程，solver的第四个参数为预处理器，在此暂时不深究预处理的细节。

```c++
solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
```

最后边完成了求解的工作，变量solution中存储了未知量的解。



### Step3::output_results

有限元法的最后一个部分便是输出计算结果以及后处理。在此没有什么后处理，不过仍然需要将计算结果写到文件中。

```c++
void Step3::output_results() const
{
```

为了将结果输出到文件中，需要有一个类来定义输出的格式，在此便是DataOut类。

```c++
DataOut<2> data_out;
```

为了将需要的数据输出，首先需要链接DoFHandler，然后便是添加所需的数据，并为其设置标签，比如将solution向量输出，并添加标签"solution"：

```c++
data_out.attach_dof_handler(dof_handler);
data_out.add_data_vector(solution, "solution");
```

为了将数据转换成各种不同的格式，需要将数据转化为一种中间格式，然后根据要求再从中间格式变为相应的输出格式：

```c++
data_out.build_patches();
```

最后，将数据根据GNUPLOT格式输出：

```c++
  std::ofstream output("solution.gpl");
  data_out.write_gnuplot(output);
}
```

### Step3::run

实现了所有的模块之后，在run()函数中将其组装起来：

```c++
void Step3::run()
{
  make_grid();
  setup_system();
  assemble_system();
  solve();
  output_results();
}
```



### main函数

在main函数中没有什么具体的实现，一般都是定义一些类，并调用相应的函数。



首先，定义deallog来输出程序的运行信息，deallog的输出信息可以写在命令行、文件等等。默认情况下deallog不会进行输出，除非显式声明。这里参数为2的意思是输出的前缀数目，比如这里可以看到DEAL:CG,这里便有两个前缀。



```c++
int main()
{
  deallog.depth_console(2);
  Step3 laplace_problem;
  laplace_problem.run();
  return 0;
}
```



# 计算结果

计算结果如下所示：

```c++
Number of active cells: 1024
Number of degrees of freedom: 1089
DEAL:cg::Starting value 0.121094
DEAL:cg::Convergence step 48 value 5.33692e-13
```



前两行是调用cout输出的结果，输出单元以及后两行是deal.II中CG求解器自动输出的日志。

除了输出的日志以外，计算的结果还输出到了solution.gpl中，通过GNUPLOT可以绘制相应的图像。

```bash
examples/step-3> gnuplot
        G N U P L O T
        Version 3.7 patchlevel 3
        last modified Thu Dec 12 13:00:00 GMT 2002
        System: Linux 2.6.11.4-21.10-default
        Copyright(C) 1986 - 1993, 1998 - 2002
        Thomas Williams, Colin Kelley and many others
        Type `help` to access the on-line reference manual
        The gnuplot FAQ is available from
        http://www.gnuplot.info/gnuplot-faq.html
        Send comments and requests for help to <info-gnuplot@dartmouth.edu>
        Send bugs, suggestions and mods to <bug-gnuplot@dartmouth.edu>
Terminal type set to 'x11'
gnuplot> set style data lines
gnuplot> splot "solution.gpl"
```

上述命令会生成交互式界面，如果需要直接输出png图像，可以通过以下命令：

```bash
gnuplot> set output "1.png"
gnuplot> set terminal png truecolor
gnuplot> set style data lines
gnuplot> splot "solution.gpl"
```

上述生成的图像堆叠在了一起，可以通过命令不显示背面的图像：

```bash
gnuplot> set hidden3d
```



<center>
<div style="display:inline-block;width:350px;">{%asset_img  1.png%}</div>
<div style="display:inline-block;margin-left:10px;width:350px">{%asset_img  2.png%}</div>
</center>



## 可能的拓展

如果需要有进一步的变化，可以尝试以下的方式：

* 改变几何与网格：在程序中通过GridGenerator::hyper_cube生成了正方形计算域。除了正方形域之外，GridGenerator函数中还有众多其它的函数，比如L形域、环形或者其它的形状。

* 改变边界条件：在示例代码中采用了ZeroFunction，将边界条件都设为0。可以采用ConstantFunction<2>代替代码中的ZerosFunction<2>，从而将边界条件设为常数。除此之外在Functions命名空间中还有许多其它函数，可以从中挑选处一个特殊的边界值。

* 修改边界条件的类型：目前在示例程序中采用的是Dirichlet条件。在边界上默认的标识数字是0，通过VectorTools::interpolate_boundary_values函数可以将所有标识数字为0的边界值设为0。那么相应的可以在一些边界点上设置不同的标识数字。比如，在调用GridGenerator::hyper_cube之后加上：

  ```c++
  triangulation.begin_active()->face(0)->set_boundary_id(1);
  ```

  在最初始时，网格没有细化，所以网格为一个正方形，将正方形的一条边的边界条件标识设为1。在这之后，网格细化之后，新的网格会继承粗网格的边界条件标识符。紧接着，需要设置边界条件的值，interpolate_boundary_values函数会对标识为0的边界设置边界值0，并且对标识为1的边界不做处理。这样做之后将会对标识为0的边界设置Dirichelet边界条件，对标识为1的边界设置Neuman边界条件（切向梯度为0）。

* 可以对标识为1的边界设置边界值，再次调用interpolate_boundary_values函数进行设置：

  ```c++
  VectorTools::interpolate_boundary_values (dof_handler,
                                            1,
                                            ConstantFunction<2>(1.),
                                            boundary_values);
  ```

  调用这个函数会将标识为1的边界的值设为1。这样在边界上是不连续的，三面为0，一面为1，中间没有过渡。

* 观察程序的收敛性：在Step-7会讨论计算误差的范数，在这里仅仅做一个简单的收敛性检查。比如采用不同的全局优化次数，计算点$(\frac{1}{3},\frac{1}{3})$的函数值，可以将求函数值的代码添加到LaplaceProblem_output_results函数中：

  ```c++
  std::cout << "Solution at (1/3,1/3): "
            << VectorTools::point_value (dof_handler, solution,
                                         Point<2>(1./3, 1./3))
            << std::endl;
  ```

  得到的一系列值为：

| of refinements | $u_h(\frac{1}{3},\frac{1}{3})$ |
| :------------: | :----------------------------: |
|       1        |            0.166667            |
|       2        |            0.227381            |
|       3        |            0.237375            |
|       4        |            0.240435            |
|       5        |            0.241140            |
|       6        |            0.241324            |
|       7        |            0.241369            |
|       8        |            0.241380            |
|       9        |            0.241383            |


  观察两次计算间的差值，大致以4倍的倍数减小，可以猜想最终的精确值约为0.241384。可以猜测这个序列的收敛速率为$O(h^2)$，理论上，收敛阶数应该为$O(h^2|\mathrm{log}h|)$,由于计算域是对称的，实际上的收敛速度会更快。

  还可以修改单元里形函数多项式的阶数，探究不同阶数的单元的计算效果。

* 查看解的平均值的收敛性：

  ```c++
  std::cout << "Mean value: "
            << VectorTools::compute_mean_value (dof_handler,
                                                QGauss<2>(fe.degree + 1),
                                                solution,
                                                0)
            << std::endl;
  ```

  计算的结果如下：



| of refinements | $\int_{\Omega}u_h(x)dx$ |
| :------------: | :---------------------: |
|       0        |       0.09375000        |
|       1        |       0.12790179        |
|       2        |       0.13733440        |
|       3        |       0.13976069        |
|       4        |       0.14037251        |
|       5        |       0.14052586        |
|       6        |       0.14056422        |
|       7        |       0.14057382        |
|       8        |       0.14057622        |
