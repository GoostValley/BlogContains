---
title: deal.II学习--Step-2
date: 2019-12-11 10:32:05
tags:
- dealii
categories:
- dealii
---

# 简介

在 Step-1 中已经介绍了如何定义网格，在这个部分将会介绍如何在网格的基础上定义自由度。这个示例程序采用线性单元，即一阶单元（$Q_1$），单元中的自由度的数目与节点的数目相同。在之后的示例程序中还会展示如何定义高阶单元，在高阶单元中，自由度与边、面或体有关，与节点数目关联不大。



在有限元法中，自由度（degree of freedom）有两层意思：

* 在有限元法中，解可以表示为形函数的线性组合，即$u_h(x)=\sum^{N-1}_{j=0} U_j \phi_j(x)$。在此，$U_j$是展开系数向量，值未知。需要组成线性方程组进行求解，因此 它们被称为“未知量”或者“自由度”。
* 有限元法的数学解释为：寻求一个解空间$\u_h \in V_h$，对所有的测试函数$\phi_h \in V_h$，满足方程$a(u_h,\phi_h)=(f,\phi_h)$。为了求解这个问题，我们需要选取空间$V_h$中的一组基，在有限元法中基要定义在网格与单元之上。从这个意义上理解，自由度即为基函数的空间$V_h$的度量。在deal.II中，类DoFHandler提供了自由度的相关功能。



在deal.II中，由于函数库已经完成了相关的实现，所以定义自由度是一件非常简单的事情。基本上，要做的事情首先是定义一个有限元类，然后通过DoFHandler::distribute_dofs函数来实现自由度的定义。DoFHandler可以解决与自由度有关的相关问题，如：“总体上有多少个自由度”，“局部单元中的下标在总体坐标系下下标是多少”。



## 稀疏性

在定义了自由度之后，便需要根据微分方程，计算相应的线性方程组，在Step-3中会讲解相应的过程。在此想要叙述的是有限元法中重要的一点：有限元法线性方程组中的矩阵是稀疏的，大多数系数为0 。



更为准确地说，稀疏矩阵每一行/列中非零项地个数是有上限的，这个上限是一个与总自由度无关的常数，一般与数值方法有关。比如求解拉普拉斯方程的5点差分法，用有限元求解，其矩阵每一行中非零项为5个。对于更加复杂的问题，如Step-22中的斯托克斯问题，每行中的非零项可以达到数百个。但是重点是这个数字与矩阵的大小无关：如果网格变得精细，每行中非零项的个数依然保持不变。



假设最后形成的线性方程组种矩阵有$N$行，每行中非零项的个数都小于一个常数，那么其空间复杂度为$O(N)$，与向量相乘的时间复杂度也为$O(N)$。因为求解线性方程组只需要固定次数的矩阵向量乘法，因此求解的时间复杂度也为$O(N)$。这是相当可观的。如果矩阵不是稀疏的，那么求解的时间复杂度为$O(N^s)$，其中$s>1$，这还是在采用了如多重网格法这样的求解器的情况下。



在有限元法中，稀疏性来源于形函数都定义在局部的单元中而非总体，而且差分操作只与相邻的单元有关。即，未知量$i$，$j$所在的单元如果没有毗邻，则在最终的矩阵中$A_{ij}=0$。



## 自由度计算

DoFHandler类默认随机处理网格上的自由度，因此，稀疏矩阵的排列也不是最优的。为了展示这一点，在代码中会输出矩阵的“稀疏模式”（即矩阵非零项的位置）。



对于大多数算法，自由度的排列是无关紧要的，如CG算法。但是也有一些算法对自由度的排布较为敏感，如SSOR、LU分解、Cholesky分解等。因此，deal.II中也包含了可以重新排列未知量的工具：DoFRenumbering。这个函数相比与随机排布有一定的优化。



# 示例程序

同上个例子相同的，需要引入以下的库：

```c++
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
```

接下来引入处理自由度的库，这也是这个例子中的核心：

```c++
#include <deal.II/dofs/dof_handler.h>
```

接下来引入与单元有关的库，在这个库种，包含了各种维度，各种阶次的拉格朗日单元的实现：

```c++
#include <deal.II/fe/fe_q.h>
```

处理自由度：

```c++
#include <deal.II/dofs/dof_tools.h>
```
对非零项进行可视化：
```c++
#include <deal.II/lac/sparse_matrix.h>
```

存储稀疏模式：
```c++
#include <deal.II/lac/dynamic_sparsity_pattern.h>
```

对自由度进行重新排列：
```c++
#include <deal.II/dofs/dof_renumbering.h>
```

文件输出：
```c++
#include <fstream>
```

命名空间：
```c++
using namespace dealii;
```



## 网格生成

在此网格生成采用的是step-1中的网格生成方法，生成了一个圆环域并进行局部优化。

```c++
void make_grid(Triangulation<2> &triangulation)
{
  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 5);
  for (unsigned int step = 0; step < 3; ++step)
    {
      for (auto &cell : triangulation.active_cell_iterators())
        for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v)
          {
            const double distance_from_center =
              center.distance(cell->vertex(v));
            if (std::fabs(distance_from_center - inner_radius) < 1e-10)
              {
                cell->set_refine_flag();
                break;
              }
          }
      triangulation.execute_coarsening_and_refinement();
    }
}
```



## 自由度处理函数DoFHandler

至此，网格生成完成，节点与单元的几何性质已经确定好了。为了能够运用数值算法，需要将节点等信息与自由度连接起来。DoFHandler就是来完成这项工作的。在处理自由度信息之前，首先需要建立有限元单元。类FE_Q可以处理拉格朗日单元，它需要一个参数表示阶数，在本程序中为1表明单元的形函数是一阶多项式。



首先需要建立有限元单元，并将其传给DoFHandler：

```c++
void distribute_dofs(DoFHandler<2> &dof_handler)
{
  const FE_Q<2> finite_element(1);
  dof_handler.distribute_dofs(finite_element);
```



网格中每个节点都会有一个形函数，比如我们要求解拉普拉斯方程，那么最终矩阵中的非零项便为每对形函数梯度的和。显而易见地，在单元相邻相对应的位置，矩阵中才会有非零元素。由于单元与节点编号的随机性，最终非零元素的排列可能不尽如人意，那么还可以进一步进行处理。



首先需要一个数据结构来存储非零元素的位置，那么便可以根据这个存储下来的“稀疏模式”来存储矩阵中的值。SparsityPattern类便是完成这项任务的，但是其缺点在于需要预设每行会出现的非零项个数。在二维情况下用DoFHandler::max_couplings_between_dofs()函数可以较为准确的估计，但是三维情况下其估计的值严重偏大，带来大量的内存浪费。为了避免这样的情况，在示例程序中采用了DynamicSparsityPattern类，其采用了不同的数据结构，可以复制到SparsityPattern类中，其初始化需要给出矩阵的大小。

```c++
DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
                                                dof_handler.n_dofs());
```

随后生成非零元素所在的位置：

```c++
DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);
```

从DynamicSparsityPattern中复制“稀疏模式”：

```c++
SparsityPattern sparsity_pattern;
sparsity_pattern.copy_from(dynamic_sparsity_pattern);
```

随后，便可以讲结果写到svg文件中：

```c++
std::ofstream out("sparsity_pattern1.svg");
sparsity_pattern.print_svg(out);
```

在生成的svg文件中，非零元素的位置是红色的小正方形。



观察输出的文件可以发现矩阵是对称的。这并不奇怪，因为我们并没有给DoFTools::make_sparsity_pattern任何生成非对称矩阵的设定。此外还可以注意到矩阵中有明显的几个区域，表明自由度的排布从稀疏的网格一直排到精细的网格。



## 自由度重新排布

在之前产生的矩阵中，非零元素离主对角线较远。对于一些算法，如LU分解或者GS迭代来说，这是极为不利的。但是这种情况可以通过重新排布未知量来优化。



$(i,j)$表示矩阵中非零项的位置，注意到单元$i$，$j$相邻时，其在矩阵中的值才会是非零的。同时，我们希望非零元素在矩阵的对角线附近。这就意味着在排布的时候，毗邻的单元$i$、$j$的序号不应该相差太大。



这个可以通过一个简单的算法来实现，首先取一个初始点，如果它的相邻的单元在它之后进行排列，便排的离它尽量近，随后对它相邻的单元进行同样的处理，以此类推。



在示例程序中采用了一个较为复杂的算法，由Cuthill和McKee提出。在代码中调用DoFRenumbering::Cuthill_McKee即可：

```c++
void renumber_dofs(DoFHandler<2> &dof_handler)
{
  DoFRenumbering::Cuthill_McKee(dof_handler);
  DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
                                                  dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);
  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dynamic_sparsity_pattern);
  std::ofstream out("sparsity_pattern2.svg");
  sparsity_pattern.print_svg(out);
}
```

那么，同样地会输出矩阵中非零元素的排列。



值得注意的是，在DoFRenumbering类中提供了很多其他的算法，实现特定的功能。比如，所有的非零元素排列在上三角或者下三角，这对与对称的结构来说是不可能的，但在一些特殊的情况下，如流体从进流边界流到出流边界的情况。

## 主函数 

最后，需要主函数来调用各个部分，形成最终的程序：

```c++
int main()
{
  Triangulation<2> triangulation;
  make_grid(triangulation);
  DoFHandler<2> dof_handler(triangulation);
  distribute_dofs(dof_handler);
  renumber_dofs(dof_handler);
}
```



# 结果

运行程序后会生成两个“稀疏模式”，svg文件可以在浏览器中打开。



在第一张图片中分为两个区域，左上角与其他部分的网格，由于细化的程度不一样，所以离对角线的距离也不一样。第二张图片在自由度重新排列之后，非零元素离对角线的距离更近。



## 可能的拓展

首先，可以尝试更改单元的阶数，比如改成3或5。

或者，如果想要观察网格细化对矩阵的影响，可以尝试

最后，还可以尝试其他的排布算法。

在可视化的方式上，可以采用gnuplot，在源代码种讲print_svg()改为print_gnuplot()。

```bash
examples/step-2> gnuplot
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
gnuplot> set style data points
gnuplot> plot "sparsity_pattern.1"
```



