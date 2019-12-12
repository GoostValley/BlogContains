---
title: deal.II学习--Step-1
date: 2019-12-03 11:10:18
tags:
- dealii
categories:
- dealii
mathjax: true
---

# 引文

deal.II 是一个开源的有限元法求解器，支持大规模并行计算，自适应网格。采用C++编写，实现优雅。其文档完整丰富，文档共有三个级别，由浅入深：

* tutorial：一系列的教学程序，共有64步教学步骤，通过tutorial的学习可以对dealii有整体的认识。
* manual：对每个类以及相应的函数的介绍。适合用于查询类与函数的具体用法。
* Modules：介绍了实现某一个功能需要用到的一系列类与函数，比如 Sparsity patterns 介绍了存储稀疏矩阵相关的内容。

这个系列将学习 deal.II 的 [tutorial](https://www.dealii.org/current/doxygen/deal.II/Tutorial.html) 部分，有关dealii的安装，如果要添加第三方库，可以参考[这篇文章](https://ghostvalley.top/2019/11/02/pi-BEM-dealii-deal2lkit的安装/)，使用candi安装；此外，也可以参考官网上的安装方法。tutorial 的 64 个例子的源代码在项目的 example 文件夹。



Step-1 的源代码在 example/step-1，编译代码只需要在命令行中输入

```bash
cmake .
make
./step-1
```

或者

```bash
cmake .
make run
```

此后项目代码的编译都是如此。

# 代码的主要内容

Step-1 的主要内容是生成网格，在这个程序中，共尝试生成了两种不同的二维网格，一种是正方形网格，一种是圆环网格，并展示了对网格进行迭代细化的过程。

首先需要加载头文件，最为重要的是 Triangulation 类，其用于生成单元，Triangulation<1>表示一维单元，以此类推，可以表示二维单元和三维单元，对于边界单元，Triangulation<1,2>表示面边界上的曲线单元，Triangulation<2,3>表示体边界上的曲面单元。

``` c++ 
# inlcude <deal.II/grid/tria.h>
```

这两个头文件实现了单元的存取与迭代：

```c++
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
```

生成一些具有标准几何形状的网格：

```c++
#include <deal.II/grid/grid_generator.h>
```

网格的格式化输出：

```c++
#include <deal.II/grid/grid_out.h>
```

C++的标准输入输出：

```c++
#include <iostream>
#include <fstream>
```

数学运算库：

```c++
#include <cmath>
```

最后加上命名空间 dealii ，确保调用的函数来自 dealii 库。

```c++
using namespace dealii;
```

## 正方形域

在接下来的代码中，将产生单位正方形域，并进行全局细化，最终生成网格。首先申明二维单元triangulation，然后设置单元为立方体形状，并对单元进行四次全局细化。最终将网格绘制成eps图像输出。

```c++
void first_grid()
{
  Triangulation<2> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(4);
  std::ofstream out("grid-1.eps");
  GridOut       grid_out;
  grid_out.write_eps(triangulation, out);
  std::cout << "Grid written to grid-1.eps" << std::endl;
}
```



从最终结果的图像中可以看出，一共有$4^4=256$个网格，在全局细化中，每个正方形单元都会变成四个小正方形单元，那么细化四次便是最终的结果。

## 圆环域

接下来介绍圆环域网格的生成过程：

```c++
void second_grid()
{
  Triangulation<2> triangulation;
  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 10);
```



同样的，先声明单元triangulation，圆环域的中心点，内径以及外径，然后用 hyper_shell 函数生成网格，在这里10指将圆环沿环向等分为10份。



一般情况下，Triangulation 假设所有边界都是直线。但是我们可以定义一些特别的几何形状，这样在网格细化的过程中新的点和单元会更加贴合实际位置。这时候需要用到manifold indicator的方式，用来标明已有特殊几何形状的单元。在这里生成圆环域所用的hyper_shell函数在源码中与SphericalManifold相配合，因此可以实现圆弧的细化。



利用GridGenerator函数（如 GridGenerator::hyper_shell 或者 GridGenerator::hyper_ball）生成网格时，每个网格都有一个默认的 flat_manifold_id, 如果网格的边界是直线，可以手动将其设为1，否则网格的边界默认为曲线。



随后对网格进行细化，首先确定细化的步数：

```c++
  for (unsigned int step = 0; step < 5; ++step)
    {
```



然后循环遍历所有网格：

```c++
for (auto it = triangulation.active_cell_iterators().begin();
     it != triangulation.active_cell_iterators().end();
     ++it)
  {
    auto cell = *it;
    // Then a miracle occurs...
  }
```

在此通过 auto 来自动设置变量类型，在变量类型非常复杂的时候auto关键字是极为有效的。active_cell_iterators的意思是激活的网格，在网格细化的过程中，一些粗网格任然保留，但是其状态不是激活的，因此这个循环遍历的是当前最细的网格。上面的写法需要声明循环的初始与结束条件，在此还可以采用更为简便的方式：

```c++
for (auto &cell : triangulation.active_cell_iterators())
  {
```

随后需要遍历网格中的所有节点，利用GeometryInfo可以得到每个网格中的节点数目。

```c++
for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v)
  {
```

在此细化的策略为：如果网格是靠近圆环内边沿的网格，即有节点在圆环内径上，那么就加密这个网格。由于计算机机内小数表示有误差，因此不能直接判断等于，而是差值绝对值小于一个极小数的形式。

```c++
      const double distance_from_center =
        center.distance(cell->vertex(v));

      if (std::fabs(distance_from_center - inner_radius) < 1e-10)
        {
          cell->set_refine_flag();
          break;
        }
    }
}
```

通过循环已经标记好了需要细化的节点，那么执行网格的稀疏与细化函数即可：

```c++
  triangulation.execute_coarsening_and_refinement();
}
```



最后，输出网格的图片。

```c++
  std::ofstream out("grid-2.eps");
  GridOut       grid_out;
  grid_out.write_eps(triangulation, out);

  std::cout << "Grid written to grid-2.eps" << std::endl;
}
```



## main函数

之前实现了两个网格的生成函数，在main函数中进行调用：

```c++
int main()
{
  first_grid();
  second_grid();
}
```

# 结果展示

运行程序后可以得到（grid-1.eps 和 grid-2.eps）。



<center>

<div style="width:300px;display:inline-block;">{%asset_img  grid-1.png%}</div>

<div style="width:300px;display:inline-block;margin-left:10px;">{%asset_img  grid-2.png%}</div>

</center>



# 拓展

对于圆环域的例子，修改单元的manifold_id，可以得到不同的结果。

```c++
void second_grid()
{
  Triangulation<2> triangulation;

  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 10);
  for (unsigned int step = 0; step < 5; ++step)
    {
      for (auto &cell : triangulation.active_cell_iterators())
        {
          cell->set_manifold_id(1);
        }
      triangulation.refine_global(1);
    }


  std::ofstream out("grid-3.eps");
  GridOut       grid_out;
  grid_out.write_eps(triangulation, out);

  std::cout << "Grid written to grid-3.eps" << std::endl;
}
```

在此将manifold_id设为1，可以看出单元在径向上有明显的分割痕迹。

<div style="width:300px; margin:auto">{% asset_img grid-3.png%}</div>

将上述代码中manifold_id设为0，则网格更加圆滑。

<div style="width:300px; margin:auto">{% asset_img grid-4.png%}</div>

但是从上述两幅图可以看出，在内径和外径上，网格一直都很光滑，这说明设置manifold_id没有影响到边界上的几何要素。