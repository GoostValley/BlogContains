---
title: hexo图片并排显示
date: 2019-12-12 09:49:47
tags:
- hexo 
categories:
- hexo 
type: "picture"
---

## 问题

在用hexo编写博客时，如果不加设置，插入两张图片时一行一个，显得非常难看。因此需要并排显示。



<div style="width:300px;margin:auto;align:center">{% asset_img 1.PNG %}</div>





## 解决方法1

hexo next主题中提供了group picture工具，具体的实现如

```markdown
{% gp 1-2 %}
{%asset_img  sparsity_pattern1.png%}{%asset_img  sparsity_pattern1.png%}
{% endgp %}
```

效果如以下所示，图片可以并排显示，但是无法紧贴在一起，博主尝试过修改group_picture.styl文件，但是没有调整出理想的效果，作罢。



<center>

<div style="margin: auto;">{% asset_img 3.PNG %}</div>

</center>



## 解决方法2

利用inline-block的方式，来生成并排的图片,图片可以居中显示，且靠的比较紧密。

```markdown
<center>
<div style="display:inline-block;">{%asset_img  sparsity_pattern1.png%}</div>
<div style="display:inline-block;margin-left:10px;">{%asset_img  sparsity_pattern2.png%}</div>
</center>

```

显示效果不错。



<div style="width:600px;margin:auto;align:center">{% asset_img 2.PNG %}</div>

