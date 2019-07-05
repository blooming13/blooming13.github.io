---
title: Matlab实验报告
date: 2019-07-01 08:06:29
tags: test
mathjax: true
---



**一、**    **实验题目**

利用中心有限差分法求解
$$
\begin{array}{l}{-u^{\prime \prime}(x)+\left(3+x^{2}\right) u(x)=f(x)} \\ {u(0)=1, u(1)=-1}\end{array}
$$
真解为                     
$$
u(x)=\cos (\pi x)
$$
空间步长h分别取为1/4，1/8，1/16，1/32，在同一个图上画出真解和不同空间步长网格上的数值解，计算数值解与真解在网格节点处的最大模误差，列表总结误差变化规律。







 二、    **实现算法**

用二阶中心差商代替二阶微商
$$
{ {\left[\frac{ { {d}^{2}}u}{d{ {x}^{2}}}\right]}_{i}}\approx \frac{u\left({ {x}_{i+1}}\right)-2u\left({ {x}_{i}}\right)+u\left({ {x}_{i-1}}\right)}{ { {h}^{2}}}
$$
原方程等价于
$$
-\frac{ { {u}_{i+1}}-2{ {u}_{i}}+{ {u}_{i-1}}}{ { {h}^{2}}}+{ {u}_{j}}={ {f}_{i}},i=1,\cdots
,n-1
$$

$$
{ {u}_{0}}=1,\quad { {u}_{n}}=\text{-}1
$$
化简得
$$
-\frac{1}{ { {h}^{2}}}{ {u}_{i-1}}+(\frac{2}{ { {h}^{2}}}+1){ {u}_{i}}-\frac{1}{ { {h}^{2}}}{ {u}_{i+1}}={ {f}_{i}},\quad i=2,\cdots ,n-2
$$

$$
{ {u}_{0}}=1,\quad { {u}_{n}}=\text{-}1
$$

由该差分格式，等价于求解线性代数方程组的三对角稀疏矩阵,即把 [0,1] 均匀分成 n 段，每段长度为 h, 用中心有限差分离散后,的矩阵向量形式为：
$$
\left[ \begin{matrix}
   1 & 0 & 0 & \cdots  & \cdots  & \cdots  & 0  \\
   { {a}_{1}} & { {d}_{1}} & { {c}_{1}} & 0 & {} & {} & \vdots   \\
   0 & { {a}_{2}} & { {d}_{2}} & { {c}_{2}} & 0 & {} & \vdots   \\
   \vdots  & \ddots  & \ddots  & \ddots  & \ddots  & \ddots  & \vdots   \\
   \vdots  & {} & 0 & { {a}_{n-2}} & { {d}_{n-2}} & { {c}_{n-2}} & 0  \\
   \vdots  & {} & {} & 0 & { {a}_{n-1}} & { {d}_{n-1}} & { {c}_{n-1}}  \\
   0 & \cdots  & \cdots  & \cdots  & 0 & 0 & 1  \\
\end{matrix} \right]\left[ \begin{matrix}
   { {u}_{0}}  \\
   { {u}_{1}}  \\
   { {u}_{2}}  \\
   \vdots   \\
   { {u}_{n-2}}  \\
   { {u}_{n-1}}  \\
   { {u}_{n}}  \\
\end{matrix} \right]=\left[ \begin{matrix}
   \alpha   \\
   { {f}_{1}}  \\
   { {f}_{2}}  \\
   \vdots   \\
   { {f}_{n-2}}  \\
   { {f}_{n-1}}  \\
   \beta   \\
\end{matrix} \right]
$$
其中
$$
{ {a}_{i}}=-\frac{1}{ { {h}^{2}}},\quad { {d}_{i}}=\frac{2}{ { {h}^{2}}}+3\text{+}{ {\text{x}}_{i}}^{2},\quad { {c}_{i}}=-\frac{1}{ { {h}^{2}}}
$$
 用A = sparse(i,j,v,m,n)命令可创建m x n阶稀疏矩阵

误差后处理: 

记误差向量
$$
e=\left(e_{1}, \cdots, e_{n-1}\right)^{T} \in R^{n-1}
$$
，其中  $ e_{i}=u\left(x_{i}\right)-u_{i}$  为真解与数值解在网格节点$x_{i}$ 处的误差.

为了度量误差，引入最大模范数

$ \|e\|_{C}=\max _{1 \leq i \leq n-1}\left|e_{i}\right| $

 

 



**三、**    **程序代码**

**1)**     **首先定义模型数据**

```matlab
function model = model_data(l, r)

%% MODEL_DATA 模型数据

% 

% Input

% -----

%  l: 区间左端点

%  r: 区间右端点

 

L = l;

R = r;

 

model = struct('init_mesh', @init_mesh, 'solution', @solution,...

    'source', @source);

 

function [X, h] = init_mesh(NS) %网格剖分

    X = linspace(L, R,NS+1)';

    h = (R - L)/NS;

end

 

function u = solution(x) %真解函数

    u = cos(pi*x);

end

 

function f = source(x) %右端项函数

    f = (3+pi*pi+x.*x).*cos(pi*x);

end

 

end
```

**2)** **编写测试框架**

```matlab
%% 测试脚本 FD1d_bvp_test.m

 

clear all

close all

 

 % 初始化相关数据

NS = [4, 8, 16, 32]; %设置剖分步长 h=1/4，h=1/8，h=1/16，h=1/32

L = 0;     %左端点

R = 1;     %右端点

 

model = model_data(L, R);

 

emax = zeros(4,1);

 

 

%% 求解并计算误差

for i = 1:4

    [uh, x] = FD1d_bvp(model, NS(i));

    [emax(i)]=FD1d_error(model.solution, uh, x);

    X{i} = x;

    U{i} = uh;

end

 

u = model.solution(X{4});

 

%% 显示真解及不同网格剖分下的数值解

%画图

plot( X{4}, u, '-k*',X{1}, U{1}, '-ro',...

     X{2},U{2}, '-gs', X{3}, U{3}, '-bd',...

    X{4}, U{4}, '-ch');

title('The solution plot');

xlabel('x');  ylabel('u');

legend('exact','h=1/4','h=1/8','h=1/16','h=1/32');

 

%% 显示误差

format shorte

disp('     emax    ');

disp([emax]);

 

%% 显示误差变化规律

ratio = emax(1:end-1)./emax(2:end);

disp('最大模误差变化规律');

disp(ratio);
```

**3)** **实现核心算法**

```matlab
function [uh, x] = FD1d_bvp(model, NS)

%% FD1d_bvp 利用中心差分格式求解两点边值问题.

%

% Input

% -----

%   model : 模型数据

%   NS    ：网格剖分段数

% Output

% ------

%   uh ：列向量，长度为 NS+1， 解向量

%   x  : 列向量，长度为 NS+1， 网格节点

%      

 

[x, h] = model.init_mesh(NS);

NV = NS + 1;

 

%

%  创建线性差分方程组系数矩阵

%

 

c1 = -1/h/h*ones(1,2*NV-4);

c2 = (2/h/h+3)*ones(1,NV-2) + (x(2:end-1).*x(2:end-1))';

 

i=[1:NV,2:NV-1,2:NV-1];

j=[1:NV,3:NV,1:NV-2];

v=[1,c2,1,c1];

A = sparse(i,j,v,NV,NV);

 

%

%  创建线性差分方程组右端项

%

 

rhs = model.source(x);

rhs(1) = model.solution(x(1));

rhs(end) = model.solution(x(end));

 

%

%  求解上述代数系统.

%

 

uh = A \ rhs;

end
```

 

**4)** **编写误差分析函数**

```matlab
function [emax] = FD1d_error(solution, uh, X)

%% FD1d_error 计算数值解与真解的最大模误差

%

% Input

% -----

%   model.solution : 真解模型数据

%   uh ：列向量，长度为 NS+1， 解向量

%   X  : 列向量，长度为 NS+1， 网格节点

% Output

% ------

%   emax ：列向量，长度为 4， 最大模误差

%      

 

NN = length(X);

h = (X(end) - X(1))/(NN -1);

u = solution(X);

ee= u - uh;

 

emax=max(abs(ee));

 

end
```



 

 

**四、** **实验结果与分析**

   真解及不同剖分步长下数值解函数的图像：

   总体图：                              

<img src=".\imgs\1.png" width = "500" height = "500" div align=center/>

 

局部图：

<img src=".\imgs\2.png" width = "500" height = "500" div align=center/>



不同剖分长度下最大模误差及误差变化规律：

<img src=".\imgs\3.png" width = "250" height = "300" div align=center/>

 

 

实验结果表明:

​    由误差范数值，函数整体、局部图可知，数值解函数逼近真解函数的程度与区间的剖分程度有关，区间剖分得越细，数值解函数与真值解函数越接近，且误差越小。由误差变化规律可看出，用二阶中心差分求解该一维问题达到了$ O\left( { {\text{h}}^{2}} \right) $的收敛阶，该实验也验证了二阶中心差分格式的稳定性与收敛性。

