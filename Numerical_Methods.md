#数值分析笔记
[TOC]

-------

## 1.误差分析
###1.1.绝对误差    
$$
E_p=|p-\tilde{p}|
$$
###1.2.相对误差    
$$
R_p=\frac{|p-\tilde{p}|}{|p|}
$$
###1.3.有效数字    
$$
\frac{|p-\tilde{p}|}{|p|}<\frac{10^{1-d}}{2}
$$

0.0123的有效位数为3位。
###1.4.截断误差
用一个基本表达式替换一个复杂表达式时引入的误差，例如用截断泰勒级数表示某个函数式。
###1.5.舍入舍去    
舍去  $0.123456\approx0.12345$

舍入  $0.123456\approx0.12346$

###1.6.精度损失
$p=0.1234567,q=0.1234568,p-q=0.0000001$精度由7位有效数字降到了1位。
典型例题，利用6位有效数字计算$f(500)$和$g(500)$，前者精度远小于后者。截断泰勒级数和***嵌套乘法***(霍纳方法）在一定程度上可以减小损失。
$$
f(x)=x(\sqrt{x+1}-\sqrt{x}),g(x)=\frac{x}{\sqrt{x+1}+\sqrt{x}}
$$
###1.7.大O表示法
即如果$|f(h)| \leq C|g(h)|$对于任意的$h \leq c$，则可以用$f(h)=O(g(h))$表示，例如：
$$
{f(x)}_ {x\rightarrow0}=1+2x+2x^2+3x^3+……\approx1+2x+O(x^2)
$$
### 1.8.误差传播
加法运算的绝对误差为每个加数绝对误差之和，乘法运算的相对误差大致为每个乘数的相对误差之和。

-------
## 2.非线性方程的解法
###2.1.解x=g(x)
####2.2.1.迭代
$$
\begin{array}{a}
p_n=g(p_{n-1})\\
P=p_\infty \\
P=g(P)\\
\end{array}
$$
上述过程称为不动点迭代，$P$称为不动点。

####2.1.2.不动点存在性
* 存在不动点$\leftarrow$对于所有的$x \in [a,b]$，$y=g(x) \in [a,b]$

* 存在唯一被动点$\leftarrow​$满足上述条件的前提下，$g’(x)​$在$(a,b)​$有定义且在$K<1​$，在$(a,b)​$上，$|g’(x)| \leq K<1​$恒成立。

  上述仅仅是充分不必要条件。

####2.1.3.不动点收敛性
1. $g,g’ \in C[a,b]$

2. $K>0$

3. $p_0 \in (a,b)$   

4. $\forall x \in [a,b],g(x) \in [a,b]$   

满足上述五点，则$[a,b]$内一定存在不动点$P(P \neq p_0)$：
* $\forall x \in [a,b],|g’(x)| \leq K<1,p_n=g(p_{n-1}) \rightarrow p_\infty=P$此时$P$点唯一，称为***吸引不动点。***
* $\forall x \in [a,b],|g’(x)|>K,p_n=g(p_{n-1}) \rightarrow {p_n}_{n \rightarrow \infty} \neq P$此时，迭代发散，不会收敛到，$P$称为排斥不动点。
* $\exists x \in [a,b],|g’(x)|=1,if\quad p_0>P \rightarrow p_\infty=P,if\quad p_0<P \rightarrow p_\infty \neq P$，当存在斜率绝对值为1时，收敛与$p_0$的初值有关。
####2.1.4.收敛结束判断依据
常用连续项之差$err=|p_n-p_{n-1}|<tol $最为迭代结束条件，但是这歌判断依据并不能保证精度，$tol$表示容错程度。
### 2.2.解f(x)=0
####2.2.1.二分法
取区间$[a_0,b_0]$和中点$c_0=\frac{a_0+b_0}{2}$，计算$f(a_0),f(b_0),f(c_0)$，如果$f(a_0)f(b_0)<0$，取$f(a_0),f(b_0)$中与$f(c_0)$同号的点，用点$c_0$取而代之，更新区间为$[a_1,b_1]$。例如，$f(b_0)f(c_0)>0 \rightarrow a_1=a_0 \quad b_1=c_0$。最终，找到无穷小的区间$[a_\infty,b_\infty]$，即得到解$r=c_\infty$
####2.2.2.二分法误差界(精度）
$$
|r-c_n| \leq \frac{b-a}{2^{n+1}}
$$
误差随迭代次数$n$的增加而减小，因此可以根据精度要求来控制迭代次数来实现，$\delta$表示要求精确度，则至少迭代次数$N$为：
$$
N=int[\frac{\ln{(b-a)}-\ln{\delta}}{\ln{2}}]
$$

####2.2.3.试值法
取中点的速度太慢，因此考虑用弦与横轴的焦点代替中点。弦的斜率为$m=\frac{f(b)-f(a)}{b-a}$，弦与横轴交点为$c=b-\frac{f(b)(b-a)}{f(b)-f(a)}$，$c$的替换思路与中值法相同。试值法最终收敛。
####2.2.4.牛顿-拉夫森法
之前的求根方法是在一个区间$[a,b]$上通过不断的缩小区间来实现的，称为***全局收敛法***。之后的方法是通过先估计一个近似根$p_0$，通过不断迭代得到最终值，称为***局部收敛法***。
将区间是为一个点$p_0$，$p_0$处切线点斜率为$m=f’(p_0)$，切线交横坐标于$p_1=p_0-\frac{p_0}{f’(p_0)}$，通过迭代最终$r=p_\infty$，此方法称为牛顿-拉夫森法，迭代公式为
$$
p_{n+1}=p_n-\frac{f(p_n)}{f’(p_n)}
$$
其收敛的充分条件为：
$$
\frac{|f(x)f’’(x)|}{|f’(x)|^2}<1
$$

####2.2.5.割线法
考虑到牛顿-拉夫森法需要求导，割线法用$m=\frac{f(p_n)-f(p_{n-1})}{p_n-p_{n-1}}$来代替斜率$m=\frac{f(b)-f(a)}{b-a}$，迭代公式更改为
$$
p{n+1}=p_n-\frac{f(p_n)}{m}=p_n-\frac{f(p_n)(p_n-p_{n-1})}{f(p_n)-f(p_{n-1})}
$$

####2.2.6.收敛速度

若对于$A \neq 0 \quad \exists R>0$有
$$
lim_{n \rightarrow \infty} \frac{|p-p_{n+1}|}{|p-p_n|^R}
$$
若$R=1$，则称为线性收敛；若$R=2$，则称为***二次收敛***。对于单根解$r$，牛顿-拉夫森法为二次收敛，对于$M$重根的解$r_1,r_2,...,r_M$，牛顿-拉夫森法为线性收敛，通过改进算法如下，可以使$M$重根二次收敛。
$$
p_{n+1}=p_n-M\frac{f(p_n)}{f’(p_n)}
$$

###2.3各种方法的循环结束条件

一般情况下，在程序中设定多个结束判断条件，例如对于$f(x)=0$的根的求解，$|r-p_n|<\delta$和$|f(p_n)|<\epsilon$的并集或者交集是两种判断方式。另外，一定会有最大的循环次数上限。

-------

##3.线性方程组
### 3.1.坐标变换

$$
R_x(\alpha)=
\begin{bmatrix}
1&0&0\\
0&\cos{\alpha}&-\sin{\alpha}\\
0&\sin{\alpha}&\cos{\alpha}\\
\end{bmatrix}

\quad
R_y(\beta)=
\begin{bmatrix}
\cos{\beta}&0&\sin{\beta}\\
0&1&0\\
-\sin{\beta}&0&\cos{\beta}\\
\end{bmatrix}
\quad

R_z(\gamma)=
\begin{bmatrix}
\cos{\gamma}&-\sin{\gamma}&0\\
\sin{\gamma}&\cos{\gamma}&0\\
0&0&1\\
\end{bmatrix}
$$

这三个矩阵分别是以角度$\alpha,\beta,\gamma​$ 绕$x,y,z​$ 轴旋转的变换矩阵。
$$
\begin{bmatrix}
x'\\y'\\z'\\
\end{bmatrix}
=R
\begin{bmatrix}
x\\y\\z\\
\end{bmatrix}
\quad
R=R_x(\alpha)\quad or \quad R_y(\beta) \quad or \quad R_z(\gamma)
$$

###3.2.三角回代法

对于三角系数矩阵可以使用***回代法***求解方程，例如下面的上三角方程组：
$$
\begin{bmatrix}
a_{11}&a_{12}&a_{13}\\
0&a_{22}&a_{23}\\
0&0&a_{33}\\
\end{bmatrix}
\begin{bmatrix}
x_{1}\\x_{2}\\x_{3}\\
\end{bmatrix}
=\begin{bmatrix}
b_{1}\\b_{2}\\b_{3}\\
\end{bmatrix}\\
\rightarrow
x_3=\frac{b_3}{a_{33}}\\
\rightarrow
x_2=\frac{b_2-a_{23}x_3}{a_{22}}\\
\rightarrow
x_1=\frac{b_1-a_{13}x_3-a_{12}x_2}{a_{11}}
$$
具体表达式写作
$$
x_k=\frac{b_k-\sum_{j=k+1}^{N}a_{kj}x_j}{a_kk}\quad k=N-1,N-2,\dots,1
$$
同理可以求解下三角方程组：
$$
\begin{bmatrix}
a_{11}&0&0\\
a_{21}&a_{22}&0\\
a_{31}&a_{32}&a_{33}\\
\end{bmatrix}
\begin{bmatrix}
x_{1}\\x_{2}\\x_{3}\\
\end{bmatrix}
=\begin{bmatrix}
b_{1}\\b_{2}\\b_{3}\\
\end{bmatrix}\\
$$

###3.3.高斯消去法

$$
\begin{bmatrix}
a_{11}&a_{12}&a_{13}\\
a_{21}&a_{22}&a_{23}\\
a_{31}&a_{32}&a_{33}\\
\end{bmatrix}
\begin{bmatrix}
x_{1}\\x_{2}\\x_{3}\\
\end{bmatrix}
=\begin{bmatrix}
b_{1}\\b_{2}\\b_{3}\\
\end{bmatrix}\\
$$

建立增广矩阵
$$
\left[
\begin{array}{ccc|c}
a_{11}^{(1)}&a_{12}^{(1)}&a_{13}^{(1)}&b_1^{(1)}\\
a_{21}^{(1)}&a_{22}^{(1)}&a_{23}^{(1)}&b_2^{(1)}\\
a_{31}^{(1)}&a_{32}^{(1)}&a_{33}^{(1)}&b_3^{(1)}\\
\end{array}
\right]
\rightarrow
\left[
\begin{array}{ccc|c}
a_{11}^{(1)}&a_{12}^{(1)}&a_{13}^{(1)}&b_1^{(1)}\\
0&a_{22}^{(2)}&a_{23}^{(2)}&b_2^{(2)}\\
0&a_{32}^{(2)}&a_{33}^{(2)}&b_3^{(2)}\\
\end{array}
\right]\\
\rightarrow
\left[
\begin{array}{ccc|c}
a_{11}^{(1)}&a_{12}^{(1)}&a_{13}^{(1)}&b_1^{(1)}\\
0&a_{22}^{(2)}&a_{23}^{(2)}&b_2^{(2)}\\
0&0&a_{33}^{(3)}&b_3^{(3)}\\
\end{array}
\right]
\rightarrow
\quad三角回代法\quad
$$

####3.3.1.主元选择

平凡选主元：如果出现$a_{kk}^{k}=0$ ，就与$a_{kk}^{k}$ 同列的之后的元素互换，如果最终找不到可以互换的，就说明方程组无解。
$$
\left[
\begin{array}{ccc|c}
(0)&3&4&5\\
1&2&3&4\\
3&4&5&6\\
\end{array}
\right]
\rightarrow
\left[
\begin{array}{ccc|c}
1&2&3&4\\
(0)&3&4&5\\
3&4&5&6\\
\end{array}
\right]
$$
***偏序选主元***：为减小误差，将元素中绝对值最大的移动到主对角线。
$$
\left[
\begin{array}{ccc|c}
1&2&3&4\\
0&3&4&5\\
(3)&4&5&6\\
\end{array}
\right]
\rightarrow
\left[
\begin{array}{ccc|c}
(3)&4&5&6\\
1&2&3&4\\
0&3&4&5\\
\end{array}
\right]
$$

####3.3.2.病态矩阵

当系数矩阵$A$ 的行列式接近于0时，这样的方程组是***病态***的。
$$
\left\{
\begin{array}{a}
x+2y-2.00=0\\
2x+3y-3.40=0\\
\end{array}    
\right.
$$
其正确解为 $\left\{\begin{array}{a}x=0.8\\y=0.6\\\end{array}    \right.$。对于$\left\{\begin{array}{a}x=1\\y=0.48\\\end{array}    \right.\quad \left\{\begin{array}{a}x=0.85\\y=0.6\\\end{array}    \right.$两组解，带入方程，结果分别为$\left\{\begin{array}{a}-0.04\\0.04\\\end{array}    \right.\quad \left\{\begin{array}{a}0.5\\1\\\end{array}    \right.$

看上去前者更接近$\left\{\begin{array}{a}0\\0\\\end{array}    \right.$ ，但是，事实是后者更精确。综上所述，对于病态方程组，直接带入判断结果准确性是无效的。

###3.3.三角分解法

一般的系数矩阵$A$可以分解$A=LU$ 的形式，$L$表示下三角矩阵，$U$表示上三角矩阵。如果不能分解，将$A$乘上置换矩阵$P$即可，即$PA=LU$，其中$P$为一个置换矩阵，其各行各列有且只有一个元素为1，其他元素为0。
$$
\begin{bmatrix}
1&0&0\\
l_{21}&1&0\\
l_{31}&l_{32}&1\\
\end{bmatrix}
\begin{bmatrix}
u_{11}&u_{12}&u_{13}\\
0&u_{22}&u_{23}\\
0&0&u_{33}\\
\end{bmatrix}
\begin{bmatrix}
x_{1}\\x_{2}\\x_{3}\\
\end{bmatrix}
=P \begin{bmatrix}
b_{1}\\b_{2}\\b_{3}\\
\end{bmatrix}\\
$$
线性方程$AX=B$可以变换如下：
$$
AX=B\rightarrow PAX=PB\rightarrow LUX=PB\\
define:UX=Y
\rightarrow
\left\{ 
\begin{array}{a}
LY=PB \rightarrow Y=L^{-1}PB\\
UX=Y \rightarrow X=U^{-1}Y\\
\end{array}
\right.
$$
其中，乘以逆矩阵的方式只是数学公式上的表示，实际就是利用回代法解决。

###3.4.迭代法

####3.4.1.雅可比迭代

对于方程组
$$
\begin{bmatrix}
a_{11}&a_{12}&a_{13}\\
a_{21}&a_{22}&a_{23}\\
a_{31}&a_{32}&a_{33}\\
\end{bmatrix}
\begin{bmatrix}
x\\y\\z
\end{bmatrix}
=
\begin{bmatrix}
b_1\\b_2\\b_3
\end{bmatrix}
$$
雅可比迭代为：
$$
\left\{
\begin{array}{a}
x_{k+1}=\frac{b_1-a_{12}y_k-a_{13}z_k}{a_{11}}\\
y_{k+1}=\frac{b_2-a_{22}y_k-a_{23}z_k}{a_{21}}\\
z_{k+1}=\frac{b_3-a_{32}y_k-a_{33}z_k}{a_{31}}\\
\end{array}
\right.
$$

####3.4.2.高斯-赛德尔迭代

高斯-赛德尔迭代为：
$$
\left\{
\begin{array}{a}
x_{k+1}=\frac{b_1-a_{12}y_k-a_{13}z_k}{a_{11}}\\
y_{k+1}=\frac{b_2-a_{22}y_{k+1}-a_{23}z_{k}}{a_{21}}\\
z_{k+1}=\frac{b_3-a_{32}y_{k+1}-a_{33}z_{k+1}}{a_{31}}\\
\end{array}
\right.
$$
当系数矩阵具有严格对角优势的时候，两种迭代方法一点是收敛的，这是一个充分不必要条件。

严格对角优势是指矩阵的对角元素的绝对值是该行唯一的最大值。即：

####3.4.3.收敛条件

当系数矩阵具有***严格对角优势***的时候，两种迭代方法一点是收敛的，这是一个充分不必要条件。

严格对角优势是指矩阵的对角元素的绝对值大于该行其他元素绝对值之和。即：
$$
|a_{kk}|>
\sum^{N}_{
\begin{array}{a}
j=1\\ j \neq k
\end{array}
}
{|a_{kj}|}
\quad 
k=1,2,\dots,N
$$

------

##4.多项式插值

###4.1.霍纳方法

所有的函数都可以用$N$阶多项式$P_N(x)$和误差项$E_N(x)$的方式表示：
$$
f(x)=P_N(x)+E_N(x)
$$
利用霍纳方法改写$y=3x^3+x^2+4x+1=((3x+1)x+4)x+1$ 。

计算流程：
$$
\begin{array}{a}
y=a_Nx^N+a_{N-1}x^{N-1}+ \dots a_1x+a_0\\
\rightarrow b_N=a_N\\
\rightarrow b_{N-1}=b_N*x+a_{N-1}\\
 \dots  \dots \\
\rightarrow b_1=b_2*x+a_1\\
\rightarrow b_0=b_1*x+a_0\\
\rightarrow y=b_0
\end{array}
$$

###4.2.泰勒多项式

泰勒多项式是根据泰勒级数的而来的，以下是$N$阶的泰勒多项式。
$$
\left\{
\begin{array}{a}
f(x)=P_N(x)+E_N(x)\\
P_N(x)=f(x_0)+f'(x_0)(x-x_0)+\frac{f''(x_0)}{2!}{(x-x_0)}^2+\dots +\frac{f^{(N)}(x_0)}{N!}{(x-x_0)}^N\\
E_N(x)=\frac{f^{(N+1)}(c)}{(N+1)!}{(x-x_0)}^{(N+1)}\quad c \in (x,x_0)or(x_0,x)\\
\end{array}
\right.
$$
由于$E_N(x)$无法实际计算而得，利用误差界估计多项式逼近时的误差。
$$
\left\{
\begin{array}{a}
|E_N(x)|=\frac{|f^{(N+1)}(c)|}{(N+1)!}{|x-x_0|}^{(N+1)}\leq \frac{MR^{N+1}}{(N+1)!}\quad c \in (x,x_0)or(x_0,x) \\
M=\max|f^{(N+1)}(c)|\\
R=|x-x_0|\\
\end{array}
\right.
$$

实际计算的时候，函数的导数是未知的，可以利用已知函数通过的点建立方程组求解系数，这种方法称为插值。

对于点$(x_0,y_0)(x_1,y_1)\cdots(x_N,y_N)$一共$N+1$个点，用$N$阶多项式$y=a_Nx^N+a_{N-1}x^{N-1}+ \dots a_1x+a_0$插值。解系数方程组。
$$
\begin{bmatrix}
x_0^N&x_0^{N-1}&\cdots& x_0&1\\
x_1^N&x_1^{N-1}&\cdots &x_1&1\\
 \vdots & \vdots &\ddots& \vdots &\vdots\\
x_{N-1}^N&x_{N-1}^{N-1}&\cdots& x_{N-1}&1\\
x_N^N&x_N^{N-1}&\cdots& x_N&1\\
\end{bmatrix}
\begin{bmatrix}
a_N \\ a_{N-1}\ \\ \vdots \\ a_1 \\a_0
\end{bmatrix}
=\begin{bmatrix}
y_N \\ y_{N-1}\ \\ \vdots \\ y_1 \\y_0
\end{bmatrix}
$$

###4.3.拉格朗日多项式

泰勒多项式要求函数的导数已知，但是我们往往已知的是函数的部分（插值）点，利用矩阵估计系数计算量又过大，因此可以改用拉格朗日插值。两点$(x_0,y_0)(x_1,y_1)$的拉格朗日插值就是直线公式。
$$
y=y_0+\frac{y_1-y_0}{x_1-x_0}(x-x_0)=y_0\frac{x-x_1}{x_0-x_1}+y_1\frac{x-x_0}{x_1-x_0}
$$
升级表达式，三点$(x_0,y_0)(x_1,y_1)(x_2,y_2)$可以表示为，
$$
y=y_0+\frac{y_1-y_0}{x_1-x_0}(x-x_0)=y_0\frac{(x-x_1)(x-x_2)}{(x_0-x_1)(x_0-x_2)}+y_1\frac{(x-x_0)(x-x_2)}{(x_1-x_0)(x_1-x_2)}+y_2\frac{(x-x_0)(x-x_1)}{(x_2-x_0)(x_2-x_1)}
$$
用公式总结拉格朗日多项式。
$$
\left\{
\begin{array}{a}
f(x)=P_N(x)+E_N(x)\\
P_N(x)=\sum_{k=0}^{N}y_kL_{N,k}(x)\quad while \quad L_{N,k}(x)=\prod_{j=0, j \neq k}^{N} \frac{x-x_j}{x_k-x_j}\\
E_N(x)=\frac{f^{(N+1)}(c)}{(N+1)!} \prod^N_{j=0}{(x-x_j)}\quad c \in (x,x_0)or(x_0,x)\\
\end{array}
\right.
$$

###4.4.牛顿多项式

利用拉格朗日多项式的缺点就是如果加入一个新点，整个插值系数都要发生变动，计算量还是较大。因此设计产生了牛顿多项式，对于新加入一个点产生的新多项式$P_{N+1}(x)$可以由原来的多项式$P_N(x)$推导，不需要重新计算原来的系数。
$$
\begin{array}{a}
P_1(x)=a_0+a_1(x-x_0)\\
P_2(x)=a_0+a_1(x-x_0)+a_2(x-x_0)(x-x_1)\\
\dots \dots \\
P_N(x)=a_0+a_1(x-x_0)+a_2(x-x_0)(x-x_1)+\dots+a_N\prod_{j=0}^{N-1}(x-x_j) \\
\rightarrow P_{N+1}=P_N(x)+a_N\prod_{j=0}^{N-1}(x-x_j) \\
\end{array}
$$
用差商法推导系数
$$
\left\{
\begin{array}{a}
f(x) \approx P_N(x)\\
f(x_0)=P_N(x_0)=a_0\\
f(x_1)=P_N(x_1)=a_0+a_1(x_1-x_0)\\
f(x_2)=P_N(x_2)=a_0+a_1(x_2-x_0)+a_2(x_2-x_0)(x_2-x_1)\\
\dots\\
\end{array}
\right.
\\ \rightarrow
\left\{
\begin{array}{a}
a_0=f(x_0)\rightarrow f[x_0]\\
a_1=\frac{f(x_1)-f(x_0)}{x_1-x_0}\rightarrow f[x_0,x_1]\\
a_2=\frac{\frac{f(x_2)-f(x_1)}{x_2-x_1}-\frac{f(x_1)-f(x_0)}{x_1-x_0}}{x_2-x_0}\rightarrow f[x_0,x_1,x_2]\\
\dots \\
f[x_{k-j},x_{k-j+1},\dots,x_k]=\frac{f[x_{k-j+1},x_{k-j+2},\dots,x_k]-f[x_{k-j},x_{k-j+1},\dots,x_{k-1}]}{x_k-x_{k-j}}
\end{array}
\right.
$$
利用差商表计算比较方便。

| $x_k$ | $f[x_k]=f(x_k)$ | $f[x_{k-1},x_k]=\frac{f[x_k]-f[x_{k-1}]}{x_k-x_{k-1}}$ | $f[,,]$                                                    | $f[,,,]$ |      |
| ----- | --------------- | ------------------------------------------------------ | ---------------------------------------------------------- | -------- | ---- |
| $x_0$ | $a_0=f[x_0]$    |                                                        |                                                            |          |      |
| $x_1$ | $f[x_1]$        | $a_1=f[x_0,x_1]=\frac{f[x_1]-f[x_0]}{x_1-x_0}$         |                                                            |          |      |
| $x_2$ | $f[x_2]$        | $f[x_1,x_2]=\frac{f[x_2]-f[x_1]}{x_2-x_1}$             | $a_2=f[x_0,x_1,x_2]=\frac{f[x_1,x_2]-f[x_0,x_1]}{x_2-x_0}$ |          |      |
| $x_3$ | $f[x_3]$        | $f[x_2,x_3]$                                           | $f[x_1,x_2,x_3]=\frac{f[x_2,x_3]-f[x_1,x_2]}{x_3-x_1}$     | $a_3$    |      |
| $x_4$ | $f[x_4]$        | $f[x_3,x_4]$                                           | $f[x_2,x_3,x_4]$                                           | $\dots$  |      |

牛顿多项式的$E_N(x)$与拉格朗日多项式相同。
$$
E_N(x)=\frac{f^{(N+1)}(c)}{(N+1)!} \prod^N_{j=0}{(x-x_j)}\quad c \in (x,x_0)or(x_0,x)\\
$$

------

##5.曲线拟合

### 5.1.最小二乘法

不同于插值，拟合不需要精确通过已知点，目的是找到使误差$E$最小的拟合函数系数。

拟合函数的误差可以用最小二乘法定义为：
$$
E=[\frac{1}{N}\sum_{n=0}^N(x_n-f(x_n))^2]^{\frac{1}{2}}
$$
或者简单地定义为：
$$
E=\sum_{n=0}^N(x_n-f(x_n))^2
$$
以线性拟合$y=Ax+B$为例介绍正归方程：
$$
\begin{array}{a}
(\sum_{k=1}^Nx_k^2)A+(\sum_{k=1}^Nx_k)B=\sum_{k=1}^Nx_ky_k\\
(\sum_{k=1}^Nx_k)A+NB=\sum_{k=1}^Ny_k\\
\end{array}
$$

### 5.2.三次样条插值

三次样条插值是插值的一种方法，对于点$(x_0,y_0)(x_1,y_1)\cdots(x_N,y_N)$一共$N+1$个点进行分段的插值。

插值定义

三次样条插值的***构造条件***可以列为：
$$
\begin{array}{a}
 1.S(x)=S_k(x)=s_{k,0}+s_{k,1}(x-x_k)+s_{k,2}{(x-x_k)}^2+s_{k,3}{(x-x_k)}^3
& x \in [x_k,x_{k+1}]\\& k=0,1,\dots,N-1 \\
2. S(x_k)=y_k & k=0,1,\dots,N \\
3.  S_k(x_{k+1})=S_{k+1}(x_{k+1}) & k=0,1,\dots,N-2  \\
4. S'_k(x_{k+1})=S'_{k+1}(x_{k+1}) & k=0,1,\dots,N-2 \\
5.  S''_k(x_{k+1})=S''_{k+1}(x_{k+1}) & k=0,1,\dots,N-2 \\
\end{array}
$$
条件1表示插值多项式为分段函数，条件2表示函数穿过插值点，条件3说明各分段间连续，与条件2 不同的是，条件3强调的是各分段函数首尾相连，条件4说明各分段间光滑，条件5定义函数上各点二阶导数也是连续的，存在曲率半径。

三次样条插值的***系数推导***如下：

根据条件5，各段的$S''_k(x)$为直线，利用拉格朗日插值：
$$
S''_k(x)=S''(x_k)\frac{x-x_{k+1}}{x_k-x_{k+1}}+S''(x_{k+1})\frac{x-x_k}{x_{k+1}-x_k}
$$
记$S''(x_k)=m_k,x_{k+1}-x_k=h_k$ :
$$
S''_k(x)=\frac{m_k}{h_k}(x_{k+1}-x)+\frac{m_{k+1}}{h_k}(x-x_k)
$$
做两次积分：
$$
S_k(x)=\frac{m_k}{6h_k}(x_{k+1}-x)^3+\frac{m_{k+1}}{6h_k}(x-x_k)^3+p_k(x_{k+1}-x)+q_k(x-x_k)
$$
带入$(x_k,y_k)(x_{k+1},y_{k+1})$得求得：
$$
S_k(x)=\frac{m_k}{6h_k}(x_{k+1}-x)^3+\frac{m_{k+1}}{6h_k}(x-x_k)^3+(\frac{y_k}{h_k}-\frac{m_kh_k}{6})(x_{k+1}-x)+(\frac{y_{k+1}}{h_k}-\frac{m_{k+1}h_k}{6})(x-x_k)
$$
求导：
$$
S'_k(x)=-\frac{m_k}{2h_k}(x_{k+1}-x)^2+\frac{m_{k+1}}{2h_k}(x-x_k)^2-(\frac{y_k}{h_k}-\frac{m_kh_k}{6})+(\frac{y_{k+1}}{h_k}-\frac{m_{k+1}h_k}{6})
$$
记$d_k=\frac{y_{k+1}-y_k}{h_k}$，用$S'_k(x),S'_{k-1}(x)$计算$S'(x_k)$：
$$
\begin{array}{a}
S'_k(x_k)=-\frac{m_k}{3}h_k-\frac{m_{k+1}}{6}h_k+d_k\\
S'_{k-1}(x_k)=\frac{m_k}{3}h_{k-1}+\frac{m_{k-1}}{6}h_{k-1}+d_{k-1}\\
\end{array}
$$
根据条件4，$S'_k(x_k)=S'_{k-1}(x_k)$可以合并两式，得到：
$$
h_{k-1}m_{k-1}+2(h_{k-1}+h_k)m_k+h_km_{k+1}=u_k
$$
其中$u_k=6(d_k-d_{k-1})$

推导***总结***：
$$
\left\{
\begin{array}{a}
m_k=S''(x_k)&unknown\\
h_k=x_{k+1}-x_k &known\\
d_k=\frac{y_{k+1}-y_k}{h_k}&known\\
u_k=6(d_k-d_{k-1})& known\\
\\
S'_k(x_k)=-\frac{m_k}{3}h_k-\frac{m_{k+1}}{6}h_k+d_k\\
S'_{k-1}(x_k)=\frac{m_k}{3}h_{k-1}+\frac{m_{k-1}}{6}h_{k-1}+d_{k-1}\\
h_{k-1}m_{k-1}+2(h_{k-1}+h_k)m_k+h_km_{k+1}=u_k\\
\end{array}
\right.
\qquad \rightarrow
m_k=S''(x_k)
$$
***三次压紧样条插值***：

已知端点约束条件$S'(x_0),S'(x_N)$
$$
\begin{array}{a}
\left\{
\begin{array}{a}
h_k=x_{k+1}-x_k \\
\rightarrow d_k=\frac{y_{k+1}-y_k}{h_k}\\
\rightarrow  u_k=6(d_k-d_{k-1})\\
\end{array}\right.
\\
\\
 \rightarrow
\left\{
\begin{array}{a}
S'(x_0)=-\frac{m_0}{3}h_0-\frac{m_{1}}{6}h_0+d_0\\
S'(x_N)=\frac{m_N}{3}h_{N-1}+\frac{m_{N-1}}{6}h_{N-1}+d_{N-1}\\
h_{k-1}m_{k-1}+2(h_{k-1}+h_k)m_k+h_km_{k+1}=u_k\\
\end{array}
\right.
\\
\\  \rightarrow
m_k=S''(x_k)\quad k=0,1,\cdots,N\\
\\  \rightarrow

\left\{
\begin{array}{a}
s_{k,0}=y_k\\
s_{k,1}=d_k-\frac{h_k(2m_k+m_{k+1})}{6}\\
s_{k,2}=\frac{m_k}{2}\\
s_{k,3}=\frac{m_{k+1}-m_{k}}{6h_k}
\end{array}\right.
\\
\end{array}
$$

***四点的三次压紧样条插值***：

$$
\begin{array}{a}

\left\{
\begin{array}{a}
(x_0,y_0)(x_1,y_1)(x_2,y_2)(x_3,y_3)\\
S'(x_0),S'(x_3)\\
\end{array}
\right.
\\
\rightarrow
\left\{
\begin{array}{a}
h_0=x_1-x_0&h_1=x_2-x_1&h_2=x_3-x_0\\
d_0=\frac{y_1-y_0}{h_0}&d_1=\frac{y_2-y_1}{h_1}&d_2=\frac{y_3-y_2}{h_2}\\
u_1=6(d_1-d_0)&u_2=6(d_2-d_1)\\
\end{array}
\right.

\\
\rightarrow
\left\{
\begin{array}{a}
S'(x_0)=-\frac{m_0}{3}h_0-\frac{m_{1}}{6}h_0+d_0\\
S'(x_N)=\frac{m_N}{3}h_{N-1}+\frac{m_{N-1}}{6}h_{N-1}+d_{N-1}\\
h_0m_0+2(h_0+h_1)m_1+h_1m_2=u_1\\
h_1m_1+2(h_1+h_2)m_2+h_2m_3=u_2\\
\end{array}
\right.
\\

\rightarrow
\left\{
\begin{array}{a}
m_0,m_1,m_2,m_3 \rightarrow \\
S_k(x)=\frac{m_k}{6h_k}(x_{k+1}-x)^3+\frac{m_{k+1}}{6h_k}(x-x_k)^3+(\frac{y_k}{h_k}-\frac{m_kh_k}{6})(x_{k+1}-x)+(\frac{y_{k+1}}{h_k}-\frac{m_{k+1}h_k}{6})(x-x_k)
\end{array}
\right.
\\
\end{array}
$$

### 5.3.三角多项式

若$f(x)$周期为$2\pi$且有$N+1$个均匀分布的点，则其$M$阶三角多项式为：
$$
\begin{array}{a}
T_M(x)=\frac{a_0}{2}+\sum_{j=0}^M(a_j\cos jx+b_j \sin jx)\\
x_j=-\pi +\frac{2j\pi}{N} & j=0,1,2,\cdots,N\\
a_j=\frac{2}{N}\sum_{k=1}^Nf(x_k)\cos jx_k & j=0,1,2,\cdots,M\\
b_j=\frac{2}{N}\sum_{k=1}^Nf(x_k)\sin jx_k & j=1,2,\cdots,M\\
\end{array}
$$
当$N\rightarrow \infty$时，多项式系数趋向于傅里叶级数的系数。

------

## 7.积分

### 7.1.闭型牛顿-科特斯公式

对函数$f(x)$的积分估计，就是对其插值多项式$P_N(x)$的积分运算。

由$f(x)=P_N(x)+E_N(x)$得到。
$$
\int_a^bf(x)dx \approx \int_a^bP_N(x)dx
$$
当等间隔插值点数不同时，计算的公式也不同，设$[a,b]$上均匀插值点数为$N+1$个，将积分区间等分为$N$份，如$N=2$，即插值点为$x_0=a,x_1=\frac{a+b}{2},x_2=b$。相邻两点的间隔用$h$表示，则插值点可以表示为$x_k=a+kh$，这样的积分方式统称为闭型牛顿-科特斯公式，对于$N=1,N=2$具体表达式如下。

***梯形公式***，积分区间选首尾两点作为插值点：
$$
\int_a^bf(x)dx \approx \frac{h}{2}(f(x_0)+f(x_1))
$$
***辛普森公式***，积分区间选首中尾三点作为插值点：
$$
\int_a^bf(x)dx \approx \frac{h}{3}(f(x_0)+4f(x_1)+f(x_2))
$$
辛普森$\frac{3}{8}$公式，积分区间选首中尾四点作为插值点：
$$
\int_a^bf(x)dx \approx \frac{3h}{8}(f(x_0)+3f(x_1)+3f(x_2)+f(x_3))
$$
布尔公式，积分区间选首中尾五点作为插值点：
$$
\int_a^bf(x)dx \approx \frac{2h}{45}(7f(x_0)+32f(x_1)+12f(x_2)+32f(x_3)+7f(x_4))
$$
当一个存在多个积分区间时，在每个积分区间上使用梯形公式或者辛普森公式，再求和，则成为组合公式。

组合梯形公式，先等分积分区间为$M$份，每一份小的区间两个插值点$x_0,x_1​$，使用梯形公式：
$$
T(f,h)=\int_a^bf(x)dx \approx \sum^M \frac{h}{2}(f(x_0)+f(x_1))
$$
组合辛普森公式，先等分积分区间为$M$份，每一份小的区间三个插值点$x_0,x_1,x_2$，使用辛普森公式：
$$
S(f,h)=\int_a^bf(x)dx \approx \sum^M \frac{h}{3}(f(x_0)+4f(x_1)+f(x_2))
$$
用$T(J)=T(f,h)$表示在积分区间$[a,b]$上进行$2^J$次组合梯形公式，共需要插值点数$2^J+1$个，且$h=\frac{b-a}{2^J}$。。同理用$S(J)=S(f,h)$表示在积分区间$[a,b]$上进行$2^J$次组合辛普森公式，共需要插值点数$2^{J+1}+1$个，且$h=\frac{b-a}{2*2^J}$。

误差分析梯形公式精度为1，辛普森公式和辛普森$\frac{3}{8}$公式精度为3，布尔公式精度为5，精度表示积分函数$f(x)$的泰勒展开阶数小于等于精度，积分不存在误差。如梯形公式求线性一次函数的积分，就不存在误差。

组合积分公式的误差为每个小区间误差之和。
$$
\begin{array}{a}
E_T(f,h)=\frac{-(b-a)f^{(2)}(c)h^2}{12}=O(h^2)\\
E_S(f,h)=\frac{-(b-a)f^{(4)}(c)h^4}{180}=O(h^4)\\
\end{array}
$$

### 7.2.龙贝格积分

用$T(J)$表示了划分区间为$2^J$个的组合梯形积分。$T(J)$内部存在迭代关系。
$$
T(J)=\frac{T(J-1)}{2}+h\sum_{k=1}^M f(x_{2k-1})
$$
$x_{2k-1}$表示$T(J-1)$到$T(J)$新插入的$2^{J-1}$个点。而$S(J)$可以用$T(J)$表示：
$$
S(J)=\frac{4T(J)-T(J-1)}{4-1}
$$
遵循规律构建龙贝格积分
$$
\begin{array}{a}
R(J,0)=T(J)=\frac{T(J-1)}{2}+h\sum_{k=1}^M f(x_{2k-1})\\
R(J,1)=S(J)=\frac{4T(J)-T(J-1)}{4-1}\\
R(J,2)=B(J)\\
\cdots\\
R(J,K)=\frac{4^KR(J,K-1)-R(J-1,K-1)}{4^K-1}\\
\end{array}
$$

$R(J,K)$的误差为$O(h^{2K+2})$。

-------

## 9.微分方程

### 9.1.欧拉方法

利用导数的性质。


$$
\begin{array}{a}
t_k=a+kh\\
f(t_k,y_k)=y'(t)|_{t=t_k}\\
y_{k+1}=y_k+hf(t_k,y_k)\\
\end{array}
$$

### 9.2.休恩方法

将欧拉方法和梯形公式组合。
$$
\begin{array}{a}
t_k=a+kh\\
f(t_k,y_k)=y'(t)|_{t=t_k}\\
p_{k+1}=y_k+hf(t_k,y_k)\\
y_{k+1}=y_k+\frac{h}{2}[f(t_k,y_k)+f(t_{k+1},p_{k+1})]\\
\end{array}
$$

### 9.3.龙格-库塔方法











