#数值分析作业

------
[TOC]

------

##work1.3.10-1		求根公式

```matlab
function work1_3_10_1()
a=1;b=-1000.001;c=1;%ax^2+bx+c=0;
[root1,root2]=root(a,b,c);
disp(root1);
disp(root2);
end

function [root1,root2]=root(a,b,c)
delta=sqrt(b.^2-4.*a.*c);
if b<=0
    root1=(-b+delta)./(2.*a);
    root2=(-2.*c)./(b-delta);
else if b>0
        root1=(-2.*c)./(b+delta);
        root2=(-b-delta)./(2.*a);
    end
end 
end
```

- 运行结果

```matlab
>> work1_3_10_1
root1 =
        1000
root2 =
   1.0000e-03
```

- 结果分析

考虑到$b$与$\delta=\sqrt{b^2-4ac}$ 做差会导致精度损失，因此求根公式在选取的时候要保证运算中的精度尽量不损失。
$$
if \quad b \leq 0 \rightarrow
\left\{ 
\begin{array}{a}
root_1=\frac{-b+\delta}{2a}\\
root_2=\frac{-2c}{b-\delta}\\
\end{array}
\right.
or\quad
if \quad b > 0 \rightarrow
\left\{ 
\begin{array}{a}
root_1=\frac{-2c}{b+\delta}\\
root_2=\frac{-b-\delta}{2a}\\
\end{array}
\right.
$$

------

## work2.1.5-1		迭代法

```matlab
function work2_1_5_1()
    close all;
    ga=@(x)(x.^5-3.*x.^3-2.*x.^2+2);
    gb=@(x)(cos(sin(x)));
    gc=@(x)(x.^2-sin(x+0.15));
    gd=@(x)(x.^(x-cos(x)));
    x=-3:0.1:3;
    figure;plot(x,feval(ga,x)); hold on;plot(x,x);grid on;
    figure;plot(x,feval(gb,x)); hold on;plot(x,x);grid on;
    figure;plot(x,feval(gc,x)); hold on;plot(x,x);grid on;
    figure;plot(x,feval(gd,x)); hold on;plot(x,x);grid on;
    max=1000;
    tol=10^-12;   
    pa0=1;%
    pb0=0;%
    pc0=0;%
    pd0=1;%
    [k1,p1,err1,P1]=fixpt(ga,pa0,tol,max);
    [k2,p2,err2,P2]=fixpt(gb,pb0,tol,max);
    [k3,p3,err3,P3]=fixpt(gc,pc0,tol,max);
    [k4,p4,err4,P4]=fixpt(gd,pd0,tol,max);
    %p2=round(p2.*10^12)./(10^12);
    p1,p2,p3,p4
    
end

function [k,p,err,P]=fixpt(g,p0,tol,max)
P(1)=p0;
for k=2:max
	 P(k)=feval(g,P(k-1));
   	 err=abs(P(k)-P(k-1));
   	 relerr=err/abs(P(k)+eps);
    p=P(k);
    if(err<tol)||(relerr<tol),break;end
end
    if k==max
    disp('maximun number of iterations exceeded')
    end 
    P=P';
end
```

- 运行结果

```matlab
work2_1_5_1
Warning: Imaginary parts of complex X and/or Y arguments ignored 
>> In work2_1_5_1 (line 11) 
maximun number of iterations exceeded
p1 =
     2
p2 =
    0.7682
p3 =
   -0.3513
p4 =
     1
```

- 结果分析

* 函数$g_a(x)$ 

设初始迭代值$p_0=0$，最终收敛于$x=2$ ，修改迭代值，发现只能收敛于$x=2$ 。

![2_1_5_1_a](C:\Users\Chemizi\Desktop\数值分析\2_1_5_1_a.png)

* 函数$g_b(x)$ 

由图案的观察可以发现只有一个交点，因此得到的迭代结果$x=0.7682$ 就是最终的不动点。

![2_1_5_1_b](C:\Users\Chemizi\Desktop\数值分析\2_1_5_1_b.png)

* 函数$g_c(x)$ 

由图可以发现，$x=-0.3513$ 并不是我们希望得到的不动点，进一步查看迭代过程，我们可以发现迭代点在$x=-0.3513$ 和$x=0.3235$ 之间跳变。

![2_1_5_1_c](C:\Users\Chemizi\Desktop\数值分析\2_1_5_1_c.png)

* 函数$g_d(x)$ 

$x=1$ 是函数的一个不动点，无法收敛到其他不动点。

![2_1_5_1_d](C:\Users\Chemizi\Desktop\数值分析\2_1_5_1_d.png)

------

## work2.2.4-2		试值法

```matlab
function work2_2_4_2()
    format long
    r=0.15;
    rho=0.710;   
    f=@(d)(pi.*(d.^3-3.*d.^2*r+4.*r.^3.*rho)/3);
    x=0:0.01:2*r;
    y=feval(f,x);
    plot(x,y);
    grid on
    a=0;
    b=0.3;
    delta=10^-8;
    epsilon=0;
    max1=10000;
    [c,err,yc]=regula(f,a,b,delta,epsilon,max1);
    c
    yc
    err
end

function [c,err,yc]=regula(f,a,b,delta,epsilon,max1)
    ya=feval(f,a);
    yb=feval(f,b);
    if ya*yb>0
        disp('Note:f(a)*f(b)>0');
        return;
    end
    for k=1:max1
        dx=yb*(b-a)/(yb-ya);
        c=b-dx;
        ac=c-a;
        yc=feval(f,c);
        if yc==0,break;
        elseif yb*yc>0
            b=c;
            yb=yc;
        else
            a=c;
            ya=yc;
        end
        dx=min(abs(dx),ac);
        if abs(dx)<delta,break,end
        if abs(yc)<epsilon,break,end
    end
    c;
    err=abs(b-a)/2;
    yc=feval(f,c);

end
```

- 运行结果

```matlab
>> work2_2_4_2
c =
   0.193193886651364
yc =
    -3.195041955020870e-12
err =
   0.001660234464923
```

- 结果分析

取中点的速度太慢，因此考虑用弦与横轴的焦点代替中点。弦的斜率为$m=\frac{f(b)-f(a)}{b-a}$，弦与横轴交点为$c=b-\frac{f(b)(b-a)}{f(b)-f(a)}$，$c$的替换思路与中值法相同。试值法最终收敛。

根据题目条件和阿基米德原理，列出方程$\frac{\pi(d^3-3d^2r+4r^3\rho)}{3}=0$，设置方程解的容许误差为$\delta=10^{-8}$，当迭代的前后两次结果之差小于$\delta$时，就可以认为，前八位小数有有效数字。使用试值法求得有效半径为$d=0.19319388m$ 。

------

## work2.3.4-1		近似根估计

```matlab
function work2_3_4_1()
    f=@(x)(1000000.*x.^3-111000.*x.^2+1110.*x-1);
    X=-2:0.001:2;%plot(X,feval(f,X));grid on;
    epsilon1=10.^-12;
    R=approot(f,X,epsilon1);
    R
    delta=10^-12;
    epsilon=0;
    max1=1000;
    [c1,err,yc]=regula(f,0.001,0.0015,delta,epsilon,max1);
    [c2,err,yc]=regula(f,0.0015,0.0095,delta,epsilon,max1);
    [c3,err,yc]=regula(f,0.0095,0.01,delta,epsilon,max1);
    [c4,err,yc]=regula(f,0.01,0.1,delta,epsilon,max1);
    [c5,err,yc]=regula(f,0.1,0.1005,delta,epsilon,max1);
    c1,c2,c3,c4,c5
end

function R=approot(f,X,epsilon)
    Y=feval(f,X);
    yrange=max(Y)-min(Y);
    epsilon2=yrange*epsilon;
    n=length(X);
    m=0;
    X(n+1)=X(n);
    Y(n+1)=Y(n);
    for k=2:n
        if Y(k-1)*Y(k)<=0
            m=m+1;
            R(m)=(X(k-1)+X(k))/2;
        end
        s=(Y(k)-Y(k-1))*(Y(k+1)*Y(k));
        if (abs(Y(k))<epsilon2)&&(s<=0),
            m=m+1;
            R(m)=X(k);
        end
    end
end

function [c,err,yc]=regula(f,a,b,delta,epsilon,max1)
    c=NaN;err=NaN;yc=NaN;
    ya=feval(f,a);
    yb=feval(f,b);
    if ya*yb>0
        disp('Note:f(a)*f(b)>0'),return,end
    for k=1:max1
        dx=yb*(b-a)/(yb-ya);
        c=b-dx;
        ac=c-a;
        yc=feval(f,c);
        if yc==0,break;
        elseif yb*yc>0
            b=c;
            yb=yc;
        else
            a=c;
            ya=yc;
        end
        dx=min(abs(dx),ac);
        if abs(dx)<delta,break,end
        if abs(yc)<epsilon,break,end
    end
    c;
    err=abs(b-a)/2;
    yc=feval(f,c);
end
```

- 运行结果

```matlab
>> work2_3_4_1
R =
   0.001000000000000   0.001500000000000   0.009500000000000   0.010000000000000   0.100000000000000   0.100500000000000
Note:f(a)*f(b)>0
c1 =
     1.000000000000000e-03
c2 =
   NaN
c3 =
   0.010000000000000
c4 =
   NaN
c5 =
   0.100000000000000
```

函数$f(x)$

![2_3_4_1](C:\Users\Chemizi\Desktop\数值分析\2_3_4_1.png)

- 结果分析

首先求解根的近似位置，得到4位小数的近似根$$R =[0.0010,0.0015,0.0095,0.0100,0.1000,0.1005]$$ ，

这6个点划分出5个区间，在5个区间中找到可能存在的唯一零点即可。最终三个解。

$$x_1=0.001,x_2=0.01,x_3=0.1$$

------

## work2.4.8-2		牛顿-拉夫森法

```matlab
function work2_4_8_2()
format long
f=inline('4800*(1-exp(-t/10))-320*t','t');
df=inline('480*exp(-t/10)-320','t');
p0=8;
delta=10.^-11;
epsilon=0;
max1=4;
result=newton(f,df,p0,delta,epsilon,max1);
%result=num2str(result,'%20.9f');
disp(result);
end

function result=newton(f,df,p0,delta,epsilon,max1)
result=zeros(max1+1,4);
index=1;
for k=1:max1+1
    p1=p0-feval(f,p0)./feval(df,p0);
    err=abs(p1-p0);
    relerr=2*err/(abs(p1)+delta);
    result(index,1)=index-1;
    result(index,2)=p0;
    result(index,3)=p1-p0;
    result(index,4)=feval(f,p0);
    index=index+1;
    y=feval(f,p1);
    p0=p1;
    if(err<delta)||(relerr<delta)||(abs(y)<epsilon)
        break;
    end
end
end
```

- 运行结果

```matlab
>> work2_4_8_2
                   0   8.000000000000000   0.797731012432170  83.220972237336355
   1.000000000000000   8.797731012432170  -0.055301598883606  -6.683697193629087
   2.000000000000000   8.742429413548564  -0.000254750135253  -0.030507523597407
   3.000000000000000   8.742174663413310  -0.000000005426140  -0.000000649778030
   4.000000000000000   8.742174657987171                   0                   0
```

- 结果分析

将区间是为一个点$p_0$，$p_0$处切线点斜率为$m=f’(p_0)$，切线交横坐标于$p_1=p_0-\frac{p_0}{f’(p_0)}$，通过迭代最终$r=p_\infty$，此方法称为牛顿-拉夫森法，迭代公式为
$$
p_{n+1}=p_n-\frac{f(p_n)}{f’(p_n)}
$$
其收敛的充分条件为：
$$
\frac{|f(x)f’’(x)|}{|f’(x)|^2}<1
$$

| $k$  | $p_k$             | $p_{k+1}-p_k$      | $f(p_k)=4800(1-e^{-\frac{p_k}{10}})-320p_k$ |
| ---- | ----------------- | ------------------ | ------------------------------------------- |
| 0    | 8.000000000000000 | 0.797731012432170  | 83.220972237336355                          |
| 1    | 8.797731012432170 | -0.055301598883606 | -6.683697193629087                          |
| 2    | 8.742429413548564 | -0.000254750135253 | -0.030507523597407                          |
| 3    | 8.742174663413310 | -0.000000005426140 | -0.000000649778030                          |
| 4    | 8.742174657987171 | 0.000000000000000  | 0.000000000000000                           |

四次迭代之后，$f(p_4) \approx 0$，并且$t \approx p_4 =8.742174657987171 $。

------

## work3.2.8-3		三维图形旋转

```matlab
function work3_2_8_3()
close all;
clear all;
a1=[0;0;0];
b1=[1;0;0];
c1=[0;1;0];
d1=[0;0;1];
C1=[a1 b1 c1 d1 b1 a1 d1 a1 c1];
plot3(C1(1,:),C1(2,:),C1(3,:));axis equal;grid on;
xlim([-1,1]);ylim([-1,1]);zlim([-1,1]);
alpha=2.7;
beta=0.15;
gamma=-1.5;
R_x=[1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)];
R_y=[cos(beta),0,sin(beta);0,1,0;-sin(beta),0,cos(beta)];
R_z=[cos(gamma),-sin(gamma),0;sin(gamma),cos(gamma),0;0,0,1];
C2=R_y*C1;
figure(1);
hold on;
plot3(C2(1,:),C2(2,:),C2(3,:));axis equal;grid on;
xlim([-1,1]);ylim([-1,1]);zlim([-1,1]);
C3=R_z*C2;
figure(1);
hold on;
plot3(C3(1,:),C3(2,:),C3(3,:));axis equal;grid on;
xlim([-1,1]);ylim([-1,1]);zlim([-1,1]);
C4=R_x*C3;
figure(1);
hold on;
plot3(C4(1,:),C4(2,:),C4(3,:));axis equal;grid on;
xlim([-1,1]);ylim([-1,1]);zlim([-1,1]);
end
```

- 运行结果

![3_2_8_3](C:\Users\Chemizi\Desktop\数值分析\3_2_8_3.png)

- 结果分析

用0，1，2，3分别表示原始图案，第一次变换，第二次变换，第三次变换的结果。

利用的是坐标旋转的变换矩阵：
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
这三个矩阵分别是以角度$\alpha,\beta,\gamma$ 绕$x,y,z$ 轴旋转的变换矩阵。

------

## work3.3.2-1		回代法

```matlab
function work3_3_2_1()
    U=zeros(10,10);
    for i=1:10
        for j=1:10
            if i<=j
                U(i,j)=cos(i*j);
            end
        end
    end
    B=tan(1:10);
    X=backsub(U,B);
    X
end

function X=backsub(A,B)
    n=length(B);
    X=zeros(n,1);
    X(n)=B(n)/A(n,n);
    for k=n-1:-1:1
        X(k)=(B(k)-A(k,k+1:n)*X(k+1:n))/A(k,k);
    end
end
```

- 运行结果

```matlab
>> work3_3_2_1
X =
 -748.9196
 -398.4748
 -207.2299
  -83.9979
   42.8938
  -75.7439
   51.4968
  -17.5075
   -0.1486
    0.7519
```

- 结果分析

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
代码使用的是上三角矩阵回代法，从最后一行依次向上求 解。对于下三角矩阵的解法如下，此方法又称前向替代法：

```matlab
function X=forsub(A,B)
    n=length(B);
    X=zeros(n,1);
    X(1)=B(1)/A(1,1);
    for k=2:1:n
        X(k)=(B(k)-A(k,1:k-1)*X(1:k-1))/A(k,k);
    end
end
```

------

##work3.4.6-2		上三角变换与回代法

```matlab
function work3_4_6_2()
    %A=[a1;a2;a3;a4;a5;a6;a7];
    %UA=B;
    X=@(x)([1,x,x.^2,x.^3,x.^4,x.^5,x.^6]);
    U=[feval(X,0);feval(X,1);feval(X,2);feval(X,3);feval(X,4);feval(X,5);feval(X,6)];
    B=[1;3;2;1;3;2;1];
    A=uptrbk(U,B);
    A
    f=@(x)(A(1)+A(2)*x+A(3)*x.^2+A(4)*x.^3+A(5)*x.^4+A(6)*x.^5+A(7)*x.^6);
    x=0:0.01:6;
   plot(x,feval(f,x));grid on;
end
function X=uptrbk(A,B)
    [N N]=size(A);
    X=zeros(N,1);
    C=zeros(1,N+1);
    Aug=[A B];
    for p=1:N-1
        [Y,j]=max(abs(Aug(p:N,p)));
        C=Aug(p,:);
        Aug(p,:)=Aug(j+p-1,:);
        Aug(j+p-1,:)=C;
        if Aug(p,p)==0
            disp('A is singular cannot solve');
            return;
        end
        for k=p+1:N
            m=Aug(k,p)/Aug(p,p);
            Aug(k,p:N+1)=Aug(k,p:N+1)-m*Aug(p,p:N+1);
        end
    end
    X=backsub(Aug(1:N,1:N),Aug(1:N,N+1));
end

function X=backsub(A,B)
    n=length(B);
    X=zeros(n,1);
    X(n)=B(n)/A(n,n);
    for k=n-1:-1:1
        X(k)=(B(k)-A(k,k+1:n)*X(k+1:n))/A(k,k);
    end
end
```

- 运行结果

```matlab
>> work3_4_6_2
A =
    1.0000
   -1.8000
   11.0250
  -10.5625
    3.9375
   -0.6375
    0.0375
```

- 结果分析

以一个三元一次方程组说明高斯消去法
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

------

## work3.5.8-1		三角分解与回代法

```matlab
function work3_5_8_1()
    A=[1,3,5,7;2,-1,3,5;0,0,2,5;-2,-6,-3,1];
    B=[1;2;3;4];
    X=lufact(A,B);
    X
 
    [L,U,P]=lu(A);
    L,U,P
    Y=forsub(L,P*B);
    X_=backsub(U,Y);
    X_
end


function X=lufact(A,B)
    [N,N]=size(A);
    X=zeros(N,1);
    Y=zeros(N,1);
    C=zeros(1,N);
    R=1:N;
    for p=1:N-1
        [max1,j]=max(abs(A(p:N,p)));
        C=A(p,:);
        A(p,:)=A(j+p-1,:);
        A(j+p-1,:)=C;
        d=R(p);
        R(p)=R(j+p-1);
        R(j+p-1)=d;
        if A(p,p)==0
            disp('A is sigular. No unique solution');
            return;
        end
        for k=p+1:N
            mult=A(k,p)/A(p,p);
            A(k,p)=mult;
            A(k,p+1:N)=A(k,p+1:N)-mult*A(p,p+1:N);
        end
    end
    Y(1)=B(R(1));
    for k=2:N
        Y(k)=B(R(k))-A(k,1:k-1)*Y(1:k-1);
    end
    X(N)=Y(N)/A(N,N);
    for k=N-1:-1:1
        X(k)=(Y(k)-A(k,k+1:N)*X(k+1:N))/A(k,k);
    end
end

function X=backsub(A,B)
    n=length(B);
    X=zeros(n,1);
    X(n)=B(n)/A(n,n);
    for k=n-1:-1:1
        X(k)=(B(k)-A(k,k+1:n)*X(k+1:n))/A(k,k);
    end
end
function X=forsub(A,B)
    n=length(B);
    X=zeros(n,1);
    X(1)=B(1)/A(1,1);
    for k=2:1:n
        X(k)=(B(k)-A(k,1:k-1)*X(1:k-1))/A(k,k);
    end
end
```

- 运行结果

```matlab
>> work3_5_8_1
X =
    1.3429
    0.6857
   -3.0000
    1.8000
L =
    1.0000         0         0         0
   -1.0000    1.0000         0         0
    0.5000   -0.5000    1.0000         0
         0         0    0.5714    1.0000
U =
    2.0000   -1.0000    3.0000    5.0000
         0   -7.0000         0    6.0000
         0         0    3.5000    7.5000
         0         0         0    0.7143
P =
     0     1     0     0
     0     0     0     1
     1     0     0     0
     0     0     1     0
X_ =
    1.3429
    0.6857
   -3.0000
    1.8000
```

- 结果分析

对于线性方程$AX=B$，可以将系数矩阵$A$直接或者间接地分解为上下两个三角矩阵之积，即$PA=LU$ ，其中$P$为一个置换矩阵，其各行各列有且只有一个元素为1，其他元素为0。则方程变换为：
$$
AX=B\rightarrow PAX=PB\rightarrow LUX=PB\\
difine:UX=Y
\rightarrow
\left\{ 
\begin{array}{a}
LY=PB \rightarrow Y=L^{-1}PB\\
UX=Y \rightarrow X=U^{-1}Y\\
\end{array}
\right.
$$
其中，乘以逆矩阵的方式只是数学公式上的表示，实际就是利用回代法解决。

------

## work3.6.5-1		雅可比迭代与高斯-赛德尔迭代

```matlab
function work3_6_5_1()
    format long 
    delta=10^-9;max1=1000;
    %1
    A1=[4,-1;1,5];B1=[15;9];P=[0;0];
    X1_jacobi=jacobi(A1,B1,P,delta,max1);
    X1_gseid=gseid(A1,B1,P,delta,max1);
    X1_jacobi,X1_gseid
    %2
    A2=[8,-3;-1,4];B2=[10;6];P=[0;0];
    X2_jacobi=jacobi(A2,B2,P,delta,max1);
    X2_gseid=gseid(A2,B2,P,delta,max1);
    X2_jacobi,X2_gseid
    %3
    A3=[-1,3;6,-2];B3=[1;2];P=[0;0];
    X3_jacobi=jacobi(A3,B3,P,delta,max1);
    X3_gseid=gseid(A3,B3,P,delta,max1);
    X3_jacobi,X3_gseid
    %4
    A4=[2,3;7,-2];B4=[1;1];P=[0;0];
    X4_jacobi=jacobi(A4,B4,P,delta,max1);
    X4_gseid=gseid(A4,B4,P,delta,max1);
    X4_jacobi,X4_gseid
    %5
    A5=[5,-1,1;2,8,-1;-1,1,4];B5=[10;11;3];P=[0;0;0];
    X5_jacobi=jacobi(A5,B5,P,delta,max1);
    X5_gseid=gseid(A5,B5,P,delta,max1);
    X5_jacobi,X5_gseid
    %6
    A6=[2,8,-1;5,-1,1;-1,1,4];B6=[11;10;3];P=[0;0;0];
    X6_jacobi=jacobi(A6,B6,P,delta,max1);
    X6_gseid=gseid(A6,B6,P,delta,max1);
    X6_jacobi,X6_gseid
    %7
    A7=[1,-5,-1;4,1,-1;2,-1,-6];B7=[-8;13;-2];P=[0;0;0];
    X7_jacobi=jacobi(A7,B7,P,delta,max1);
    X7_gseid=gseid(A7,B7,P,delta,max1);
    X7_jacobi,X7_gseid
    %8
    A8=[4,1,-1;1,-5,-1;2,-1,-6];B8=[13;-8;-2];P=[0;0;0];
    X8_jacobi=jacobi(A8,B8,P,delta,max1);
    X8_gseid=gseid(A8,B8,P,delta,max1);
    X8_jacobi,X8_gseid
end


function X=jacobi(A,B,P,delta,max1)
    N=length(B);
    X=zeros(1,N);
    for k=1:max1
        for j=1:N
            X(j)=(B(j)-A(j,[1:j-1,j+1:N])*P([1:j-1,j+1:N]))/A(j,j);
        end
        err=abs(norm(X'-P));
        relerr=err/(norm(X)+eps);
        P=X';
        if (err<delta)||(relerr<delta)
            return;
        end
        if(sum(isinf(X)))
            return;
        end
    end
    X=X';
end


function X=gseid(A,B,P,delta,max1)
    N=length(B);
    for k=1:max1
        for j=1:N
            if j==1
                X(1)=(B(1)-A(1,2:N)*P(2:N))/A(1,1);
            elseif j==N
                X(N)=(B(N)-A(N,1:N-1)*(X(1:N-1))')/A(N,N);
            else
                X(j)=(B(j)-A(j,1:j-1)*X(1:j-1)'-A(j,j+1:N)*P(j+1:N))/A(j,j);
            end
        end
    err=abs(norm(X'-P));
    relerr=err/(norm(X)+eps);
    P=X';
    if (err<delta)||(relerr<delta)
        return;
    end
    if(sum(isinf(X)))
        return;
    end
    end
    X=X';
end
```

- 运行结果

```matlab
>> work3_6_5_1
X1_jacobi =
   4.000000000195312   0.999999999375000
X1_gseid =
   3.999999999990235   1.000000000001953
X2_jacobi =
   1.999999999580432   1.999999999720288
X2_gseid =
   1.999999999960666   1.999999999990166
X3_jacobi =
  -Inf  -Inf
X3_gseid =
  -Inf  -Inf
X4_jacobi =
  1.0e+307 *
   5.073446820241841                -Inf
X4_gseid =
  1.0e+307 *
   5.073446820241841                 Inf
X5_jacobi =
   1.999999999370842   1.000000000327143   1.000000000225465
X5_gseid =
   2.000000000155297   0.999999999920040   1.000000000058814
X6_jacobi =
   Inf   Inf   Inf
X6_gseid =
   Inf   Inf   NaN
X7_jacobi =
   Inf  -Inf   Inf
X7_gseid =
  1.0e+308 *
   1.051782244356424                -Inf                 Inf
X8_jacobi =
   2.999999999384407   2.000000000054841   0.999999999300220
X8_gseid =
   2.999999999801222   2.000000000048810   0.999999999925606
```

- 结果分析

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

## work4.1.3-1		泰勒级数

```matlab
function work4_1_3_1()
    format long 
    x=linspace(-1,1,10);
%     x=-1:1;
    y=sin(x);
    P_5= x-x.^3/factorial(3)+x.^5/factorial(5);
    P_7= x-x.^3/factorial(3)+x.^5/factorial(5)-x.^7/factorial(7);
    P_9= x-x.^3/factorial(3)+x.^5/factorial(5)-x.^7/factorial(7)+x.^9/factorial(9);
    plot(x,y);grid on;
    hold on;plot(x,P_5);
    hold on;plot(x,P_7);
    hold on;plot(x,P_9);
    legend('y','P_5','P_7','P_9');
    D=[x' y' P_5' P_7' P_9'];
    D
end
```

- 运行结果

```matlab
>> work4_1_3_1
D =
  -1.000000000000000  -0.841470984807897  -0.841666666666667  -0.841468253968254  -0.841471009700176
  -0.777777777777778  -0.701697876146735  -0.701731753854144  -0.701697590682923  -0.701697877719170
  -0.555555555555556  -0.527415385771866  -0.527418612790507  -0.527415371918140  -0.527415385810769
  -0.333333333333333  -0.327194696796152  -0.327194787379973  -0.327194696656288  -0.327194696796294
  -0.111111111111111  -0.110882628509953  -0.110882628551429  -0.110882628509946  -0.110882628509953
   0.111111111111111   0.110882628509953   0.110882628551429   0.110882628509946   0.110882628509953
   0.333333333333333   0.327194696796152   0.327194787379973   0.327194696656288   0.327194696796294
   0.555555555555556   0.527415385771866   0.527418612790507   0.527415371918140   0.527415385810769
   0.777777777777778   0.701697876146735   0.701731753854144   0.701697590682923   0.701697877719170
   1.000000000000000   0.841470984807897   0.841666666666667   0.841468253968254   0.841471009700176
```

- 结果分析

四个泰勒多项式的图像基本重合。

![4_1_3_1_a](C:\Users\Chemizi\Desktop\数值分析\4_1_3_1_a.png)

放大细节可以看到，泰勒多项式次数越高，误差越小。

![4_1_3_1_b](C:\Users\Chemizi\Desktop\数值分析\4_1_3_1_b.png)

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

------

## work4.2.2-1		霍纳方法

```matlab
function work4_2_2_1()
    P=[1,2,4];%y=x^2+2x+4
    y=horner(P,2);
    y_deri=horner_deri(P,2);
    y_inte=horner_inte(P,2);
    
    y,y_deri,y_inte
end
function y=horner(P,c)
    n=length(P);
    a=P;
    b(1)=a(1);
    for k=2:n;
        b(k)=a(k)+c*b(k-1);
    end
    y=b(n);
end

function y=horner_deri(P,c)
    n=length(P);
    a=P.*[n-1:-1:0];
    a=a(1:end-1);
    n=n-1;
    b(1)=a(1);
    for k=2:n;
        b(k)=a(k)+c*b(k-1);
    end
    y=b(n);
end

function y=horner_inte(P,c)
    n=length(P);
    a=P./[n:-1:1];
    a=[a,0];
    n=n+1;
    b(1)=a(1);
    for k=2:n;
        b(k)=a(k)+c*b(k-1);
    end
    y=b(n);
end
```

运行结果

```matlab
>> work4_2_2_1
y =
    12
y_deri =
     6
y_inte =
  14.666666666666666
```

- 结果分析

使用的多项式为$y=x^2+2x+4$ ，利用霍纳方法可以改写为$y=x(1x+2)+4$ ， 至于$y'=2x+2$ 和$\int y dx=\frac{1}{3}x^3+x^2+4x=((\frac{1}{3}x+1)x+4)x$ ，运算方法同理。

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

------

## work4.3.5-2		拉格朗日插值法

```matlab
function work4_3_5_2()
    X=1:6;
    Y=[66,66,65,64,63,63];
    [C,L]=lagran(X,Y);%C=[Cn,Cn-1,...,C0];
    C,L
    x=1:0.1:6;
    y=polyval(C,x);
    close all;
    plot(x,y);grid on;
    hold on;scatter(X,Y);
    ave=(horner_inte(C,6)-horner_inte(C,1))/5;
    ave
end

function [C,L]=lagran(X,Y)
    w=length(X);
    n=w-1;
    L=zeros(w,w);
    for k=1:n+1
        V=1;
        for j=1:n+1
            if k~=j
                V=conv(V,poly(X(j)))/(X(k)-X(j));
            end
        end
        L(k,:)=V;
    end
    C=Y*L;
end

function y=horner_inte(P,c)
    n=length(P);
    a=P./[n:-1:1];
    a=[a,0];
    n=n+1;
    b(1)=a(1);
    for k=2:n;
        b(k)=a(k)+c*b(k-1);
    end
    y=b(n);
end
```

- 运行结果

```matlab
>> work4_3_5_2
C =
   0.016666666666667  -0.291666666666657   2.000000000000000  -6.708333333333030   9.983333333333576  61.000000000000000
L =
  -0.008333333333333   0.166666666666667  -1.291666666666667   4.833333333333333  -8.699999999999999   6.000000000000000
   0.041666666666667  -0.791666666666667   5.708333333333334 -19.208333333333332  29.250000000000000 -15.000000000000000
  -0.083333333333333   1.500000000000000 -10.083333333333334  31.000000000000000 -42.333333333333336  20.000000000000000
   0.083333333333333  -1.416666666666667   8.916666666666666 -25.583333333333332  33.000000000000000 -15.000000000000000
  -0.041666666666667   0.666666666666667  -3.958333333333333  10.833333333333332 -13.499999999999998   6.000000000000000
   0.008333333333333  -0.125000000000000   0.708333333333333  -1.875000000000000   2.283333333333334  -1.000000000000000
ave =
  64.500000000009422

```

6个测量点和插值曲线

![4_3_5_2](C:\Users\Chemizi\Desktop\数值分析\4_3_5_2.png)



- 结果分析

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

------

## work4.4.4-1		牛顿插值法

```matlab
function work4_4_4_1()
    clear all;
    close all;
    %1
    X=[1,2,2.5];
    Y=X+2./X;
    [C,D]=newpoly(X,Y);
    x=0:0.1:3;
    y=polyval(C,x);
    figure(1);
    plot(x,y);grid on;
    hold on;scatter(X,Y);
    y_15=polyval(C,1.5);
    y_12=polyval(C,1.2);
    y_15,y_12
    %2
    clear all;
    X=[0.5,1,2,2.5];
    Y=X+2./X;
    [C,D]=newpoly(X,Y);
    x=0:0.1:3;
    y=polyval(C,x);
    figure(2);
    plot(x,y);grid on;
    hold on;scatter(X,Y);
    y_15=polyval(C,1.5);
    y_12=polyval(C,1.2);
    y_15,y_12
   
end

function [C,D]=newpoly(X,Y)
    n=length(X);
    D=zeros(n,n);
    D(:,1)=Y';
    for j=2:n
        for k=j:n
            D(k,j)=(D(k,j-1)-D(k-1,j-1))/(X(k)-X(k-j+1));
        end
    end
    C=D(n,n);
    for k=(n-1):-1:1
        C=conv(C,poly(X(k)));
        m=length(C);
        C(m)=C(m)+D(k,k);
    end
end
```

- 运行结果

```matlab
>> work4_4_4_1
y_15 =
   2.900000000000000
y_12 =
   2.936000000000000
y_15 =
   2.700000000000000
y_12 =
   2.769600000000001
```

3点的牛顿多项式插值

![4_4_4_1_a](C:\Users\Chemizi\Desktop\数值分析\4_4_4_1_a.png)

4点的牛顿多项式插值

![4_4_4_1_b](C:\Users\Chemizi\Desktop\数值分析\4_4_4_1_b.png)

- 结果分析

利用拉格朗日多项式的缺点就是如果加入一个新点，整个插值系数都要发生变动，计算量还是较大。因此设计产生了牛顿多项式，对于新加入一个点产生的新多项式$P_{N+1}(x)$可以由原来的多项式$P_N(x)$推导，不需要重新计算原来的系数。
$$
\begin{array}{a}
P_1(x)=a_0+a_1(x-x_0)\\
P_2(x)=a_0+a_1(x-x_0)+a_2(x-x_0)(x-x_1)\\
\dots \dots \\
P_N(x)=a_0+a_1(x-x_0)+a_2(x-x_0)(x-x_1)+\dots+a_N\prod_{j=0}^N(x-x_j) \\
\rightarrow P_{N+1}=P_N(x)+a_N\prod_{j=0}^N(x-x_j) \\
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
随着插值点的增多，多项式曲线越接近与原来函数的样子。

------

## work5.1.4-1		线性拟合

```matlab
function work5_1_4_1()
    %a
    x_a=0.2:0.2:1.0;
    F_a=[3.6,7.3,10.9,14.5,18.2];
    [A_a,B_a]=lsline(x_a,F_a);
    A_a,B_a
    %b
    x_b=0.2:0.2:1.0;
    F_b=[5.3,10.6,15.9,21.2,26.4];
    [A_b,B_b]=lsline(x_b,F_b);
    A_b,B_b
end

function [A,B]=lsline(X,Y)
    xmean=mean(X);
    ymean=mean(Y);
    sumx2=(X-xmean)*(X-xmean)';
    sumxy=(Y-ymean)*(X-xmean)';
    A=sumxy/sumx2;
    B=ymean-A*xmean;
end
```

- 运行结果

```matlab
>> work5_1_4_1
A_a =
  18.199999999999999
B_a =
  -0.020000000000000
A_b =
  26.399999999999995
B_b =
   0.040000000000004
```

- 结果分析

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
利用直线拟合得到$F_a=18.20x_a-0.02$ ，$F_b=26.40x_b+0.04$ 。

-------

## work5.2.9-1		最小二乘曲线

```matlab
function work5_2_9_1()
   close all;
   clear all;

    x=1:24;
    y=[58,58,58,58,57,57,57,58,60,64,67,68,66,66,65,64,63,63,62,61,60,60,59,58];
    E=@(u)(sum((u(1).*cos(u(2).*x)+u(3).*sin(u(4).*x)+u(5)-y).^2));
    [u,err]=fmincon(E,[1,1,-3,1,63]);
    u,err
    figure;
    y_=@(x)(u(1).*cos(u(2).*x)+u(3).*sin(u(4).*x)+u(5));
    x_=1:0.01:24;
    plot(x_,feval(y_,x_));grid on;
    hold on;scatter(x,y);
    
end
```

- 运行结果

```matlab
>> work5_2_9_1
u =
  -1.443066176033485   0.166409623692637  -3.881169108052652   0.351366363915423  61.576146393718510
err =
  33.224846943051176
```

![5_2_9_1](C:\Users\Chemizi\Desktop\数值分析\5_2_9_1.png)

- 结果分析

利用$y=A\cos(Bx)+C\sin(Dx)+E$ 来拟合温度变化点，初始的参数对于拟合的结果有很大的影响。

------

## work5.3.8-1		三次样条插值

```matlab
function work5_3_8_1()
    t=0:2:8;
    d=[0,40,160,300,480];
    dx0=0;
    dxn=98;

    S=csfit(t,d,dx0,dxn);
    S
    t_=0:0.01:8;
    d_=zeros(1,length(t_));
    for i=1:length(t_)
        if t_(i)<=2
            d_(i)=polyval(S(1,:),t_(i));
        elseif t_(i)<=4
            d_(i)=polyval(S(2,:),t_(i)-2);
        elseif t_(i)<=6
            d_(i)=polyval(S(3,:),t_(i)-4);
        else
            d_(i)=polyval(S(4,:),t_(i)-6);
        end
    end
    plot(t_,d_);
    hold on;scatter(t,d);grid on;
end
function S=csfit(X,Y,dx0,dxn)
    N=length(X)-1;
    H=diff(X);
    D=diff(Y)./H;
    A=H(2:N-1);
    B=2*(H(1:N-1)+H(2:N));
    C=H(2:N);
    U=6*diff(D);
    
    B(1)=B(1)-H(1)/2;
    U(1)=U(1)-3*(D(1)-dx0);
    B(N-1)=B(N-1)-H(N)/2;
    U(N-1)=U(N-1)-3*(dxn-D(N));
    
    for k=2:N-1
        temp=A(k-1)/B(k-1);
        B(k)=B(k)-temp*C(k-1);
        U(k)=U(k)-temp*U(k-1);
    end
    M(N)=U(N-1)/B(N-1);
    for k=N-2:-1:1
        M(k+1)=(U(k)-C(k)*M(k+2))/B(k);
    end
    M(1)=3*(D(1)-dx0)/H(1)-M(2)/2;
    M(N+1)=3*(dxn-D(N))/H(N)-M(N)/2;
    for k=0:N-1
        S(k+1,1)=(M(k+2)-M(k+1))/(6*H(k+1));
        S(k+1,2)=M(k+1)/2;
        S(k+1,3)=D(k+1)-H(k+1)*(2*M(k+1)+M(k+2))/6;
        S(k+1,4)=Y(k+1);
    end
end
    

```

- 运行结果

```matlab
>> work5_3_8_1
S =
    0.8125    8.3750         0         0
   -2.4375   13.2500   43.2500   40.0000
    1.4375   -1.3750   67.0000  160.0000
   -0.8125    7.2500   78.7500  300.0000
```

三次样条插值结果

![5_3_8_1](C:\Users\Chemizi\Desktop\数值分析\5_3_8_1.png)

- 结果分析

得到的矩阵$S$为三次样条插值系数，每一行表示一个区间内三次多项式的四个系数。需要注意的是，该系数并非直接的多项式，例如第二个区间内的曲线为

$$S_1(x)=-2.4375(x-2)^3+13.2500(x-2)^2+43.2500(x-2)+40.0000$$ 

三次样条插值的条件可以列为：
$$
\begin{array}{a}
1. \quad S(x)=S_k(x)=s_{k,0}+s_{k,1}(x-x_k)+s_{k,2}{(x-x_k)}^2+s_{k,3}{(x-x_k)}^3 \\
\qquad \qquad \qquad \qquad \qquad \qquad  x \in [x_k,x_{k+1}]  \quad k=0,1,\dots,N-1 \\
2. \quad S(x_k)=y_k \qquad  \qquad \qquad k=0,1,\dots,N \\
3. \quad S_k(x_{k+1})=S_{k+1}(x_{k+1}) \quad k=0,1,\dots,N-2  \\
4. \quad S'_k(x_{k+1})=S'_{k+1}(x_{k+1}) \quad k=0,1,\dots,N-2 \\
5. \quad S''_k(x_{k+1})=S''_{k+1}(x_{k+1}) \quad k=0,1,\dots,N-2 \\
\end{array}
$$
条件1表示插值多项式为分段函数，条件2表示函数穿过插值点，条件3说明各分段间连续，与条件2 不同的是，条件3强调的是各分段函数首尾相连，条件4说明各分段间光滑，条件5定义函数上各点二阶导数也是连续的，存在曲率半径。

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

-------

## work5.4.3-2		三角多项式

```matlab
function work5_4_3_2()
    clear all;
    close all;
    N=12;
    k=0:N;
    x_12=-pi+2*k*pi/N;
%     x_12=linspace(-pi,pi,N);
    y_12=0.5.*x_12;
    [A_12,B_12]=tpcoeff(x_12,y_12,5);
    x_12_=-pi:0.01:pi;
    y_12_=tp(A_12,B_12,x_12_,5);
    plot(x_12_,y_12_);
    hold on;scatter(x_12,y_12);
    grid on;
    
    clear all;
    N=60;
    k=0:N;
    x_60=-pi+2*k*pi/N;
%     x_60=linspace(-pi,pi,N);
    y_60=0.5.*x_60;
    [A_60,B_60]=tpcoeff(x_60,y_60,5);
    A_60,B_60
    x_60_=-pi:0.01:pi;
    y_60_=tp(A_60,B_60,x_60_,5);
    plot(x_60_,y_60_);
    %hold on;scatter(x_60,y_60);
    
    clear all;
    N=360;
    k=0:N;
    x_360=-pi+2*k*pi/N;
%     x_360=linspace(-pi,pi,N);
    y_360=0.5.*x_360;
    [A_360,B_360]=tpcoeff(x_360,y_360,5);
    A_360,B_360
    x_360_=-pi:0.01:pi;
    y_360_=tp(A_360,B_360,x_360_,5);
    hold on;plot(x_360_,y_360_);
    grid on
    legend('12','12','60','360');
end

function [A,B]=tpcoeff(X,Y,M)
    N=length(X)-1;
    max1=fix((N-1)/2);
    if M>max1
        M=max1;
    end
    A=zeros(1,M+1);
    B=zeros(1,M+1);
    Yends=(Y(1)+Y(N+1))/2;
    Y(1)=Yends;
    Y(N+1)=Yends;
    A(1)=sum(Y);
    for j=1:M
        A(j+1)=cos(j*X)*Y';
        B(j+1)=sin(j*X)*Y';
    end
    A=2*A/N;
    B=2*B/N;
    A(1)=A(1)/2;
end

function z=tp(A,B,x,M)
    z=A(1);
    for j=1:M
        z=z+A(j+1)*cos(j*x)+B(j+1)*sin(j*x);
    end
end
```

- 运行结果

```matlab
>> work5_4_3_2
A_60 =
   1.0e-15 *
   0.062912638062092  -0.022204460492503   0.229446091755866  -0.296059473233375   0.105471187339390  -0.229446091755866
B_60 =
                   0   0.999085980671827  -0.498170957882692   0.330587256251576  -0.246333856496581   0.195409723331371
A_360 =
   1.0e-16 *
  -0.074014868308344   0.148029736616688  -0.592118946466750   0.641462192005646  -0.160365548001411  -0.444089209850063
B_360 =
                   0   0.999974615086140  -0.499949229398986   0.333257175498468  -0.249898452610956   0.199873059962485
```

12点，60点，360点的三角多项式

![5_4_3_2](C:\Users\Chemizi\Desktop\数值分析\5_4_3_2.png)

- 结果分析

随着点数的增加，拟合效果并没有进一步提升了，5阶多项式的级数变化越来越小，最终近似于傅里叶级数$1.0,-0.5,0.33,-0.25,0.2$。三角多项式的形式如下：
$$
P(x)=\frac{a_0}{2}+\sum^M_{j=1}(a_j\cos{jx}+b_j\sin{jx})
$$

------

## work7.2.3-2		组合辛普森积分

```matlab
function work7_2_3_2()
    format long
    f_a=@(x)(x.^3);f_a_deri=@(x)(3.*x.^2);L_a=@(x)(2*pi*feval(f_a,x)*sqrt(1+(feval(f_a_deri,x)).^2));
    f_b=@(x)(sin(x));f_b_deri=@(x)(cos(x));L_b=@(x)(2*pi*feval(f_b,x)*sqrt(1+(feval(f_b_deri,x)).^2));
    f_c=@(x)(exp(-x));f_c_deri=@(x)(-exp(-x));L_c=@(x)(2*pi*feval(f_c,x)*sqrt(1+(feval(f_c_deri,x)).^2));
    
    delta=10.^-11;
    a=0;b=1;
    M=fix(abs(b-a)/(2*(delta*180/abs(b-a)).^0.25))+1;
    s_a=simprl(L_a,a,b,M);
    s_a
   
    
    delta=10.^-11;
    a=0;b=pi/4;
    M=fix(abs(b-a)/(2*(delta*180/abs(b-a)).^0.25))+1;
    s_b=simprl(L_b,a,b,M);
    s_b
    
    delta=10.^-11;
    a=0;b=1;
    M=fix(abs(b-a)/(2*(delta*180/abs(b-a)).^0.25))+1;
    s_c=simprl(L_c,a,b,M);
    s_c

end
function s=simprl(f,a,b,M)
    h=(b-a)/(2*M);
    s1=0;
    s2=0;
    for k=1:M
        x=a+h*(2*k-1);
        s1=s1+feval(f,x);
    end
    for k=1:(M-1)
        x=a+h*2*k;
        s2=s2+feval(f,x);
    end
    s=h*(feval(f,a)+feval(f,b)+4*s1+2*s2)/3;
end
```

- 运行结果

```matlab
>> work7_2_3_2
s_a =
   3.563121862823188
s_b =
   2.422428051161270
s_c =
   4.849218491907431
```

- 结果分析

组合辛普森公式即在各个区间上使用辛普森公式，再求和。

辛普森公式为
$$
\int_{x_0}^{x_2}f(x)dx \approx \frac{h}{3}(f_0+4f_1+f_2)
$$
根据公式$S=2\pi{\int^b_a}f(x) \sqrt{1+(f'(x))^2}dx$ ，求得$$S_a=3.563121862823188,S_b=2.422428051161270,S_c=4.849218491907431$$

------

## work7.3.2-2		龙贝格积分

```matlab
function work7_3_2_2()
    tol=10.-10;n=5;
    f=@(x)(sqrt(4*x-x.^2));
    a=0;b=2;
    [R_a,quad,err,h]=romber(f,a,b,n,tol);
    R_a
    
    tol=10.-10;n=5;
    f=@(x)(4/(1+x.^2));
    a=0;b=1;
    [R_b,quad,err,h]=romber(f,a,b,n,tol);
    R_b
    
end
function [R,quad,err,h]=romber(f,a,b,n,tol)
    M=1;
    h=b-a;
    err=1;
    J=0;
    R=zeros(4,4);
    R(1,1)=h*(f(a)+f(b))/2;

    while((err>tol)&&(J<n))||(J<4)
        J=J+1;
        h=h/2;
    s=0;
    for p=1:M
    	x=a+h*(2*p-1);
        s=s+f(x);
    end
    R(J+1,1)=R(J,1)/2+h*s;
    M=2*M;
    for K=1:J
        R(J+1,K+1)=R(J+1,K)+(R(J+1,K)-R(J,K))/(4^K-1);
    end
    err=abs(R(J,J)-R(J+1,K+1));
    end
    quad=R(J+1,J+1);
end

```

- 运行结果

```matlab
>> work7_3_2_2
R_a =
    2.0000         0         0         0         0         0
    2.7321    2.9761         0         0         0         0
    2.9957    3.0836    3.0908         0         0         0
    3.0898    3.1212    3.1237    3.1242         0         0
    3.1233    3.1344    3.1353    3.1355    3.1355         0
    3.1351    3.1391    3.1394    3.1394    3.1394    3.1394
R_b =
    3.0000         0         0         0         0         0
    3.1000    3.1333         0         0         0         0
    3.1312    3.1416    3.1416    3.1416    3.1416         0
    3.1414    3.1416    3.1416    3.1416    3.1416    3.1416
```

- 结果分析

可以明显看到，第二种的龙贝格序列积分速度明显快于第一种，开方会明显降低速度。

龙贝格积分是在连续积分的背景下，改进提出的，每一次改进后的误差项为$O(h^{2N+2})$。

未改进公式（递归梯形公式）$\rightarrow$第一次改进公式（递归辛普森公式）$\rightarrow$第二次改进公式（递归布尔公式）$\rightarrow \dots$
$$
\begin{array}{a} 
T(0)=(h/2)(f(a)+f(b))\\
R(J,0)=T(J)=\frac{T(J-1)}{2}+h\sum_{k=1}^{M}f(x_{2k-1}) \quad while \{h=\frac{b-a}{2^J}  \quad x_k=a+kh\}\\
R(J,1)=S(J)=\frac{4T(J)-T(J-1)}{3}\\
R(J,2)=B(J)=\frac{16S(J)-S(J-1)}{15}\\
R(J,K)=\frac{4^KR(J,K-1)-R(J-1,K-1)}{4^K-1} \quad J \geq K\\
\end{array}
$$

------

## work9.2.4-8		欧拉方法

```matlab
function work9_2_4_8()
    f=@(t,v)(32-0.032*v.^1.5);
    a=0;
    b=6;
    ya=0;
    M=(b-a)/0.05;
    E=euler(f,a,b,ya,M);
    E
end

function E=euler(f,a,b,ya,M)
    h=(b-a)/M;
    T=zeros(1,M+1);
    Y=zeros(1,M+1);
    T=a:h:b;
    Y(1)=ya;
    for j=1:M
        Y(j+1)=Y(j)+h*feval(f,T(j),Y(j));
    end
    E=[T' Y'];
end
```

- 运行结果

```matlab
>> work9_2_4_8
E =
         0         0
    0.0500    1.6000
    0.1000    3.1968
    0.1500    4.7876
    0.2000    6.3709
    0.2500    7.9451
    0.3000    9.5093
    0.3500   11.0624
    0.4000   12.6035
    0.4500   14.1319
    0.5000   15.6469
    0.5500   17.1479
    0.6000   18.6343
    0.6500   20.1056
    0.7000   21.5613
    0.7500   23.0011
    0.8000   24.4246
    0.8500   25.8315
    0.9000   27.2214
    0.9500   28.5942
    1.0000   29.9496
    1.0500   31.2873
    1.1000   32.6073
    1.1500   33.9094
    1.2000   35.1934
    1.2500   36.4594
    1.3000   37.7072
    1.3500   38.9367
    1.4000   40.1480
    1.4500   41.3409
    1.5000   42.5156
    1.5500   43.6721
    1.6000   44.8103
    1.6500   45.9304
    1.7000   47.0323
    1.7500   48.1163
    1.8000   49.1822
    1.8500   50.2304
    1.9000   51.2608
    1.9500   52.2736
    2.0000   53.2688
    2.0500   54.2468
    2.1000   55.2075
    2.1500   56.1512
    2.2000   57.0780
    2.2500   57.9880
    2.3000   58.8815
    2.3500   59.7586
    2.4000   60.6195
    2.4500   61.4643
    2.5000   62.2933
    2.5500   63.1066
    2.6000   63.9045
    2.6500   64.6872
    2.7000   65.4547
    2.7500   66.2074
    2.8000   66.9455
    2.8500   67.6691
    2.9000   68.3785
    2.9500   69.0738
    3.0000   69.7552
    3.0500   70.4231
    3.1000   71.0775
    3.1500   71.7188
    3.2000   72.3470
    3.2500   72.9624
    3.3000   73.5652
    3.3500   74.1557
    3.4000   74.7339
    3.4500   75.3002
    3.5000   75.8548
    3.5500   76.3977
    3.6000   76.9293
    3.6500   77.4497
    3.7000   77.9591
    3.7500   78.4578
    3.8000   78.9459
    3.8500   79.4236
    3.9000   79.8911
    3.9500   80.3485
    4.0000   80.7962
    4.0500   81.2342
    4.1000   81.6627
    4.1500   82.0820
    4.2000   82.4921
    4.2500   82.8933
    4.3000   83.2858
    4.3500   83.6697
    4.4000   84.0451
    4.4500   84.4124
    4.5000   84.7715
    4.5500   85.1227
    4.6000   85.4661
    4.6500   85.8019
    4.7000   86.1303
    4.7500   86.4513
    4.8000   86.7652
    4.8500   87.0721
    4.9000   87.3721
    4.9500   87.6654
    5.0000   87.9521
    5.0500   88.2324
    5.1000   88.5063
    5.1500   88.7741
    5.2000   89.0358
    5.2500   89.2916
    5.3000   89.5416
    5.3500   89.7859
    5.4000   90.0247
    5.4500   90.2580
    5.5000   90.4860
    5.5500   90.7088
    5.6000   90.9265
    5.6500   91.1393
    5.7000   91.3472
    5.7500   91.5503
    5.8000   91.7487
    5.8500   91.9426
    5.9000   92.1320
    5.9500   92.3171
    6.0000   92.4979
```

- 结果分析

对于方程$v'=32-0.032v^{\frac{3}{2}}$ 利用欧拉方法$y_{k+1}=y_k+hf(t_k,y_k)$ 迭代计算，其中$h=\frac{b-a}{M}$，$M$表示划分的区间个数，$t_k=a+kh$。至于$f(t_k,y_k)=v'=32-0.032v^{\frac{3}{2}}$。

------

## work9.3.3-6		休恩方法

```matlab
function work9_3_3_6()
format short
    f=@(t,v)(-32-0.1*v);
    a=0;
    b=30;
    ya=160;
    M=(b-a)/0.5;
    H=heun(f,a,b,ya,M);
    H
    scatter(H(:,1),H(:,2));
    t=0:0.01:30;
    v=480*exp(-t/10)-320;
    hold on;plot(t,v);
    legend('Heun','Real');
    grid on;
end

function H=heun(f,a,b,ya,M)
    h=(b-a)/M;
    T=zeros(1,M+1);
    Y=zeros(1,M+1);
    T=a:h:b;
    Y(1)=ya;
    for j=1:M
        k1=feval(f,T(j),Y(j));
        k2=feval(f,T(j+1),Y(j)+h*k1);
        Y(j+1)=Y(j)+(h/2)*(k1+k2);
    end
    H=[T' Y'];
end
```

- 运行结果

```matlab
>> work9_3_3_6
H =
         0  160.0000
    0.5000  136.6000
    1.0000  114.3407
    1.5000   93.1666
    2.0000   73.0248
    2.5000   53.8648
    3.0000   35.6389
    3.5000   18.3015
    4.0000    1.8093
    4.5000  -13.8789
    5.0000  -28.8023
    5.5000  -42.9982
    6.0000  -56.5020
    6.5000  -69.3476
    7.0000  -81.5669
    7.5000  -93.1905
    8.0000 -104.2474
    8.5000 -114.7654
    9.0000 -124.7706
    9.5000 -134.2880
   10.0000 -143.3415
   10.5000 -151.9536
   11.0000 -160.1458
   11.5000 -167.9387
   12.0000 -175.3517
   12.5000 -182.4033
   13.0000 -189.1112
   13.5000 -195.4920
   14.0000 -201.5617
   14.5000 -207.3356
   15.0000 -212.8280
   15.5000 -218.0526
   16.0000 -223.0226
   16.5000 -227.7502
   17.0000 -232.2474
   17.5000 -236.5253
   18.0000 -240.5947
   18.5000 -244.4657
   19.0000 -248.1480
   19.5000 -251.6508
   20.0000 -254.9828
   20.5000 -258.1524
   21.0000 -261.1675
   21.5000 -264.0356
   22.0000 -266.7638
   22.5000 -269.3591
   23.0000 -271.8278
   23.5000 -274.1762
   24.0000 -276.4101
   24.5000 -278.5352
   25.0000 -280.5566
   25.5000 -282.4794
   26.0000 -284.3086
   26.5000 -286.0485
   27.0000 -287.7037
   27.5000 -289.2781
   28.0000 -290.7758
   28.5000 -292.2005
   29.0000 -293.5557
   29.5000 -294.8449
   30.0000 -296.0712
```

![9_3_3_6](C:\Users\Chemizi\Desktop\数值分析\9_3_3_6.png)

- 结果分析

休恩方法是对于欧拉方法的改进，基本步骤为：

$$p_{k+1}=y_k+hf(t_k,y_k)\rightarrow y_{k+1}=y_k+\frac{h}{2}(f(t_k,y_k)+f(t_{k+1},p_{k+1}))$$

------

## work9.4.2-1		三次泰勒方法

```matlab
function work9_4_2_1()
    Dy=@(t,y)(t.^2+y.^2);%y'=t^2+y^2;
    D2y=@(t,y)(2*t+2*y*feval(Dy,t,y));%y''=2t+2y*y';
    D3y=@(t,y)(2+2*(feval(Dy,t,y)).^2+2*y*feval(D2y,t,y));
    df={Dy,D2y,D3y};
    a=0;
    b=0.8;
    ya=0;
    h1=0.05;h2=0.025;h3=0.0125;h4=0.00625;
    T3_1=taylor(df,a,b,ya,(b-a)/h1);
%     subplot(221);
    hold on;scatter(T3_1(:,1),T3_1(:,2),'*');grid on;
    T3_2=taylor(df,a,b,ya,(b-a)/h2);
%     subplot(222);
    hold on;scatter(T3_2(:,1),T3_2(:,2),'o');grid on;
    T3_3=taylor(df,a,b,ya,(b-a)/h3);
%     subplot(223);
    hold on;scatter(T3_3(:,1),T3_3(:,2),'^');grid on;
    T3_4=taylor(df,a,b,ya,(b-a)/h4);
%     subplot(224);
    hold on;scatter(T3_4(:,1),T3_4(:,2),'+');grid on;
    legend('0.05','0.025','0..125','0.00625');
end

function T3=taylor(df,a,b,ya,M)
    h=(b-a)/M;
    T=zeros(1,M+1);
    Y=zeros(1,M+1);
    T=a:h:b;
    Y(1)=ya;

    for j=1:M
    	D(1)=feval(df{1},T(j),Y(j));%df=[y' y'' y'''];
        D(2)=feval(df{2},T(j),Y(j));
        D(3)=feval(df{3},T(j),Y(j));
        Y(j+1)=Y(j)+h*(D(1)+h*(D(2)/2+h*(D(3)/6)));
    end
    T3=[T' Y'];
end
```

- 运行结果

![9_4_2_1](C:\Users\Chemizi\Desktop\数值分析\9_4_2_1.png)

- 结果分析

三阶泰勒方法为：

$$y_{k+1}=y_k+d_1h+\frac{d_2h^2}{2!}+\frac{d_3h^3}{3!}$$

其中，$d_j=y^{(j)}(t_k)\quad j=1,2,3$

------

## work9.5.6-6		龙格-库塔方法

```matlab
function work9_5_6_6()
    f=@(t,y)(9*t*exp(3*t));
    a=0;b=3;
    ya=0;
    h=0.1;
    M=(b-a)/h;
    tol=10.^-7;
    R=rkf45(f,a,b,ya,M,tol);
    scatter(R(:,1),R(:,2),'x');
    t=0:0.01:2.5;
    y=3.*t.*exp(3.*t)-exp(3.*t)+1;
    hold on;
    plot(t,y);grid on;
    
    
end

function R=rk4(f,a,b,ya,M)

h=(b-a)/M;
T=zeros(1,M+1);
Y=zeros(1,M+1);
T=a:h:b;
Y(1)=ya;
for j=1:M
   k1=h*f(T(j),Y(j));
   k2=h*f(T(j)+h/2,Y(j)+k1/2);
   k3=h*f(T(j)+h/2,Y(j)+k2/2);
   k4=h*f(T(j)+h,Y(j)+k3);
   Y(j+1)=Y(j)+(k1+2*k2+2*k3+k4)/6;
end

R=[T' Y'];
end

function R = rkf45(f,a,b,ya,m,tol)

a2 = 1/4; b2 = 1/4; a3 = 3/8; b3 = 3/32; c3 = 9/32; a4 = 12/13;
b4 = 1932/2197; c4 = -7200/2197; d4 = 7296/2197; a5 = 1;
b5 = 439/216; c5 = -8; d5 = 3680/513; e5 = -845/4104; a6 = 1/2;
b6 = -8/27; c6 = 2; d6 = -3544/2565; e6 = 1859/4104; f6 = -11/40;
r1 = 1/360; r3 = -128/4275; r4 = -2197/75240; r5 = 1/50;
r6 = 2/55; n1 = 25/216; n3 = 1408/2565; n4 = 2197/4104; n5 = -1/5;
big = 1e15;
h = (b-a)/m;
hmin = h/64;
hmax = 64*h;
max1 = 200;
Y(1) = ya;
T(1) = a;
j = 1;
tj = T(1);
br = b - 0.00001*abs(b);
while (T(j)<b),
  if ((T(j)+h)>br), h = b - T(j); end
  tj = T(j);
  yj = Y(j);
  y1 = yj;
  k1 = h*f(tj,y1);
  y2 = yj+b2*k1;                         if big<abs(y2) break, end
  k2 = h*f(tj+a2*h,y2);
  y3 = yj+b3*k1+c3*k2;                   if big<abs(y3) break, end
  k3 = h*f(tj+a3*h,y3);
  y4 = yj+b4*k1+c4*k2+d4*k3;             if big<abs(y4) break, end
  k4 = h*f(tj+a4*h,y4);
  y5 = yj+b5*k1+c5*k2+d5*k3+e5*k4;       if big<abs(y5) break, end
  k5 = h*f(tj+a5*h,y5);
  y6 = yj+b6*k1+c6*k2+d6*k3+e6*k4+f6*k5; if big<abs(y6) break, end
  k6 = h*f(tj+a6*h,y6);
  err = abs(r1*k1+r3*k3+r4*k4+r5*k5+r6*k6);
  ynew = yj+n1*k1+n3*k3+n4*k4+n5*k5;
  if ((err<tol)|(h<2*hmin)),
    Y(j+1) = ynew;
    if ((tj+h)>br),
      T(j+1) = b;
    else
      T(j+1) = tj + h;
    end
    j = j+1;
    tj = T(j); 
  end
  if (err==0),
    s = 0;
  else
    s = 0.84*(tol*h/err)^(0.25);
  end
  if ((s<0.75)&(h>2*hmin)), h = h/2; end
  if ((s>1.50)&(2*h<hmax)), h = 2*h; end
  if ((big<abs(Y(j)))|(max1==j)), break, end
  mend = j;
  if (b>T(j)),
    m = j+1;
  else
    m = j;
  end
end
R=[T' Y'];
end
```

- 运行结果

![9_5_6_6](C:\Users\Chemizi\Desktop\数值分析\9_5_6_6.png)

------

## work10.1.7-1		双曲型方程

```matlab
function [U ]= work10_1_7_1()
    f=@(x)(sin(pi.*x));
    g=@(x)(0);
    a=1;b=1;c=0;
    h=0.1;k=0.1;
    n=a/h;m=b/k;
    U=finedif(f,g,a,b,c,n,m);
%     contour(U);
    x=linspace(0,a,n);
    t=linspace(0,b,m);
    [X,T]=meshgrid(x,t);
    surf(X,T,U);
end 

function U = finedif(f,g,a,b,c,n,m)
    h = a/(n-1);
    k = b/(m-1);
    r = c*k/h;
    r2=r^2;
    r22=r^2/2;
    s1 = 1 - r^2;
    s2 = 2 - 2*r^2;
    U = zeros(n,m);

    %Comput first and second rows
    for i=2:n-1
        U(i,1)=f(h*(i-1));
        U(i,2)=s1*f(h*(i-1))+k*g(h*(i-1)) ...
          +r22*(f(h*i)+f(h*(i-2)));
    end
    
    %Compute remaining rows of U
    for j=3:m,
        for i=2:(n-1),
            U(i,j) = s2*U(i,j-1)+r2*(U(i-1,j-1)+U(i+1,j-1))-U(i,j-2);
        end
    end
    U=U';
end
```

- 运行结果

![10_1_7_1](C:\Users\Chemizi\Desktop\数值分析\10_1_7_1.png)

- 结果分析

利用差分法求解波动方程$u_{tt}(x,t)=c^2u_{xx}(x,t)$，需要的初始条件为求解区间$R=\{ (x,t)|0 \leq x \leq a, 0 \leq t \leq b\}$，边界条件$u(0,t)=0,u(a,t)=0$，已知$u(x,0)=f(x),u_t(x,0)=g(x)$。

利用下列关系式求解：
$$
\begin{array}{a}
u_{i,j+1}=u_{i+1,j}+u_{i-1,j}-u_{i,j-1}\\
u_{i,1}=u(x_i,0)\\
u_{i,2}=u(x_i,k)\\
i=1,2,\dots,n\\
\end{array}
$$

------

## work10.2.5-1		抛物型方程

```matlab
function work10_2_5_1()
    f=@(x)(sin(pi.*x)+sin(2*pi.*x));
   c1=1;c2=1;a=1;b=0.1;r=1;
    h=0.1;k=0.01;c=sqrt(r*h^2/k);
    n=a/h;m=b/k;
    U=crnich(f,c1,c2,a,b,c,n,m);
%     contour(U);
    x=linspace(0,a,n);
    t=linspace(0,b,m);
    [X,T]=meshgrid(x,t);
    surf(X,T,U);
end 
function U=crnich(f,c1,c2,a,b,c,n,m)

h=a/(n-1);
k=b/(m-1);
r=c^2*k/h^2;
s1=2+2/r;
s2=2/r-2;
U=zeros(n,m);

%Boundary conditions

U(1,1:m)=c1;
U(n,1:m)=c2;

%Generate first row

U(2:n-1,1)=f(h:h:(n-2)*h)';

%Form the diagonal and off-diagonal elemnts of A and 
%the constant vector B and solve tridiagonal system AX=B

Vd(1,1:n)=s1*ones(1,n);
Vd(1)=1;
Vd(n)=1;
Va=-ones(1,n-1);
Va(n-1)=0;
Vc=-ones(1,n-1);
Vc(1)=0;
Vb(1)=c1;
Vb(n)=c2;
for j=2:m
   for i=2:n-1
      Vb(i)=U(i-1,j-1)+U(i+1,j-1)+s2*U(i,j-1);
   end
   X=trisys(Va,Vd,Vc,Vb);
   U(1:n,j)=X';
end

U=U';
end

function X=trisys(A,D,C,B)


N=length(B);
for k=2:N
   mult=A(k-1)/D(k-1);
   D(k)=D(k)-mult*C(k-1);
   B(k)=B(k)-mult*B(k-1);
end

X(N)=B(N)/D(N);

for k= N-1:-1:1
   X(k)=(B(k)-C(k)*X(k+1))/D(k);
end
end

```

- 运行结果

![10_2_5_1](C:\Users\Chemizi\Desktop\数值分析\10_2_5_1.png)

------

## work10.3.8-1		椭圆型方程



```matlab
function work10_3_8_1()
    f1=@(x)(x.^4);
    f2=@(x)(x.^4-13.5.*x.^2+5.0625);
    f3=@(y)(y.^4);
    f4=@(y)(y.^4-13.5.*y.^2+5.0625);
    
   a=1.5;b=1.5;
    h=0.5;
    tol=10^-9;
    max1=100000;
    U=dirich(f1,f2,f3,f4,a,b,h,tol,max1);
%     contour(U);
    x=0:h:a;
    y=0:h:b;
    [X,Y]=meshgrid(x,y);
    surf(X,Y,U);
    
    
    hold on;
    a=1.5;b=1.5;
    x_=0:0.1:a;
    y_=0:0.1:b;
    [X_,Y_]=meshgrid(x_,y_);
    U_=X_.^4-6*X_.^2.*Y_.^2+Y_.^4;
     h=surf(X_,Y_,U_);
     set(h,'edgecolor','none');
%      shading interp
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
end 

function U=dirich(f1,f2,f3,f4,a,b,h,tol,max1)



n=fix(a/h)+1;
m=fix(b/h)+1;
ave=(a*(f1(0)+f2(0)) ...
   +b*(f3(0)+f4(0)))/(2*a+2*b);
U=ave*ones(n,m);

%Boundary conditions

U(1,1:m)=f3(0:h:(m-1)*h)';
U(n,1:m)=f4(0:h:(m-1)*h)';
U(1:n,1)=f1(0:h:(n-1)*h);
U(1:n,m)=f2(0:h:(n-1)*h);
U(1,1)=(U(1,2)+U(2,1))/2;
U(1,m)=(U(1,m-1)+U(2,m))/2;
U(n,1)=(U(n-1,1)+U(n,2))/2;
U(n,m)=(U(n-1,m)+U(n,m-1))/2;

%SOR parameter

w=4/(2+sqrt(4-(cos(pi/(n-1))+cos(pi/(m-1)))^2));

%Refine approximations and sweep operator throughout the grid

err=1;
cnt=0;
while((err>tol)&(cnt<=max1))
   err=0;
   for j=2:m-1
      for i=2:n-1
         relx=w*(U(i,j+1)+U(i,j-1)+U(i+1,j)+U(i-1,j)-4*U(i,j))/4;
         U(i,j)=U(i,j)+relx;
         if (err<=abs(relx))
           err=abs(relx);
         end
      end   
   end
cnt=cnt+1;
end

% U=flipud(U');

end
```

- 运行结果

![10_3_8_1](C:\Users\Chemizi\Desktop\数值分析\10_3_8_1.png)

- 结果分析

带网格的为估计的曲面形状，不带网格的为准确函数曲面。

------

