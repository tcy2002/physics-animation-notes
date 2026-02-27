# 基于力的约束

以半径为1的圆环上运动的小球为例，小球的位置为$x$，其满足的约束方程为：

$$
C(x)=\frac{1}{2}x \cdot x-\frac{1}{2}=0
$$

几何含义为小球到环心的距离始终为1( $|x|=1$ )，乘 $\frac{1}{2}$ 可以消掉一二阶导的系数。一二阶导分别为：

$$
\dot C(x)=x\cdot \dot x=0
$$

$$
\ddot C(x)=x\cdot \ddot x+\dot x \cdot \dot x=0
$$

要使小球满足约束，且不往破坏约束条件的情况下演化，需要通过一定方式使得上述三个方程全部成立。

## 约束模型

对于一般约束 $C$ ， $\frac{\partial C}{\partial x}$ 表示约束 $C$ 的梯度，那么：

$$
\dot C(x)=\frac{\partial C}{\partial x}\cdot \dot x=0\
$$

$$
\ddot C(x)=\frac{\partial \dot C}{\partial x}\cdot \dot x+\frac{\partial C}{\partial x}\cdot \ddot x=0
$$

$\ddot x=\frac{f+\hat f}{m}$ ， $f$ 为外力， $\hat f$ 为约束力，代入得到

$$
\frac{\partial C}{\partial x}\cdot \frac{\hat f}{m}=-\frac{\partial C}{\partial x}\cdot \frac{f}{m}-\frac{\partial \dot C}{\partial x}\cdot \dot x
$$

另外还有能量约束，即约束力不能做功造成能量衰减，那么需要保证 $\hat f$ 与速度 $\dot x$ 的方向垂直，即： $\hat f \cdot \dot x=0$ ；又根据约束一阶导满足的等式 $\frac{\partial C}{\partial x}\cdot \dot x=0$ ，说明 $\hat f$ 的方向应该与 $\frac{\partial C}{\partial x}$ 一致，
可以理解为 $\hat f$ 必须指向约束的梯度方向，因此

$$
\hat f=\lambda\frac{\partial C}{\partial x}
$$

$$
\begin{aligned}
\frac{\partial C}{\partial x}\cdot (\frac{1}{m} \frac{\partial C}{\partial x} \lambda)=-\frac{\partial C}{\partial x}\cdot \frac{f}{m}-\frac{\partial \dot C}{\partial x}\cdot \dot x 
\end{aligned}
\tag 1
$$

然后可以得到 $\lambda$ ，进而可以得到约束力 $\hat f$ ，这就是基于力的约束模型。

在式(1)的基础上也可以增加一种反馈机制，即： $\ddot C=-k_sC-k_d\dot C$ ，第一项 $-k_sC$ 和弹簧约束类似，当约束发生偏离时，添加一个反向的反馈促使约束复原；第二项 $-k_d\dot C$ 则可以使约束偏离的速度的不至于增大。

推广到复杂场景，存在多个约束 $C_1,C_2,...,C_m$ ，雅可比矩阵 $J\in R^{m\times 3n}$ ，逆质量矩阵 $W\in R^{3n\times 3n}$ ，位置矢量 $q=[x_1^T,x_2^T,...,x_n^T]^T$ ，外力矢量 $Q=[f_1^T,f_2^T,...,f_n^T]^T$ 

若不是粒子而是刚体，除位置外还需考虑旋转朝向，则多3个旋转自由度， $J\in R^{m\times 6n}$ ， $W\in R^{6n\times 6n}$ ， $q=[x_1^T,\alpha_1^T,x_2^T,\alpha_2^T...,x_n^T,\alpha_n^T]$ ， $Q=[f_1^T,t_1^T,f_2^T,t_2^T...,f_n^T,t_n^T]$ ， $\dot \alpha=\omega$ ， $t$ 为外力矩， $W$ 相应位置填充惯性张量矩阵的逆。

式(1)扩展为：

$$
JWJ^T\lambda=-JWQ-\dot J\dot q
$$

对于经典距离约束 $C(x_1,x_2)=\frac{1}{2}(x_1-x_2)\cdot(x_1-x_2)-\frac{1}{2}d^2=0$ ，雅可比矩阵为

$$
\begin{aligned}
J&=[(\frac{\partial C}{\partial x_1})^T,(\frac{\partial C}{\partial x_2})^T] \\
 &=[x_1-x_2,x_2-x_1] \\
 &=\left[
  \begin{matrix}
  x_{1x} - x_{2x} & x_{2x} - x_{1x} \\
  x_{1y} - x_{2y} & x_{2y} - x_{1y} \\
  x_{1z} - x_{2z} & x_{2z} - x_{1z}
  \end{matrix}
 \right]
\end{aligned}
$$

对于一般的不考虑外力作用的场景而言，有：

$$
JWJ^T\lambda=-\dot J\dot q_1-k_sC-k_d\dot C
$$

这种基于力的方法的思路是：根据以上公式求出约束力，然后综合约束力和外力，通过积分方法来求运动过程。
