## 约束方程
### 例子
以半径为1的圆环上运动的小球为例，小球的位置为$x$，其满足的约束方程为：
$$C(x)=\frac{1}{2}x \cdot x-\frac{1}{2}=0$$
几何含义为小球到环心的距离始终为1($|x|=1$)，乘$\frac{1}{2}$可以消掉一二阶导的系数。一二阶导分别为：
$$\dot C(x)=x\cdot \dot x=0$$
$$\ddot C(x)=x\cdot \ddot x+\dot x \cdot \dot x=0$$
要使小球满足约束，且不往破坏约束条件的情况下演化，需要通过一定方式使得上述三个方程全部成立。

### 传统约束模型（基于力的约束）
对于一般约束$C$，$\frac{\partial C}{\partial x}$表示约束$C$的梯度，那么：

$$\dot C(x)=\frac{\partial C}{\partial x}\cdot \dot x=0\$$
$$\ddot C(x)=\frac{\partial \dot C}{\partial x}\cdot \dot x+\frac{\partial C}{\partial x}\cdot \ddot x=0$$

$\ddot x=\frac{f+\hat f}{m}$，$f$为外力，$\hat f$为约束力，代入得到
$$\frac{\partial C}{\partial x}\cdot \hat f==-\frac{\partial C}{\partial x}\cdot f-m\frac{\partial \dot C}{\partial x}\cdot \dot x\tag{1}$$

另外还有能量约束，即约束力不能做功造成能量衰减，那么需要保证$\hat f$与速度$\dot x$的方向垂直，即：$\hat f \cdot \dot x=0$；又根据约束一阶导满足的等式$\frac{\partial C}{\partial x}\cdot \dot x=0$，说明$\hat f$的方向应该与$\frac{\partial C}{\partial x}$一致，
可以理解为$\hat f$必须指向约束的梯度方向，因此
$$\hat f=\lambda\frac{\partial C}{\partial x}$$
然后可以得到$\lambda$，进而可以得到约束力$\hat f$，这就是经典的基于力的约束模型。

在式(1)的基础上也可以增加一种反馈机制，即：$\ddot C=-k_sC-k_d\dot C$，第一项$-k_sC$和弹簧约束类似，当约束发生偏离时，添加一个反向的反馈促使约束复原；第二项$-k_d\dot C$则可以使约束偏离的速度的不至于增大。

推广到复杂场景，存在多个约束$C_1,C_2,...,C_m$，雅可比矩阵$J\in R^{m\times 3n}$，逆质量矩阵$W\in R^{3n\times 3n}$，位置矢量$q=[x_1^T,x_2^T,...,x_n^T]^T$，外力矢量$Q=[f_1^T,f_2^T,...,f_n^T]^T$ (若不是粒子，需要考虑刚体的旋转朝向，则多3个旋转自由度，$J\in R^{m\times 6n}$，$W\in R^{6n\times 6n}$，$q=[x_1^T,\alpha_1^T,x_2^T,\alpha_2^T...,x_n^T,\alpha_n^T]$，$Q=[f_1^T,t_1^T,f_2^T,t_2^T...,f_n^T,t_n^T]$，$\dot \alpha=\omega$，$t$为外力矩，$W$相应位置填充惯性张量矩阵的转置)

式(1)扩展为：
$$JW\hat Q=-JWQ-\dot J\dot q$$
其中
$$\hat Q=J^T\lambda$$

对于经典距离约束$C(x_1,x_2)=\frac{1}{2}(x_1-x_2)\cdot(x_1-x_2)-\frac{1}{2}d^2$，雅可比矩阵为
$$J=[(\frac{\partial C}{\partial x_1})^T,(\frac{\partial C}{\partial x_2})^T]$$

对于一般的不考虑外力作用的场景而言，有：
$$JWJ^T\lambda=-\dot J\dot q_1\tag 2-k_sC-k_d\dot C$$

### Sequential Impulse约束求解框架（面向刚体）
对于式(2)，正常解矩阵方程的方法不太高效，由此产生了Sequential Impulse方法，逐个求解每个约束$C_1,C_2,...,C_m$，相当于逐个对约束所作用的两个刚体施加冲量，类似Gauss-Seidel-type迭代。由于刚体的碰撞和摩擦约束是非线性的，存在LCP（线性互补问题），所以有一些特殊处理。

求解过程：
- 计算没有约束情况下的下一帧速度$\dot q_{n+1}^*$：
  $$\dot q_{n+1}^*=\dot q_1+\Delta tWQ$$
- 对每个约束$C_k$，逐个求解出$\lambda_k$，进而得到$J_k^T\lambda_k$
- 逐个向$\dot q_{n+1}^*$施加约束，每一步得到的$\dot q_{n+1}^{(k)}$作为下一步的输入：
  $$\dot q_{n+1}^{(k+1)}=\dot q_{n+1}^{(k)}+WJ_k^T\lambda_k$$
- 最后更新位置$q_{n+1}=q_n+\Delta t\dot q_n$

### 矩阵求解方法
常见的矩阵线性求解方法有：
- 高斯消元法：通过行变换将矩阵化为上三角矩阵，然后回代求解；
- LU分解：将矩阵分解为下三角矩阵$L$和上三角矩阵$U$，然后通过前向后向替换求解；
- Cholesky分解：适用与对称正定矩阵，将矩阵分解为$LL^T$，其中$L$是下三角矩阵；
- QR分解：将矩阵分解为正交矩阵$Q$和上三角矩阵$R$，然后通过回代求解。

还有一些迭代方法：

- 雅可比迭代：对于矩阵方程$Ax=b$，其中$A$非奇异且$a_{ii}\neq 0$，将A分解为：$A=D+L+U$，其中$D$为包含元素$a_{ii}$的对角矩阵，$L$为包含元素$a_{ij},i<j$的严格下三角矩阵，$U$为包含元素$a_{ij},i>j$的严格上三角矩阵，迭代过程为：
  $$x^{(k+1)}=-D^{-1}(L+U)x^{(k)}+D^{-1}b$$
  直观理解：
  $$\begin{cases}
  8x_1-3x_2+2x_3=20 \\
  4x_1+11x_2-x_3=33 \\
  6x_1+3x_2+13x_3=36
  \end{cases}\Rightarrow
  \begin{cases}
  x_1=(20+3x_2-2x_3)/8 \\
  x_2=(33-4x_1+x_3)/11 \\
  x_3=(36-6x_1-3x_2)/12
  \end{cases}$$
  使用后一个方程组进行迭代，精确解为$x=[3,2,1]^T$，选择$x^{(0)}=[0,0,0]^T$，迭代5次有$x=[2.999843,2.000072,1.000061]^T$
- Gauss-Seidel迭代：雅可比迭代的改进，在计算$x_i^{(k+1)}$时，使用已经更新的$x_j^{(k+1)},j<i$，而不是直接使用旧的$x^{(k)}$，即：
  $$\begin{cases}
  x_1^{(k+1)}=(20+3x_2^{(k)}-2x_3^{(k)})/8 \\
  x_2^{(k+1)}=(33-4x_1^{(k+1)}+x_3^{(k)})/11 \\
  x_3^{(k+1)}=(36-6x_1^{(k+1)}-3x_2^{(k+1)})/12
  \end{cases} \tag 3$$
  即，迭代方程修改为：
  $$x^{(k+1)}=(D-L)^{-1}Ux^{(k)}+(D-L)^{-1}b$$
  高斯-赛德尔迭代收敛更快，当系数矩阵$A$严格对角占优或对称正定时，迭代必定收敛。
- Gauss-Seidel-type迭代：观察方程组(3)，每次$x_n^{(k+1)}$是独立更新的，并且后面会用到当前求出的新值，这便是Gauss-Seidel-type迭代的基本特征：独立更新每个约束。在PBD的论文中，采用这种Gauss-Seidel-type来处理位置约束，并且将约束变化投影到梯度方向上（比如两点距离约束的梯度方向指向对面的点），称为Projected Gauss Seidel

Eigen库的线性求解器
- PartialPivLU：直接求解，部分选主元，适用于一般矩阵（选主元：指选择列中绝对值最大的元素作为主元，并交换行到对角线上）
- FUllPivLU：直接求解，完全选主元，适用于一般矩阵，可计算秩和逆
- HouseholderQR：直接求解，基于Householder变换的QR分解，适用于一般矩阵
- ColPivHouseholderQR：直接求解，列选主元的QR分解，适用于亏秩矩阵
- FullPivHouseholderQR：直接求解，完全选主元的QR分解，适用于亏秩矩阵
- LLT：直接求解，Cholesky分解，适用于对称正定矩阵
- LDLT：直接求解，LDLT分解，适用于对称半正定矩阵
- ConjugateGradient：迭代求解，共轭梯度法，适用与对称正定矩阵
- BiCGSTAB：迭代求解，双共轭梯度稳定法，适用于一般矩阵
- SparseLU：迭代求解，稀疏LU分解，适用于稀疏矩阵
- SparseQR：迭代求解，稀疏QR分解，适用于稀疏矩阵
- BDCSVD：最小二乘求解，分治法奇异值分解，适用于一般矩阵和最小二乘问题
- JacobiSVD：最小二乘求解，Jacobi奇异值分解，适用于小规模矩阵的最小二乘问题
- 一般而言，小型密集矩阵用PartialPivLU或HouseholderQR，对称正定矩阵用LLT或LDLT，大规模稀疏矩阵用ConjugateGradient或BiCGSTAB，最小二乘问题用BDCSVD或JacobiSVD
