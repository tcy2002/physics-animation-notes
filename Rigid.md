# 刚体动力学（非粒子）
## 刚体属性 
相比于普通粒子或质点约束，刚体有体积和取向，除了3个位置自由度$x,y,z$外，还多了3个旋转自由度$w_x,w_y,w_z$. 另外，为描述刚体的有体积质量，需要用到惯性张量矩阵$I\in R^{3\times 3}$，$I$对角线的值即刚体相对于坐标轴的转动惯量，刚体相对于任意轴$n$的转动惯量为$n^TIn$.
一般形状相对于质心的惯性张量矩阵为：
$$
\begin{matrix}
\int y^2+z^2dm & -\int xydm & -\int xzdm\\
-\int xydm & \int x^2+z^2dm & -\int yzdm\\
-\int xzdm & -\int yzdm & \int x^2+y^2dm 
\end{matrix}
$$
对于中心对称的几何体，除对角线外，矩阵的其他位置均为0. 一些常见形状相对于质心的惯性张量矩阵如下（质量为1）
- 盒体：$extent = (a, b, c)$
  $$
  \begin{matrix}
    \frac{b^2+c^2}{12} & 0 & 0\\
    0 & \frac{a^2+c^2}{12} & 0\\
    0 & 0 & \frac{a^2+b^2}{12}
  \end{matrix}
  $$
- 圆柱体：$radius = r$，$height = h$
  $$
  \begin{matrix}
    \frac{3r^2+h^2}{12} & 0 & 0\\
    0 & \frac{r^2}{2} & 0\\
    0 & 0 & \frac{3r^2+h^2}{12}
  \end{matrix}
  $$
- 球体：$radius = r$
  $$
  \begin{matrix}
    \frac{2}{5}r^2 & 0 & 0\\
    0 & \frac{2}{5}r^2 & 0\\
    0 & 0 & \frac{2}{5}r^2
  \end{matrix}
  $$
- 椭球体：$extent = (a, b, c)$
  $$
  \begin{matrix}
    \frac{1}{5}(b^2+c^2) & 0 & 0\\
    0 & \frac{1}{5}(a^2+c^2) & 0\\
    0 & 0 & \frac{1}{5}(a^2+b^2)
  \end{matrix}
  $$
- 半球体：$radius = r$（底面平行于$xOz$，质心到底面距离为$\frac{3}{8}r$）
  $$
  \begin{matrix}
    \frac{2}{5}r^2-(\frac{3}{8}r)^2 & 0 & 0\\
    0 & \frac{2}{5}r^2 & 0\\
    0 & 0 & \frac{2}{5}r^2-(\frac{3}{8}r)^2
  \end{matrix}
  $$
- 胶囊体，游戏物理中经常用到，因为碰撞检测很方便：$radius = r$，$height = h$（根据半球体和圆柱体叠加推导）
  $$
  \begin{matrix}
    I_{xx,semi-sphere}+I_{xx,cylinder} & 0 & 0\\
    0 & I_{yy,semi-sphere}+I_{yy,cylinder} & 0\\
    0 & 0 & I_{xx,semi-sphere}+I_{xx,cylinder}
  \end{matrix}
  $$
  $$I_{xx,semi-sphere}=I_{zz,semi-sphere}=(\frac{2}{5}r^2-(\frac{3}{8}r)^2+(\frac{1}{2}h+\frac{3}{8}r)^2)\cdot \frac{4r}{4r+3h}$$
  $$I_{xx,cylinder}=I_{zz,cylinder}=\frac{3r^2+h^2}{12}\cdot \frac{3h}{4r+3h}$$
  $$I_{yy,semi-sphere}=\frac{2}{5}r^2\cdot \frac{4r}{4r+3h}$$
  $$I_{yy,cylinder}=\frac{r^2}{2}\cdot \frac{3h}{4r+3h}$$
- 圆锥：$height = h$, $radius = r$
  $$
  \begin{matrix}
    \frac{3r^2+2h^2}{20} & 0 & 0\\
    0 & \frac{3}{10}r^2 & 0\\
    0 & 0 & \frac{3r^2+2h^2}{20}
  \end{matrix}
  $$
- 甜甜圈：环半径为$R$，截面半径为$r$
  $$
  \begin{matrix}
    \frac{5}{8}r^2+\frac{1}{2}R^2 & 0 & 0\\
    0 & \frac{3}{4}r^2+R^2 & 0\\
    0 & 0 & \frac{5}{8}r^2+\frac{1}{2}R^2
  \end{matrix}
  $$
- 多面体网格：推导比较复杂，[参考链接](https://www.geometrictools.com/Documentation/PolyhedralMassProperties.pdf)，[代码实现](https://github.com/tcy2002/Physics-Engine/blob/main/src/phys/shape/convex_mesh_shape.cpp)，也适用于非凸刚体。其实一般不需要精确计算多面体网格的惯性张量矩阵，在bullet、react physics中，均使用AABB包围盒来估计惯性张量；react physics精确计算了胶囊体的惯性张量，但bullet也直接用AABB包围盒来估计。

刚体满足的动力学条件：
$$F=M\dot v, T=I\dot \omega$$
为方便统一，在数值求解器中，通常使用整体形式：
$$Q=M\dot u$$
这里的$M\in R^{6\times 6}$，$v\in R^{6\times 1}$：
$$M=\left[
  \begin{matrix}
  \frac{1}{m} & 0 & 0 & 0 & 0 & 0\\
  0 & \frac{1}{m} & 0 & 0 & 0 & 0\\
  0 & 0 & \frac{1}{m} &  & 0 & 0\\
  0 & 0 & 0 & I_{xx} & I_{xy} & I_{xz}\\
  0 & 0 & 0 & I_{xy} & I_{yy} & I_{yz}\\
  0 & 0 & 0 & I_{xz} & I_{yz} & I_{zz}
  \end{matrix}
\right], u=\left[
  \begin{matrix}
  \dot v_x\\
  \dot v_y\\
  \dot v_z\\
  \dot \omega_x\\
  \dot \omega_y\\
  \dot \omega_z\\
  \end{matrix}
\right]$$
对于有$n$个刚体的物理场景，存在$6n$个自由度，可以用$6n$维的的稀疏矩阵来表示。

## 刚体约束动力学
### naive：基于冲量（速度）的约束

### LCP：线性互补条件&库伦摩擦模型

### 无矩阵的Sequential Impulse+Gauss Seidel迭代

### 一些高精度数值方法（优化问题）
- Rigid-IPC
- ABD
- Primal-Dual

## 主流刚体引擎架构（碰撞解算和响应是核心）
### 宽域碰撞

### 窄域碰撞

### 约束求解

### 引擎的运行逻辑：每一帧的执行顺序

### 内存池、多线程

## 拓展：充分发挥硬件性能
### simd

### 定点数

### gpu
