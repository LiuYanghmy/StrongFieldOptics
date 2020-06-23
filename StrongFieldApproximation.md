<font size = 3 color = grey>
强场近似模型(KFR-Approximation; L.V.Keldysh. SOVIET PHYSICS JETP, 1965[1]; F.H.M.Faisal. J. Phys. B, 1973[2]; H.R.Reiss. PHYS. REV. A, 1980[3])

以下笔记的内容主要基于Reiss那篇以及另一篇ATI(Above-Thershold-Ionization)综述文章：AB0VE-THRESHOLD IONIZATION: FROM CLASSICA L FEATURES TO QUANTUM EFFECTS. ADVANCES IN ATOMIC, MOLECULAR, AND OPTICAL PHYSICS, VOL. 48; pages:35-98[4]。

关于SFA中的规范问题，参看这篇：PHYSICAL REVIEW A 72, 023415 (2005)[5]。

跃迁振幅的求解：鞍点积分法（由于指数项快速振荡）

鞍点的物理意义在[4]中有所体现（与费恩曼路径积分的联系）。

<!------TODO：详解鞍点法(Arfken. Mathematical Methods for Physicists; Courant & Hilbert. Methods of Mathematical Physics; Riley, Hobson, and Bence. Mathematical Methods for Physics and Engineering)----->
</font>


---------------------
**TODO LIST**
- [x] 直接电子跃迁振幅
- [ ] 散射电子跃迁振幅
- [ ] 鞍点法 
- [x] 规范问题
- [ ] 失效情形
- [ ] 库仑修正
- [ ] Dyson Equation（时间演化算符的积分方程）
- [ ] Volkov解

<!------
详解鞍点法：Arfken. Mathematical Methods for Physicists; Courant & Hilbert. Methods of Mathematical Physics; Riley, Hobson, and Bence. Mathematical Methods for Physics and Engineering
库仑修正：https://doi.org/10.1080/09500340802161881
失效情形：static-field limit in the total (energy-angular integrated) ionization rate(在以上网址的文章中提及); Bashkansky, M.; Bucksbaum, P.H.; Schumacher, D.W. Phys. Rev. Lett. 1988, 60, 2458–2461 ; Rudenko, A.; Zrost, K.; Ergler, Th.; Voitkiv, A.B.; Najjari, B.; de Jesus, V.L.B.; Feuerstein, B.; Schro¨ ter, C.D.; Moshammer, R.; Ullrich, J. J. Phys. B: At. Mol. Opt. Phys. 2005, 38, L191–L198
----->
--------------------

## 引入
哈密顿量<a href="https://www.codecogs.com/eqnedit.php?latex=H=T&plus;V_F&plus;V_B" target="_blank"><img src="https://latex.codecogs.com/png.latex?H=T&plus;V_F&plus;V_B" title="H=T+V_F+V_B" /></a>，其中<a href="https://www.codecogs.com/eqnedit.php?latex=V_B" target="_blank"><img src="https://latex.codecogs.com/png.latex?V_B" title="V_B" /></a>为由核引入的束缚势，<a href="https://www.codecogs.com/eqnedit.php?latex=V_F" target="_blank"><img src="https://latex.codecogs.com/png.latex?V_F" title="V_F" /></a>表示外场的作用。

电离之前，系统初态为<a href="https://www.codecogs.com/eqnedit.php?latex=|\psi_i\rangle" target="_blank"><img src="https://latex.codecogs.com/png.latex?|\psi_i\rangle" title="|\psi_i\rangle" /></a>；电离之后，系统处于末态<a href="https://www.codecogs.com/eqnedit.php?latex=|\varPsi_f\rangle" target="_blank"><img src="https://latex.codecogs.com/png.latex?|\varPsi_f\rangle" title="|\varPsi_f\rangle" /></a>。

跃迁振幅为<a href="https://www.codecogs.com/eqnedit.php?latex=S_{fi}=\lim_{t\rightarrow-\infty}\langle\varPsi_f|\psi_i\rangle" target="_blank"><img src="https://latex.codecogs.com/png.latex?S_{fi}=\lim_{t\rightarrow-\infty}\langle\varPsi_f|\psi_i\rangle" title="S_{fi}=\lim_{t\rightarrow-\infty}\langle\varPsi_f|\psi_i\rangle" /></a>。

为求解跃迁振幅，引入两个态(<a href="https://www.codecogs.com/eqnedit.php?latex=\hbar=1" target="_blank"><img src="https://latex.codecogs.com/png.latex?\hbar=1" title="\hbar=1" /></a>)：

1. <a href="https://www.codecogs.com/eqnedit.php?latex=(i\partial_t-T-V_B)\psi_0(t)=0" target="_blank"><img src="https://latex.codecogs.com/png.latex?(i\partial_t-T-V_B)\psi_0(t)=0" title="(i\partial_t-T-V_B)\psi_0(t)=0" /></a>，即忽略外场时系统的状态；
2. <a href="https://www.codecogs.com/eqnedit.php?latex=(i\partial_t-T-V_F)\psi_f(t)=0" target="_blank"><img src="https://latex.codecogs.com/png.latex?(i\partial_t-T-V_F)\psi_f(t)=0" title="(i\partial_t-T-V_F)\psi_f(t)=0" /></a>，即仅受外场作用时电子的状态<a href="https://www.codecogs.com/eqnedit.php?latex=\psi^{Vv}_p(t)" target="_blank"><img src="https://latex.codecogs.com/png.latex?\psi^{Vv}_p(t)" title="\psi^{Vv}_p(t)" /></a>（Volkov态）。

**忽略束缚势对末态的影响**，Reiss[3]将跃迁振幅简化为

<a href="https://www.codecogs.com/eqnedit.php?latex=(S-1)_{fi}=-i\int^{(-\infty,\infty)}&space;d\tau\langle\psi_p^{Vv}(\tau)|V_F(\tau)|\psi_i(\tau)\rangle" target="_blank"><img src="https://latex.codecogs.com/png.latex?(S-1)_{fi}=-i\int^{(-\infty,\infty)}&space;d\tau\langle\psi_p^{Vv}(\tau)|V_F(\tau)|\psi_i(\tau)\rangle" title="(S-1)_{fi}=-i\int^{(-\infty,\infty)} d\tau\langle\psi_p^{Vv}(\tau)|V_F(\tau)|\psi_i(\tau)\rangle" /></a>

该近似的有效性参看Reiss[3]；其中最为关键的是被忽略的部分比上简化表达式的相对大小：

<a href="https://www.codecogs.com/eqnedit.php?latex=\sum_j(-i\int&space;d\tau\langle\psi_{f}|V_B|\psi_{fj}\rangle)\ll1" target="_blank"><img src="https://latex.codecogs.com/png.latex?\sum_j(-i\int&space;d\tau\langle\psi_{f}|V_B|\psi_{fj}\rangle)\ll1" title="\sum_j(-i\int d\tau\langle\psi_{f}|V_B|\psi_{fj}\rangle)\ll1" /></a>

其中j为所有可能的末态指标。

在有限范围势场<a href="https://www.codecogs.com/eqnedit.php?latex=V_B" target="_blank"><img src="https://latex.codecogs.com/png.latex?V_B" title="V_B" /></a>的情况下（如：负离子电离，因而核呈中性），上式可以粗略估计在<a href="https://www.codecogs.com/eqnedit.php?latex=\sum_j&space;10^{-12}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\sum_j&space;10^{-12}" title="\sum_j 10^{-12}" /></a>的量级。

**在电离之前，忽略外场的作用**，则初态可以用<a href="https://www.codecogs.com/eqnedit.php?latex=\psi_0(\tau)" target="_blank"><img src="https://latex.codecogs.com/png.latex?\psi_0(\tau)" title="\psi_0(\tau)" /></a>代替；上式即描述直接电离的SFA跃迁振幅表达式。

## 另一种导出方法
<font size = 3 color = grey>
Reiss[3]详细讨论了采用SFA带来的误差，但不免有些繁琐。

下面的方法使用时间演化算符，个人觉得推导过程简单直观，易于理解。
</font>

跃迁振幅<a href="https://www.codecogs.com/eqnedit.php?latex=M_p=\lim_{t\rightarrow\infty,t'\rightarrow-\infty}\langle\psi_p(t)|U(t,t')|\psi_0(t')\rangle" target="_blank"><img src="https://latex.codecogs.com/png.latex?M_p=\lim_{t\rightarrow\infty,t'\rightarrow-\infty}\langle\psi_p(t)|U(t,t')|\psi_0(t')\rangle" title="M_p=\lim_{t\rightarrow\infty,t'\rightarrow-\infty}\langle\psi_p(t)|U(t,t')|\psi_0(t')\rangle" /></a>（与上一节的S矩阵元<a href="https://www.codecogs.com/eqnedit.php?latex=S_{fi}" target="_blank"><img src="https://latex.codecogs.com/png.latex?S_{fi}" title="S_{fi}" /></a>相同，只是显式地标出末态的动量p）

这里引入了时间演化算符，因此将哈密顿量写为含时和不含时的两部分：

<a href="https://www.codecogs.com/eqnedit.php?latex=H=H_0&plus;H_I" target="_blank"><img src="https://latex.codecogs.com/png.latex?H=H_0&plus;H_I" title="H=H_0+H_I" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=U,&space;U_0" target="_blank"><img src="https://latex.codecogs.com/png.latex?U,&space;U_0" title="U, U_0" /></a>分别对应<a href="https://www.codecogs.com/eqnedit.php?latex=H,&space;H_0" target="_blank"><img src="https://latex.codecogs.com/png.latex?H,&space;H_0" title="H, H_0" /></a>的时间演化算符，它们满足积分方程

<a href="https://www.codecogs.com/eqnedit.php?latex=U(t,&space;t_0)&space;=&space;U_0(t,t_0)&space;-&space;i\int_{t_0}^t&space;d\tau&space;U(t,\tau)H_I(\tau)U_0(\tau,t_0)" target="_blank"><img src="https://latex.codecogs.com/png.latex?U(t,&space;t_0)&space;=&space;U_0(t,t_0)&space;-&space;i\int_{t_0}^t&space;d\tau&space;U(t,\tau)H_I(\tau)U_0(\tau,t_0)" title="U(t, t_0) = U_0(t,t_0) - i\int_{t_0}^t d\tau U(t,\tau)H_I(\tau)U_0(\tau,t_0)" /></a>

代入跃迁振幅 ，并注意束缚态<a href="https://www.codecogs.com/eqnedit.php?latex=\psi_0(t)" target="_blank"><img src="https://latex.codecogs.com/png.latex?\psi_0(t)" title="\psi_0(t)" /></a>与散射态<a href="https://www.codecogs.com/eqnedit.php?latex=\psi_p(t)" target="_blank"><img src="https://latex.codecogs.com/png.latex?\psi_p(t)" title="\psi_p(t)" /></a>正交（电离前后波函数没有重叠），得到

<a href="https://www.codecogs.com/eqnedit.php?latex=M_p=-i\lim_{t\rightarrow\infty}\int_{-\infty}^t&space;d\tau\langle\psi_p(t)|U(t,\tau)H_I|\psi_0(\tau)\rangle" target="_blank"><img src="https://latex.codecogs.com/png.latex?M_p=-i\lim_{t\rightarrow\infty}\int_{-\infty}^t&space;d\tau\langle\psi_p(t)|U(t,\tau)H_I|\psi_0(\tau)\rangle" title="M_p=-i\lim_{t\rightarrow\infty}\int_{-\infty}^t d\tau\langle\psi_p(t)|U(t,\tau)H_I|\psi_0(\tau)\rangle" /></a>

**引入强场近似**，即将末态代以Volkov态，并认为：<a href="https://www.codecogs.com/eqnedit.php?latex=\psi_p^{Vv}(t)U(t,t_0)\rightarrow\psi_p^{Vv}(t_0)" target="_blank"><img src="https://latex.codecogs.com/png.latex?\psi_p^{Vv}(t)U(t,t_0)\rightarrow\psi_p^{Vv}(t_0)" title="\psi_p^{Vv}(t)U(t,t_0)\rightarrow\psi_p^{Vv}(t_0)" /></a>，得到

<a href="https://www.codecogs.com/eqnedit.php?latex=M_p=-i\int_{-\infty}^\infty&space;d\tau\langle\psi_p^{Vv}(\tau)|H_I(\tau)|\psi_0(\tau)\rangle" target="_blank"><img src="https://latex.codecogs.com/png.latex?M_p=-i\int_{-\infty}^\infty&space;d\tau\langle\psi_p^{Vv}(\tau)|H_I(\tau)|\psi_0(\tau)\rangle" title="M_p=-i\int_{-\infty}^\infty d\tau\langle\psi_p^{Vv}(\tau)|H_I(\tau)|\psi_0(\tau)\rangle" /></a>

上式是直接电离的跃迁振幅。
## 鞍点法(saddle point method)
以长度规范为例，<a href="https://www.codecogs.com/eqnedit.php?latex=H_I=-erE(t)" target="_blank"><img src="https://latex.codecogs.com/png.latex?H_I=-erE(t)" title="H_I=-erE(t)" /></a>

Volkov波函数形式为<a href="https://www.codecogs.com/eqnedit.php?latex=|\psi_p^{Vv}(t)\rangle&space;=&space;|p-eA(t)\rangle&space;e^{-iS_p(t)}" target="_blank"><img src="https://latex.codecogs.com/png.latex?|\psi_p^{Vv}(t)\rangle&space;=&space;|p-eA(t)\rangle&space;e^{-iS_p(t)}" title="|\psi_p^{Vv}(t)\rangle = |p-eA(t)\rangle e^{-iS_p(t)}" /></a>，其中<a href="https://www.codecogs.com/eqnedit.php?latex=S_p(t)=\frac{1}{2m}\int^t&space;d\tau[p-eA(\tau)]^2" target="_blank"><img src="https://latex.codecogs.com/png.latex?S_p(t)=\frac{1}{2m}\int^t&space;d\tau[p-eA(\tau)]^2" title="S_p(t)=\frac{1}{2m}\int^t d\tau[p-eA(\tau)]^2" /></a>；

初态<a href="https://www.codecogs.com/eqnedit.php?latex=|\psi_0(t)\rangle=e^{iE_{IP}t}|\psi_0\rangle" target="_blank"><img src="https://latex.codecogs.com/png.latex?|\psi_0(t)\rangle=e^{iE_{IP}t}|\psi_0\rangle" title="|\psi_0(t)\rangle=e^{iE_{IP}t}|\psi_0\rangle" /></a>。

<a href="https://www.codecogs.com/eqnedit.php?latex=M_p=-i\int_{-\infty}^\infty&space;d\tau\langle&space;p-eA(\tau)|H_I(\tau)|\psi_0\rangle&space;e^{i[S_p(\tau)&plus;E_{IP}\tau]}" target="_blank"><img src="https://latex.codecogs.com/png.latex?M_p=-i\int_{-\infty}^\infty&space;d\tau\langle&space;p-eA(\tau)|H_I(\tau)|\psi_0\rangle&space;e^{i[S_p(\tau)&plus;E_{IP}\tau]}" title="M_p=-i\int_{-\infty}^\infty d\tau\langle p-eA(\tau)|H_I(\tau)|\psi_0\rangle e^{i[S_p(\tau)+E_{IP}\tau]}" /></a>

对于强外场，指数因子<a href="https://www.codecogs.com/eqnedit.php?latex=g(t)=i[S_p(\tau)&plus;E_{IP}\tau]" target="_blank"><img src="https://latex.codecogs.com/png.latex?g(t)=i[S_p(\tau)&plus;E_{IP}\tau]" title="g(t)=i[S_p(\tau)+E_{IP}\tau]" /></a>沿积分路径快速变化，因此可以采用鞍点近似

<a href="https://www.codecogs.com/eqnedit.php?latex=M_p=-i\sum_{t_s}e^{i\theta}\sqrt{\frac{2\pi}{g''(t_s)}}\langle&space;p-eA(t_s)|H_I(t_s)|\psi_0\rangle&space;e^{g(t_s)}" target="_blank"><img src="https://latex.codecogs.com/png.latex?M_p=-i\sum_{t_s}e^{i\theta}\sqrt{\frac{2\pi}{g''(t_s)}}\langle&space;p-eA(t_s)|H_I(t_s)|\psi_0\rangle&space;e^{g(t_s)}" title="M_p=-i\sum_{t_s}e^{i\theta}\sqrt{\frac{2\pi}{g''(t_s)}}\langle p-eA(t_s)|H_I(t_s)|\psi_0\rangle e^{g(t_s)}" /></a>

其中鞍点<a href="https://www.codecogs.com/eqnedit.php?latex=t_s" target="_blank"><img src="https://latex.codecogs.com/png.latex?t_s" title="t_s" /></a>由方程<a href="https://www.codecogs.com/eqnedit.php?latex=g'(t_s)=0\Rightarrow&space;E_{IP}&plus;\frac{1}{2}[p-eA(t_s)]^2=0" target="_blank"><img src="https://latex.codecogs.com/png.latex?g'(t_s)=0\Rightarrow&space;E_{IP}&plus;\frac{1}{2}[p-eA(t_s)]^2=0" title="g'(t_s)=0\Rightarrow E_{IP}+\frac{1}{2}[p-eA(t_s)]^2=0" /></a>确定，最速下降方向由<a href="https://www.codecogs.com/eqnedit.php?latex=\theta&space;=&space;-\frac{\arg(g''(t_s))}{2}&plus;\frac{\pi}{2}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\theta&space;=&space;-\frac{\arg(g''(t_s))}{2}&plus;\frac{\pi}{2}" title="\theta = -\frac{\arg(g''(t_s))}{2}+\frac{\pi}{2}" /></a>确定。

## 规范
<font size = 3 color = grey>
虽然强场近似在一定程度上取得了很好的效果，但忽略束缚势仍然导致了一些问题

其一就是在强场近似下，不再具有规范不变性。
</font>

不同规范下，相互作用能和Volkov波函数的表达式有所不同

- 长度规范：<a href="https://www.codecogs.com/eqnedit.php?latex=H_I=-erE(t)\qquad\&space;|\psi_p^{Vv}(t)\rangle=|p\rangle&space;e^{-iS_p(t)}" target="_blank"><img src="https://latex.codecogs.com/png.latex?H_I=-erE(t)\qquad\&space;|\psi_p^{Vv}(t)\rangle=|p\rangle&space;e^{-iS_p(t)}" title="H_I=-erE(t)\qquad\ |\psi_p^{Vv}(t)\rangle=|p\rangle e^{-iS_p(t)}" /></a>
- 速度规范：<a href="https://www.codecogs.com/eqnedit.php?latex=H_I=-\frac{e}{m}pA(t)\qquad|\psi_p^{Vv}(t)\rangle=|p-eA(t)\rangle&space;e^{-iS_p(t)}" target="_blank"><img src="https://latex.codecogs.com/png.latex?H_I=-\frac{e}{m}pA(t)\qquad|\psi_p^{Vv}(t)\rangle=|p-eA(t)\rangle&space;e^{-iS_p(t)}" title="H_I=-\frac{e}{m}pA(t)\qquad|\psi_p^{Vv}(t)\rangle=|p-eA(t)\rangle e^{-iS_p(t)}" /></a>

强场近似的跃迁振幅可以通过分部积分(?)写成以下形式

<a href="https://www.codecogs.com/eqnedit.php?latex=M_p=-i\int_{-\infty}^\infty&space;d\tau\langle\psi_p^{Vv}(\tau)|V_B(r)|\psi_0(\tau)\rangle" target="_blank"><img src="https://latex.codecogs.com/png.latex?M_p=-i\int_{-\infty}^\infty&space;d\tau\langle\psi_p^{Vv}(\tau)|V_B(r)|\psi_0(\tau)\rangle" title="M_p=-i\int_{-\infty}^\infty d\tau\langle\psi_p^{Vv}(\tau)|V_B(r)|\psi_0(\tau)\rangle" /></a>

这样，不同规范下就只有Volkov波函数的差别；采用鞍点法，这一差别将体现在指数前因子中：

- 长度规范：<a href="https://www.codecogs.com/eqnedit.php?latex=\langle&space;p-eA(t_s)|V_B(r)|\psi_0\rangle" target="_blank"><img src="https://latex.codecogs.com/png.latex?\langle&space;p-eA(t_s)|V_B(r)|\psi_0\rangle" title="\langle p-eA(t_s)|V_B(r)|\psi_0\rangle" /></a>

- 速度规范：<a href="https://www.codecogs.com/eqnedit.php?latex=\langle&space;p|V_B(r)|\psi_0\rangle" target="_blank"><img src="https://latex.codecogs.com/png.latex?\langle&space;p|V_B(r)|\psi_0\rangle" title="\langle p|V_B(r)|\psi_0\rangle" /></a>

[5]中给出了一个例子：外场A为单色线偏场，在一个周期内，鞍点方程有两个解；对于长度规范，<a href="https://www.codecogs.com/eqnedit.php?latex=p-eA(t_s)" target="_blank"><img src="https://latex.codecogs.com/png.latex?p-eA(t_s)" title="p-eA(t_s)" /></a>在这两个解处沿外场方向的分量符号相反，取决于初态<a href="https://www.codecogs.com/eqnedit.php?latex=\psi_0" target="_blank"><img src="https://latex.codecogs.com/png.latex?\psi_0" title="\psi_0" /></a>的宇称，这两个鞍点的指数前因子将同号或反号，它们的和将相长或相消；而对于速度规范，不同鞍点的和总是相长的。

对于有限范围势场，采用长度规范可能会更好。
## 库仑修正
<!---------------[S.V. Popruzhenko & D. Bauer. Journal of Modern Optics 55:16(2008); pages: 2573-2589 ](https://doi.org/10.1080/09500340802161881)--------->


## 附录A Dyson Integral equation

哈密顿量<a href="https://www.codecogs.com/eqnedit.php?latex=H=T&plus;V(r)&plus;H_I(t)" target="_blank"><img src="https://latex.codecogs.com/png.latex?H=T&plus;V(r)&plus;H_I(t)" title="H=T+V(r)+H_I(t)" /></a>，对应的时间演化算符为<a href="https://www.codecogs.com/eqnedit.php?latex=U" target="_blank"><img src="https://latex.codecogs.com/png.latex?U" title="U" /></a>

哈密顿量<a href="https://www.codecogs.com/eqnedit.php?latex=H_a=T&plus;V(r)" target="_blank"><img src="https://latex.codecogs.com/png.latex?H_a=T&plus;V(r)" title="H_a=T+V(r)" /></a>，对应的时间演化算符为<a href="https://www.codecogs.com/eqnedit.php?latex=U_a" target="_blank"><img src="https://latex.codecogs.com/png.latex?U_a" title="U_a" /></a>

哈密顿量<a href="https://www.codecogs.com/eqnedit.php?latex=H_f=T&plus;H_I(t)" target="_blank"><img src="https://latex.codecogs.com/png.latex?H_f=T&plus;H_I(t)" title="H_f=T+H_I(t)" /></a>，对应的时间演化算符为<a href="https://www.codecogs.com/eqnedit.php?latex=U_f" target="_blank"><img src="https://latex.codecogs.com/png.latex?U_f" title="U_f" /></a>

Dyson Equations 

<a href="https://www.codecogs.com/eqnedit.php?latex=U(t,t')=U_a(t,t')-i\int_{t'}^tU(t,\tau)H_I(\tau)U_a(\tau,t')" target="_blank"><img src="https://latex.codecogs.com/png.latex?U(t,t')=U_a(t,t')-i\int_{t'}^tU(t,\tau)H_I(\tau)U_a(\tau,t')" title="U(t,t')=U_a(t,t')-i\int_{t'}^tU(t,\tau)H_I(\tau)U_a(\tau,t')" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=U(t,t')=U_f(t,t')-i\int_{t'}^tU_f(t,\tau)V(r)U(\tau,t')" target="_blank"><img src="https://latex.codecogs.com/png.latex?U(t,t')=U_f(t,t')-i\int_{t'}^tU_f(t,\tau)V(r)U(\tau,t')" title="U(t,t')=U_f(t,t')-i\int_{t'}^tU_f(t,\tau)V(r)U(\tau,t')" /></a>

## 附录B Volkov态

## 附录C L-guage & V-guage
电场作用下的哈密顿量<a href="https://www.codecogs.com/eqnedit.php?latex=H=\frac{[p-eA(t)]^2}{2m}&plus;V(r)=\frac{p^2}{2m}&plus;V(r)-\frac{e}{m}pA(t)&plus;\frac{e^2}{2m}A^2(t)" target="_blank"><img src="https://latex.codecogs.com/png.latex?H=\frac{[p-eA(t)]^2}{2m}&plus;V(r)=\frac{p^2}{2m}&plus;V(r)-\frac{e}{m}pA(t)&plus;\frac{e^2}{2m}A^2(t)" title="H=\frac{[p-eA(t)]^2}{2m}+V(r)=\frac{p^2}{2m}+V(r)-\frac{e}{m}pA(t)+\frac{e^2}{2m}A^2(t)" /></a>

薛定谔方程<a href="https://www.codecogs.com/eqnedit.php?latex=i\partial_t\psi&space;=&space;H\psi" target="_blank"><img src="https://latex.codecogs.com/png.latex?i\partial_t\psi&space;=&space;H\psi" title="i\partial_t\psi = H\psi" /></a>

插入相位因子<a href="https://www.codecogs.com/eqnedit.php?latex=U=e^{i(\frac{e^2}{2m}\int^t&space;d\tau&space;A^2(\tau))}" target="_blank"><img src="https://latex.codecogs.com/png.latex?U=e^{i(\frac{e^2}{2m}\int^t&space;d\tau&space;A^2(\tau))}" title="U=e^{i(\frac{e^2}{2m}\int^t d\tau A^2(\tau))}" /></a>，<a href="https://www.codecogs.com/eqnedit.php?latex=i\partial_t&space;U^*U\psi&space;=&space;H&space;U^*U\psi" target="_blank"><img src="https://latex.codecogs.com/png.latex?i\partial_t&space;U^*U\psi&space;=&space;H&space;U^*U\psi" title="i\partial_t U^*U\psi = H U^*U\psi" /></a>，并视<a href="https://www.codecogs.com/eqnedit.php?latex=U\psi" target="_blank"><img src="https://latex.codecogs.com/png.latex?U\psi" title="U\psi" /></a>为新的波函数<a href="https://www.codecogs.com/eqnedit.php?latex=\psi'" target="_blank"><img src="https://latex.codecogs.com/png.latex?\psi'" title="\psi'" /></a>。整理得到<a href="https://www.codecogs.com/eqnedit.php?latex=\psi'" target="_blank"><img src="https://latex.codecogs.com/png.latex?\psi'" title="\psi'" /></a>满足的薛定谔方程

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{e^2}{2m}A^2(t)U^*&space;i\partial_t\psi'&space;&plus;&space;U^*&space;i\partial_t\psi'&space;=&space;HU^*&space;\psi'\Rightarrow&space;i\partial_t\psi'&space;=&space;U(H-\frac{e^2}{2m}A^2(t))U^*\psi'" target="_blank"><img src="https://latex.codecogs.com/png.latex?\frac{e^2}{2m}A^2(t)U^*&space;i\partial_t\psi'&space;&plus;&space;U^*&space;i\partial_t\psi'&space;=&space;HU^*&space;\psi'\Rightarrow&space;i\partial_t\psi'&space;=&space;U(H-\frac{e^2}{2m}A^2(t))U^*\psi'" title="\frac{e^2}{2m}A^2(t)U^* i\partial_t\psi' + U^* i\partial_t\psi' = HU^* \psi'\Rightarrow i\partial_t\psi' = U(H-\frac{e^2}{2m}A^2(t))U^*\psi'" /></a>

变换后，相互作用项<a href="https://www.codecogs.com/eqnedit.php?latex=H_I=-\frac{e}{m}pA(t)" target="_blank"><img src="https://latex.codecogs.com/png.latex?H_I=-\frac{e}{m}pA(t)" title="H_I=-\frac{e}{m}pA(t)" /></a>。

而若是插入相位因子<a href="https://www.codecogs.com/eqnedit.php?latex=U=e^{i(-erA(t)&plus;\frac{e^2}{2m}\int^t&space;d\tau&space;A^2(\tau))}" target="_blank"><img src="https://latex.codecogs.com/png.latex?U=e^{i(-erA(t)&plus;\frac{e^2}{2m}\int^t&space;d\tau&space;A^2(\tau))}" title="U=e^{i(-erA(t)+\frac{e^2}{2m}\int^t d\tau A^2(\tau))}" /></a>，同样地，可以得到变换后的相互作用项

<a href="https://www.codecogs.com/eqnedit.php?latex=H_I=-erE(t)" target="_blank"><img src="https://latex.codecogs.com/png.latex?H_I=-erE(t)" title="H_I=-erE(t)" /></a>

注意两种变换的波函数仅相差一个相位因子<a href="https://www.codecogs.com/eqnedit.php?latex=e^{ierA(t)}" target="_blank"><img src="https://latex.codecogs.com/png.latex?e^{ierA(t)}" title="e^{ierA(t)}" /></a>。
