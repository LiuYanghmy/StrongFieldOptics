<font size = 3 color = grey>
强场近似模型(also called KFR-Approximation; L.V.Keldysh. SOVIET PHYSICS JETP, 1965[1]; F.H.M.Faisal. J. Phys. B, 1973[2]; H.R.Reiss. PHYS. REV. A, 1980[3])

以下笔记的内容主要基于Reiss那篇以及另一篇ATI(Above-Thershold-Ionization)综述文章：AB0VE-THRESHOLD IONIZATION: FROM CLASSICA L FEATURES TO QUANTUM EFFECTS. ADVANCES IN ATOMIC, MOLECULAR, AND OPTICAL PHYSICS, VOL. 48; pages:35-98[4]。

关于SFA中的规范问题，参看这篇：PHYSICAL REVIEW A 72, 023415 (2005)[5]。

跃迁振幅的求解：鞍点积分法（由于指数项快速振荡）

鞍点的物理意义在[4]中有所体现（与费恩曼路径积分的联系）。

<!------TODO：详解鞍点法(Arfken. Mathematical Methods for Physicists; Courant & Hilbert. Methods of Mathematical Physics; Riley, Hobson, and Bence. Mathematical Methods for Physics and Engineering)----->
</font>


---------------------
**TODO LIST**
- [x] 直接电子跃迁振幅
- [x] 散射电子跃迁振幅
- [ ] 鞍点法 
- [x] 规范问题
- [ ] 失效情形
- [x] 库仑修正

<!------
详解鞍点法：Arfken. Mathematical Methods for Physicists; Courant & Hilbert. Methods of Mathematical Physics; Riley, Hobson, and Bence. Mathematical Methods for Physics and Engineering
库仑修正：https://doi.org/10.1080/09500340802161881
失效情形：static-field limit in the total (energy-angular integrated) ionization rate(在以上网址的文章中提及); Bashkansky, M.; Bucksbaum, P.H.; Schumacher, D.W. Phys. Rev. Lett. 1988, 60, 2458–2461 ; Rudenko, A.; Zrost, K.; Ergler, Th.; Voitkiv, A.B.; Najjari, B.; de Jesus, V.L.B.; Feuerstein, B.; Schro¨ ter, C.D.; Moshammer, R.; Ullrich, J. J. Phys. B: At. Mol. Opt. Phys. 2005, 38, L191–L198
----->
--------------------

## 引入

考虑电子在外场作用下发生电离的过程。通常，单光子的能量不足以使电子跃迁至连续态，但在强外场下，电子所处的势场受到了明显的改变，这使得电子吸收多光子从而电离成为可能。

哈密顿量$H=T+V_F+V_B$，其中$V_B$为由核引入的束缚势，$V_F$表示外场的作用。

电离之前，系统初态为$|\psi_i\rangle$；电离之后，系统处于末态$|\varPsi_f\rangle$。

跃迁振幅为$S_{fi}=\lim_{t\rightarrow-\infty}\langle\varPsi_f|\psi_i\rangle$。

为求解跃迁振幅，引入两个态($\hbar=1$)：

1. $(i\partial_t-T-V_B)\psi_0(t)=0$，即忽略外场时系统的状态；
2. $(i\partial_t-T-V_F)\psi_f(t)=0$，即仅受外场作用时电子的状态$\psi^{Vv}_p(t)$（Volkov态）。

**忽略束缚势$V_B$对末态的影响**，Reiss[3]将跃迁振幅简化为

$(S-1)_{fi}=-i\int^{(-\infty,\infty)} d\tau\langle\psi_p^{Vv}(\tau)|V_F(\tau)|\psi_i(\tau)\rangle$

该近似的有效性参看Reiss[3]；其中最为关键的是被忽略的部分比上简化表达式的相对大小：

$\sum_j(-i\int d\tau\langle\psi_{f}|V_B|\psi_{fj}\rangle)\ll1$

其中$j$为所有可能的末态指标。

在有限范围势场$V_B$的情况下（如：负离子电离，因而核呈中性），上式可以粗略估计在$\sum_j 10^{-12}$的量级。

**在电离之前，忽略外场的作用**，则初态可以用$\psi_0(\tau)$代替；此即描述直接电离的SFA跃迁振幅表达式。

-------------------------------------

在单色场的情况下，对$A(t)$的积分将贡献一个正弦型函数，表达式（代入初末态具体表达式后的跃迁振幅）中形如$e^{i\sin\omega t}$的项可以用贝塞尔函数展开，剩余的指数因子在积分后变为$\delta$函数，具体来说，具有以下形式

$S_{fi}\sim s_{fi}\sum_n \alpha_nJ_n\delta(\frac{p^2}{2m}-n\omega+z\omega+E_{IP})$

其中$E_{IP}$是电离能，$z=\frac{e^2A^2}{4m\omega}$是自由电子在外场中的最小平均能量。

因此，对光电子跃迁概率有主要贡献的是具有动能$\frac{p^2}{2m}=n\omega-z\omega-E_{IP}$的这些电子。这一结论具有很清晰的物理图景：吸收n个光子并克服电离能电离，再扣除与电场的作用能后余下的能量即光电子动能。

------------------------------------

实际情形下，该近似受到更多的限制，例如：

1. 样品损耗的问题。在总跃迁概率为$W$下（简便起见，认为它不随时间变化），$N_0$的原子中有$N=N_0e^{-iWt}$的原子剩余，则光电子比率为$Y=\frac{N_0-N}{N}=1-e^{-iWt}$；为了使$Y$足够小，从而$Y\approx Wt$以便实验观测（脉冲时间t总是非常小的，因此$Y$随t的变化难以准确测量），需要有$Wt\ll1$；

2. 采用散射矩阵的描述是否合适。实际过程中，电子从$t=0$开始演化，到$t$时刻结束演化，这时$t\rightarrow\infty$的条件不再满足；为了使有限演化时间带来的影响足够小，且最终跃迁振幅仍有类似上述$\delta$函数的行为，Reiss[3]给出了以下限制：$\omega t\gg1$。


### 另一种导出方法
<font size = 3 color = grey>
Reiss[3]详细讨论了采用SFA带来的误差。

下面的方法本质上并无不同，但使用时间演化算符表示，个人觉得物理图像简单直观，易于理解。
</font>

跃迁振幅$M_p=\lim_{t\rightarrow\infty,t'\rightarrow-\infty}\langle\psi_p(t)|U(t,t')|\psi_0(t')\rangle$（与上一节的S矩阵元$S_{fi}$相同，只是显式地标出末态的动量$p$）

哈密顿量：$H=H_0+H_I$，其中$H_0=T+V_B$，$H_I=V_F$表示与外场的相互作用。

$U, U_0$分别对应$H, H_0$的时间演化算符，它们满足积分方程

$U(t, t_0) = U_0(t,t_0) - i\int_{t_0}^t d\tau U(t,\tau)H_I(\tau)U_0(\tau,t_0)$

代入跃迁振幅 ，并注意束缚态$\psi_0(t)$与散射态$\psi_p(t)$正交（电离前后波函数没有重叠），得到

$M_p=-i\lim_{t\rightarrow\infty}\int_{-\infty}^t d\tau\langle\psi_p(t)|U(t,\tau)H_I|\psi_0(\tau)\rangle$

**引入强场近似**，即将末态代以Volkov态，并有：$\psi_p^{Vv}(t)U(t,t_0)\rightarrow\psi_p^{Vv}(t_0)$，得到

$M_p=-i\int_{-\infty}^\infty d\tau\langle\psi_p^{Vv}(\tau)|H_I(\tau)|\psi_0(\tau)\rangle$

上式是直接电离的跃迁振幅。

这一表达式的物理图像还是比较清晰的：电子在$\tau$时刻电离，在外场作用下演化，直到末态，末态在具有动量$p$的Volkov态的投影对应被积函数；对所有时刻电离的电子求和，得到（具有动量$p$的）总跃迁振幅。

将跃迁振幅平方，即可得到动量谱；对动量谱积分可导出能量分布或角度分布。
### 散射电子
散射电子在$t_0$时刻电离，随后在外场作用下运动，并在$t_1$时刻返回束缚势，发生散射。

上一节中，对代入积分方程后得到的表达式，如果再次代入积分方程（形式略有不同，见附录B），将得到两项

$M_p^1=-i\lim_{t\rightarrow\infty}\int_{-\infty}^t dt_0\langle\psi_p(t)|U_f(t,t_0)H_I|\psi_0(t_0)\rangle$

$M_p^2=-\lim_{t\rightarrow\infty}\int_{t_0}^{t}dt_1\int_{-\infty}^tdt_0\langle\psi_p(t)|U_f(t,t_1)V(r)U(t_1,t_0)H_I(t_0)|\psi_0(t_0)\rangle$

其中第一项即直接电离，第二项为经过散射后的电离。上一节中显然忽略了散射。

这两项可以通过分部积分<font color = red>(?)</font>合并为一项（据[4]）

$M_p=-i\int_{-\infty}^\infty dt_1\int_{-\infty}^{t_1}dt_0\langle\psi_p^{Vv}(t_1)|V(r)U(t_1,t_0)V(r)|\psi_0(t_0)\rangle$

引入强场近似：将末态代以Volkov态$\psi_p\rightarrow \psi_p^{Vv}$，并且时间演化算符$U\rightarrow U_f=\int d^3k|\psi_k^{Vv}(t)\rangle\langle\psi_k^{Vv}(t')|$；散射电子的跃迁振幅

$M_p^2=-\int_{t_0}^\infty dt_1\int_{-\infty}^\infty dt_0\int d^3k \langle\psi_p^{Vv}(t_1)|V(r)|\psi_k^{Vv}(t_1)\rangle\langle\psi_k^{Vv}(t_0)|H_I(t_0)|\psi_0(t_0)\rangle$
## 鞍点法(saddle point method)
以长度规范为例，$H_I=-erE(t)$

Volkov波函数形式为$|\psi_p^{Vv}(t)\rangle = |p-eA(t)\rangle e^{-iS_p(t)}$，其中$S_p(t)=\frac{1}{2m}\int^t d\tau[p-eA(\tau)]^2$；

初态$|\psi_0(t)\rangle=e^{iE_{IP}t}|\psi_0\rangle$。

$M_p=-i\int_{-\infty}^\infty d\tau\langle p-eA(\tau)|H_I(\tau)|\psi_0\rangle e^{i[S_p(\tau)+E_{IP}\tau]}$

对于强外场，指数因子$g(t)=i[S_p(\tau)+E_{IP}\tau]$沿积分路径快速变化，因此可以采用鞍点近似

$M_p=-i\sum_{t_s}e^{i\theta}\sqrt{\frac{2\pi}{|g''(t_s)|}}\langle p-eA(t_s)|H_I(t_s)|\psi_0\rangle e^{g(t_s)}$

其中鞍点$t_s$由方程$g'(t_s)=0\Rightarrow E_{IP}+\frac{1}{2}[p-eA(t_s)]^2=0$确定，最速下降方向由$\theta = -\frac{\arg(g''(t_s))}{2}+\frac{\pi}{2}$确定。

--------------------------------
<!--------------表达式未经验证地引用于[4]---------->

对于散射电子，其指数因子$g(t_1,t_0,\mathbf{k})=i\\{-\frac{1}{2m}\int_{t_1}d\tau[\mathbf{p}-e\mathbf{A}(\tau)]^2-\frac{1}{2m}\int_{t_0}^{t_1}d\tau[\mathbf{k}-e\mathbf{A}(\tau)]^2+t_0E_{IP}\\}$（使用加粗的$\bf k$以强调$S_p$具有5个独立变量，因而其鞍点方程有5个，后记为$q_j$）

指数前因子$m_p(\mathbf q)=\langle \mathbf p-e\mathbf A(t_1)|V(r)|\mathbf k-e\mathbf A(t_1)\rangle\langle \mathbf k-e\mathbf A(t_0)|H_I(t_0)| \psi_0(t_0)\rangle$

鞍点方程

$[\mathbf k-e\mathbf A(t_0)]^2=-2mE_{IP}$

$(t_1-t_0)\mathbf k=\int_{t_0}^{t_1}d\tau e\mathbf A(\tau)$

$[\mathbf p-e\mathbf A(t_1)]^2=[\mathbf k-e\mathbf A(t_1)]^2$

鞍点近似表达式

$M_p\sim \sum_s m_p(q)(\frac{(2\pi i)^5}{det(\partial^2 g(q)/\partial q_j\partial q_l)})^{1/2}e^{g(q)}$

----------------------------------------

### 鞍点与复轨道
$M_p=\sum_s \frac{m_p(t_s)}{\sqrt{|S''(t_s)|}}e^{iS(t_s)},\qquad S(t)=S_p(t)+E_{IP}t$

简便起见，令$m=1$，则在形式上有$v(t)=p-eA(t)$。

-------------------------------

经典的牛顿方程：$\ddot r = \dot v = -eE(t)$

给出初始条件和边界条件：$v^2(t_s)=-2E_{IP},v(\infty)=p$

可以解出轨道方程（从$t_s$时刻开始演化，直到末态速度稳定在$p$）

$r_p(t)=C+p(t-t_s)+G(t)-G(t_s)$

其中$G(t)=\int^td\tau eA(\tau)$。

特别要注意的是，这条轨道的**初始时刻和初始速度都是复数**，$v(t_s)$甚至是纯虚的，因此$r_p(t)$也是复数。这样的轨道被称为“复轨道”(complex trajectory)。

但如果$t$为实数，则$v(t)=p-eA(t)$也为实数；进一步地，如果选取合适的积分常数$C=p\textnormal Im(t_s)+\textnormal Im(G(t_s))$，$r_p(t)$也为实数。

考虑这样的演化过程：时间从$t_s$沿虚轴演化，演化至实轴$t_0=\textnormal Re(t_s)$时，沿实轴演化。则电子从$t_0$时刻“出现”在电场中，随后沿着实轨道$r_p(t)$运动。

我们将**每个鞍点与这样的一条复轨道联系起来**，并赋予了它物理意义，尽管在ATI这样的量子体系中提及经典的轨道概念不甚合适。

---------------------------------------------

沿着该复轨道写出作用量

$W_p(t_s)=\int_{t_s}^\infty dt(v^2/2-erE(t)-E_{IP})+v(t_s)r(t_s)-pr(\infty)$

且$S(t_s)=\int^{t_s}dt(v^2/2+E_{IP}+\frac{d}{dt}(vr)-\frac{d}{dt}(vr))=\int^{t_s}(-v^2/2-\dot vr+E_{IP})+v(t_s)r(t_s)-v(\infty)r(\infty)=W_p(t_s)$

也就是说，跃迁振幅可以是所有这些复轨道作用量的“和”：$M_p=\sum_s \frac{m_p(t_s)}{\sqrt{|W_p''(t_s)|}}e^{iW_p(t_s)}$

这一点与费恩曼路径积分的思想一致。

----------------------------------------------

对于散射电子，也可以写出这样的复轨道，但不同的是，对应每个鞍点$q_s$，存在两条复轨道（电离／散射），它们的起始时间分别为$t_{0s}, t_{1s}$，末态动量分别为$k, p$。
## 规范
<font size = 3 color = grey>
虽然强场近似在一定程度上取得了很好的效果，但忽略束缚势仍然导致了一些问题

其一就是在强场近似下，不再具有规范不变性。
</font>

不同规范下，相互作用能和Volkov波函数的表达式有所不同

- 长度规范：$H_I=-erE(t)\qquad\ |\psi_p^{Vv}(t)\rangle=|p-eA(t)\rangle e^{-iS_p(t)}$
- 速度规范：$H_I=-\frac{e}{m}pA(t)\qquad|\psi_p^{Vv}(t)\rangle=|p\rangle e^{-iS_p(t)}$

强场近似的跃迁振幅可以通过分部积分<font color = red>(?)</font>写成以下形式

$M_p=-i\int_{-\infty}^\infty d\tau\langle\psi_p^{Vv}(\tau)|V_B(r)|\psi_0(\tau)\rangle$

这样，不同规范下就只有Volkov波函数的差别；采用鞍点法，这一差别将体现在指数前因子中：

- 长度规范：$\langle p-eA(t_s)|V_B(r)|\psi_0\rangle$

- 速度规范：$\langle p|V_B(r)|\psi_0\rangle$

[5]中给出了一个例子：外场$A$为单色线偏场，在一个周期内，鞍点方程有两个解；对于长度规范，$p-eA(t_s)$在这两个解处沿外场方向的分量符号相反，取决于初态$\psi_0$的宇称，这两个鞍点的指数前因子将同号或反号，它们的和将相长或相消；而对于速度规范，不同鞍点的和总是相长的。

对于有限范围势场，采用长度规范可能会更好。
## 库仑修正
<font color = grey size = 3>

人们对SFA理论的库仑修正做出了许多尝试，以下只是其中一种。

</font>

[S.V. Popruzhenko & D. Bauer. Journal of Modern Optics 55:16(2008); pages: 2573-2589](https://doi.org/10.1080/09500340802161881)[6]

既然跃迁振幅与复轨道存在联系，可以通过修正复轨道的作用量来修正跃迁振幅。

记原复轨道为$r_0(p,t_s)$，引入库仑作用后的复轨道为$r(p,t_s)=r_0(p,t_s)+r_1(p,t_s)$；库仑作用将导致两项修正

$W_C^1(p,t)=-\int_{t_s(p)}^\infty U_c(r)dt =Ze\int_{t_s(p)}^\infty\frac{dt}{|r(p,t)|}$

$W_C^2(p,t)=\int_{t_s(p)}^\infty [v_0\cdot v_1-eE(t)\cdot r_1]dt=\int_{t_s(p)}^\infty(\dot r_0\cdot\dot r_1+\ddot r_0\cdot r_1)dt=v_0(\infty)r_1(\infty)-v_0(t_s)r_1(t_s)$

其中$r_1$满足牛顿方程$\ddot r_1=\frac{-eZ(r_0+r_1)}{|r_0+r_1|^3}$

这里将库仑作用视为微扰。第一项为库仑作用导致的作用量，第二项为$r_0\rightarrow r_0+r_1$代入$W_p$后多出的项。

但一般来说，当$t\rightarrow t_s$时，电子离核很近，库仑作用相当明显，此时的修正不能简单地视为微扰，比如$W_C^1$中库仑势的发散（在这时，电子可能表现出一些奇异的行为，比如再次与核散射）；[6]中认为对于初态为s态的情形，$W_C^1$和$W_C^2$在$t\in [t_s,t_0]$的部分对总作用量没有影响，即我们不需要考虑$t\rightarrow t_s$时的发散问题。

----------------------------------

此外，由于末态$v(\infty)=v_0(\infty)+v_1(\infty)\neq p$，这样的轨道$r_0(p,t)$不是我们要修正的轨道。假设另一条轨道$\widetilde r_0(p,t)=r_0(\widetilde p,t)$，它的初始时刻为$\widetilde t_s(p)=t_s(\widetilde p)$，且满足条件$\widetilde v_0(\widetilde t_s)=-2E_{IP},\widetilde v_0(\infty)=\widetilde p$；经过修正后，$\widetilde v_0(\infty)+v_1(\infty)=p$。

确定轨道$\widetilde r_0$的步骤如下：
1. 按SFA求解原轨道$r_0(p,t)$；
2. 求解$r_1$，得出$v_1(\infty)$，进而得出$\widetilde p=p-v_1(\infty)$，关于$r_1$的初始条件为$r(t_0)=r_0(t_0),v(t_0)=v_0(t_0)$（我们仅需要$r_1$的实数部分，见下）；
3. 按SFA求解关于$\widetilde p$的轨道$r_0(\widetilde p,t)$，鞍点$\widetilde t_s(p)$，修正后的轨道$r=\widetilde r_0 +r_1$；
4. 但此时$r_1$不满足$\ddot r_1=\frac{-eZr}{|r|^3}$，因此需要以$\widetilde r_0$作为新的$r_0$进行迭代，直到$\ddot r_1\approx \frac{-eZr}{|r|^3}$...

解出$\widetilde r_0$和$r$后，就可以对作用量进行修正，但需要注意几点：

1. $t\in [\widetilde t_s,\widetilde t_0]$时，由于$W_C^1$和$W_C^2$没有贡献，这部分的作用量应沿着轨道$\widetilde r_0$，为$W(\widetilde p,\widetilde t_s)$；
2. $t\in[\widetilde t_0,\infty]$时，$W_C^2$已经包含在沿轨道$r$的作用量$W(p,\widetilde t_s)$中，此外还要加上库仑作用项$W_C^1$。

---------------------------------------------------

以上修正仅对于s态的电离过程具有规范不变性；对于其它态，仍然存在规范问题，但比传统SFA“峰－谷”这样极端的差别要小一些。

由于引入了库仑修正，SFA不再局限于短程势，对处理中性原子或分子电离更加有效。
## 附录A Dyson Integral equation

哈密顿量$H=T+V(r)+H_I(t)$，对应的时间演化算符为$U$

哈密顿量$H_a=T+V(r)$，对应的时间演化算符为$U_a$

哈密顿量$H_f=T+H_I(t)$，对应的时间演化算符为$U_f$

薛定谔方程：

1. $(i\partial_t - T - V - H_I)\psi(t) = 0$
2. $(i\partial_t - T - V)\psi_a(t) = 0$
3. $(i\partial _t - T - H_I)\psi_f(t) = 0$

---------------------------------------------

对于方程1，给出一个格林算符表示的解 

$|\psi(t)\rangle = |\psi_f(t)\rangle + \int dt_0G_a(t,t_0)H_I(t_0)|\psi(t_0)\rangle$

其中$G_a(t,t_0)=-i\sum_j|\psi_{aj}(t)\rangle\langle\psi_{aj}(t_0)|$满足与方程2对应的方程：

$(i\partial_t - T - V)G_a(t,t_0) = \delta(t-t_0)$

上面的解可以表示为时间演化算符的形式

$U(t,t')=U_a(t,t')-i\int_{t'}^td\tau U_a(t,\tau)H_I(\tau)U(\tau,t')$

上式中的$U_a$和$U$可以交换，通过取复共轭并交换$t$和$t'$；可以理解为反向演化。

---------------------------------------------

关于时间演化算符的Dyson方程 ：

$U(t,t')=U_a(t,t')-i\int_{t'}^td\tau U(t,\tau)H_I(\tau)U_a(\tau,t')$

$U(t,t')=U_f(t,t')-i\int_{t'}^td\tau U_f(t,\tau)V(r)U(\tau,t')$

## 附录B L-guage & V-guage
电场作用下的哈密顿量$H=\frac{[p-eA(t)]^2}{2m}+V(r)=\frac{p^2}{2m}+V(r)-\frac{e}{m}pA(t)+\frac{e^2}{2m}A^2(t)$

薛定谔方程$i\partial_t\psi = H\psi$

插入相位因子$U=e^{i(\frac{e^2}{2m}\int^t d\tau A^2(\tau))}$，$i\partial_t U^*U\psi = H U^*U\psi$，并视$U\psi$为新的波函数$\psi'$。整理得到$\psi'$满足的薛定谔方程

$\frac{e^2}{2m}A^2(t)U^* \psi' + U^* i\partial_t\psi' = HU^* \psi'\Rightarrow i\partial_t\psi' = U(H-\frac{e^2}{2m}A^2(t))U^*\psi'$

变换后，相互作用项$H_I=-\frac{e}{m}pA(t)$。

而若是插入相位因子$U=e^{i(-erA(t)+\frac{e^2}{2m}\int^t d\tau A^2(\tau))}$，同样地，可以得到变换后的相互作用项

$H_I=-erE(t)$

注意两种变换的波函数仅相差一个相位因子$e^{ierA(t)}$。
