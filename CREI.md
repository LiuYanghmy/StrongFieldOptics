Charge-Resonance-Enhanced Ionization(CREI)

-----------------------------

TODO LIST

- [x] 概述：CREI
- [ ] Charge-Resonant(CR) States
- [ ] Multiple Ionization Bursts
- [ ] momentum gates

-----------------------------

[1] T. Zuo and A. D. Bandrauk, PHYS. REV. A, 52(4), 1995
[2] Tamar Seideman, * M. Yu. Ivanov, and P. B. Corkum, PHYS. REV. LETT. 75, 2819, 1995 
[3] Norio Takemoto* and Andreas Becker, PHYS. REV. A, 84, 023401, 2011

## 概述

双原子分子离子（例：$\mathtt H_2^+$）电子在线偏场下电离，电离率随核间距$R$变化（见[1] fig.1）；随核间距增大，电离率增大，随后出现峰值，进而减小，并趋于单原子电离率。

-----------------------------------------

在电场＋库仑势的作用下，电子的概率密度分布并非均分到两个核，而是随电场在两个核间振荡，这一性质与$\mathtt H_2^+$的基态$|g\rangle$和第一激发态$|u\rangle$（被称为"charge-resonant states (pair)"）有关；由此产生的偶极矩$d\sim R/2$。

记两个核的连线方向为z轴，且原点设置在中点。由于外场作用，势垒在z轴方向呈现“左低右高”的结构，并被中间势垒隔开，在核附近形成一高一低两个势阱（见[2] fig. 1）；两势阱的能量差约为$ER$。

处于高势阱中的电子具有更高的能量，容易越过中间的势垒发生电离；而若是此时电子主要集中在高势阱中，就可以获得较大的电离率。高低势阱的位置随电场周期性交替，而电荷分布也随电场同步振荡。

- $R$较小时，电子的不对称分布并不明显，处于高势阱中的电子不算太多，因此电离率不大；

- $R$很大时，两个核几乎可以看作是独立的原子，中间势垒很高，即使电子主要分布在高势阱中，仍然难以电离，极端情况下趋于原子电离率；

- 只有在不太小和不太大的核间距$R$下，电离率才会得到显著增强。

## 电子分布

使用$|R\rangle=\frac{1}{\sqrt{2}}(|g\rangle+|u\rangle)$表示电子处于右侧（$z>0$）的状态；使用$|L\rangle=\frac{1}{\sqrt{2}}(|g\rangle-|u\rangle)$表示电子处于左侧（$z<0$）的状态(?)。

$P_L=|\langle L|\varPsi\rangle|^2$，$P_R=|\langle R|\varPsi\rangle|^2$分别为电子处于左侧和右侧的概率密度，态$\varPsi$可以通过求解含时薛定谔方程得到，或表示为Floquet态的叠加：$|\varPsi\rangle=c_1e^{-i\epsilon_1\tau/\hbar}|\psi_1^F\rangle+c_2e^{-i\epsilon_2\tau/\hbar}|\psi_2^F\rangle$（[3] sec. 3）

$P_L\approx|c_1+c_2|^2/2+\Delta P$

$P_R\approx|c_1-c_2|^2/2-\Delta P$

$P_L,P_R$在一个脉冲周期内的变化见[3] fig. 1，$\Delta P$随时间和外场强度的变化见[3] fig. 4。

（[3] fig. 1）当电场指向$-z$方向时，右侧的势垒高，此时恰好$P_R>0.5$，电子主要占据右侧的核；并且在半个周期内，占据数出现了多个峰值，这意味着在半个周期内可以出现多次增强电离（Multiple Ionization Bursts,MIBs）。这一结论在[3] fig. 4中也可以看到。

-------------------------------------

$P_L,P_R$的第一项$|c_1\pm c_2|^2/2$不随时间改变，这意味着除非电子呆在其中一个Floquet态，否则电子密度总是不对称分布的；在$R$比较小时，振荡的密度分布（对应一个偶极矩$d\sim R/2$）比较小，此时，这种静态的不对称分布可能也会导致增强电离（[1] fig. 1中靠左的峰）。