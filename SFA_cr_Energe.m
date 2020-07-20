clear
c=2.99792e8;
h=6.62607e-34;
hbar=1.05457e-34;
e=1.60217e-19;
me=9.10938e-31;
epsilon=8.85419e-12;
%basic translation of atomic unit to SI
L0= 5.29177e-11;
T0= 2.41888e-17; 
En0= 4.35974e-18;
Mom0=1.99285e-24;  
F0= 8.23872e-8;  
E0=5.1422e11; 
V0= 2.18769e6; 
%%
lambda=800e-9;
omega=2*pi*c/lambda;    %SI
period=lambda/c;    %SI
omega0=2*pi*c/V0/(lambda/L0);   %atomic unit
omega02=2*omega0;
period0=period/T0;  %atomic unit
lambda0=lambda/L0;  %atomic unit

global phi;
global rotation;
global ratio;
global eta;
global I;
phi = 0;
rotation = -1;
ratio = 4;
eta = 1;
nI = 100;

scale = 2;
tprecision = 2;
Ip=0.57;
l=1;
R=2*scale;
nrho=100*scale;
ntheta=90*tprecision;
krho=linspace(0,R,nrho);
ktheta=linspace(0,2*pi,ntheta);
drho=krho(2)-krho(1);
dtheta=ktheta(2)-ktheta(1);
dR=((krho+drho).^2-krho.^2).*dtheta./2;
M0=zeros(ntheta,nrho);
M1=zeros(ntheta,nrho);
SS=zeros(ntheta,nrho);
[Krho,Ktheta]=meshgrid(krho,ktheta);
Prho=linspace(0,R,nrho);
Irate=linspace(0,1,nI);
[PP,II]=meshgrid(Prho,Irate);
P=zeros(nrho,nI);
options = optimset('TolX',1e-15,'display','off'); 

%--- I0=1e14
% t1 = [132.181871631921 + 6.96141169070234i,186.411467849943 - 0.0895005032587793i];
% t2 = [132.199257498573 + 6.53324238638841i,186.327367475254 - 0.107186789909792i]; 
% td = 128.706819580963 + 6.77890181189868i; 
%------------
%--- I0=0.6e14
t1 = [132.13066220896783 + 8.166373781758889i,186.68205624791742 - 0.01268137029923173i];
t2 = [132.13066220896783 + 8.16637378175889i,186.68205624791744 - 0.01268137029923066i]; 
td = 128.70681962907867 + 7.97812527129341i; 
tsd = zeros(nI,1);
tss1 = zeros(nI,2);
tss2 = zeros(nI,2);

nImin = nI;
 nT = 30;
tis=zeros(ntheta,nrho,nT*3);
trs=zeros(ntheta,nrho,nT*3);
tis1=zeros(ntheta,nrho,nT*3);
SS0=zeros(ntheta,nrho,nT*3);
SS1=zeros(ntheta,nrho,nT*3);
for rate = nI:-1:nImin
I=0.6e14*rate/nI;
I2=ratio*I;

AmpE=2744*sqrt(I)/E0;
AmpE2=2744*sqrt(I2)/E0;
T=3*period0;

Ex=@(t) AmpE*sin(omega0*t)+AmpE2*sin(2*omega0*t+phi);
Ey=@(t) eta*(AmpE*cos(omega0*t)+rotation*AmpE2*cos(2*omega0*t+phi));
Ax=@(t) AmpE/omega0.*cos(omega0*t)+AmpE2/(2*omega0).*cos(2*omega0*t+phi);
Ay=@(t) -eta*(AmpE/omega0.*sin(omega0*t)+rotation*AmpE2/(2*omega0).*sin(2*omega0*t+phi));
Xx=@(t) AmpE/omega0/omega0.*sin(omega0*t)+AmpE2/(2*omega0)/(2*omega0).*sin(2*omega0*t+phi);
Xy=@(t) eta*(AmpE/omega0/omega0.*cos(omega0*t)+rotation*AmpE2/(2*omega0)/(2*omega0).*cos(2*omega0*t+phi));

Up=((AmpE/omega0)^2+(AmpE2/(2*omega0))^2);

r2x = @(t) (2*Up)^0.5.*cos(omega0*t);
r2y = @(t) (2*Up)^0.5.*sin(omega0*t);
r5x = @(t) (5*Up)^0.5.*cos(omega0*t);
r5y = @(t) (5*Up)^0.5.*sin(omega0*t);
r10x = @(t) (10*Up)^0.5.*cos(omega0*t);
r10y = @(t) (10*Up)^0.5.*sin(omega0*t);

if 0
    t = linspace(0,T,1000);
    figure
    plot(Ax(t),Ay(t),'green')
    figure
    plot(Xx(t),Xy(t),'red')
    figure
    plot(Ex(t),Ey(t),'blue')
end

 global px py;

 tmp = [0,0];
 tmpd = 0;
 DT = [period0/3,period0/3];

 %------singlecycle
for n = 0:2
 TT=2*period0;
 jd = mod(92*tprecision+n*ntheta/3,ntheta);
 j00 = mod(32*tprecision+n*ntheta/3,ntheta);
 jnn = mod(32*tprecision-1+n*ntheta/3,ntheta);
%% direct
 i0 = 1;
 j0 = mod(61*tprecision+n*ntheta/3, ntheta);
  if n == 0 %更新每次迭代的t0
       px=krho(i0)*cos(ktheta(j0));
       py=krho(i0)*sin(ktheta(j0));
       f=@(t) (px+Ax(t)).^2+(py+Ay(t)).^2+2*Ip;
       [tmpd,~]=fsolve(f,td,options);
  td = tmpd;
  tsd(rate,1) = td;
 end
  t0 = td - n*period0/3;
 for ii=i0
    ii;
    for jj=j0
       px=krho(ii)*cos(ktheta(jj));
       py=krho(ii)*sin(ktheta(jj));
      f=@(t) (px+Ax(t)).^2+(py+Ay(t)).^2+2*Ip;
      [ts,~]=fsolve(f,t0,options);
      tis1(jj,ii,n+1)=ts;
    end
 end

  for ii=i0+1:nrho
    ii;
    for jj=j0
        px=krho(ii)*cos(ktheta(jj));
        py=krho(ii)*sin(ktheta(jj));
      f=@(t) (px+Ax(t)).^2+(py+Ay(t)).^2+2*Ip;
      [ts,~]=fsolve(f,tis1(jj,ii-1,n+1),options);
      tis1(jj,ii,n+1)=ts;
    end
  end
  for ii=i0-1:-1:1
    ii;
    for jj=j0
        px=krho(ii)*cos(ktheta(jj));
        py=krho(ii)*sin(ktheta(jj));
      f=@(t) (px+Ax(t)).^2+(py+Ay(t)).^2+2*Ip;
      [ts,~]=fsolve(f,tis1(jj,ii+1,n+1),options); 
      tis1(jj,ii,n+1)=ts;
    end
  end 
 for ii=1:nrho
    ii;
    for jj=j0+1:ntheta
        px=krho(ii)*cos(ktheta(mod(jj-1,ntheta)+1));
        py=krho(ii)*sin(ktheta(mod(jj-1,ntheta)+1));
      f=@(t) (px+Ax(t)).^2+(py+Ay(t)).^2+2*Ip;
      [ts,~]=fsolve(f,tis1(mod(jj-1-1,ntheta)+1,ii,n+1),options); 
      tis1(mod(jj-1+ntheta,ntheta)+1,ii,n+1)=ts;     
    end
end
for ii=1:nrho
    ii;
    for jj=j0-1:-1:1
        px=krho(ii)*cos(ktheta(mod(jj+ntheta-1,ntheta)+1));
        py=krho(ii)*sin(ktheta(mod(jj+ntheta-1,ntheta)+1));
      f=@(t) (px+Ax(t)).^2+(py+Ay(t)).^2+2*Ip;
      [ts,~]=fsolve(f,tis1(mod(jj+ntheta,ntheta)+1,ii,n+1),options); 
      tis1(mod(jj-1+ntheta,ntheta)+1,ii,n+1)=ts;        
    end
end


%% scaterring
%% theta < jd
 i0 = 1;
 j0 = mod(61*tprecision+n*ntheta/3, ntheta);
%       t0 = [131.995701216215 - n*period0/3 + 6.49625241138080i, ...
%     189.020656896083 - n*period0/3 - 0.164141711565128i];  
 if n == 0 %更新每次迭代的t0
       px=krho(i0)*cos(ktheta(j0));
       py=krho(i0)*sin(ktheta(j0));
  [tmp,~]=fsolve(@myfun1D,t1,options);
  t1 = tmp;
  tss1(rate,1)=t1(1);
  tss1(rate,2)=t1(2);
 end
  t0 = t1 - n.*DT;
 for ii=i0
    ii;
    for jj=j0
       px=krho(ii)*cos(ktheta(jj));
       py=krho(ii)*sin(ktheta(jj));
      [ts,~]=fsolve(@myfun1D,t0,options);
      tis(jj,ii,n+1)=ts(1);
      trs(jj,ii,n+1)=ts(2);
    end
 end

  for ii=i0+1:nrho
    ii;
    for jj=j0
        px=krho(ii)*cos(ktheta(jj));
        py=krho(ii)*sin(ktheta(jj));
      [ts,~]=fsolve(@myfun1D,[tis(jj,ii-1,n+1);trs(jj,ii-1,n+1)],options);
      tis(jj,ii,n+1)=ts(1);
      trs(jj,ii,n+1)=ts(2);
    end
  end
  for ii=i0-1:-1:1
    ii;
    for jj=j0
        px=krho(ii)*cos(ktheta(jj));
        py=krho(ii)*sin(ktheta(jj));
      [ts,~]=fsolve(@myfun1D,[tis(jj,ii+1,n+1);trs(jj,ii+1,n+1)],options);
      tis(jj,ii,n+1)=ts(1);
      trs(jj,ii,n+1)=ts(2);
    end
  end 
 for ii=1:nrho
    ii;
    for jj=j0+1:j0+mod(jd-j0+ntheta,ntheta)-1
        px=krho(ii)*cos(ktheta(mod(jj-1,ntheta)+1));
        py=krho(ii)*sin(ktheta(mod(jj-1,ntheta)+1));
      [ts,~]=fsolve(@myfun1D,[tis(mod(jj-1-1,ntheta)+1,ii,n+1);trs(mod(jj-1-1,ntheta)+1,ii,n+1)],options);
      tis(mod(jj-1+ntheta,ntheta)+1,ii,n+1)=ts(1);
      trs(mod(jj-1+ntheta,ntheta)+1,ii,n+1)=ts(2); 
    end
end
for ii=1:nrho
    ii;
    for jj=j0-1:-1:j0-mod(j0-j00+ntheta,ntheta)
        px=krho(ii)*cos(ktheta(mod(jj+ntheta-1,ntheta)+1));
        py=krho(ii)*sin(ktheta(mod(jj+ntheta-1,ntheta)+1));
      [ts,~]=fsolve(@myfun1D,[ tis(mod(jj+ntheta,ntheta)+1,ii,n+1); trs(mod(jj+ntheta,ntheta)+1,ii,n+1)],options);
      tis(mod(jj-1+ntheta,ntheta)+1,ii,n+1)=ts(1);
      trs(mod(jj-1+ntheta,ntheta)+1,ii,n+1)=ts(2); 
    end
end
%% theta > jd
 i0 = 1;
 j0 = mod(10*tprecision+n*ntheta/3, ntheta);
%       t0 = [133.305261871640 - n*period0/3 + 6.87432663079718i, ...
%     178.499402603998 - n*period0/3 + 0.706015929240085i];
 if n == 0
       px=krho(i0)*cos(ktheta(j0));
       py=krho(i0)*sin(ktheta(j0)); 
  [tmp,~]=fsolve(@myfun1D,t2,options);
  t2 = tmp;   
  tss2(rate,1)=t2(1);
  tss2(rate,2)=t2(2);
 end
  t0 = t2 - n.*DT;
 for ii=i0
    ii;
    for jj=j0
       px=krho(ii)*cos(ktheta(jj));
       py=krho(ii)*sin(ktheta(jj));
      [ts,~]=fsolve(@myfun1D,t0,options);
      tis(jj,ii,n+1)=ts(1);
      trs(jj,ii,n+1)=ts(2);
    end
 end

  for ii=i0+1:nrho
    ii;
    for jj=j0
        px=krho(ii)*cos(ktheta(jj));
        py=krho(ii)*sin(ktheta(jj));
      [ts,~]=fsolve(@myfun1D,[tis(jj,ii-1,n+1);trs(jj,ii-1,n+1)],options);
      tis(jj,ii,n+1)=ts(1);
      trs(jj,ii,n+1)=ts(2);
    end
  end
  for ii=i0-1:-1:1
    ii;
    for jj=j0
        px=krho(ii)*cos(ktheta(jj));
        py=krho(ii)*sin(ktheta(jj));
      [ts,~]=fsolve(@myfun1D,[tis(jj,ii+1,n+1);trs(jj,ii+1,n+1)],options);
      tis(jj,ii,n+1)=ts(1);
      trs(jj,ii,n+1)=ts(2);
    end
  end 
 for ii=1:nrho
    ii;
    for jj=j0+1:j0+mod(jnn-j0+ntheta,ntheta)
        px=krho(ii)*cos(ktheta(mod(jj-1,ntheta)+1));
        py=krho(ii)*sin(ktheta(mod(jj-1,ntheta)+1));
      [ts,~]=fsolve(@myfun1D,[tis(mod(jj-1-1,ntheta)+1,ii,n+1);trs(mod(jj-1-1,ntheta)+1,ii,n+1)],options);
      tis(mod(jj-1+ntheta,ntheta)+1,ii,n+1)=ts(1);
      trs(mod(jj-1+ntheta,ntheta)+1,ii,n+1)=ts(2);      
    end
end
for ii=1:nrho
    ii;
    for jj=j0-1:-1:j0-mod(j0-jd+ntheta,ntheta)
        px=krho(ii)*cos(ktheta(mod(jj+ntheta-1,ntheta)+1));
        py=krho(ii)*sin(ktheta(mod(jj+ntheta-1,ntheta)+1));
      [ts,~]=fsolve(@myfun1D,[ tis(mod(jj+ntheta,ntheta)+1,ii,n+1); trs(mod(jj+ntheta,ntheta)+1,ii,n+1)],options);
      tis(mod(jj-1+ntheta,ntheta)+1,ii,n+1)=ts(1);
      trs(mod(jj-1+ntheta,ntheta)+1,ii,n+1)=ts(2); 
    end
end
end
%-----------multicycle
for kk = 4:nT*3
    tis1(:,:,kk)=tis1(:,:,kk-3)+period0;
    tis(:,:,kk)=tis(:,:,kk-3)+period0;
    trs(:,:,kk)=trs(:,:,kk-3)+period0;
end
%% calcu.
 for ii=1:nrho
    ii;
    for jj=1:ntheta
            px=krho(ii)*cos(ktheta(jj));
            py=krho(ii)*sin(ktheta(jj));
        for kk = 1:nT*3
            A1=AmpE/omega0;
            A2=AmpE2/(2*omega0);
         %direct
            ti=tis1(jj,ii,kk);
%             ss=quad(@(t) (py+Ay(t)).^2./2+(px+Ax(t)).^2./2+Ip,TT,tis1(jj,ii,kk));
            ss=((px^2+py^2)/2+Ip)*(ti-TT) + px*Xx(ti)-px*Xx(TT)+py*Xy(ti)-py*Xy(TT) + ...
              ((A1^2+A2^2)*(ti-TT))/2 + A1*A2*(sin(3*omega0*ti)-sin(3*omega0*TT))/(3*omega0);% eta=1,rotation=-1,phi=0
            SS0(jj,ii,kk)= exp(1i*ss);
         %scattering
            t_m1=tis(jj,ii,kk);
            t_m2=trs(jj,ii,kk);
            kxx=1/(t_m1-t_m2)*(Xx(t_m2)-Xx(t_m1));
            kyy=1/(t_m1-t_m2)*(Xy(t_m2)-Xy(t_m1));
%             Sr_t=-(integral(@(y) ((kxx+Ax(y)).^2+(kyy+Ay(y)).^2)./2,t_m1,real(t_m1)))+Ip*(t_m1-TT);
%             Sr_c=-(integral(@(y) ((kxx+Ax(y)).^2+(kyy+Ay(y)).^2)./2,real(t_m1),real(t_m2)));
%             Sr_r=-(integral(@(y) ((kxx+Ax(y)).^2+(kyy+Ay(y)).^2)./2,real(t_m2),t_m2))+(integral(@(y) ((px+Ax(y)).^2./2+(py+Ay(y)).^2./2),real(t_m2),t_m2));
%             Sr_l=-(integral(@(y) (px+Ax(y)).^2./2+(py+Ay(y)).^2./2,real(t_m2),TT));

            Sr_t=((kxx^2+kyy^2)/2)*(t_m1-real(t_m1))+Ip*(t_m1-TT) + kxx*Xx(t_m1)-kxx*Xx(real(t_m1))+kyy*Xy(t_m1)-kyy*Xy(real(t_m1)) + ...
                 ((A1^2+A2^2)*(t_m1-real(t_m1)))/2 + A1*A2*(sin(3*omega0*t_m1)-sin(3*omega0*real(t_m1)))/(3*omega0);
            Sr_c=((kxx^2+kyy^2)/2)*(real(t_m1)-real(t_m2)) + kxx*Xx(real(t_m1))-kxx*Xx(real(t_m2))+kyy*Xy(real(t_m1))-kyy*Xy(real(t_m2)) + ...
                 ((A1^2+A2^2)*(real(t_m1)-real(t_m2)))/2 + A1*A2*(sin(3*omega0*real(t_m1))-sin(3*omega0*real(t_m2)))/(3*omega0);
            Sr_r=((kxx^2+kyy^2)/2)*(real(t_m2)-t_m2) + kxx*Xx(real(t_m2))-kxx*Xx(t_m2)+kyy*Xy(real(t_m2))-kyy*Xy(t_m2) + ...
                 ((A1^2+A2^2)*(real(t_m2)-t_m2))/2 + A1*A2*(sin(3*omega0*real(t_m2))-sin(3*omega0*t_m2))/(3*omega0) + ...
                 ((px^2+py^2)/2)*(t_m2-real(t_m2)) + px*Xx(t_m2)-px*Xx(real(t_m2))+py*Xy(t_m2)-py*Xy(real(t_m2)) + ...
                 ((A1^2+A2^2)*(t_m2-real(t_m2)))/2 + A1*A2*(sin(3*omega0*t_m2)-sin(3*omega0*real(t_m2)))/(3*omega0);
            Sr_l=((px^2+py^2)/2)*(real(t_m2)-TT) + px*Xx(real(t_m2))-px*Xx(TT)+py*Xy(real(t_m2))-py*Xy(TT) + ...
                 ((A1^2+A2^2)*(real(t_m2)-TT))/2 + A1*A2*(sin(3*omega0*real(t_m2))-sin(3*omega0*TT))/(3*omega0);
            Sr=Sr_t+Sr_c+Sr_r+Sr_l;
            SS1(jj,ii,kk)= exp(1i*Sr);
        end
    end
 end 
%  SS = SS1(:,:,1)+SS1(:,:,2)+SS1(:,:,3);
%  M1=abs(SS).^2; 
%  SS = SS0(:,:,1)+SS0(:,:,2)+SS0(:,:,3);
%  M0=abs(SS).^2; 
%  M=M1+M0;
 SS = sum(SS1,3)+sum(SS0,3);
 M=abs(SS).^2;
 P(:,rate)=sum(M,1);
end
save('data','P','tsd','tss1','tss2');
%  figure('Color',[1 1 1]);
%  X=krho'*cos(ktheta);
%  Y=krho'*sin(ktheta);
%  h=pcolor(X,Y,M');
% h=pcolor(Ktheta,Krho,M1);
% h=pcolor(II(nImin:nI,:)',PP(nImin:nI,:)',P(nImin:nI,:));
%  set(h, 'LineStyle','none');


  %%
 function f=myfun1D(t)
c=2.99792e8;
h=6.62607e-34;
hbar=1.05457e-34;
e=1.60217e-19;
me=9.10938e-31;
epsilon=8.85419e-12;
%basic translation of atomic unit to SI
L0= 5.29177e-11;
T0= 2.41888e-17; 
En0= 4.35974e-18;
Mom0=1.99285e-24;  
F0= 8.23872e-8;  
E0=5.1422e11; 
V0= 2.18769e6; 
%%
lambda=800e-9;
omega=2*pi*c/lambda;    %SI
period=lambda/c;    %SI
omega0=2*pi*c/V0/(lambda/L0);   %atomic unit
omega02=2*omega0;
period0=period/T0;  %atomic unit
lambda0=lambda/L0;  %atomic unit
global phi;
global rotation;
global ratio;
global eta;
global I;
I2=ratio*I;
AmpE=2744*sqrt(I)/E0;
AmpE2=2744*sqrt(I2)/E0;
Ip=0.57;
global px py;
Ax=@(t) AmpE/omega0.*cos(omega0*t)+AmpE2/(2*omega0).*cos(2*omega0*t+phi);
Ay=@(t) -eta*(AmpE/omega0.*sin(omega0*t)+rotation*AmpE2/(2*omega0).*sin(2*omega0*t+phi));
Xx=@(t) AmpE/omega0/omega0.*sin(omega0*t)+AmpE2/(2*omega0)/(2*omega0).*sin(2*omega0*t+phi);
Xy=@(t) eta*(AmpE/omega0/omega0.*cos(omega0*t)+rotation*AmpE2/(2*omega0)/(2*omega0).*cos(2*omega0*t+phi));
 f(1)=((Xx(t(2))-Xx(t(1)))+((t(1)-t(2))).*Ax(t(1))).^2+((Xy(t(2))-Xy(t(1)))+((t(1)-t(2))).*Ay(t(1))).^2+2*Ip*((t(1)-t(2))).^2;
 f(2)=(px+Ax(t(2))).^2*((t(1)-t(2))).^2+(py+Ay(t(2))).^2*((t(1)-t(2))).^2-(((Xx(t(2))-Xx(t(1)))+((t(1)-t(2)))*Ax(t(2))).^2+((Xy(t(2))-Xy(t(1)))+((t(1)-t(2)))*Ay(t(2))).^2);
 end
%%
%由前向散射确定初值后在动量谱上扫一遍之后在动量谱另一半得到的是背向散射,不需要在另算背向散射