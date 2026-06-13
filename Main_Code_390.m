%---------------------------------------
%    DNS data for Re=(u_c\delta)/nu=4400
%----------------------------------------
%
clear
close all
clc

ny=257; %-- Number of grid points in y-wall normal direction (Chebychev modes)
nx=256; %-- Number of grid points in x-wall parallel direction (Fourier modes)
nz=256; %-- Number of grid points in z-wall parallel direction (Fourier modes)

load('Dataset_390.mat')

%==========================================================================
Re=4400;
Ret=390;
uc=mean(mean(utest(:,:,round((ny/2))+1),1),2);
nu=uc/Re;
uta=Ret*nu;
yplus=flip(y*(uta/nu));
%==========================================================================
%                         1.wall distance
%==========================================================================
Y=flip(y);
for j=1:256
    for i=1:129
        d(i,j)=Y(129)-Y(i);
    end
end
for j=1:256
    for i=129:257
        d(i,j)=Y(i)-Y(129);
    end
end
f1=d;
%==========================================================================
%                             2.Uave
%==========================================================================
Uave=reshape(mean(utest,1),[256,257]);
f2=Uave';
%==========================================================================
%                             3.Vave
%==========================================================================
Vave=reshape(mean(vtest,1),[256,257]);
f3=Vave';
%==========================================================================
%                             4.Wave
%==========================================================================
Wave=reshape(mean(wtest,1),[256,257]);
f4=Wave';
%==========================================================================
%                       --fluctuation velocities--
%==========================================================================
U_ave=reshape(mean(mean(utest,1),2),[257,1]);
V_ave=reshape(mean(mean(vtest,1),2),[257,1]);
W_ave=reshape(mean(mean(wtest,1),2),[257,1]);
for i=1:256
    for k=1:256
        for j=1:257
            uprime(i,k,j)=utest(i,k,j)-U_ave(j);
            vprime(i,k,j)=vtest(i,k,j)-V_ave(j);
            wprime(i,k,j)=wtest(i,k,j)-W_ave(j);
        end
    end
end
%=============================F5:TKE=======================================
u1=reshape(mean(power(uprime,2),1),[256,257]);
u2=reshape(mean(power(vprime,2),1),[256,257]);
u3=reshape(mean(power(wprime,2),1),[256,257]);
K=0.5.*(u1+u2+u3);
f5=K';
%=============================F6:du_dy=====================================
clear('Uave')
for i=1:256
    Uave=reshape(mean(utest(i,:,:),2),[257,1]);
    for j=1:255
        alpha(j)=(Y(j+2)-Y(j+1))/(Y(j+1)-Y(j));
        A1(j)=-Uave(j+2);
        A2(j)=((1+alpha(j)).^2).*(Uave(j+1));
        A3(j)=-alpha(j).*(alpha(j)+2).*Uave(j);
        A4(j)=alpha(j).*(alpha(j)+1).*(Y(j+1)-Y(j));
        du_dy(i,j)=(A1(j)+A2(j)+A3(j))./A4(j);
    end
    clear('Uave', 'A1' , 'A2' , 'A3' , 'A4')
end
Uave=reshape(mean(utest(:,end,:),1),[257,1]);
du_dy(:,256)=(Uave(257)-Uave(255))./((Y(257)-Y(255)));
du_dy(:,257)=(Uave(257)-Uave(256))./(Y(257)-Y(256));
f6=du_dy';
%===========================F7:dissipation Rate============================
for i=1:256
    for k=1:256
        for j=1:255
            alpha(j)=(y(j+2)-y(j+1))/(y(j+1)-y(j));
            Fr(j)=alpha(j)*((1+alpha(j))^2)*(y(j+1)-y(j));
            ey1(i,k,j)=(-alpha(j).*(alpha(j)+2)*uprime(i,k,j)+(1+alpha(j))^2*uprime(i,k,j+1)-uprime(i,k,j+2))./Fr(j);
            ey2(i,k,j)=(-alpha(j).*(alpha(j)+2)*vprime(i,k,j)+(1+alpha(j))^2*vprime(i,k,j+1)-vprime(i,k,j+2))./Fr(j);
            ey3(i,k,j)=(-alpha(j).*(alpha(j)+2)*wprime(i,k,j)+(1+alpha(j))^2*wprime(i,k,j+1)-wprime(i,k,j+2))./Fr(j);
        end
    end
end
ey1(:,:,256)=(uprime(:,:,256)-uprime(:,:,255))./(2*(y(256)-y(255)));
ey2(:,:,256)=(vprime(:,:,256)-vprime(:,:,255))./(2*(y(256)-y(255)));
ey3(:,:,256)=(wprime(:,:,256)-wprime(:,:,255))./(2*(y(256)-y(255)));
ey1(:,:,257)=(uprime(:,:,257)-uprime(:,:,256))./(2*(y(257)-y(256)));
ey2(:,:,257)=(vprime(:,:,257)-vprime(:,:,256))./(2*(y(257)-y(256)));
ey3(:,:,257)=(wprime(:,:,257)-wprime(:,:,256))./(2*(y(257)-y(256)));
ey1=reshape(mean(ey1.^2,1),[256,257]);
ey2=reshape(mean(ey2.^2,1),[256,257]);
ey3=reshape(mean(ey3.^2,1),[256,257]);

for i=1:256
    for j=1:257
        for k=2:255
            ez1(i,k,j)=((uprime(i,k+1,j)-uprime(i,k-1,j)))./2;
            ez2(i,k,j)=((vprime(i,k+1,j)-vprime(i,k-1,j)))./2;
            ez3(i,k,j)=((wprime(i,k+1,j)-wprime(i,k-1,j)))./2;
            
        end
    end
end
ez1(:,256,:)=(uprime(:,256,:)-uprime(:,254,:))./2;
ez1(:,1,:)=ez1(:,2,:);
ez2(:,256,:)=(vprime(:,256,:)-vprime(:,254,:))./2;
ez2(:,1,:)=ez2(:,2,:);
ez3(:,256,:)=(wprime(:,256,:)-wprime(:,254,:))./2;
ez3(:,1,:)=ez3(:,2,:);
ez1=reshape(mean(ez1.^2,1),[256,257]);
ez2=reshape(mean(ez2.^2,1),[256,257]);
ez3=reshape(mean(ez3.^2,1),[256,257]);

for j=1:257
    for k=1:256
        for i=2:255   
            ex1(i,k,j)=((uprime(i+1,k,j)-uprime(i-1,k,j)))./2;
            ex2(i,k,j)=((vprime(i+1,k,j)-vprime(i-1,k,j)))./2;
            ex3(i,k,j)=((wprime(i+1,k,j)-wprime(i-1,k,j)))./2;
        end
    end
end
ex1(256,:,:)=((uprime(256,:,:)-uprime(255,:,:)));
ex2(256,:,:)=((vprime(256,:,:)-vprime(255,:,:)));
ex3(256,:,:)=(wprime(256,:,:)-wprime(255,:,:));
ex1(1,:,:)=ex1(2,:,:);
ex2(1,:,:)=ex2(2,:,:);
ex3(1,:,:)=ex3(2,:,:);
ex1=reshape(mean(ex1.^2,1),[256,257]);
ex2=reshape(mean(ex2.^2,1),[256,257]);
ex3=reshape(mean(ex3.^2,1),[256,257]);
epsilon=nu.*(ex1+ex2+ex3+ey1+ey2+ey3+ez1+ez2+ez3);
f7=epsilon';
%============================F8:Cross Correlation==========================
for j=1:257
        for i=2:255   
            ccx1(i,j)=((K(i+1,j)-K(i-1,j)))./2;
            ccx2(i,j)=((epsilon(i+1,j)-epsilon(i-1,j)))./2;
        end
end
ccx1(256,:)=(K(256,:)-K(255,:));
ccx2(256,:)=(epsilon(256,:)-epsilon(256,:));
cc1=ccx1.*ccx2;

for i=1:256
        for j=1:255
            alpha0(j)=(y(j+2)-y(j+1))/(y(j+1)-y(j));
            Fr0(j)=alpha(j)*((1+alpha(j))^2)*(y(j+1)-y(j));
            ccy1(i,j)=(-alpha0(j).*(alpha0(j)+2)*K(i,j)+(1+alpha0(j))^2*K(i,j+1)-K(i,j+2))./Fr0(j);
            ccy2(i,j)=(-alpha0(j).*(alpha0(j)+2)*epsilon(i,j)+(1+alpha0(j))^2*epsilon(i,j+1)-epsilon(i,j+2))./Fr0(j);
        end
end
ccy1(:,256)=(K(:,257)-K(:,255))./(y(257)-y(255));
ccy2(:,256)=(epsilon(:,257)-epsilon(:,255))./(y(257)-y(255));

ccy1(:,257)=(K(:,257)-K(:,256))./(2*(y(257)-y(256)));
ccy2(:,257)=(epsilon(:,257)-epsilon(:,256))./(y(257)-y(256));

cc2=ccy1.*ccy2;
f8=(cc1+cc2)';
%===================F9_output1 : Eddy Viscosity============================
%Global Eddy Viscosity
R1=reshape(mean(mean(uprime.*vprime,1),2),[257,1]);
R2=reshape(mean(du_dy,1),[257,1]);
E1=R1./R2;
for i=1:256
    for j=1:257
        GEV(i,j)=E1(j);
    end
end
f9=GEV';
%===================F9_output1 : DNS Stress================================
uprime_vprime=reshape(mean(uprime.*vprime,1),[256,257]);
f10=uprime_vprime';
%========================[Import to Python]================================
save("data_properties_390.mat" , 'f1' ,'f2' ,'f3' , 'f4' , 'f5', 'f6' , 'f7', 'f8' , 'f9' , 'f10')
%==========================================================================




