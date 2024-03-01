r0=zeros(21,1);
rhs=zeros(21,1);
r0(1,1)=0.05;

for j=2:21
    r0(j,1)=r0(j-1,1)+0.01;
end

k1=1.99925e01;
k2=-25.5267;
k3=465251;


for j=1:21
    rhs(j,1)=((k1)*(exp(k2*r0(j,1))))+k3;
end

plot(r0,rhs,"LineWidth",3)