function [MassMat]= LocMassMatrix ( eID, msh )
%function that return the local mass matrix
LocMaMat=zeros(2,2); %Initialise the local mass matrix
GQ=CreateGQScheme(2); %Use a GQ scheme of 2

LocMaMat(1,1)=(GQ.gsw(1)*msh.elem(eID).J*((1-GQ.xipts(1))/2)^2)+(GQ.gsw(2)*msh.elem(eID).J*((1-GQ.xipts(2))/2)^2); %Replace the integral with the GQ summation function of different gaussian weight and point values
LocMaMat(1,2)=(msh.elem(eID).J*((1-GQ.xipts(1))/2)*((1+GQ.xipts(1))/2))*GQ.gsw(1)+(msh.elem(eID).J*((1-GQ.xipts(2))/2)*((1+GQ.xipts(2))/2))*GQ.gsw(2); %Which true for this sceme only since gsw(1)=gsw(2)=1
LocMaMat(2,1)=(msh.elem(eID).J*((1-GQ.xipts(1))/2)*((1+GQ.xipts(1))/2))*GQ.gsw(1)+(msh.elem(eID).J*((1-GQ.xipts(2))/2)*((1+GQ.xipts(2))/2))*GQ.gsw(2);
LocMaMat(2,2)=(GQ.gsw(1)*msh.elem(eID).J*((1-GQ.xipts(1))/2)^2)+(GQ.gsw(1)*msh.elem(eID).J*((1-GQ.xipts(2))/2)^2);
MassMat=LocMaMat;
end