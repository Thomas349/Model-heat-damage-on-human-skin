function [ SourceVec ] =LocSourceVector(F,eID,msh)
%function that creates the Local Source Vector
GQ=CreateGQScheme(2); %Integration is done by using the GQ Scheme

Flocal=zeros(2,1); %Initialize the local source Vector

Flocal(1,1)=(GQ.gsw(1)*F*msh.elem(eID).J*((1-GQ.xipts(1))/2)); %Assign the first integral by using the GQ summation formula
Flocal(2,1)=(GQ.gsw(2)*F*msh.elem(eID).J*((1-GQ.xipts(2))/2)); %Assign the second integral by using the GQ summation formula
SourceVec=Flocal;
end