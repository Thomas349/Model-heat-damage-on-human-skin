function [ StifMat ] = LocStiffnessMatrix(D,L,eID, msh)
%function that calculates the local stiffness matrix K for our DE

gq = CreateGQScheme(1); %Call this function in order to create the appropriate Gaussian Quadrature.
% N=1 for a polynomial of degree one

LocElMat = zeros (2); % Define a local 2x2 matrix to represent the local element matrix
LocElMat(1,1) = gq.gsw(1)*D/(4*(msh.elem(eID).J)); % Substitude the 1st row 1st column with value Int00. Instead of using the integral the GQScheme summation formula is used.
%the corresponding Gauss point and Gauss weight(defined in "CreateGQScheme") is taken in to account. Thus for our case there is no î, but only weighting which is equal to gq.gsw (1) for a first degree polynomial
LocElMat(1,2) = -gq.gsw(1)*D/(4*(msh.elem(eID).J));% Substitude the 1st row 2nd column with value Int01
LocElMat(2,1) = -gq.gsw(1)*D/(4*(msh.elem(eID).J));% Substitude the 2nd row 1st column with value Int10
LocElMat(2,2) = gq.gsw(1)*D/(4*(msh.elem(eID).J)); % Substitude the 2nd row 2nd column with value Int11

GQ = CreateGQScheme(2); %Call function to create Gaussian Quadrature terms but for a second degree polynomial this time
LinearReact = zeros (2,2); % Define local 2x2 matrix to represent the local reaction matrix
LinearReact(1,1) = (GQ.gsw(1)*L*msh.elem(eID).J*((1-GQ.xipts(1))/2)^2)+(GQ.gsw(2)*L*msh.elem(eID).J*((1-GQ.xipts(2))/2)^2); %This time the integral is of second order and thus the summation formula of GQScheme will be used. The function in used within the summation terms is the basis function. However, position of î are replaced with corresponding GQ terms for a second order polynomia, i.e. GQ.xipt(1), GQ.xipt (2).
LinearReact(1,2) = (L*msh.elem(eID).J*((1-GQ.xipts(1))/2)*((1+GQ.xipts(1))/2))*GQ.gsw(1)+(L*msh.elem(eID).J*((1-GQ.xipts(2))/2)*((1+GQ.xipts(2))/2))*GQ.gsw(2);% Substitude the 1st row 2nd column with value Int01. This is true since weights are equal, i.e. GQ.gsw(1)=GQ.gsw(2)=1
LinearReact(2,1) = (L*msh.elem(eID).J*((1-GQ.xipts(1))/2)*((1+GQ.xipts(1))/2))*GQ.gsw(1)+(L*msh.elem(eID).J*((1-GQ.xipts(2))/2)*((1+GQ.xipts(2))/2))*GQ.gsw(2);% Substitude the 2nd row 1st column with value Int10
LinearReact(2,2) = (GQ.gsw(1)*L*msh.elem(eID).J*((1-GQ.xipts(1))/2)^2)+(GQ.gsw(2)*L*msh.elem(eID).J*((1-GQ.xipts(2))/2)^2);% Substitude the 2nd row 2nd column with value Int11

StifMat = LocElMat-LinearReact;% The local stiffnes matrix is given by the difference of the local elemental matrix and linear reaction matrix 
end