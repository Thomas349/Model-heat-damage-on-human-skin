function [TestSolution, ErrorF]= TransientSolution(N)
%function that solves the simple equation with D=1, L=0, f=0
msh = OneDimLinearMeshGen(0,1,10); %Create a simple mesh os 10 element equally spaced
theta = 0.5;
dt = 0.01; %Define a time step
Ne = 10; %Number of element
D = 1; %Diffusion coefficient
L = 0;

GM = zeros(Ne+1); %Initialise Global Matrix
M = zeros(Ne+1); %Initialise Global Mass Matrix
K = zeros(Ne+1); %Initialise Global Stiffness Matrix
GV = zeros(Ne+1,1); %Initialise Global Vector

Ccurrent = zeros (Ne+1,1); % Define a solution vector for this timestep
Cnext = zeros (Ne+1,1); %Define a solution vector for next timestep

Solution = zeros(Ne+1,N+1); %Initialise the matrix the matrix where solutions will be presented
Solution(:,1) = Ccurrent; %Set First column of solution matrix equal to boundary conditions (=Ccurrent)
ErrorF=zeros(12,101);
jdx = 2; %Set a parameter that will change for each iteration of time loop in order to be used for storing purposes

for j = dt:dt:N*dt %Loop over time with timestep dt. 0 is not considered as at 0 boundary condition are applied
    
    for i = 1:Ne %Loop over elements

        Kloc1 = LocStiffnessMatrix(D,L,i,msh); %Calculate the local stiffnes matrix for the element examined
        Mloc1 = LocMassMatrix (i, msh); %Calculate the local mass matrix for the element examined
        if i==1 %for only the first element
            Kloc2 = LocStiffnessMatrix(D,L,i+1,msh); %Calculate the local Stiffness matrix for second element of that examined 
            K(1,1:2) = Kloc1(1,1:2); %Place at the first row and first two columns of Global stiffness matrix the first row from the elemental stiffnes matrix

            K(i+1,i) = Kloc1(2,1); %Place at the second row first column of Global stiffness the first value from the second row of local stiffnes
            K(i+1,i+1) = Kloc1(2,2)+Kloc2(1,1);%At the second row second column add the second column second row value from local stiffnes matrix of this element with the same of the next element to be examined
            K(i+1,i+2) = Kloc2(1,2);%At second row third column of Global stiffness, place the first row last column value of local stiffnes

            Mloc2 = LocMassMatrix(i+1,msh); %Calculate the local mass matrix for second element of the one being examined
            M(1,1:2) = Mloc1(1,1:2); %Place at the first row and first two columns of Global Mass matrix the first row of the local mass matrix

            M(i+1,i) = Mloc1(2,1);%Place at second row first column of Global mass the first value from second row of local mass
            M(i+1,i+1) = Mloc1(2,2)+Mloc2(1,1); %At the second row second column add the second column row value from local mass matrix of this element with the same of the next element to be examined
            M(i+1,i+2) = Mloc2(1,2);%At second row third column of Global Mass, place the first row last column value of local mass

        elseif i==Ne %for final element
            K(Ne+1,Ne:Ne+1) = Kloc1(2,:);%Place in the Global stiffness matrix at last row last two columns, the last row from the local stiffness matrix
            M(Ne+1,Ne:Ne+1) = Mloc1(2,:);%Place in the Global mass matrix at last row last two column, the last row from the local mass matrix
        
        else %for any other element inbetween
            Kloc2 = LocStiffnessMatrix(D,L,i+1,msh);%Calculate the local stiffnes matrix for the next element of that examined
            K(i+1,i) = Kloc1(2,1); %Place in the Global stiffness matrix in the correct position the details for the element examined from its local stiffness matrix
            K(i+1,i+1) = Kloc1(2,2)+Kloc2(1,1); %Add at the correct position in the Global stiffness matrix, the second row second column value of the local stiffness matrix with the first row first column of the local stiffness of the next element to be examined
            K(i+1,i+2) = Kloc2(1,2); %Place in the Global stiffness matrix in the correct position the details for the next element to be examined from the local stiffness matrix

            Mloc2 = LocMassMatrix(i+1,msh);%Calculate the local mass matrix for the next element of that examined
            M(i+1,i) = Mloc1(2,1); %Place in the Global Mass matrix in the correct position the details for the element examined from its local mass matrix
            M(i+1,i+1) = Mloc1(2,2)+Mloc2(1,1); %Add at the correct position in the Global mass matrix, the second row second column value of the local mass matrix with the first row first column of the local mass of the next element to be examined
            M(i+1,i+2) = Mloc2(1,2);%Place in the Global mass matrix in the correct position the details for the next element to be examined from the local mass matrix
        end
    end

GM = M+theta*dt*K; %Calculate the Global Matrix according to equation that theory provides
GM(1,1:2) = [1 0]; %Prepare first row of Global Matrix for Dirichlet Boundary Conditions
GM(Ne+1,Ne:Ne+1) = [0 1]; %Prepare last row of Global Matrix for Dirichlet Boundary Conditions
GV = (M-(1-theta)*dt*K)*Ccurrent; %Calculate the Global Vector Accordind to equation that theory provides. Ccurrent represents the solution of previous time iteration apart from first time where it consists of BC's
GV(1,1) = 0; %Set Dirichlet Boundary Conditions in the first row of Global Vector
GV(end,1) = 1;%Set Dirichlet Boundary Conditions in the last row of Global Vector

Cnext = GM\GV %Inverst Global Matrix to find Global solution vector for this time step
Solution(:,jdx) = Cnext; %Store the solution as a vector in the a unique column in solution matrix
Ccurrent = Cnext %Set Condition for Ccurrent to be followed for the next iteration.

for i=1:Ne %Loop over time
        X0=TransientAnalyticSoln(msh.nvec(1,i),j); %create the local analytical solution for the examinable element at the beginning nodal position
        X1=TransientAnalyticSoln(msh.nvec(1,i+1),j); %create the local analytical solution for the examinable element at the final nodal position
        C0=Cnext(i,1); %Assign the local numerical solution for the examinable element at the beginning nodal position
        C1=Cnext(i+1,1); %Assign the local numerical solution for the examinable element at the final nodal position

        GQ=CreateGQScheme(2); %Create a GQ scheme of order 2, N=2

        Cex(1)=X0*(1-GQ.xipts(1))/2+X1*(1+GQ.xipts(1))/2; %Create the first value of the first term of the integrational sum
        C(1)=  C0*(1-GQ.xipts(1))/2+C1*(1+GQ.xipts(1))/2; %Create the second value of the first term of the integrational sum

        Cex(2)=X0*(1-GQ.xipts(2))/2+X1*(1+GQ.xipts(2))/2; %Create the first value of the second term of the integrational sum
        C(2)=  C0*(1-GQ.xipts(2))/2+C1*(1+GQ.xipts(2))/2; %Create the second value of the second term of the integrational sum

        E=GQ.gsw(1)*(Cex(1)-C(1))^2*msh.elem(i).J+GQ.gsw(2)*(Cex(2)-C(2))^2*msh.elem(i).J; %Evaluate the integral by using the GQ summation method
        Error(i,jdx-1)=E; %Store the Error at a big matrix. Each cell of the matrix will contain an error for one specific element in spacetime. Each column will represent the error of different time increments since the formula is within the time loop. Row at each column will represent the error along the space domain.
    
end
Error(12,jdx-1)=sqrt(sum(Error(:,jdx-1)));%Within the matrix create an additional row that will contain the L2 norm value of each time increment. The Error for each time increment is calculated by taking the sum of each elemental error in the space domain

GM = zeros(Ne+1); %Re-initialise Global Matrix
M = zeros(Ne+1); %Re-initialise Mass Matrix
K = zeros(Ne+1); %Re-initialise Stiffnes Matrix
GV = zeros(Ne+1,1); %Re-initialise Global Vector

jdx = jdx+1; %Set the next condition to this paramaeter to be used for storing purposes
end
TestSolution = Solution; %Assign the function's output to the matrix containing the at each column the corresponding Global vector solution
ErrorF=Error; %Assign the error as an output of the function

%  x=msh.nvec;
%  z1=Solution(:,6);
%  z2=Solution(:,11);
%  z3=Solution(:,31);
%  z4=Solution(:,end);
%  plot(x,z1,x,z2,x,z3,x,z4);
%  xlabel('x');
%  ylabel('c(x,t)');
%  title('Transient Solution for different time incidents');
%  legend('t=0.05','t=0.1','t=0.3','t=1');

     idx=1;
     for i=0:dt:N*dt
     P(idx,1)=TransientAnalyticSoln(0.8,i);
     idx=idx+1;
     end
   
     C=Solution(9,:);
 t=0:dt:N*dt
 PLER=ErrorF(end,:);
  figure  
     plot(t,C,t,P);
     xlabel('time/ (s)');
     ylabel('c(x,t)');
     title('Analytical and Numerical Solution for x=0.8 and time domain [0-1]s');
     legend('Numerical','Analytical');

  hold on
  plot (t, PLER);
  legend('Backward Euler','Crank-Nicolson');
  xlabel('Time (sec)');
  ylabel('L2-Norm');
  title('L2-Norm error distribution in time domain');
  Area=trapz(t,PLER)

end