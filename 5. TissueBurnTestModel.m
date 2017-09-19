function [Temperature, Gama]= TissueBurnTestModel(N, Thealthy)
msh = OneDimLinearMeshGen(0,0.01,18); %generate a mesh of 18 elements of equal length
theta = 1; %Backward Euler
dt=0.005; %Set time step
Ne=18; %number of elements
L=0;
f=0;

GM=zeros(Ne+1);
M=zeros(Ne+1);
K=zeros(Ne+1);
GV=zeros(Ne+1,1);

Ccurrent = zeros (Ne+1,1);
Ccurrent(1,1)=Thealthy; %The temperature causing the burn is set here
Ccurrent(2:end,1)=37;   %The temperature along the thickness of the skin remains constant at 37 degrees
Cnext = zeros (Ne+1,1);
Solution=zeros(Ne+1,N+1);
Solution(:,1)=zeros(Ne+1,1);

jdx=2;

for j=dt:dt:N*dt
    
    for i=1:Ne
        
        Mloc1=LocMassMatrix (i, msh);
        
        if i==1 %for the first element
            D=1/158400; %epidermis diffusion coefficient is set
            Kloc1 = LocStiffnessMatrix(D,L,i,msh);%Calculate the local stiffnes matrix for the element examined
            Kloc2=LocStiffnessMatrix(D,L,i+1,msh);%Calculate the local Stiffness matrix for second element of that examined
            K(1,1:2)=Kloc1(1,1:2);%Place at the first row and first two columns of Global stiffness matrix the first row from the elemental stiffnes matrix

            K(i+1,i)=Kloc1(2,1);%Place at the second row first column of Global stiffness the first value from the second row of local stiffnes
            K(i+1,i+1)=Kloc1(2,2)+Kloc2(1,1);%At the second row second column add the second column second row value from local stiffnes matrix of this element with the same of the next element to be examined
            K(i+1,i+2)=Kloc2(1,2);%At second row third column of Global stiffness, place the first row last column value of local stiffnes
            
            Mloc2=LocMassMatrix(i+1,msh);%Calculate the local mass matrix for second element of the one being examined
            M(1,1:2)=Mloc1(1,1:2);%Place at the first row and first two columns of Global Mass matrix the first row of the local mass matrix

            M(i+1,i)=Mloc1(2,1);%Place at second row first column of Global mass the first value from second row of local mass
            M(i+1,i+1)=Mloc1(2,2)+Mloc2(1,1);%At the second row second column add the second column row value from local mass matrix of this element with the same of the next element to be examined
            M(i+1,i+2)=Mloc2(1,2);%At second row third column of Global Mass, place the first row last column value of local mass

         elseif i==2%for second tile only
            D=1/158400; %Set diffusion epidermis coefficient
            Kloc1 = LocStiffnessMatrix(D,L,i,msh);
            Kloc2=LocStiffnessMatrix(D,L,i+1,msh);
            K(i+1,i)=Kloc1(2,1);
            K(i+1,i+1)=Kloc1(2,2)+Kloc2(1,1);
            K(i+1,i+2)=Kloc2(1,2);

            Mloc2=LocMassMatrix(i+1,msh);
            M(i+1,i)=Mloc1(2,1);
            M(i+1,i+1)=Mloc1(2,2)+Mloc2(1,1);
            M(i+1,i+2)=Mloc2(1,2);
         
        elseif i==3%For third tile
            
            Kloc1 = LocStiffnessMatrix((1/158400),L,i,msh); %input D for Epidermis
            Kloc2=LocStiffnessMatrix((1/99000),L,i+1,msh); %input D for Dermis
            K(i+1,i)=Kloc1(2,1);
            K(i+1,i+1)=Kloc1(2,2)+Kloc2(1,1);
            K(i+1,i+2)=Kloc2(1,2);

            Mloc2=LocMassMatrix(i+1,msh);
            M(i+1,i)=Mloc1(2,1);
            M(i+1,i+1)=Mloc1(2,2)+Mloc2(1,1);
            M(i+1,i+2)=Mloc2(1,2);

        elseif 3<i && i<9 % For 4th untill 8th tile
            D=1/99000; % Set D for Dermis layer
            Kloc1 = LocStiffnessMatrix(D,L,i,msh);
            Kloc2=LocStiffnessMatrix(D,L,i+1,msh);
            K(i+1,i)=Kloc1(2,1);
            K(i+1,i+1)=Kloc1(2,2)+Kloc2(1,1);
            K(i+1,i+2)=Kloc2(1,2);

            Mloc2=LocMassMatrix(i+1,msh);
            M(i+1,i)=Mloc1(2,1);
            M(i+1,i+1)=Mloc1(2,2)+Mloc2(1,1);
            M(i+1,i+2)=Mloc2(1,2);
     
        elseif i==9 %for 9th tile
            Kloc1 = LocStiffnessMatrix((1/99000),L,i,msh); %set D for Dermis layer
            Kloc2=LocStiffnessMatrix((1/198000),L,i+1,msh); %set D for Sub-cutaneous layer
            K(i+1,i)=Kloc1(2,1);
            K(i+1,i+1)=Kloc1(2,2)+Kloc2(1,1);
            K(i+1,i+2)=Kloc2(1,2);

            Mloc2=LocMassMatrix(i+1,msh);
            M(i+1,i)=Mloc1(2,1);
            M(i+1,i+1)=Mloc1(2,2)+Mloc2(1,1);
            M(i+1,i+2)=Mloc2(1,2);

        elseif 9<i && i<Ne %for the 9th tile until the pre-last
            D=1/198000; %Set D for sub-cutaneous           
            Kloc1 = LocStiffnessMatrix(D,L,i,msh);
            Kloc2=LocStiffnessMatrix(D,L,i+1,msh);
            K(i+1,i)=Kloc1(2,1);
            K(i+1,i+1)=Kloc1(2,2)+Kloc2(1,1);
            K(i+1,i+2)=Kloc2(1,2);

            Mloc2=LocMassMatrix(i+1,msh);
            M(i+1,i)=Mloc1(2,1);
            M(i+1,i+1)=Mloc1(2,2)+Mloc2(1,1);
            M(i+1,i+2)=Mloc2(1,2);

        elseif i==Ne % for last tile
            D=1/198000;
            Kloc1 = LocStiffnessMatrix(D,L,i,msh);
            K(Ne+1,Ne:Ne+1)=Kloc1(2,:);
            M(Ne+1,Ne:Ne+1)=Mloc1(2,:);
        end
    end
GM=M+theta*dt*K;
GM(1,1:2)=[1 0]; % prepare first row of Global Matrix for Dirichlet BC
GM(Ne+1,Ne:Ne+1)=[0 1]; %prepare last row of Global Matrix for Dirichlet BC
GV=(M-(1-theta)*dt*K)*Ccurrent; %Calculate the Global Vector Accordind to equation that theory provides. Ccurrent represents the solution of previous time iteration apart from first time where it consists of BC's. Equation reduces to this form since there is no source effect.
GV(1,1)=Thealthy; %Set Dirichlet BC
GV(end,1)=37; %Set Dirichlet BC

Cnext=GM\GV;

Solution(:,jdx)=Cnext;

Ccurrent=Cnext;

GM=zeros(Ne+1);
M=zeros(Ne+1);
K=zeros(Ne+1);
GV=zeros(Ne+1,1);

jdx=jdx+1;
end
Temperature=Solution; %Assign the function's output to the matrix containing at each column the corresponding Global vector solution 
Temperature(1,1)=Thealthy; %set the initial conditions at first column of the matrix containing the global solution vector
Temperature(2:end,1)=37; %set the initial conditions at first column of the matrix containing the global solution vector

Tgama = Temperature (4, :); %Create an array containing all the temperatures specifically for mesh node of x=E
Tgama = Tgama(Tgama>44); %Eliminate from the row any value that is less than 44 Celsius. As temperature increases in a linear rate, once 44 is exceed at the array all remaining value will be greater than 44.
Gama = trapz(2*10^98*exp(-12017./Tgama))*dt; %Formula for calculating the gama integral by using the trapezoidal rule.

%   PlotTemp=zeros(Ne+1,101); %Create a smaller matrix than the functions output in order to be plotted
%     for i=1:101  %loop of the numer of time instances that will be presented
%         j=(i*100)-99; % Create a parameter that will extract from the Temperature matrix the Global solution vector corresponding to the time instance selected
%         PlotTemp(:,i)=Temperature(:,j); %fill the new matrix to be plotted with specific solutions from the Temperature matrix
%     end
%    
%   %Assign the correct values to the axis
%   x=0:dt*100:N*dt; %Set the time domain for wich we are plotting the results (every 0.5sec)
%   y=msh.nvec; %Set the space domain for which we are plotting the results (skin thickness)
%   [t,X]=meshgrid(x,y); %create a mesh grid for the new space-time relation representing the actual solution of the equation
%   surf(t,X,PlotTemp) %Create a 3D plot of the spatial Temperature distribution along with corrected time assigned
%   title('Spatial temperature distribution in tissue for t=0-50sec');
%   ylabel('Skin Thickness (x)');
%   xlabel('Time (sec)');
%   zlabel('Temperature (degrees)');
%   c=colorbar;
%   c.Label.String='Temperature in Celsius';
end
