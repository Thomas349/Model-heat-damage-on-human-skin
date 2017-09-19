function [Temperature, Gama]= TissueBurnFullModel(N,Ne,Texposed)
    msh = OneDimLinearMeshGen(0,0.01,Ne);
    theta = 1;
    dt=0.005;
    L=0;

    GM=zeros(Ne+1);
    M=zeros(Ne+1);
    K=zeros(Ne+1);
    GV=zeros(Ne+1,1);

    Ccurrent = zeros (Ne+1,1);
    Ccurrent(1,1)=Texposed;
    Ccurrent(2:end,1)=37;
    Cnext=zeros(Ne+1,1);
    Fcurrent=zeros(Ne+1,1); %Create a source vector that will represent the previous time iterations vector
    Fnext=zeros(Ne+1,1); % Create a source vector that will represent the current time examined vector
    Solution=zeros(Ne+1,N+1);
    Solution(:,1)=Ccurrent;

    jdx=2;

    for j=dt:dt:N*dt

        for i=1:Ne

            if i==(1/18)*Ne %Generalise code by allowing a multiple of 18 element to be used
                D=1/158400;
                f=0;
                Kloc1 = LocStiffnessMatrix(D,L,i,msh);
                Kloc2=LocStiffnessMatrix(D,L,i+1,msh);
                K(1,1:2)=Kloc1(1,1:2);

                K(i+1,i)=Kloc1(2,1);
                K(i+1,i+1)=Kloc1(2,2)+Kloc2(1,1);
                K(i+1,i+2)=Kloc2(1,2);

                Mloc1=LocMassMatrix (i, msh);
                Mloc2=LocMassMatrix(i+1,msh);
                M(1,1:2)=Mloc1(1,1:2);

                M(i+1,i)=Mloc1(2,1);
                M(i+1,i+1)=Mloc1(2,2)+Mloc2(1,1);
                M(i+1,i+2)=Mloc2(1,2);

                Floc1=SourcelocalMatrix (f,i,msh);
                Floc2=SourcelocalMatrix(f,i+1,msh);
                Fnext(1,1)= Floc1(1,1); %Assign at the first row value of the global source vector the first row of the the local source vector for the first element
                Fnext(2,1)=Floc1(2,1)+Floc2(1,1); % Assign at the second row of the global source vector the second row of the local source vector for the first element and the first row value from the local source vector of the next element

             elseif i==(2/18)*Ne
                D=1/158400;
                f=0;
                Kloc1 = LocStiffnessMatrix(D,L,i,msh);
                Kloc2=LocStiffnessMatrix(D,L,i+1,msh);
                K(i+1,i)=Kloc1(2,1);
                K(i+1,i+1)=Kloc1(2,2)+Kloc2(1,1);
                K(i+1,i+2)=Kloc2(1,2);

                Mloc1=LocMassMatrix (i, msh);
                Mloc2=LocMassMatrix(i+1,msh);
                M(i+1,i)=Mloc1(2,1);
                M(i+1,i+1)=Mloc1(2,2)+Mloc2(1,1);
                M(i+1,i+2)=Mloc2(1,2);

                Floc1=SourcelocalMatrix (f,i,msh);
                Floc2=SourcelocalMatrix(f,i+1,msh);
                Fnext(i+1,1)=Floc1(2,1)+Floc2(1,1); %assign at each row of the Global source vector the sum of the second value of the local source vector for the element examined and the first row value for the local source vector of the next element to be examined 

            elseif i==(3/18)*Ne
                f=0;
                Kloc1 = LocStiffnessMatrix((1/158400),L,i,msh);
                Kloc2=LocStiffnessMatrix((1/99000),L,i+1,msh);
                K(i+1,i)=Kloc1(2,1);
                K(i+1,i+1)=Kloc1(2,2)+Kloc2(1,1);
                K(i+1,i+2)=Kloc2(1,2);

                Mloc1=LocMassMatrix (i, msh);
                Mloc2=LocMassMatrix(i+1,msh);
                M(i+1,i)=Mloc1(2,1);
                M(i+1,i+1)=Mloc1(2,2)+Mloc2(1,1);
                M(i+1,i+2)=Mloc2(1,2);

                Floc1=SourcelocalMatrix (f,i,msh);
                Floc2=SourcelocalMatrix(((0.0375*1060*3770)/(1200*3300))*37,i+1,msh);
                Fnext(i+1,1)=Floc1(2,1)+Floc2(1,1);
            elseif (3/18)*Ne<i && i<(9/18)*Ne
                D=1/99000;
                L=-(0.0375*1060*3770)/(1200*3300);
                f=-L*37; %Set the value of f for this layer
                Kloc1 = LocStiffnessMatrix(D,L,i,msh);
                Kloc2=LocStiffnessMatrix(D,L,i+1,msh);
                K(i+1,i)=Kloc1(2,1);
                K(i+1,i+1)=Kloc1(2,2)+Kloc2(1,1);
                K(i+1,i+2)=Kloc2(1,2);

                Mloc1=LocMassMatrix (i, msh);
                Mloc2=LocMassMatrix(i+1,msh);
                M(i+1,i)=Mloc1(2,1);
                M(i+1,i+1)=Mloc1(2,2)+Mloc2(1,1);
                M(i+1,i+2)=Mloc2(1,2);

                Floc1=SourcelocalMatrix (f,i,msh);
                Floc2=SourcelocalMatrix(f,i+1,msh);
                Fnext(i+1,1)=Floc1(2,1)+Floc2(1,1);
            elseif i==(9/18)*Ne
                L=-(0.0375*1060*3770)/(1200*3300);
                f=-L*37;%Set value of f for this layer
                Kloc1 = LocStiffnessMatrix((1/99000),L,i,msh);
                Kloc2=LocStiffnessMatrix((1/198000),L,i+1,msh);
                K(i+1,i)=Kloc1(2,1);
                K(i+1,i+1)=Kloc1(2,2)+Kloc2(1,1);
                K(i+1,i+2)=Kloc2(1,2);

                Mloc1=LocMassMatrix (i, msh);
                Mloc2=LocMassMatrix(i+1,msh);
                M(i+1,i)=Mloc1(2,1);
                M(i+1,i+1)=Mloc1(2,2)+Mloc2(1,1);
                M(i+1,i+2)=Mloc2(1,2);

                Floc1=SourcelocalMatrix (f,i,msh);
                Floc2=SourcelocalMatrix(f,i+1,msh);
                Fnext(i+1,1)=Floc1(2,1)+Floc2(1,1);
            elseif (9/18)*Ne<i && i<(Ne/18)*Ne
                D=1/198000;
                f=-L*37; %set value of L for this layer
                L=(0.0375*1060*3770)/(1200*3300);
                Kloc1 = LocStiffnessMatrix(D,L,i,msh);
                Kloc2=LocStiffnessMatrix(D,L,i+1,msh);
                K(i+1,i)=Kloc1(2,1);
                K(i+1,i+1)=Kloc1(2,2)+Kloc2(1,1);
                K(i+1,i+2)=Kloc2(1,2);

                Mloc1=LocMassMatrix (i, msh);
                Mloc2=LocMassMatrix(i+1,msh);
                M(i+1,i)=Mloc1(2,1);
                M(i+1,i+1)=Mloc1(2,2)+Mloc2(1,1);
                M(i+1,i+2)=Mloc2(1,2);

                Floc1=SourcelocalMatrix (f,i,msh);
                Floc2=SourcelocalMatrix(f,i+1,msh);
                Fnext(i+1,1)=Floc1(2,1)+Floc2(1,1);
            elseif i==(Ne/18)*Ne
                D=1/198000;
                L=-(0.0375*1060*3770)/(1200*3300);
                f=-L*37;
                Kloc1 = LocStiffnessMatrix(D,L,i,msh);
                K(Ne+1,Ne:Ne+1)=Kloc1(2,:);

                Mloc1=LocMassMatrix (i, msh);
                M(Ne+1,Ne:Ne+1)=Mloc1(2,:);

                Floc1=SourcelocalMatrix(f,i,msh);
                Fnext(Ne+1,1)=Floc1(2,1);
            end
        end
    GM=M+theta*dt*K;
    GM(1,1:2)=[1 0];
    GM(Ne+1,Ne:Ne+1)=[0 1];
    GV=(M-(1-theta)*dt*K)*Ccurrent+dt*theta*Fnext+dt*(1-theta)*Fcurrent; %Use the full expansion of the equation now that reaction and source are present
    GV(1,1)=Texposed;
    GV(end,1)=37;

    Cnext=GM\GV;

    Solution(:,jdx)=Cnext;

    Ccurrent=Cnext;
    Fcurrent=Fnext;
    GM=zeros(Ne+1);
    M=zeros(Ne+1);
    K=zeros(Ne+1);
    GV=zeros(Ne+1,1);
    Fnext=zeros(Ne+1,1);
    Cnext=zeros(Ne+1,1);

    jdx=jdx+1;
    end

    Temperature=Solution; % Assign output of the function

    Tgama = Temperature (4, :);
    Tgama = Tgama(Tgama>44);
    Gama = trapz(2*10^98*exp(-12017./Tgama))*dt; %Measure Gama factor for x=E

      PlotTemp=zeros(Ne+1,101); %Plot selected data from the solution matrix
       for i=1:101 
           j=(i*100)-99;
           PlotTemp(:,i)=Temperature(:,j);
       end
     x=0:dt*100:N*dt; %Set the time domain for wich we are plotting the results (every 0.5sec)
     y=msh.nvec; %Set the space domain for which we are plotting the results (skin thickness)
     [t,X]=meshgrid(x,y); %create a mesh grid for the new space-time relation representing the actual solution of the equation
     surf(t,X,PlotTemp) %Create a 3D plot of the spatial Temperature distribution along with corrected time assigned
     title('Spatial temperature distribution in tissue for t=0-50sec');
     ylabel('Skin Thickness (x)');
     xlabel('Time (sec)');
     zlabel('Temperature (degrees)');
     c=colorbar;
     c.Label.String='Temperature in Celsius';
end