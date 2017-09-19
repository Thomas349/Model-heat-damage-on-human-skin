function [ gq ] = CreateGQScheme(N)
%CreateGQScheme Creates Gaussian Quadrature Scheme of polynomials of orderN<4
%   Creates and initialises a data structure 
    gq.npts = N;
    if (N > 0) && (N < 4)
        %order of quadrature scheme i.e. %number of Gauss points
        gq.gsw = zeros(N,1); %structural array of Gauss weights
        gq.xipts = zeros(N,1); %Structural array of Gauss points
        switch N
            case 1 %For First Degree Polynomial
              gq.gsw(1) = 2; %weight; 
              gq.xipts(1) = 0; %value of î;
            case 2 %2nd degree
              gq.gsw(1) = 1;  
              gq.gsw(2) = 1; 
              gq.xipts(1) = -(1/3)^0.5;
              gq.xipts(2) = (1/3)^0.5;
            case 3 %3rd degree
              gq.gsw(1) = 5/9;  
              gq.gsw(2) = 8/9;
              gq.gsw(3) = 5/9;
              gq.xipts(1) = -(3/5)^0.5;
              gq.xipts(2) = 0;
              gq.xipts(3) = (3/5)^0.5;
        end        
    else % If input is larger than 3
      fprintf('Invalid number of Gauss points specified');
    end 
end

