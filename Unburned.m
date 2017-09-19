
for Texposed=60:-0.5:44 %Set The Dirichlet BC to be input at the model and loop in order to change temperature in a descending linear order by -0.5 degrees. In order to decrease the computation time it was found that the BC is close to 50, thus the iteration can start from 60:-0.5:44.
     [Temperature, Gama]=TissueBurnFullModel(10000,Texposed); %Run the model for the Dirichlet BC stated above
     if Gama < 1 %when Gama becomes less than one, second degree burn does not occur
        fprintf('Temperature at x=0 must be %3.1f Celsius \n', Texposed) %The BC condition for which safety is ensured is displayed
        break %the Loop breaks and the programm ends when the BC is found
     end
end