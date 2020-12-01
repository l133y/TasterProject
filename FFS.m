tic %Begin timing the run time

%Rate of transitions between I and S
gamma = 1;
beta = 1.2;

%Total rate of transitions
totalrate = 0;


uniform = 0;    %Uniform random variable
counter = 0;    %Counter variable
FromA = 1;      %Whether the first stage simulation has been back to A since the previous crossing

factor = 3;                             %Scaling the size of the system
N = 100;                                %Total number of individuals
N=N*factor;                             %Scaling up the number of individuals
lambdaA = ceil(N*(1 - (gamma/beta)));   %The interfaces at the end of the steady state A
I =  lambdaA + 3*factor ;               %The initial number of infected
t = 0;                                  %Starting from t=0
flux = 0;                               %Initialising the flux from A to the first interface
lambda = [lambdaA:-3*factor:0];         %The interfaces between A and B
if lambda[length(lambda)] ~= 0;         
    lambda = [lambda,0];                %Make sure that 0 is included as an interface
end


%First stage of Forward Flux Sampling, requires we cross the first
%interface a certain number of times before moving onto the second stage
while (counter < 50)
    totalrate = (((beta/N)*I*(N-I))+(gamma*I));             %Total rate of transitions
    D = exprnd(1/totalrate);                                %Time to next transition
    t = t + D;                                              %Total time 
    uniform = rand;                                         %Generates a uniform random variable
    if uniform >= ((((beta/N)*I*(N-I)))/totalrate)          %Decides the next transition
        I = I - 1;                                          %Transition from I to S
    else
        I = I + 1;                                          %Transition from S to I
    end
    if (I == lambdaA - 1)                                   %Crossing of the first interface
        if (FromA == 1)
            counter = counter + 1;                          %Increase the counter
            FromA = 0;                                      %Must return to A again
        end
    end
    if (I == lambdaA)                                       %Simulation has returned to A
        FromA = 1;                                          %Update the status of the simulation
    end
    if (I == 0)                                             %If transitioned to B
        counter = 50;                                       %Set the condition to end the loop
        disp("I = 0")                                       %Display that the first stage has
        disp("------------------")                          %failed
        flux = 1                                            %Set the condition to end
    end
    
end

%Second stage of Forward Flux Sampling
if (flux == 0)                                              %If the first stage has succeeded
    flux = (counter/t);                                     %Define the flux
    prob = zeros(1,length(lambda)-1);                       %The probabilities of reaching 
                                                            %the next interface

    for j = 1:length(lambda)-1                              %Loop over the interfaces
        I = lambda(j);                                      %Initial value of the starting interface
        successes = 0;                                      %Successful trial runs
        failures = 0;                                       %Unsuccessful trial runs
       while (successes < 100)                              %Until a certain number of successes
            totalrate = (((beta/N)*I*(N-I))+(gamma*I));     %Set the total rate
            D = exprnd(1/totalrate);                        %Time between transitions
            t = t + D;                                      %Total time
            uniform = rand;                                 %Uniform random variable
            if (uniform < ((((beta/N)*I*(N-I)))/totalrate)) %Decides the next transition
                I = I + 1;                                  %Transition from S to I
            else
                I = I - 1;                                  %Transition from I to S
            end
            if (I == lambda(j+1))                           %If we reach the next interface
                successes = successes + 1;                  %Successful trial run
                I = lambda(j);                              %Restart from the previous interface
            elseif (I >= lambdaA)                           %If we return to A
                I = lambda(j);                              %Restart from the previous interface
                failures = failures + 1;                    %Unsuccessful trial run
            end
        end
        prob(j) = successes/(successes+failures);           %Probability of reaching the next interface
    end

    disp(prob)                          %Display the probabilities
    fluxAB = flux * prod(prob)          %Flux from A to B
    extinctionT = 1/fluxAB              %Extinction time

end
disp(lambda)    %Display the interfaces

toc %Finish timing the run time
