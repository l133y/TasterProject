tic
repeats = 1000;

gamma = 1;
beta = 1.2;
totalrate = 0;
uniform = 0;
extinction = zeros(repeats,1);
counter = 0;

factor = 3;
N = 100;
N=N*factor;
lambdaA = ceil(N*(1 - (gamma/beta)));
I =  lambdaA + 3*factor ;
t = 0;
flux = 0;
lambda = [lambdaA:-3*factor:0];




while (counter < 50)
    %disp(I)
    totalrate = (((beta/N)*I*(N-I))+(gamma*I));
    D = exprnd(1/totalrate);
    t = t + D;
    uniform = rand;
    if or(uniform >= ((((beta/N)*I*(N-I)))/totalrate),I==N)
        I = I - 1;
    else
        I = I + 1;
    end
    if (I == lambdaA - 1)
        %Only when going the right way!!!
        counter = counter + 1;
    end
    if (I == 0)
        counter = 50;
        disp("I = 0")
        disp("------------------")
        flux = 1
    end
    
end

if (flux == 0)
    flux = (counter/t);


    if lambda[length(lambda)] ~= 0;
        lambda = [lambda,0];
    end
    prob = zeros(1,length(lambda)-1);

    for j = 1:length(lambda)-1
        I = lambda(j);
        successes = 0;
        failures = 0;
       while (successes < 100)
            totalrate = (((beta/N)*I*(N-I))+(gamma*I));
            D = exprnd(1/totalrate);
            t = t + D;
            uniform = rand;
            if (uniform < ((((beta/N)*I*(N-I)))/totalrate))
                I = I + 1;
            else
                I = I - 1;
            end
            if (I == lambda(j+1))
                successes = successes + 1;
                I = lambda(j);
            elseif (I >= lambdaA)
                I = lambda(j);
                failures = failures + 1;
            end
        end
        prob(j) = successes/(successes+failures);
    end

    disp(prob)
    fluxAB = flux * prod(prob)
    extinctionT = 1/fluxAB

end
disp(lambda)

toc