function [astar, hist_f] = binary_montecarlo(func, n, niter)
    % [best,hist_f] = binary_montecarlo(func,n,iters)
    % Performs a binary Monte Carlo search. Given objective func,
    % bitstring length n, and number of iterations, this algorithm
    % will try to find the bitstring that maximizes func.
    hist_f = zeros(1,niter); 
    astar = rand(n,1) > 0.5;    
    f_best = func(astar);   
    hist_f(1) = f_best;    
    for (i=2:niter)              
        x = rand(n,1) > 0.5;       
        f_x = func(x);          
        if (f_x >= f_best)        
            astar = x;
            f_best = f_x;
        end
        hist_f(i) = f_best;
    %    plot(hist_f)
   %     drawnow()
    end
end

