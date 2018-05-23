function [parent, hist_f] = binary_evolution(func, n, niter)
    p = 1/n;
    hist_f = zeros(1,niter);
    parent = rand(n,1) > 0.5;
    f_parent = func(parent);
    hist_f(1) = f_parent;
    for i = 2:niter
        bits = (rand(n,1) < p);
        x = mod (parent + bits, 2);
        f_x = func (x); 
        if (f_x >= f_parent)
            parent = x;
            f_parent = f_x;
        end
        hist_f(i) = f_parent;
    end
end