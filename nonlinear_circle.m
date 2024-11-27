function z=nonlinear_circle(r,nu)
    if r>0.5
        z=0.1*r*(r-1)*(r-nu);
    else 
        z=10 * r*(r-1)*(r-nu);
    end
end