function z=f_repulsion(r) 
    cutoff=0.2;
    U=8.0;
    k=2.;
    z=0;
    if r<cutoff
        z= U/k * exp(-r/k); 
    end
end 