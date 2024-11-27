function z=energy_angular(b1,b2,a1,a2,c1,c2)
    dot_pcelluct = (b1-a1)*(c1-a1)+(b2-a2)*(c2-a2);
    ab_length = sqrt((b1-a1)^2+(b2-a2)^2);
    ac_length = sqrt((c1-a1)^2+(c2-a2)^2);
    if (ab_length*ac_length>0)
        z= dot_pcelluct / (ab_length*ac_length) + 1;
    else
        z=0;
    end
end