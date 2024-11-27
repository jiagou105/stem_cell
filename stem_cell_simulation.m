% first clean figures and memory from all previous calculations
clc;
clf;
clear;

% 1. Make proliferation rate dependent on consumption of EV
% 2. EVs are secreted along entire cell, not just nodes 
% 3. Cell-cell adhesions only initiated from cell edges (When cells interact they do not align)
% 4. A cell can become circular randomly, no need to be isolated. No vesicle production when circular
% 5. Add a crowding number to determine whether division happens 
% [* optional] Maybe, use a diffusion PDE for EVs, instead of individual tracking.


% Thoughts:
% 1. what is the comparison of cell size with domain size???
% 2. add a volume to cells. Cells can deform when squeezed, round when
% isolated
% 3. add repulsion and attraction to cells to add volume 
% 4. add mechanisms for polarity determination so the cell shape does have
% to be restricted by the connected springs 


% nohup matlab -nodisplay -nosplash <stem_cell_simulation.m >output.txt &
% https://stackoverflow.com/questions/3000724/running-matlab-in-the-background
% need functions separately defined in files 

number_of_cells = 40; % initial number if cells 
number_of_nodes = 7; % ODD NUMBER ; WE HAVE NUCLEUS IN THE MIDDLE
i_nucleus = int16((number_of_nodes - 1)/2)+1; % WE HAVE NUCLEUS IN THE MIDDLE OF EACH CELL
number_of_vesicles = 100;
D_v = 0.001;%0.025; % spatial diffusion of vesicles
D_circ = 0.0125; % spatial diffusion of circular cells 
switch_biology = 1;

% shape cell parameters 
k_spring = 0.2;
k_angular =0.1;
k_filopodium = 0.125; 
k_adhesion = 0.5;

k_repulsion_circular = 0.001;

L_box = 10;
L_0_eq = 2.0/number_of_nodes; % equal length of each spring
L_vesicles = 0.5*L_0_eq; 
R_bound = (sqrt(2)+0.2) * L_box;
max_cell_length = number_of_nodes * L_0_eq * 1.25; 


dt = 0.01;
frame_frequency = 200;
frame_frequency_0 = 200;%frame_frequency; 
nT = 200*frame_frequency_0+1;%90001;%10001;
dx = 0.00002;
n_repulsion = 4;  %% meaning of this number???


%filopodia parameters 
max_number_of_filopodia = 40;
filopodium_duration_time = 1 * frame_frequency * dt; % last one output frame??
filopodium_average_length = 0.35;
filopodium_birth_rate = switch_biology * 2.0;
q_not_edge = 0.05; % probability of having a filopodium on a side 

%cell division parameters
division_rate = switch_biology * 0.02;%switch_biology * 0.1;
division_age = 5.0; 

%adhesion parameters
max_number_of_adhesions = 40;
adhesion_max_length = 0.2;
adhesion_equilibrium_length = 0.3;
adhesion_regular_time = 400 * frame_frequency * dt;
adhesion_birth_rate = switch_biology * 100.8;

% vesicle parameters
vesicle_production_rate = switch_biology*5.5;
vesicle_life_expectancy = 2*frame_frequency/10 *dt; 


how_far_from_center = 1.0; % to initiate rods closer to the center of the domain
L_0 = L_0_eq;
magnitude_of_randomness_in_initial_conditions = L_0/2.5;

circular_threshold = 0.5;
overlap_positions = 0;
rate_of_circularity = switch_biology * 0.00001;

ev_avglife = frame_frequency/10 *dt;
ev_totallife = [];
ev_num = [];

% initializing or cells 
for nr=1:number_of_cells
    
    phi_initial = 2*pi*rand(); 

    
    
    %x_initial = L_box/2 + how_far_from_center*L_box/2 * (2*rand()-1);
    %y_initial = L_box/2 + how_far_from_center*L_box/2 * (2*rand()-1);
    %TEST
    if nr==1
        phi_initial = phi_initial;
        x_initial = L_box/2;
        y_initial = L_box/2; 
    end
    if nr==2
        phi_initial = 2*pi*rand();
        x_initial = L_box/2.0 - 0.5*(double(i_nucleus)) * L_0;
        y_initial = L_box/2.0 - (double(i_nucleus)-1) * L_0; 
    end
    if nr>2
        x_initial = L_box/2 + how_far_from_center*L_box/2 * (2*rand()-1);
        y_initial = L_box/2 + how_far_from_center*L_box/2 * (2*rand()-1);
    end
        
    for i=1:number_of_nodes
         x_cell (i,nr) = x_initial + (i-1)*L_0*cos(phi_initial)+ magnitude_of_randomness_in_initial_conditions *(2*rand()-1); % 7 rows n columns
         y_cell (i,nr) = y_initial + (i-1)*L_0*sin(phi_initial)+ magnitude_of_randomness_in_initial_conditions *(2*rand()-1);     
    end 
    
    % properties relating to the division
    cell_birth_time(nr) = 0.0 + rand()*100*dt;
    
    % two states circular and elongataed
    how_circular_is_the_cell(nr) = 0;
    is_the_cell_circular(nr) = 0;

    % family members 
    who_is_my_brother(nr)=0;
    time_i_was_born(nr)=0;
    
    % on overlapping 
    for i=1:number_of_nodes
        is_the_node_overlapping(nr,i)=0;
        overlapping_other_cell(nr,i)=0;
    end
    
    % filopodia properties
    for j=1:max_number_of_filopodia
        is_filopodium_exists(nr,j)=0;
        birth_time_of_filopodium(nr,j)=0;
        location_of_filopodium(nr,j)=0; % node's number
        tip_x_filopodium(nr,j)=0; 
        tip_y_filopodium(nr,j)=0;
        activation_of_filopodium(nr,j)=0;
    end
    
    % adhesion properties 
    for j=1:max_number_of_adhesions
        is_adhesion_exists(nr,j)=0;
        birth_time_of_adhesion(nr,j)=0;
        location_of_adhesion(nr,j)=0;
        the_other_cell_number(nr,j)=0;
        the_other_cell_node(nr,j)=0;
        %tip_adhesion_location(nr,j)=42; 
    end
    
    for i=1:number_of_nodes
        number_of_outcoming_adhesions(nr,i)=0;
    end

    % cell length 
    cellLength(nr)=0;
    for i=1:number_of_nodes-1
        cellLength(nr)=cellLength(nr)+sqrt((x_cell(i+1,nr)-x_cell(i,nr))^2+(y_cell(i+1,nr)-y_cell(i,nr))^2);
    end
    equilibrium_length(nr) = L_0_eq;

    % absored vesicle
    ev_totallife(nr) = 0;
    ev_num(nr) = 0;
end




% initiation of vesicles 
for nv=1:number_of_vesicles 
    xv(nv) = L_box * rand();
    yv(nv) = L_box * rand();
    vn(nv) = 0; % vesicle number - the number (label) of the cell that produced this vesicles
    agev(nv) = 0;
end



% checking initial configuration
figure(1)
clf();
for nr=2:number_of_cells
    plot(x_cell(:,nr),y_cell(:,nr),'black'); hold on;
    plot(x_cell(:,nr),y_cell(:,nr),'red.'); hold on;
end

grid off
daspect([1 1 1]);
axis([0 L_box 0 L_box])


% dynamics

for nt=1:nT
    % keep current values in back up ("old") variables 
    x_cell_old = x_cell;
    y_cell_old = y_cell;
    xv_old = xv; % extracellular vesicles 
    yv_old = yv; 
    
    number_of_new_cells = 0;
    overlap_positions = 0;
    
    for nr=1:number_of_cells % for each cell 
        
        % Value of L_0 depends on age
        cell_age = nt*dt - cell_birth_time(nr); % linear with simulation time ??
        %L_0 = L_0_eq /2 *(1+ sigma_f(cell_age / (division_age))); 
        
        if how_circular_is_the_cell(nr)<circular_threshold
            equilibrium_length(nr) = min(L_0_eq, equilibrium_length(nr)+cell_age/(dt*20*frame_frequency)*L_0_eq/2);
            %equilibrium_length(nr) = min(L_0_eq, L_0_eq/2);
        end
       

        if how_circular_is_the_cell(nr)>circular_threshold
            equilibrium_length(nr) = L_0_eq / 5;
        end
       
        L_0 = equilibrium_length(nr); 

        % Computing overlapping 
        for i=1:number_of_nodes
            is_the_node_overlapping(nr,i)=0;
            for nr2=1:number_of_cells
                 if (nr2~=nr)
                     for i2=1:(number_of_nodes-1) % why minus one?
                        distance = sqrt((x_cell_old(i,nr)-x_cell_old(i2,nr2))^2+(y_cell_old(i,nr)-y_cell_old(i2,nr2))^2);
                        % resolve overlaps
                         if (mod(nt,1)==0) && (distance<1.0) && (i~=i_nucleus) % 
                             is_there_an_overlap = 0;
                             if i>i_nucleus 
                                 i_adj = i-1;
                             end

                             if i<i_nucleus
                                 i_adj = i+1; 
                             end

                             % geometric check on overlaps
                             xA = x_cell_old(i,nr);
                             yA = y_cell_old(i,nr);

                             xB = x_cell_old(i_adj,nr);
                             yB = y_cell_old(i_adj,nr);

                             xC = x_cell_old(i2,nr2);
                             yC = y_cell_old(i2,nr2);

                             xD = x_cell_old(i2+1,nr2);
                             yD = y_cell_old(i2+1,nr2);
                             
                             % what is the algorithm here
                             detA = (xB-xA)*(yC-yD) - (xC-xD)*(yB-yA);
                             Delta_1 = (xC-xA)*(yC-yD) - (xC-xD)*(yC-yA);
                             Delta_2 = (xB-xA)*(yC-yA) - (xC-xA)*(yB-yA);

                             if detA ~= 0 
                                 lambda_1 = Delta_1 / detA;
                                 lambda_2 = Delta_2 / detA; 
                                 if (lambda_1>=0)&&(lambda_1<=1)&&(lambda_2>=0)&&(lambda_2<=1)
                                     is_there_an_overlap = 1;
                                     is_the_node_overlapping(nr,i)=1;
                                     overlapping_other_cell(nr,i)=nr2; 
                                     overlap_positions = overlap_positions+1;
                                     x_overlap(overlap_positions)=xA+lambda_1*(xB-xA);
                                     y_overlap(overlap_positions)=yA+lambda_1*(yB-yA);
                                 end
                             end
                         end
                     end
                 end
            end
        end
        
        for i=1:number_of_nodes
            %spring forces
            f_spring_x(i) = 0;
            f_spring_y(i) = 0;
            
            
            
            if i>1
                distance = sqrt((x_cell_old(i,nr)-x_cell_old(i-1,nr))^2+(y_cell_old(i,nr)-y_cell_old(i-1,nr))^2); 
                if (distance>0)
                    if is_the_node_overlapping(nr,i) == 0
                        f_spring_x(i) = f_spring_x(i)-k_spring*(distance - L_0)*(x_cell_old(i,nr)-x_cell_old(i-1,nr))/distance;
                        f_spring_y(i) = f_spring_y(i)-k_spring*(distance - L_0)*(y_cell_old(i,nr)-y_cell_old(i-1,nr))/distance;
                    end
                    if is_the_node_overlapping(nr,i) == 1 && (i>i_nucleus)
                        i_adj = i - 1;
                        xA = x_cell_old(i,nr);
                        yA = y_cell_old(i,nr);
                        xB = x_cell_old(i_adj,nr);
                        yB = y_cell_old(i_adj,nr);

                        f_spring_x(i) = f_spring_x(i) + 5* k_spring *(xB-xA); % strengthen spring, shrink spring to resolve overlap???
                        f_spring_y(i) = f_spring_y(i) + 5* k_spring *(yB-yA);
                    end
                end
            end
                 
            if i<number_of_nodes
               distance = sqrt((x_cell_old(i,nr)-x_cell_old(i+1,nr))^2+(y_cell_old(i,nr)-y_cell_old(i+1,nr))^2); 
               if (distance>0)
                   if is_the_node_overlapping(nr,i) == 0
                        f_spring_x(i) = f_spring_x(i)-k_spring*(distance - L_0)*(x_cell_old(i,nr)-x_cell_old(i+1,nr))/distance;
                        f_spring_y(i) = f_spring_y(i)-k_spring*(distance - L_0)*(y_cell_old(i,nr)-y_cell_old(i+1,nr))/distance;
                   end
                   if is_the_node_overlapping(nr,i) && (i<i_nucleus)
                        i_adj = i + 1;
                        xA = x_cell_old(i,nr);
                        yA = y_cell_old(i,nr);
                        xB = x_cell_old(i_adj,nr);
                        yB = y_cell_old(i_adj,nr);
                        
                        f_spring_x(i) = f_spring_x(i) + 5*k_spring * (xB-xA);
                        f_spring_y(i) = f_spring_y(i) + 5*k_spring * (yB-yA);
                   end
               end
            end
        end
    
        %bending forces
        for i=1:number_of_nodes
            f_bending_x(i)=0;
            f_bending_y(i)=0;
        
            % angle phi(i-1)
            if i>2
                im1=i-1;
                im2=i-2;
        
                e = energy_angular(x_cell_old(im2,nr),y_cell_old(im2,nr),x_cell_old(im1,nr),y_cell_old(im1,nr),x_cell_old(i,nr),y_cell_old(i,nr));
                e_dx = energy_angular(x_cell_old(im2,nr),y_cell_old(im2,nr),x_cell_old(im1,nr),y_cell_old(im1,nr),x_cell_old(i,nr)+dx,y_cell_old(i,nr));
                e_dy = energy_angular(x_cell_old(im2,nr),y_cell_old(im2,nr),x_cell_old(im1,nr),y_cell_old(im1,nr),x_cell_old(i,nr),y_cell_old(i,nr)+dx);
        
                f_bending_x(i) = f_bending_x(i) - k_angular*(e_dx-e)/dx;
                f_bending_y(i) = f_bending_y(i) - k_angular*(e_dy-e)/dx;
            end
        
            % angle phi(i)
            if (i>1)&&(i<number_of_nodes)
                im1=i-1;
                ip1=i+1;
        
                e0 = energy_angular(x_cell_old(im1,nr),y_cell_old(im1,nr),x_cell_old(i,nr),y_cell_old(i,nr),x_cell_old(ip1,nr),y_cell_old(ip1,nr));
                e_dx = energy_angular(x_cell_old(im1,nr),y_cell_old(im1,nr),x_cell_old(i,nr)+dx,y_cell_old(i,nr),x_cell_old(ip1,nr),y_cell_old(ip1,nr));
                e_dy = energy_angular(x_cell_old(im1,nr),y_cell_old(im1,nr),x_cell_old(i,nr),y_cell_old(i,nr)+dx,x_cell_old(ip1,nr),y_cell_old(ip1,nr));
        
                f_bending_x(i) = f_bending_x(i) - k_angular*(e_dx-e0)/dx;
                f_bending_y(i) = f_bending_y(i) - k_angular*(e_dy-e0)/dx;
            end
        
            % angle phi(i+1)
            if (i<number_of_nodes-1)
                ip1=i+1;
                ip2=i+2;
       
                e0 = energy_angular(x_cell_old(i,nr),y_cell_old(i,nr),x_cell_old(ip1,nr),y_cell_old(ip1,nr),x_cell_old(ip2,nr),y_cell_old(ip2,nr));
                e_dx = energy_angular(x_cell_old(i,nr)+dx,y_cell_old(i,nr),x_cell_old(ip1,nr),y_cell_old(ip1,nr),x_cell_old(ip2,nr),y_cell_old(ip2,nr));
                e_dy = energy_angular(x_cell_old(i,nr),y_cell_old(i,nr)+dx,x_cell_old(ip1,nr),y_cell_old(ip1,nr),x_cell_old(ip2,nr),y_cell_old(ip2,nr));
        
                f_bending_x(i) = f_bending_x(i) - k_angular*(e_dx-e0)/dx;
                f_bending_y(i) = f_bending_y(i) - k_angular*(e_dy-e0)/dx; 
            end
        end
        
      
        %steric interaction - repulsion
        for i=1:number_of_nodes
            f_repulsion_x(i) = 0;
            f_repulsion_y(i) = 0;
            for nr2=1:number_of_cells
                 if (nr2~=nr)&&(who_is_my_brother(nr)~=nr2)
                     for i2=1:(number_of_nodes-1)
                         for k2=1:(n_repulsion+1) %% meaning of this algorithm ??? push the whole line segment instant of two points???
                             x_star = x_cell_old(i2,nr2) + (k2-1)/n_repulsion * (x_cell_old(i2+1,nr2)-x_cell_old(i2,nr2));
                             y_star = y_cell_old(i2,nr2) + (k2-1)/n_repulsion * (y_cell_old(i2+1,nr2)-y_cell_old(i2,nr2));
                             distance = sqrt((x_cell_old(i,nr)-x_star)^2+(y_cell_old(i,nr)-y_star)^2);
                             if (distance>0)&&(is_the_node_overlapping(nr,i)==0)
                                 coef_repulsion = 1.0 - how_circular_is_the_cell(nr2)*(1-k_repulsion_circular);
                                 f_repulsion_x(i) = f_repulsion_x(i)+coef_repulsion*f_repulsion(distance)/n_repulsion*(x_cell_old(i,nr)-x_star)/distance;
                                 f_repulsion_y(i) = f_repulsion_y(i)+coef_repulsion*f_repulsion(distance)/n_repulsion*(y_cell_old(i,nr)-y_star)/distance;
                             end
                         end
                     end
                 end
            end
        end
        

    
    
        % filopodia propulsion
        % update filopodia
        % removal filopodia  
        for j=1:max_number_of_filopodia
            % deal with existing filopodium
            if is_filopodium_exists(nr,j)==1
                l = location_of_filopodium(nr,j);
            
                current_length = sqrt((tip_x_filopodium(nr,j)-x_cell_old(l,nr))^2+ (tip_y_filopodium(nr,j)-y_cell_old(l,nr))^2);
                
                % if filopodium is old, then kill it
                if (nt*dt - birth_time_of_filopodium(nr,j))>filopodium_duration_time
                    is_filopodium_exists(nr,j)=0;
                end
                
                % if filopodium is too long, then kill it
                if is_filopodium_exists(nr,j)==1 && current_length>filopodium_average_length
                    is_filopodium_exists(nr,j)=0;
                end
            end
            
            % deal with non-existing filopodium
            if is_filopodium_exists(nr,j)==0 && rand()<filopodium_birth_rate*dt
                is_filopodium_exists(nr,j)=1;
                birth_time_of_filopodium(nr,j)=nt*dt;
                
                % assign new filopodium to a node
                randloc = rand();
                if randloc<0.5 - q_not_edge/2.0 % filopodium emerges at front with probability 1/3
                    l = 1;
                else
                    if randloc>0.5 + q_not_edge/2.0 % filopodium emerges at front with probability 2/3
                        l=number_of_nodes;
                    else 
                        l=randi([2,number_of_nodes-1]);
                    end
                end
                location_of_filopodium(nr,j)=l;
                
                if (l==1) % if the filopodium is at front
                    phi_filopodium = - pi/8 * (2*rand()-1) + atan2(y_cell_old(1,nr)-y_cell_old(2,nr),x_cell_old(1,nr)-x_cell_old(2,nr));
                    filopodium_length = (1+ 0.25*(2*rand()-1)) * filopodium_average_length; 
                    tip_x_filopodium(nr,j)=x_cell_old(1,nr)+filopodium_length*cos(phi_filopodium); 
                    tip_y_filopodium(nr,j)=y_cell_old(1,nr)+filopodium_length*sin(phi_filopodium);
                end
                
                if (l==number_of_nodes) % if the filopodium as at back
                    n=number_of_nodes;
                    phi_filopodium = - pi/8 * (2*rand()-1) + atan2(y_cell_old(n,nr)-y_cell_old(n-1,nr),x_cell_old(n,nr)-x_cell_old(n-1,nr));
                    filopodium_length = (1+ 0.25*(2*rand()-1)) * filopodium_average_length; 
                    tip_x_filopodium(nr,j)=x_cell_old(n,nr)+filopodium_length*cos(phi_filopodium); 
                    tip_y_filopodium(nr,j)=y_cell_old(n,nr)+filopodium_length*sin(phi_filopodium);
                end
                
                if (l>1)&&(l<number_of_nodes) % if the filopodium is somewhere in the middle
                    right_or_left=randi([1,2]);
                    if right_or_left == 1
                        phi_filopodium = - pi/16 * (2*rand()-1) + 0.5*(atan2(y_cell_old(l+1,nr)-y_cell_old(l,nr),x_cell_old(l+1,nr)-x_cell_old(l,nr))+atan2(y_cell_old(l-1,nr)-y_cell_old(l,nr),x_cell_old(l-1,nr)-x_cell_old(l,nr)));
                        filopodium_length = (1+ 0.25*(2*rand()-1)) * filopodium_average_length; 
                        tip_x_filopodium(nr,j)=x_cell_old(l,nr)+filopodium_length*cos(phi_filopodium); 
                        tip_y_filopodium(nr,j)=y_cell_old(l,nr)+filopodium_length*sin(phi_filopodium);
                    else
                        phi_filopodium = pi - pi/16 * (2*rand()-1) + 0.5*(atan2(y_cell_old(l+1,nr)-y_cell_old(l,nr),x_cell_old(l+1,nr)-x_cell_old(l,nr))+atan2(y_cell_old(l-1,nr)-y_cell_old(l,nr),x_cell_old(l-1,nr)-x_cell_old(l,nr)));
                        filopodium_length = (1+ 0.25*(2*rand()-1)) * filopodium_average_length; 
                        tip_x_filopodium(nr,j)=x_cell_old(l,nr)+filopodium_length*cos(phi_filopodium); 
                        tip_y_filopodium(nr,j)=y_cell_old(l,nr)+filopodium_length*sin(phi_filopodium);
                    end
                end
                
                % activation level of filopodium (updated only once during birth)
                activation_of_filopodium(nr,j)=0;
                for nv=1:number_of_vesicles
                     dist = sqrt((tip_x_filopodium(nr,j) - xv(nv))^2+(tip_y_filopodium(nr,j) - yv(nv))^2);
                     if (dist < L_vesicles)&&(vn(nv)~=nr)&&(agev(nv)<vesicle_life_expectancy)&&(vn(nv)~=0)
                          activation_of_filopodium(nr,j)=activation_of_filopodium(nr,j)+1;
                          % if a vesicle is too close, then it is consumed
                          % (deleted from consideration)
                          if (dist < L_vesicles/2.0)
                            vn(nv)=0;
                            agev(nv)=vesicle_life_expectancy+1.0; % still count age to identify those vesicles 
                            ev_totallife(nr) = ev_totallife(nr) + ev_avglife;
                            ev_num(nr) = ev_num(nr) + 1;
                          end
                     end 
                end 
                
            end
            
           
        end
        
        %compute force coming from filopodia
        for i=1:number_of_nodes 
            force_filopodium_x(i)=0;
            force_filopodium_y(i)=0;
            for j=1:max_number_of_filopodia
                if is_filopodium_exists(nr,j)&&(location_of_filopodium(nr,j)==i)
                    if (i>1)&&(i<number_of_nodes)
                        geom_coef_filopodia = 0.1;
                    else 
                        geom_coef_filopodia = 1.0;
                    end
                    nx_aux = tip_x_filopodium(nr,j)-x_cell(i,nr);
                    ny_aux = tip_y_filopodium(nr,j)-y_cell(i,nr);

                    length_aux = sqrt(nx_aux^2+ny_aux^2);

                    if length_aux>0
                        nx_aux = nx_aux / length_aux;
                        ny_aux = ny_aux / length_aux;
                    end


                    force_filopodium_x(i) = force_filopodium_x(i) + geom_coef_filopodia*k_filopodium * sigma_f(activation_of_filopodium(nr,j)+0.1) * nx_aux;
                    force_filopodium_y(i) = force_filopodium_y(i) + geom_coef_filopodia*k_filopodium * sigma_f(activation_of_filopodium(nr,j)+0.1) * ny_aux;
                end
            end
        end
        
        
        %computation of cell-cell adhesion
        
        % adhesion properties 
        for j=1:max_number_of_adhesions
            % deal with existing adhesion
            if is_adhesion_exists(nr,j)==1
                l = location_of_adhesion(nr,j); % node of the current cell 
                nr2 = the_other_cell_number(nr,j); 
                n = the_other_cell_node(nr,j);
                
                tip_x_adhesion = x_cell_old(n,nr2);   
                tip_y_adhesion = y_cell_old(n,nr2);
                
                current_length = sqrt((tip_x_adhesion-x_cell_old(l,nr))^2+ (tip_y_adhesion-y_cell_old(l,nr))^2);
                
                % if adhesion is old, then kill it
                if (nt*dt - birth_time_of_adhesion(nr,j))>adhesion_regular_time
                    is_adhesion_exists(nr,j)=0;
                    number_of_outcoming_adhesions(nr,l)=number_of_outcoming_adhesions(nr,l)-1;
                end
                
                % if adhesion is too long, then kill it
                if is_adhesion_exists(nr,j)==1 && current_length>adhesion_max_length
                    is_adhesion_exists(nr,j)=0;
                    number_of_outcoming_adhesions(nr,l)=number_of_outcoming_adhesions(nr,l)-1;
                end
            end
            
            
            % deal with non-existing adhesion
            if is_adhesion_exists(nr,j)==0 && rand()<adhesion_birth_rate*dt
                
                % assign new adhesion to a node
                % adhesion can only appear at ends 
                randloc = rand();
                if randloc<1/2 % adhesion emerges at front with probability 1/3
                    l = 1;
                else
                    %if randloc>=1/2 % adhesion emerges at front with probability 2/3
                    l=number_of_nodes;
                    %else 
                    %    l=randi([2,number_of_nodes-1]);
                    %end
                end
                location_of_adhesion(nr,j)=l;
                
                % search of other cell springs nearby
                for nr2 = 1:number_of_cells 
                    if (nr2~=nr)
                        % change this to nucleus to nucleus distance ???
                        % distC = sqrt((x_cell_old(l,nr)-x_cell_old(i_nucleus,nr2))^2+(y_cell_old(l,nr)-y_cell_old(i_nucleus,nr2))^2);
                        distC = sqrt((x_cell_old(i_nucleus,nr)-x_cell_old(i_nucleus,nr2))^2+(y_cell_old(i_nucleus,nr)-y_cell_old(i_nucleus,nr2))^2);
                        if distC<5.0*double(i_nucleus)*L_0_eq % NOT TO CONSIDER CELLS THAT ARE TOO FAR
                            %firstcheck = 1
                            for k=1:number_of_nodes 
                                Ax = x_cell_old(k,nr2); % the other cell node is A
                                Ay = y_cell_old(k,nr2);
                                Cx = x_cell_old(l,nr); % the current cell node is C
                                Cy = y_cell_old(l,nr);
                                
                                distM = sqrt((Cx-Ax)^2+(Cy-Ay)^2);

                                % We also need to distinguish between edge and
                                % internal nodes. Edge node prefer to adhere
                                % orthogonally (removed)
                                
                                orientation_check = 1;

                                if (l==1)||(l==number_of_nodes)
                                    l_adjacent = 2; 
                                    if (l==number_of_nodes)
                                        l_adjacent = number_of_nodes-1; 
                                    end
                                    % Ex = x_cell_old(l_adjacent,nr);
                                    % Ey = y_cell_old(l_adjacent,nr);
                                    % 
                                    % dot_EC_CA = (Cx-Ex)*(Ax-Cx)+ (Cy-Ey)*(Ay-Cy);
                                    % norm_EC = sqrt((Cx-Ex)^2 + (Cy-Ey)^2);
                                    % norm_CA = sqrt((Ax-Cx)^2 + (Ay-Cy)^2);
                                    % cos_EC_CA = dot_EC_CA / (norm_EC * norm_CA);

                                    % if (cos_EC_CA)^2 > (cos(3*pi/8))^2  % why????
                                    %     orientation_check =0;
                                    % end

                                end

                                if (distM < adhesion_max_length)&&(orientation_check)&&(number_of_outcoming_adhesions(nr,l)<2) % why 2???
                                    % HERE IS AN ESSENTIAL FLAW: IF WE HAVE
                                    % MULTIPLE CANDIDATES TO ADHERE, WE CHOOSE
                                    % THE ONE WITH THE HIGHER LABEL - NEED TO
                                    % CORRECT IN FUTURE
                                    is_adhesion_exists(nr,j)=1;
                                    birth_time_of_adhesion(nr,j)=nt*dt;
                                    location_of_adhesion(nr,j)=l;
                                    the_other_cell_number(nr,j)=nr2;
                                    the_other_cell_node(nr,j)=k; 
                                    number_of_outcoming_adhesions(nr,l)=number_of_outcoming_adhesions(nr,l)+1;
                                end
                            end % k
                        end %if disct<i_nucleus*L_0_eq
                    end %  if (nr2~=nr)
                end % nr2
            end % if is_adhesion_exists(nr,j)==0 && rand()<adhesion_birth_rate*dt
       end %j
        
        % 
        for i=1:number_of_nodes 
            force_adhesion_x(i)=0;
            force_adhesion_y(i)=0;
            
            %compute force coming from outcoming adhesion
            for j=1:max_number_of_adhesions
                if is_adhesion_exists(nr,j)&&(location_of_adhesion(nr,j)==i)
                    nr2 = the_other_cell_number(nr,j); 
                    n = the_other_cell_node(nr,j);
                                   
                    tip_x_adhesion = x_cell_old(n,nr2);   
                    tip_y_adhesion = y_cell_old(n,nr2);
                
                    current_length = max(0.0000001,sqrt((tip_x_adhesion-x_cell_old(i,nr))^2+ (tip_y_adhesion-y_cell_old(i,nr))^2));
                    nx = (tip_x_adhesion-x_cell_old(i,nr))/current_length;
                    ny = (tip_y_adhesion-y_cell_old(i,nr))/current_length;
                    
                    force_adhesion_x(i) = force_adhesion_x(i) + k_adhesion*(current_length-adhesion_equilibrium_length)*nx;
                    force_adhesion_y(i) = force_adhesion_y(i) + k_adhesion*(current_length-adhesion_equilibrium_length)*ny;
                end
            end
            
            %compute force coming from incoming adhesion
            % meaning of coming and incoming adhesions???
            for nr2=1:number_of_cells
                if nr2~=nr
                    for j=1:max_number_of_adhesions
                        if is_adhesion_exists(nr2,j)&&(the_other_cell_number(nr2,j)==nr)&&(the_other_cell_node(nr2,j)==i)
                            n = location_of_adhesion(nr2,j);
                            
                            tip_x_adhesion = x_cell_old(n,nr2);   
                            tip_y_adhesion = y_cell_old(n,nr2);
                
                            current_length = max(0.0000001,sqrt((tip_x_adhesion-x_cell_old(i,nr))^2+ (tip_y_adhesion-y_cell_old(i,nr))^2));
                            nx = (tip_x_adhesion-x_cell_old(i,nr))/current_length;
                            ny = (tip_y_adhesion-y_cell_old(i,nr))/current_length;
                    
                            force_adhesion_x(i) = force_adhesion_x(i) + k_adhesion*(current_length-adhesion_equilibrium_length)*nx;
                            force_adhesion_y(i) = force_adhesion_y(i) + k_adhesion*(current_length-adhesion_equilibrium_length)*ny;
                        end
                    end
                end
            end
            
        end
        
        % compute bounding force 
        % bound back from boundary ???
        for i=1:number_of_nodes 
            how_far_from_center_of_box = sqrt((x_cell_old(i,nr)-L_box/2.0)^2+(y_cell_old(i,nr)-L_box/2.0)^2);
            force_bound_x(i) =0.0;
            force_bound_y(i) =0.0;

            if how_far_from_center_of_box > R_bound
                  nx = - (x_cell_old(i,nr)-L_box/2.0)/how_far_from_center_of_box;
                  ny = - (y_cell_old(i,nr)-L_box/2.0)/how_far_from_center_of_box;
                  force_bound_x(i) = k_spring *(how_far_from_center_of_box-R_bound)/R_bound *nx;
                  force_bound_y(i) = k_spring *(how_far_from_center_of_box-R_bound)/R_bound *ny;
            end

        end
        
        %update how_circular_is_the_cell
        if is_the_cell_circular(nr)==1
            how_circular_is_the_cell(nr)=1;
        else 
            how_circular_is_the_cell(nr)=max(0.0,how_circular_is_the_cell(nr)-0.1*dt);
        end
        
        %update x_cell and y_cell
        for i=1:number_of_nodes
            x_cell(i,nr)=x_cell_old(i,nr) + min(0.05,(1.0-how_circular_is_the_cell(nr)) * dt * (f_spring_x(i)+f_bending_x(i)+f_repulsion_x(i)+force_filopodium_x(i)+force_adhesion_x(i)+force_bound_x(i)));
            y_cell(i,nr)=y_cell_old(i,nr) + min(0.05,(1.0-how_circular_is_the_cell(nr)) * dt * (f_spring_y(i)+f_bending_y(i)+f_repulsion_y(i)+force_filopodium_y(i)+force_adhesion_y(i)+force_bound_y(i)));
        end
         
         % cell length 
        cellLength(nr)=0;
        for i=1:number_of_nodes-1
            cellLength(nr)=cellLength(nr)+sqrt((x_cell(i+1,nr)-x_cell(i,nr))^2+(y_cell(i+1,nr)-y_cell(i,nr))^2);
        end
        % sometimes displacement of nodes are very huge, to prevent it I
        % use the condition below
        if cellLength(nr) > 1.4 * number_of_nodes * L_0_eq
            for i=1:number_of_nodes
                x_cell(i,nr)=x_cell_old(i,nr)+ min(0.1,(1.0-how_circular_is_the_cell(nr)) * dt *f_spring_x(i))  ; % why add it again???
                y_cell(i,nr)=y_cell_old(i,nr)+ min(0.1,(1.0-how_circular_is_the_cell(nr)) * dt *f_spring_y(i))  ;
            end
        end
         
        
        %if cell is too short it is circular
        length_of_the_current_circle = 0.0;
        if is_the_cell_circular(nr)==0
            for i=1:number_of_nodes-1
                length_of_the_current_circle = length_of_the_current_circle + sqrt((x_cell(i+1,nr)-x_cell(i,nr))^2+(y_cell(i+1,nr)-y_cell(i,nr))^2);
            end
            if length_of_the_current_circle < L_0_eq / 3.0
                is_the_cell_circular(nr) = 0; % change to circular
            else 
                is_the_cell_circular(nr) = 0;
            end
        end
            
        %if cell is not very active it may become circular
        if is_the_cell_circular(nr)==0 
            if (sum(activation_of_filopodium(nr,:)) == 0)&&(rand()<rate_of_circularity)
                is_the_cell_circular(nr)=0; % change to circular 
                phi_cell2 = atan2(y_cell(number_of_nodes,nr)-y_cell(i_nucleus,nr),x_cell(number_of_nodes,nr)-x_cell(i_nucleus,nr));
                x_center = x_cell(i_nucleus,nr);
                y_center = y_cell(i_nucleus,nr);
                for i=1:number_of_nodes
                    x_cell(i,nr)=x_center + (double(i)-double(i_nucleus))*L_0_eq/10.0*cos(phi_cell2); % why this treatment?
                    y_cell(i,nr)=y_center + (double(i)-double(i_nucleus))*L_0_eq/10.0*sin(phi_cell2);
                end
            end
        end
        
        % cell isolation condition???
        %cell stops to be not circular if it feels the other cell
        % if is_the_cell_circular(nr)==1
        %     total_repulsion = 0; 
        %     for i=1:number_of_nodes 
        %         total_repulsion = total_repulsion + sqrt(f_repulsion_x(i)^2+f_repulsion_y(i)^2);
        %     end
        %     if total_repulsion > 0 
        %         is_the_cell_circular(nr) = 0;
        %     end
        % end
            
        
        % update x_cell and _y_cell if cell is circular
        random_displacement_x = 2*rand()-1;
        random_displacement_y = 2*rand()-1;
        if is_the_cell_circular(nr)==1
           for i=1:number_of_nodes
                x_cell(i,nr)=x_cell(i,nr) + sqrt(2*D_circ*dt)*random_displacement_x;
                y_cell(i,nr)=y_cell(i,nr) + sqrt(2*D_circ*dt)*random_displacement_y;
           end
        end
        
        % update total life of evs absorted by cells 
        ev_totallife(nr) = max(0,ev_totallife(nr)-ev_num(nr)*dt);
        
        % tempDistVec = sqrt((x_cell(i_nucleus,nr)-x_cell(i_nucleus,:)).^2 + (y_cell(i_nucleus,nr)-y_cell(i_nucleus,:)).^2);
        neighNum = 0;
        for i=1:number_of_cells
            if i~=nr
                tempDistVec = pdist2([x_cell(:,nr),y_cell(:,nr)],[x_cell(:,i),y_cell(:,i)]);
                if min(min(tempDistVec))<0.8
                    neighNum = neighNum + 1;
                end
            end
        end
        % neighNum = sum(tempDistVec<0.8);
        %cell division
        ttlife = ev_totallife(nr);
        if (rand()<(division_rate+ttlife/(ttlife+10))*dt)&&(nt*dt-cell_birth_time(nr)>division_age) && neighNum<3
            number_of_new_cells = number_of_new_cells + 1;
            
            for i=1:number_of_nodes 
                x_helper(i)=x_cell(i,nr);
                y_helper(i)=y_cell(i,nr);
            end
            
           
            
            % creation of new cell (first) at slot i  !!! I ASSUME THAT THE NUMBER
            % OF NODES IS ODD!!!
            
            for i=i_nucleus:(number_of_nodes-1)
                x_cell(2*(i-i_nucleus)+1,nr)=x_helper(i); % why this assignment ??? does not it stay the same???
                y_cell(2*(i-i_nucleus)+1,nr)=y_helper(i);
                x_cell(2*(i-i_nucleus)+2,nr)=0.5*(x_helper(i+1)+x_helper(i));
                y_cell(2*(i-i_nucleus)+2,nr)=0.5*(y_helper(i+1)+y_helper(i));
            end
            x_cell(number_of_nodes,nr)=x_helper(number_of_nodes);
            y_cell(number_of_nodes,nr)=y_helper(number_of_nodes);
            
            % properties relating to the division
            cell_birth_time(nr) = nt*dt;
            how_circular_is_the_cell(nr) = 0;
            is_the_cell_circular(nr) = is_the_cell_circular(nr);
    
            % on overlapping 
            for i=1:number_of_nodes
                is_the_node_overlapping(nr,i)=0;
                overlapping_other_cell(nr,i)=0;
            end
            

    
            % filopodia properties
            for j=1:max_number_of_filopodia
                is_filopodium_exists(nr,j)=0;
                birth_time_of_filopodium(nr,j)=0;
                location_of_filopodium(nr,j)=0; % node's number
                tip_x_filopodium(nr,j)=0; 
                tip_y_filopodium(nr,j)=0;
                activation_of_filopodium(nr,j)=0;
            end
            
            % adhesion properties 
            for j=1:max_number_of_adhesions
                is_adhesion_exists(nr,j)=0;
                birth_time_of_adhesion(nr,j)=0;
                location_of_adhesion(nr,j)=0;
                the_other_cell_number(nr,j)=0;
                the_other_cell_node(nr,j)=0;
            end
            
            for i=1:number_of_nodes
                number_of_outcoming_adhesions(nr,i)=0;
            end


            % cell length 
            cellLength(nr)=0;
            for i=1:number_of_nodes-1
                cellLength(nr)=cellLength(nr)+sqrt((x_cell(i+1,nr)-x_cell(i,nr))^2+(y_cell(i+1,nr)-y_cell(i,nr))^2);
            end
            
            equilibrium_length(nr) = equilibrium_length(nr)/2;
            
            ev_totallife(nr) = 0;
            ev_num(nr) = 0;
            % creation of new cell (second) at new slot
            nr2 = number_of_cells + number_of_new_cells;
            who_is_my_brother(nr)=nr2;
            time_i_was_born(nr)=nt*dt;
            for i=1:i_nucleus-1
                x_cell(2*(i-1)+1,nr2)=x_helper(i);
                y_cell(2*(i-1)+1,nr2)=y_helper(i);
                x_cell(2*(i-1)+2,nr2)=0.5*(x_helper(i+1)+x_helper(i));
                y_cell(2*(i-1)+2,nr2)=0.5*(y_helper(i+1)+y_helper(i));
            end
            x_cell(number_of_nodes,nr2)=x_helper(i_nucleus);
            y_cell(number_of_nodes,nr2)=y_helper(i_nucleus);
            
            % properties relating to the division
            cell_birth_time(nr2) = nt*dt;
            how_circular_is_the_cell(nr2) = 0;
            is_the_cell_circular(nr2) = is_the_cell_circular(nr);
    
            % on overlapping 
            for i=1:number_of_nodes
                is_the_node_overlapping(nr2,i)=0;
                overlapping_other_cell(nr2,i)=0;
            end
            
            who_is_my_brother(nr2)=nr;
            time_i_was_born(nr2)=nt*dt;
    
            % filopodia properties
            for j=1:max_number_of_filopodia
                is_filopodium_exists(nr2,j)=0;
                birth_time_of_filopodium(nr2,j)=0;
                location_of_filopodium(nr2,j)=0; % node's number
                tip_x_filopodium(nr2,j)=0; 
                tip_y_filopodium(nr2,j)=0;
                activation_of_filopodium(nr2,j)=0;
            end
            
            % adhesion properties 
            for j=1:max_number_of_adhesions
                is_adhesion_exists(nr2,j)=0;
                birth_time_of_adhesion(nr2,j)=0;
                location_of_adhesion(nr2,j)=0;
                the_other_cell_number(nr2,j)=0;
                the_other_cell_node(nr2,j)=0;
            end
            
            for i=1:number_of_nodes
                number_of_outcoming_adhesions(nr2,i)=0;
            end

            % cell length 
            cellLength(nr2)=0;
            for i=1:number_of_nodes-1
                cellLength(nr2)=cellLength(nr2)+sqrt((x_cell(i+1,nr2)-x_cell(i,nr2))^2+(y_cell(i+1,nr2)-y_cell(i,nr2))^2);
            end
            if is_the_cell_circular(nr2)==0  
                equilibrium_length(nr2) = equilibrium_length(nr);
            else 
                equilibrium_length(nr2) = L_0_eq/5;
            end

            % 
            ev_totallife(nr2) = 0;
            ev_num(nr2) = 0;
            
        end % if (rand()<division_rate*dt)&&(nt*dt-cell_birth_time(nr)>division_age)
        
        
    end
    
    number_of_cells = number_of_cells + number_of_new_cells;
    
    for nv=1:number_of_vesicles
        xv(nv) = xv(nv) + sqrt (2*D_v*dt)*(2*rand()-1);
        yv(nv) = yv(nv) + sqrt (2*D_v*dt)*(2*rand()-1);
        agev(nv) = agev(nv) + dt;
    end
    
   % Production of new vecicles
   number_of_new_vesicles = 0;
   for nr=1:number_of_cells
       if rand()<vesicle_production_rate*dt
           l=randi([1,number_of_nodes-1]);
           lt = rand(); % a value between 0 and 1
           have_i_recorded_the_vesicle = 0;
           ves_number = 0;
           vesind = find(agev>vesicle_life_expectancy,1);
           if ~isempty(vesind)
               xv(vesind) = x_cell(l,nr) + lt*(x_cell(l+1,nr)-x_cell(l,nr));
               yv(vesind) = y_cell(l,nr) + lt*(y_cell(l+1,nr)-y_cell(l,nr));
               vn(vesind) = nr;
               agev(vesind) = 0;
           else
               number_of_new_vesicles = number_of_new_vesicles + 1;
               xv(nv+number_of_new_vesicles) = x_cell(l,nr) + lt*(x_cell(l+1,nr)-x_cell(l,nr));
               yv(nv+number_of_new_vesicles) = y_cell(l,nr) + lt*(y_cell(l+1,nr)-y_cell(l,nr));
               vn(nv+number_of_new_vesicles)=nr;
               agev(nv+number_of_new_vesicles)=0;
           end
           % while (have_i_recorded_the_vesicle == 0)&&(ves_number<nv)
           %      ves_number = ves_number + 1;
           %      if (agev(ves_number)>vesicle_life_expectancy)
           %          xv(ves_number) = x_cell(l,nr) + 0.01*(2*rand()-1);
           %          yv(ves_number) = y_cell(l,nr) + 0.01*(2*rand()-1);
           %          vn(ves_number) = nr;
           %          agev(ves_number) = 0;
           %          have_i_recorded_the_vesicle = 1;
           %      end
           % end
           % if have_i_recorded_the_vesicle == 0
           %     number_of_new_vesicles = number_of_new_vesicles + 1;
           %     xv(nv+number_of_new_vesicles) = x_cell(l,nr) + 0.01*(2*rand()-1);
           %     yv(nv+number_of_new_vesicles) = y_cell(l,nr) + 0.01*(2*rand()-1);
           %     vn(nv+number_of_new_vesicles)=nr;
           %     agev(nv+number_of_new_vesicles)=0;
           % end
       end
   end
   number_of_vesicles = number_of_vesicles + number_of_new_vesicles;    
   
   
    
   % Drawing cells
   
   % f1=figure(2);
   % %scrsz = get(groot,'ScreenSize'); maxscrsz=min(scrsz(3),scrsz(4));
   % %set(f1,'Position',[scrsz(3)/3 0 maxscrsz maxscrsz],'Color','w') 
   % set(gcf,'Position',[100 100 1300 500])
   % if (mod(nt,frame_frequency_0)==1)
   %      clf;
   %      ik=int16((nt+1)/frame_frequency_0);
   % 
   %      subplot(1,3,1)
   % 
   %      % drawing cells,
   %     for nr=1:number_of_cells
   %          if is_the_cell_circular(nr)==0
   %              plot(x_cell(:,nr),y_cell(:,nr),'color',[0 0.5 0],'LineWidth',2); hold on;         
   %          else
   %              plot(x_cell(i_nucleus,nr),y_cell(i_nucleus,nr),'green.','MarkerSize',20); hold on;   
   %          end
   %     end
   % 
   % 
   % 
   % 
   % 
   %      ch=sprintf("Time %d, Cell (%d)",ik,number_of_cells);
   %      title(ch);
   % 
   %      grid off
   %      daspect([1 1 1]);
   %      axis([0 L_box 0 L_box]);
   % 
   % 
   %      subplot(1,3,2)
   % 
   %      % drawing vesicles
   %      % for nv=1:number_of_vesicles
   %      %      if (agev(nv)<vesicle_life_expectancy)
   %      %         if vn(nv)==1
   %      %             plot(xv(nv),yv(nv),'blue.'); hold on;
   %      %         else if vn(nv)==2
   %      %                plot(xv(nv),yv(nv),'blue.'); hold on;
   %      %             else 
   %      %                 plot(xv(nv),yv(nv),'blue.'); hold on;
   %      %             end
   %      %         end
   %      %      end
   %      %  end
   % 
   % 
   % 
   % 
   %      % drawing cells, their filopodia and adhesions
   %      for nr=1:number_of_cells
   %          plot(x_cell(:,nr),y_cell(:,nr),'green','LineWidth',2); hold on;
   %          plot(x_cell(:,nr),y_cell(:,nr),'black.','MarkerSize',5); hold on;
   %          if is_the_cell_circular(nr)==1
   %              plot(x_cell(i_nucleus,nr),y_cell(i_nucleus,nr),'black.','MarkerSize',20); hold on; 
   %          end
   %          % drawing filopodia
   %          for j=1:max_number_of_filopodia
   %              if is_filopodium_exists(nr,j)==1
   %                  l=location_of_filopodium(nr,j);
   %                  if activation_of_filopodium(nr,j)>0
   %                      plot([x_cell(l,nr) tip_x_filopodium(nr,j)],[y_cell(l,nr) tip_y_filopodium(nr,j)],'red'); hold on;
   %                  else
   %                      plot([x_cell(l,nr) tip_x_filopodium(nr,j)],[y_cell(l,nr) tip_y_filopodium(nr,j)],'black'); hold on;
   %                  end
   %              end
   %          end            
   %      end
   % 
   % 
   % 
   %      ch=sprintf("Filopodia (%d) & Vesicles (%d)",sum(sum(is_filopodium_exists)),number_of_vesicles);
   %      title(ch);
   % 
   %      grid off
   %      daspect([1 1 1]);
   %      axis([0 L_box 0 L_box]);
   % 
   % 
   %      subplot(1,3,3)
   % 
   % 
   %      % drawing cells, their adhesions
   %      for nr=1:number_of_cells
   %          plot(x_cell(:,nr),y_cell(:,nr),'green','LineWidth',2); hold on;
   %          plot(x_cell(:,nr),y_cell(:,nr),'black.','MarkerSize',5); hold on;
   % 
   %          if is_the_cell_circular(nr)==1
   %              plot(x_cell(i_nucleus,nr),y_cell(i_nucleus,nr),'black.','MarkerSize',20); hold on; 
   %          end
   % 
   %          % drawing adhesions
   %          for j=1:max_number_of_adhesions
   %              if is_adhesion_exists(nr,j)==1
   %                  l=location_of_adhesion(nr,j);
   %                  nr2 = the_other_cell_number(nr,j); 
   %                  n = the_other_cell_node(nr,j);
   % 
   %                  tip_x_adhesion = x_cell(n,nr2);   
   %                  tip_y_adhesion = y_cell(n,nr2);
   % 
   %                  distA = sqrt((tip_x_adhesion-x_cell(l,nr))^2+(tip_y_adhesion-y_cell(l,nr))^2);
   % 
   %                  if distA<adhesion_max_length
   %                      plot([x_cell(l,nr) tip_x_adhesion],[y_cell(l,nr) tip_y_adhesion],'red'); hold on;
   %                  end
   % 
   %              end
   %          end
   % 
   %      end
   % 
   %      %drawing overlaps 
   %      if overlap_positions>0
   %          plot(x_overlap([1:overlap_positions]),y_overlap([1:overlap_positions]),'blue*','MarkerSize',12); hold on;
   %      end
   % 
   %      ch=sprintf("Adhesions (%d), Overlaps (%d)",sum(sum(is_adhesion_exists)),overlap_positions);
   %      title(ch);
   % 
   % 
   %      grid off
   %      daspect([1 1 1]);
   %      axis([0 L_box 0 L_box]);
   % 
   % 
   % 
   % 
   %      ch=sprintf("figs/%d.png",ik);
   %      saveas(gcf,ch);
   % end
   if (mod(nt,frame_frequency) == 0) 
        save(['./data/testStem1_t',num2str(nt),'.mat']);
    end
end

% drawing entire tissue
% figure(3)
% clf;
% % drawing cells,
% for nr=1:number_of_cells
%     if is_the_cell_circular(nr)==0
%         plot(x_cell(:,nr),y_cell(:,nr),'black','LineWidth',1); hold on;         
%     else
%         plot(x_cell(i_nucleus,nr),y_cell(i_nucleus,nr),'black.','MarkerSize',20); hold on;   
%     end    
% end

            
        
       
        
% ch=sprintf("Time %d, Cell (%d)",ik,number_of_cells);
% title(ch);
% 
% grid off
% daspect([1 1 1]);

% exit;