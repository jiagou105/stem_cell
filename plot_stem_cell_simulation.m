% load the output of the stem cell simulation results and plot
clearvars;
timeinv = 200;
f1 = figure(1);
count = 1;
for nt=timeinv:timeinv:30400
    filename = ['./data/testStem1_t',num2str(nt),'.mat'];
    load(filename)
    scrsz = get(groot,'ScreenSize'); maxscrsz=min(scrsz(3),scrsz(4));
   set(f1,'Position',[scrsz(3)/3 0 maxscrsz maxscrsz],'Color','w') 
   set(gcf,'Position',[100 100 1300 500])

    clf;
    ik=int16((nt+1)/frame_frequency_0);

    subplot(1,3,1)

    % drawing cells,
   for nr=1:number_of_cells
        if is_the_cell_circular(nr)==0
            plot(x_cell(:,nr),y_cell(:,nr),'color',[0 0.5 0],'LineWidth',2); hold on;         
        else
            plot(x_cell(i_nucleus,nr),y_cell(i_nucleus,nr),'green.','MarkerSize',20); hold on;   
        end
   end

    ch=sprintf("Time %d, Cell (%d)",ik,number_of_cells);
    title(ch);

    grid off
    daspect([1 1 1]);
    axis([0 L_box 0 L_box]);


    subplot(1,3,2)

    % drawing vesicles
    % for nv=1:number_of_vesicles
    %      if (agev(nv)<vesicle_life_expectancy)
    %         if vn(nv)==1
    %             plot(xv(nv),yv(nv),'blue.'); hold on;
    %         else if vn(nv)==2
    %                plot(xv(nv),yv(nv),'blue.'); hold on;
    %             else 
    %                 plot(xv(nv),yv(nv),'blue.'); hold on;
    %             end
    %         end
    %      end
    %  end




    % drawing cells, their filopodia and adhesions
    for nr=1:number_of_cells
        plot(x_cell(:,nr),y_cell(:,nr),'green','LineWidth',2); hold on;
        plot(x_cell(:,nr),y_cell(:,nr),'black.','MarkerSize',5); hold on;
        if is_the_cell_circular(nr)==1
            plot(x_cell(i_nucleus,nr),y_cell(i_nucleus,nr),'black.','MarkerSize',20); hold on; 
        end
        % drawing filopodia
        for j=1:max_number_of_filopodia
            if is_filopodium_exists(nr,j)==1
                l=location_of_filopodium(nr,j);
                if activation_of_filopodium(nr,j)>0
                    plot([x_cell(l,nr) tip_x_filopodium(nr,j)],[y_cell(l,nr) tip_y_filopodium(nr,j)],'red'); hold on;
                else
                    plot([x_cell(l,nr) tip_x_filopodium(nr,j)],[y_cell(l,nr) tip_y_filopodium(nr,j)],'black'); hold on;
                end
            end
        end            
    end



    ch=sprintf("Filopodia (%d) & Vesicles (%d)",sum(sum(is_filopodium_exists)),number_of_vesicles);
    title(ch);

    grid off
    daspect([1 1 1]);
    axis([0 L_box 0 L_box]);


    subplot(1,3,3)


    % drawing cells, their adhesions
    for nr=1:number_of_cells
        plot(x_cell(:,nr),y_cell(:,nr),'green','LineWidth',2); hold on;
        plot(x_cell(:,nr),y_cell(:,nr),'black.','MarkerSize',5); hold on;

        if is_the_cell_circular(nr)==1
            plot(x_cell(i_nucleus,nr),y_cell(i_nucleus,nr),'black.','MarkerSize',20); hold on; 
        end

        % drawing adhesions
        for j=1:max_number_of_adhesions
            if is_adhesion_exists(nr,j)==1
                l=location_of_adhesion(nr,j);
                nr2 = the_other_cell_number(nr,j); 
                n = the_other_cell_node(nr,j);

                tip_x_adhesion = x_cell(n,nr2);   
                tip_y_adhesion = y_cell(n,nr2);

                distA = sqrt((tip_x_adhesion-x_cell(l,nr))^2+(tip_y_adhesion-y_cell(l,nr))^2);

                if distA<adhesion_max_length
                    plot([x_cell(l,nr) tip_x_adhesion],[y_cell(l,nr) tip_y_adhesion],'red'); hold on;
                end

            end
        end

    end

    %drawing overlaps 
    if overlap_positions>0
        plot(x_overlap([1:overlap_positions]),y_overlap([1:overlap_positions]),'blue*','MarkerSize',12); hold on;
    end

    ch=sprintf("Adhesions (%d), Overlaps (%d)",sum(sum(is_adhesion_exists)),overlap_positions);
    title(ch);


    grid off
    daspect([1 1 1]);
    axis([0 L_box 0 L_box]);

    F(count) = getframe(gcf);
    drawnow;
    count = count+1;
    % ch=sprintf("figs/%d.png",ik);
    % saveas(gcf,ch);
end

%%
 writerObj = VideoWriter('myVideo11.avi');
 writerObj.FrameRate = 10;
 % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);