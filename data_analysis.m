% analyze simulation data obtained from the stem cell model 
% 
clearvars;
timeinv = 200;
f1 = figure(1);
count = 1;
timeArray = timeinv:timeinv:40000;
timeStepNum = length(timeArray);
cellNumAll = zeros(1,timeStepNum);
evNumAll = zeros(1,timeStepNum);
cellNucDistMean = zeros(1,timeStepNum);
cellNucDistStd = zeros(1,timeStepNum);
cellAllPtDistMean = zeros(1,timeStepNum);
cellAllPtDistStd = zeros(1,timeStepNum);
orderQAllT = zeros(1,timeStepNum);
distRef = [0.5,1,3,6,10];
oriDist = zeros(5,timeStepNum);
rAll = 1:15;
cRAllT = zeros(timeStepNum,length(rAll));

cirRad = 0.5;
% avg distance between cell
% orientation between cells
for indi=1:timeStepNum
    nt = timeArray(indi);
    filename = ['./data/testStem1_t',num2str(nt),'.mat'];
    load(filename)
    cellNumAll(indi)=number_of_cells;
    evNumAll(indi) = number_of_vesicles;
    
    cellNucDistMinI = zeros(1,number_of_cells);
    cellAllPtDistMinI = ones(1,number_of_cells)*1e4;
    oriVecI = zeros(3,number_of_cells);
    orderQ = zeros(2,2);
    for indj = 1:number_of_cells
        oriVecI(1,indj) = x_cell(i_nucleus+1,indj)-x_cell(i_nucleus-1,indj);
        oriVecI(2,indj) = y_cell(i_nucleus+1,indj)-y_cell(i_nucleus-1,indj);
        oriVecI(3,indj) = atan2(oriVecI(2,indj),oriVecI(1,indj));
    end
    % temp 
    tempDistMatrix = pdist2([x_cell(i_nucleus,:).',y_cell(i_nucleus,:).'],[x_cell(i_nucleus,:).',y_cell(i_nucleus,:).']);
    for indk=1:length(rAll)
        r = rAll(indk);
        cR = 0;
        rNum = 0;
        for indm = 1:number_of_cells
            cellId = find(tempDistMatrix(indm,indm+1:end)<=r)+indm; % column ids
            if length(cellId)>=1
                phii = oriVecI(3,indm);
                phij = oriVecI(3,cellId); % vector 
                cR = cR + sum(cos(2*(phii-phij)));
                rNum = rNum+length(cellId);
            end
        end
        cRAllT(indi,indk) = cR/rNum;
    end
    for indj = 1:number_of_cells
        tempDistVec = sqrt((x_cell(i_nucleus,indj)-x_cell(i_nucleus,:)).^2 + (y_cell(i_nucleus,indj)-y_cell(i_nucleus,:)).^2);
        tempDistVecSort = sort(tempDistVec);
        cellNucDistMinI(indj) = tempDistVecSort(2); % 2nd smallest distance; the 1st one is the indj cell itself
        for indk=1:number_of_cells
            if indk~=indj
                tempPtDistVec = pdist2([x_cell(:,indj),y_cell(:,indj)],[x_cell(:,indk),y_cell(:,indk)]);
                cellAllPtDistMinI(indj) = min(cellAllPtDistMinI(indj),min(min(tempPtDistVec)));
            end
        end

        neighborInd = find(tempPtDistVec<cirRad); % need to exclude the cell itself 
        phii = oriVecI(3,indj);
        orderQ = orderQ + [cos(phii)^2-1/2,cos(phii)*sin(phii);cos(phii)*sin(phii), sin(phii)^2-1/2];
    end
    orderQ = orderQ/number_of_cells;
    lambda = 1/2*trace(orderQ)+sqrt(trace(orderQ)^2-4*det(orderQ));
    orderQAllT(indi) = lambda;
    cellNucDistMean(indi) = mean(cellNucDistMinI);
    cellNucDistStd(indi) = std(cellNucDistMinI);
    cellAllPtDistMean(indi) = mean(cellAllPtDistMinI);
    cellAllPtDistStd(indi) = std(cellAllPtDistMinI);
end


%%

for indi=1:timeStepNum
    nt = timeArray(indi);
    filename = ['./data/testStem1_t',num2str(nt),'.mat'];
    load(filename)
    groupNum = zeros(1,number_of_cells);
    numGroups = 0;
    for indj = 1:number_of_cells
        if groupNum(indj) == 0
            numGroups = numGroups + 1;
            groupNum(indj) = numGroups;
        end
        for indj1 = indj+1:number_of_cells
            if distMetric(x_cell(:,indj),y_cell(:,indj),x_cell(:,indj1),y_cell(:,indj1))==true
                groupNum(indj1) = groupNum(indj);
            end
        end
    end
    numEle = zeros(1,numGroups);
    for indj=1:numGroups
        numEle(indj) = sum(groupNum==indj);
    end
    maxGSize = max(numEle);
    numGk = zeros(1,maxGSize);
    for indj=1:maxGSize
        numGk(indj) = sum(numEle==indj);
    end
    G(indi) = 0;
    for k=1:maxGSize
        G(indi) = G(indi)+k^2*numGk(k);
    end
    G(indi) = G(indi)/number_of_cells^2;

   figure(1)
   clf;
   for nr=1:number_of_cells
        % if groupNum(nr) = 
            if is_the_cell_circular(nr)==0
                plot(x_cell(:,nr),y_cell(:,nr),'Color',[groupNum(nr) 0 0]/(10+groupNum(nr)),'LineWidth',2); hold on;         
            else
                plot(x_cell(i_nucleus,nr),y_cell(i_nucleus,nr),'MarkerSize',20); hold on;   
            end
        % end
   end

    ch=sprintf("Time %d, Cell (%d)",indi,number_of_cells);
    title(ch);

    grid off
    daspect([1 1 1]);
    axis([0 L_box 0 L_box]);



end
%%
figure(2)
plot(1:timeStepNum, G, "Color",'r');


%%
figure(3)
subplot(1,2,1)
plot(timeArray,cellNumAll,'LineWidth',1.5);
title('Cell Numbers')
xlabel('Time')
subplot(1,2,2)
plot(timeArray,evNumAll,'LineWidth',1.5);
title('EV Numbers')
xlabel('Time')
fn=sprintf("./figs/cell_ev_num.png");
saveas(gcf,fn);

%%
figure(4)
subplot(2,1,1)
plot(timeArray,cellNucDistMean,'LineWidth',1.5);
hold on;
plot(timeArray,cellNucDistMean+cellNucDistStd,'r--','LineWidth',1);
plot(timeArray,cellNucDistMean-cellNucDistStd,'r--','LineWidth',1)
hold off;
subplot(2,1,2)
plot(timeArray,cellAllPtDistMean,'LineWidth',1.5);
hold on;
plot(timeArray,cellAllPtDistMean+cellAllPtDistStd,'r--','LineWidth',1);
plot(timeArray,cellAllPtDistMean-cellAllPtDistStd,'r--','LineWidth',1)
hold off;
title('Min Distance Between Cellls')
xlabel('Time')
fn=sprintf("./figs/cell_min_dist.png");
saveas(gcf,fn);
% axis([0 max(timeArray) 0 0.01])

%% 
figure(5)
plot(timeArray,orderQAllT,'LineWidth',2);

%% 
figure(6)
% distance r increment = 100

for indi = 1:timeStepNum
    plot(rAll,cRAllT(indi,:),'LineWidth',2);
    axis([0 15 0 1])
    F(count) = getframe(gcf);
    drawnow;
    count = count+1;
    % ch=sprintf("figs/%d.png",ik);
    % saveas(gcf,ch);
end

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