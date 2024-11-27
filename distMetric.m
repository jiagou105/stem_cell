function flag = distMetric(xcell1,ycell1,xcell2,ycell2)

flag = false;
minDist = 100;
for indi = 1:length(xcell1)
    for indj = 1:length(xcell2)
        dist = sqrt((xcell1(indi)-xcell2(indj))^2+(ycell1(indi)-ycell2(indj))^2);
        if dist<minDist
            minDist = dist;
        end
    end
end
        % distxy = [sqrt((xcell1(1)-xcell2(1))^2+(ycell1(1)-ycell2(1))^2);
        %   sqrt((xcell1(1)-xcell2(end))^2+(ycell1(1)-ycell2(end))^2);
        %   sqrt((xcell1(end)-xcell2(1))^2+(ycell1(end)-ycell2(1))^2);
        %   sqrt((xcell1(end)-xcell2(end))^2+(ycell1(end)-ycell2(end))^2)];
if minDist<0.2
    flag = true;
end



