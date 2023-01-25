function nth = plot_IRF(Y,row,col,names)

% Author: Alvin Lo Hei Chun
% INPUT: 
%   1) Y: The time series which has dimension T x 3 x N, N is then number
%   of variables
%   2) row: Number of rows for each subplot figure 
%   3) col: number of column for each subplot figure
% 
% There will be N/(row*col) figures 

[T,~,N] = size(Y);
F = N/(row*col);

for i= 1:F
    figure; 
    for j = 1:row*col
        subplot(row,col,j);
        hold on
        plot(1:T,Y(:,2,(i-1)*row*col + j),'k','LineWidth',2);
        plot(1:T,Y(:,1,(i-1)*row*col + j),'--k');
        plot(1:T,Y(:,3,(i-1)*row*col + j),'--k');
        line(xlim,[0,0]);
        if nargin == 3
            title(['Variable',int2str((i-1)*row*col+j)]);
        else
            title(names((i-1)*row*col+j));
        end
    end
end

end
    