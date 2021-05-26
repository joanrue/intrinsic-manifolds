%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code by Selen Atasoy
%
% to make the figures fancy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h] = fancy_figure(h, opt)

% if nargin<3
%     font_type = 'Times New Roman';
% end;
% 
% if nargin<2
%     font_size = 30;
% end;

if nargin>0
    figure(h);
else
    h = figure;
end;

if nargin<2
    opt = [];
end;

%%
if ~isfield(opt, 'FontName')
    opt.FontName = 'Arial';
end;
if ~isfield(opt, 'FontSize')
    opt.FontSize = 15;
end;
if ~isfield(opt, 'ShowBox')
    opt.ShowBox = 'on';
end;
if ~isfield(opt, 'TickDir')
    opt.TickDir = 'in';
end;
if ~isfield(opt, 'TickLength')
    opt.TickLength = [0.02 0.02];
end;
% if ~isfield(opt, 'XTick')
%     opt.XTick = 'on';
% end;
% if ~isfield(opt, 'YTick')
%     opt.YTick = [];
% end;
% if ~isfield(opt, 'ZTick')
%     opt.ZTick = [];
% end
if ~isfield(opt, 'XMinorTick')
    opt.XMinorTick = 'on';
end;
if ~isfield(opt, 'YMinorTick')
    opt.YMinorTick = 'on';
end;
if ~isfield(opt, 'ZMinorTick')
    opt.ZMinorTick = 'on';
end;
if ~isfield(opt, 'XMinorGrid')
    opt.XMinorGrid = 'off';
end;
if ~isfield(opt, 'YMinorGrid')
    opt.YMinorGrid = 'off';
end;
if ~isfield(opt, 'ZMinorGrid')
    opt.ZMinorGrid = 'off';
end;
if ~isfield(opt, 'AxisColor')
    opt.AxisColor = [0 0 0];
end;
if ~isfield(opt, 'AxisLineWidth')
    opt.AxisLineWidth = 1;
end;

if ~isfield(opt, 'LegendBox')
    opt.LegendBox = 'on';
end;
        %%

ch=get(gcf,'children');
get(get(gcf,'children'),'type');

for i=1:length(ch)
        
    isAxes = strcmp(get(get(gcf,'children'),'type'), 'axes');
    isColorBar = strcmp(get(ch(i),'Tag'), 'Colorbar');
    isLegend = strcmp(get(ch(i),'Tag'), 'legend');
    if isAxes(i) 
            set( ch(i),          ...
                'FontName'    , opt.FontName, ...
                'FontSize'    , opt.FontSize,...
                'FontUnits','points',... 
                'FontWeight','normal',... 
                'FontSize',9,...
                'Box'         , opt.ShowBox     , ...
                'Color'       , 'none',...
                'TickDir'     , opt.TickDir     , ...
                'TickLength'  , opt.TickLength , ...
                'XMinorTick'  , opt.XMinorTick, ...
                'YMinorTick'  , opt.YMinorTick, ...
                'ZMinorTick'  , opt.ZMinorTick, ...
                'XMinorGrid'  , opt.XMinorGrid, ...
                'YMinorGrid'  , opt.YMinorGrid, ...
                'ZMinorGrid'  , opt.ZMinorGrid, ...
                'XColor'      , opt.AxisColor, ...
                'YColor'      , opt.AxisColor, ...
                'ZColor'      , opt.AxisColor, ...
                'LineWidth'   , opt.AxisLineWidth);

            if isfield(opt, 'XTick')
                 set( ch(i), 'XTick', opt.XTick);
            end;
            if isfield(opt, 'YTick')
                 set( ch(i), 'YTick', opt.YTick);
            end;
            if isfield(opt, 'ZTick')
                 set( ch(i), 'ZTick', opt.ZTick);
            end;
        
        set(ch(i), 'FontName', opt.FontName, 'FontSize', opt.FontSize); 

        h_title=get(ch(i),'Title');
        set(h_title,  'FontName', opt.FontName, 'FontSize', opt.FontSize); 

        h_title=get(ch(i),'Xlabel');
        set(h_title,  'FontName', opt.FontName, 'FontSize', opt.FontSize); 

        h_title=get(ch(i),'Ylabel');
        set(h_title,  'FontName', opt.FontName, 'FontSize', opt.FontSize); 

        if isColorBar
            set(ch(i),  'FontName', opt.FontName, 'FontSize', opt.FontSize,  ...
                'YMinorTick','off', ...
                'Box', 'on'); 
        end;
        
        if strcmp(ch(i).XScale, 'log') || strcmp(ch(i).XScale, 'log')
            set(ch(i),'YMinorTick','on', 'XminorTick', 'on');
        end;
    else
        if isLegend
            set(ch(i),  'FontName', opt.FontName, 'FontSize', round(opt.FontSize/5*3), 'box', opt.LegendBox);
        end;
    end;
end;
