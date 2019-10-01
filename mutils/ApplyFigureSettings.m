function ApplyFigureSettings( hfig )
%ApplyFigureSettings This function applies default figure settings for
%articles.
%   The function applies:
%      - Font Settings = Times New Roman
%      - Interpreter = latex

% Authors: Jan Kral <kral.j@lit.cz>
% Date: 20.1.2017

    set(findall(hfig,'-property','FontName'),'FontName','Times');
    set(findall(hfig,'-property','Interpreter'),'Interpreter','latex');
    set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
end

