% % % % % %  set(findall(gcf,'-property','FontSize'),'FontSize',20)
 set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex')
% % % % % % %  set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
grid minor

set(gcf,'Renderer','painters');
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',20,'interpreter','latex')
set(findall(figureHandle,'type','axes'),'Color','w','fontSize',20)
set(findall(figureHandle,'type','axes'),'TickLabelInterpreter','latex')
