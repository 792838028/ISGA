clc;clear;close all
box_pp = 1; 
N=50; 
lb=-100; 
ub=100;
for dim =  10 
    T=10000*dim;
for F = 1:12 

 fobj = @(x) cec22_func(x',F);

    for j = 1:30
[Top_Score_1(j,:),Top_Position_1,Convergence_curve_1(j,:)]=ISGA(N,T,lb,ub,dim,fobj);
    end 
box_plot = [Top_Score_1']; 

    figure;
    if box_pp == 1 
  mycolor =[0, 0.75, 0.75;...  
            0, 0, 0;...        
            0.75, 0.75, 0;...  
            0.75, 0, 0.75;...  
            0, 0, 1;...        
            0, 0.5, 0;...      
            0.93, 0.69, 0.13;...
            1, 0.5, 0;...      
            1, 0, 0];    

        box_figure = boxplot(box_plot','color',[0 0 0],'Symbol','o');
        set(box_figure,'Linewidth',1.2);
        boxobj = findobj(gca,'Tag','Box');
        for op = 1
            patch(get(boxobj(op),'XData'),get(boxobj(op),'YData'),mycolor(10-op,:),'FaceAlpha',0.5,...
                'LineWidth',0.7);
        end
        set(gca,'XTickLabel',{'ISGA'});
        title(['CEC2022-F',num2str(F),' (Dim=', num2str(dim),')'])
        hold on
    end 

colors = {
    [0, 0.75, 0.75],  
    [0, 0, 0],        
    [0.75, 0.75, 0],       
    [0.75, 0, 0.75],  
    [0, 0, 1],       
    [0, 0.5, 0],      
    [0.93, 0.69, 0.13],
    [1, 0.5, 0],        
    [1, 0, 0],    
    [0, 1, 0],        
    [1, 0.5, 0],      
    [0.75, 0.75, 0]   
};

figure;
semilogy(mean(Convergence_curve_1), '-', 'Color', colors{9}, 'LineWidth', 1.5);
hold on;
legend('ISGA','Interpreter','none');
title(['CEC2022-F',num2str(F),' (Dim=', num2str(dim),')'])
xlabel('FEs#');
ylabel('Mean Fitness Value');
axis tight;
grid on;
box on;
hold off;


disp('-------------------------------------');
display(['CEC2022-F',num2str(F),' Dim: ', num2str(dim),' Np: ', num2str(N), ' FEs: ', num2str(T)])
disp(['Ave: ', num2str(mean(Top_Score_1)), '，Std: ', num2str(std(Top_Score_1)), '，Best: ', num2str(min(Top_Score_1))]);
end
end
