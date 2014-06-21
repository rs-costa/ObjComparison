function generateResults(glob,solution,tit,constit,fname,task)
    %% bar plots
    f=figure('visible','off');
    if task==1 || task==2
        noOfCons=length(solution);
        noOfObj=length(solution{1});
        for j=1:noOfObj
            subplot(2,noOfObj/2,j);
            for i=1:noOfCons
                ylim([.5 noOfCons+.5])
                xlim([0 2])
                if ~isfield(solution{i}{j},'minE')
                continue;
                end
                if abs(solution{i}{j}.minE - solution{i}{j}.maxE) > .001
                    barh(i,solution{i}{j}.maxE,'BaseValue',solution{i}{j}.minE,'BarWidth',0.4); 
                else
                    plot(solution{i}{j}.minE, i,'-or')            
                end
                hold on;
                line([0 solution{i}{j}.minE],[i i],'LineStyle',':');
                hold on;
            end
            barh(subplot(2,noOfObj/2,j),0,0,'BaseValue',0);

            title(tit{j},'FontSize', 10,'FontWeight','b','Color','b');
            hold off;
            if(j==1 || (j-noOfObj/2 >= 1 && j-noOfObj/2 <= 1.5))
                ytick = 1:noOfCons;
                set(gca,'ytick',ytick);
                yticklabel=constit;
                set(gca,'yticklabel',yticklabel)
            else
                set(gca,'ytick',[]);
            end
        end
    elseif task==4
        f=figure('visible','off');
        for i=1:length(solution)
            silico_flux = solution{i}.x(glob.indices);
            silico_flux=silico_flux/100;
            vivo_flux=glob.exp_flux/100;
            l=min([vivo_flux;silico_flux]);
            u=max([vivo_flux;silico_flux]);

            correlation=corr(vivo_flux, silico_flux);
            correlation = floor(correlation * 100) / 100;
            subplot(3,3,i)
            plot(vivo_flux, silico_flux, 'o'); 
            title({tit{i};strcat('R=',num2str(correlation))});
            hold on;
            xlim([l u])
            ylim([l u])
            plot([l u], [l u], 'k:'); 
            hold off;
        end    
        
    else
%         shown=20;
        shown=length(solution);
        for j=1:shown

            xlim([.5 shown+.5])
            ylim([0 1])

            subplot(2,1,1)
            plot(j,solution(j),'-*r');grid on;
            hold on;
            line([j j],[0 solution(j)],'LineStyle',':');
            hold on;
        end
        title('Euclidean distances to in vivo fluxes','FontSize', 10,'FontWeight','b','Color','b');
        hold off;
        xtick=1:shown;
        set(gca,'xtick',xtick);
        xticklabel=tit(1:shown);
        set(gca,'xticklabel',xticklabel)
        h = get(gca,'xlabel');
        xlabelstring = get(h,'string');
        xlabelposition = get(h,'position');
        % ## construct position of new xtick labels
        yposition = xlabelposition(2);
        yposition = repmat(yposition,length(xtick),1);
        % ## disable current xtick labels
        set(gca,'xtick',[]);
        % ## set up new xtick labels and rotate
        hnew = text(xtick, yposition, xticklabel);
        set(hnew,'rotation',90,'horizontalalignment','right');

    end

    saveas(f,strcat(fname,'.fig'));
    save(strcat(fname,'.mat'),'solution','tit');
end
