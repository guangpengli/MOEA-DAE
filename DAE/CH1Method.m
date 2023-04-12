function Population=CH1Method(Population,GlobalN)
    cons=Population.cons;
    Isde=ISDE(Population);
    
    %对个体违约度进行初始化
    [SN,SM]=size(cons);
    szmax=max(cons,[],1);
    cons=(cons)./(repmat(szmax,SN,1));
    v=sum(max(0,cons),2)/SM;
    %将ISDE值与违约度v构成双目标
    PopObj=[Isde,v];
    [Next,~,~]=EnvironmentalSelection(PopObj,GlobalN);
    Population=Population(Next);
end


function Isde=ISDE(Population)
    %将产生的一个子代p和其父代种群Population合并，求其ISDE

    PopObj=Population.objs;
    %种群的大小
    [N,M]=size(PopObj);           %N为种群大小
    zmin=min(PopObj,[],1);
    zmax=max(PopObj,[],1);
    PopObj = (PopObj-repmat(zmin,N,1))./(repmat(zmax-zmin,N,1)); %规范化
    fpr=(mean(PopObj,2));  
    [~,rank]=sort(fpr);
    
    Isde=zeros(N,1);
    for j = 2 : N

        SFunctionValue = max(PopObj(rank(1:j-1),:),repmat(PopObj(rank(j),:),(j-1),1));
        
        Distance = inf(1,j-1);
            
        for i = 1 : (j-1)
            Distance(i) = norm(SFunctionValue(i,:)-PopObj(rank(j),:))/M;
        end
           
        Distance = min(Distance);
        
        Isde(rank(j)) = exp(-Distance);
    end 
    
end