function Population=CH1Method(Population,GlobalN)
    cons=Population.cons;
    Isde=ISDE(Population);
    
    %�Ը���ΥԼ�Ƚ��г�ʼ��
    [SN,SM]=size(cons);
    szmax=max(cons,[],1);
    cons=(cons)./(repmat(szmax,SN,1));
    v=sum(max(0,cons),2)/SM;
    %��ISDEֵ��ΥԼ��v����˫Ŀ��
    PopObj=[Isde,v];
    [Next,~,~]=EnvironmentalSelection(PopObj,GlobalN);
    Population=Population(Next);
end


function Isde=ISDE(Population)
    %��������һ���Ӵ�p���丸����ȺPopulation�ϲ�������ISDE

    PopObj=Population.objs;
    %��Ⱥ�Ĵ�С
    [N,M]=size(PopObj);           %NΪ��Ⱥ��С
    zmin=min(PopObj,[],1);
    zmax=max(PopObj,[],1);
    PopObj = (PopObj-repmat(zmin,N,1))./(repmat(zmax-zmin,N,1)); %�淶��
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