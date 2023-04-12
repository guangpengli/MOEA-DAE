function DAE(Global)
%<algorithm><D>
%MOEA/D-epsilonDAE

    %% 初始化种群和参考向量
    Population=Global.Initialization();
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    T  = ceil(Global.N/10);
    nr = 2;
    %% Detect the neighbours of each solution
    %计算向量间的距离
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,2:T);
    Z=min(Population.objs,[],1);
    
    %% epsilon技术初始化
    cp=3;
    Tc=0.8*Global.maxgen;
    scons_0=sum(max(0,Population.cons),2);
    epsilon_0=max(scons_0);
    Gif=0;
    epsilon_pre=zeros(1,Global.maxgen);
    
    delta=0.9;
    %外部档集
    arch = archive(Population,Global.N);
    
    %% 判断初始种群中的可行解所占比例是否超过0.95，超过则之后启动逃离策略
    feaP=scons_0==0;
    feaP=sum(feaP)/Global.N;
    flag=false;
    if feaP>=0.95
        flag=true;
        %创建矩阵，保留种群每一代的违约度和目标值
        CP=zeros(1,Global.maxgen);
        FP=zeros(1,Global.maxgen);
    end
    
    %逃离策略标志
    CH1=false;
    maxconCH1=0;
    
    while Global.NotTermination(Population)
        if flag==true
            scons_1=sum(max(0,Population.cons),2);
            CP(Global.gen)=sum(scons_1);
            fp=sum(Population.objs,2);
            FP(Global.gen)=sum(fp);
            if(Global.gen>10)
                ROC=(abs(CP(Global.gen)-CP(Global.gen-10)))/CP(Global.gen);
                ROF=(abs(FP(Global.gen)-FP(Global.gen-10)))/FP(Global.gen);
                if ROF<=0.00001
                    CH1=true;
                    T=Global.gen;
                end
            end
        end
     
        %条件1 flag=true 
        %条件2 收敛到当前可行域的PF，然后执行双目标优化法
        if flag==true&&CH1==true
            scons2=sum(max(0,Population.cons),2);
            maxcons2=max(scons2);
            if maxcons2>maxconCH1
                maxconCH1=maxcons2;
            end
            %采用逃离策略CH1---双目标法
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population.objs,Global.N);
            MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
            Offspring  = GA(Population(MatingPool));
            Population=[Population,Offspring];
            Population=CH1Method(Population,Global.N);
            T1=Global.gen-T;
            %判断结束条件 运行T代后，违约度变化率，
            if (T1>20&&ROC<0.00001)||(Global.gen>Tc)
                flag=false;
                Gif=Global.gen;
                epsilon_pre(Global.gen-1)=maxconCH1;
            end
        else
    
            scons3=sum(max(0,Population.cons),2);
            maxv=max(scons3);
            minv=min(scons3);

            [epsilon_k,epsilon_0,Gif,epsilon_pre]=up_eplison2(epsilon_0,epsilon_pre,Global.gen,Tc,cp,maxv,minv,Gif,Global.M);

            for i = 1 : Global.N
                % Choose the parents
                if rand < delta
                    P = B(i,randperm(end));
                else
                    P = randperm(Global.N);
                end

                % Generate an offspring
                Offspring = DE(Population(i),Population(P(1)),Population(P(2)));

                CV0=sum(max(0,Offspring.con));
                CVP=sum(max(0,Population(P).cons),2);

                % Update the ideal point
                Z = min(Z,Offspring.obj);

                 % Update the solutions in P by PBI
                normW   = sqrt(sum(W(P,:).^2,2));
                normP   = sqrt(sum((Population(P).objs-repmat(Z,length(P),1)).^2,2));
                normO   = sqrt(sum((Offspring.obj-Z).^2,2));
                CosineP = sum((Population(P).objs-repmat(Z,length(P),1)).*W(P,:),2)./normW./normP;
                CosineO = sum(repmat(Offspring.obj-Z,length(P),1).*W(P,:),2)./normW./normO;
                g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
                np=size(Population(P),2);
                %采用eplison约束处理技术
                next=CH0Method(g_old,g_new,np,CV0,CVP,epsilon_k,nr);
                Population(P(next))=Offspring;

            end
       
        end
        
         % Output the non-dominated and feasible solutions.
        arch = archive([arch,Population],Global.N);
        if Global.evaluated >= Global.evaluation
            Population = arch;
        end
        
    end
end

%本文所用
function [result,epsilon_0,Gif,epsilon_pre]=up_eplison2(epsilon_0,epsilon_pre,gen,Tc,cp,maxv,minv,Gif,m)
    if gen<=Tc
        if gen==1
           result=epsilon_0*((1-((gen-Gif)/Tc))^cp);
           epsilon_pre(gen)=result; 
        else
            if epsilon_pre(gen-1)<=minv
                epsilon_0=maxv/m;
                Gif=gen;
                result=epsilon_0*((1-((gen-Gif)/Tc))^cp);
                epsilon_pre(gen)=result;
            else
                result=epsilon_0*((1-((gen-Gif)/Tc))^cp);
                epsilon_pre(gen)=result;
            end
        end
        
    else
        result=0;
    end

end


function result=up_eplison(epsilon_0,gen,Tc,cp)
    if gen<=Tc
        result=epsilon_0*((1-(gen/Tc))^cp);
    else
        result=0;
    end
    
end
function result=up_eplison1(epsilon_0,gen,Tc,cp,maxv)
    if gen<=Tc
        if epsilon_0==0
            result=maxv;
        else
            result=epsilon_0*((1-(gen/Tc))^cp);
        end
    else
        result=0;
    end
end
function [result,epsilon_0,Gif,epsilon_pre]=up_eplison3(epsilon_0,epsilon_pre,gen,Tc,cp,maxv,minv,Gif,m)
    if gen<=Tc
        if gen==1
           result=epsilon_0*((1-((gen-Gif)/Tc))^cp);
           epsilon_pre(gen)=result; 
        else
            if epsilon_pre(gen-1)<=minv
                epsilon_0=maxv/m;
%                 epsilon_0=maxv;
                Gif=gen;
                result=epsilon_0*((1-((gen-Gif)/Tc))^cp);
                epsilon_pre(gen)=result;
            elseif epsilon_pre(gen-1)>maxv
                epsilon_0=maxv;
                Gif=gen;
                result=epsilon_0*((1-((gen-Gif)/Tc))^cp);
                epsilon_pre(gen)=result;
            else

                result=epsilon_0*((1-((gen-Gif)/Tc))^cp);
                epsilon_pre(gen)=result;
            end
        end
        
    else
        result=0;
    end

end
function result=up_eplison4(epsilon_pre,gen,Tc,rf,maxc)
    if gen<Tc
        if rf<0.95
            result=(1-0.1)*epsilon_pre;
        else
            result=(1+0.1)*maxc;
        end
    else
        result=0;
    end
end













