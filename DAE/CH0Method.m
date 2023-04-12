function next=EnvironSelection(g_old,g_new,n,CV0,CVP,epsilon_k,nr)
    next=false(1,n);
    for i=1:n
        if sum(next)<=nr
            if CV0<=epsilon_k&&CVP(i)<epsilon_k
                if g_old(i)>=g_new(i)
                    next(i)=true;
                end
            else
                if CV0==CVP(i)
                    if g_old(i)>=g_new(i)
                        next(i)=true;
                    end
                else
                    if CV0<CVP(i)
                        next(i)=true;
                    end
                end
                    
            end
        
        end
    end
end