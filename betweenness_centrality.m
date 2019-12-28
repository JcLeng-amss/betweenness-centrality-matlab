function Cb=betweenness_centrality(A)
p=size(A,1);
V=1:p;
Cb=zeros(1,p);
for s=V
    S=[];
    P=cell(p,1);
    sigma=zeros(1,p);
    sigma(s)=1;
    d=-ones(1,p);
    d(s)=0;
    Q=[];
    Q=[Q,s];
    while ~isempty(Q)
        v=Q(1);
        Q(1)=[];
        S=[S,v];
        neighbors_of_v=find(A(v,:)~=0);
        for w=neighbors_of_v
            if d(w)<0
                Q=[Q,w];
                d(w)=d(v)++1;
            end
            if d(w)==d(v)+1
                sigma(w)=sigma(w)+sigma(v);
                P{w}=[P{w},v];
            end
        end
    end
    delta=zeros(1,p);
    while ~isempty(S)
        w=S(end);
        S(end)=[];
        for v=P{w}
            delta(v)=delta(v)+sigma(v)/sigma(w)*(1+delta(w));
        end
        if w~=s
            Cb(w)=Cb(w)+delta(w);
        end
    end
end
end