% Computes the dual graph of a Delaunay triangulation
% in 3d

% Update July 29: Reprogrammed for speed
function [dual,relations]=dualize3d(Triangulation)
n=size(Triangulation,1);
% Graph: Just list of edges
% Would also like: The type of the edge. Maybe separate function?
%dual=[];
% Each tetrahedron has 4 faces, so at most 4 neighbors.
% By the handshake lemma the dual contains at most 2n edges
dual=zeros(2*n,2);
relations=zeros(2*n,3);
%relations=[];
%triangle_relations=[];
j=0;
for i=1:n
    
    % Search for adjacent tetrahedra
    v1=Triangulation(i,1);v2=Triangulation(i,2);
    v3=Triangulation(i,3);v4=Triangulation(i,4);
    index123=find(and(any(Triangulation==v1,2),and(any(Triangulation==v2,2),any(Triangulation==v3,2))));
    index124=find(and(any(Triangulation==v1,2),and(any(Triangulation==v2,2),any(Triangulation==v4,2))));
    index134=find(and(any(Triangulation==v1,2),and(any(Triangulation==v3,2),any(Triangulation==v4,2))));
    index234=find(and(any(Triangulation==v2,2),and(any(Triangulation==v3,2),any(Triangulation==v4,2))));
    
    % Add connections
    if length(index123)==2 && min(index123)==i
    %if length(index123)==2 & sum(and(any(dual==v1,2),and(any(dual==v2,2),any(dual==v3,2))))==0
        j=j+1;
        dual(j,:)=index123';
        %dual=[dual;index123'];
        relations(j,:)=intersect(Triangulation(index123(1),:),Triangulation(index123(2),:));
        %relations=[relations; intersect(Triangulation(index123(1),:),Triangulation(index123(2),:))];
        % Todo: Check this is not needed right?
        % triangle_relations=[triangle_relations; index123'];
    end
    
    if length(index124)==2 && min(index124)==i
    %if length(index124)==2 & sum(and(any(dual==v1,2),and(any(dual==v2,2),any(dual==v4,2))))==0
        j=j+1;
        dual(j,:)=index124';
        %dual=[dual;index123'];
        relations(j,:)=intersect(Triangulation(index124(1),:),Triangulation(index124(2),:));
    
        %dual=[dual;index124'];
        %relations=[relations; intersect(Triangulation(index124(1),:),Triangulation(index124(2),:))];
    end
    
    if length(index134)==2 && min(index134)==i
    %if length(index134)==2 & sum(and(any(dual==v1,2),and(any(dual==v3,2),any(dual==v4,2))))==0
        j=j+1;
        dual(j,:)=index134';
        %dual=[dual;index123'];
        relations(j,:)=intersect(Triangulation(index134(1),:),Triangulation(index134(2),:));
    
        %dual=[dual;index134'];
        %relations=[relations; intersect(Triangulation(index134(1),:),Triangulation(index134(2),:))];
    end
    
    
    if length(index234)==2 && min(index234)==i
    %if length(index234)==2 & sum(and(any(dual==v2,2),and(any(dual==v3,2),any(dual==v4,2))))==0
        j=j+1;
        dual(j,:)=index234';
        %dual=[dual;index123'];
        relations(j,:)=intersect(Triangulation(index234(1),:),Triangulation(index234(2),:));
    
        %dual=[dual;index234'];
        %relations=[relations; intersect(Triangulation(index234(1),:),Triangulation(index234(2),:))];
    end
    
end
dual=dual(1:j,:);
relations=relations(1:j,:);
end