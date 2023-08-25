function [vertices,Ed,Tr,Tet,flag]=Compute_superlevelset(threshold,Q,X,Y,Z,draw)
flag=0;
XMAX=max(max(max(X)))+0.01;
YMAX=max(max(max(Y)))+0.01;
ZMAX=max(max(max(Z)))+0.01;
XMIN=min(min(min(X)))-0.01;
YMIN=min(min(min(Y)))-0.01;
ZMIN=min(min(min(Z)))-0.01;

% If threshold is more than maximum, we have an empty set and we
% just return the flag
if threshold>max(max(max(Q)))
    flag=-1;
    vertices=[];Ed=[];Tr=[];Tet=[];
    return
end

% If threshold is less than the minimum, return a rectangle
if threshold <min(min(min(Q)))
    flag=1;
    vertices=[XMIN YMIN ZMIN; XMIN YMIN ZMAX;...
    XMIN YMAX ZMIN;XMIN YMAX ZMAX; XMAX YMIN ZMIN; XMAX YMIN ZMAX;...
    XMAX YMAX ZMIN; XMAX YMAX ZMAX];
    Ed=[5 1;5 7;5 8;1 7;1 8;7 8;5 6;1 6;6 8;1 3;3 7;3 8;1 4;3 4;4 8;...
        1 2;6 2;2 8;2 4];
    Tr=[5 1 7;5 1 8;5 7 8;1 7 8;5 1 6;5 6 8;1 6 8;1 3 7; 1 3 8;3 7 8;1 3 4;...
        1 4 8;3 4 8; 1 6 2;1 2 8;6 2 8;2 4 8;1 2 4];
    Tet=[5 1 7 8;5 1 6 8;1 3 7 8;1 3 4 8;1 6 2 8;1 2 4 8];
    return
end

levelset=isosurface(X,Y,Z,Q,threshold);

supplement=[levelset.vertices; XMIN YMIN ZMIN; XMIN YMIN ZMAX;...
    XMIN YMAX ZMIN;XMIN YMAX ZMAX; XMAX YMIN ZMIN; XMAX YMIN ZMAX;...
    XMAX YMAX ZMIN; XMAX YMAX ZMAX];

delT=delaunay(supplement);
[dual,adjacencies]=dualize3d(delT);



boundary_weight=size(delT,1)+1;

%% testing. Improved version

weights=zeros(1,size(adjacencies,1));
for i=1:size(adjacencies,1)
   triangle=adjacencies(i,:);
   v1=triangle(1);v2=triangle(2);v3=triangle(3);
   test=isempty(find(and(any(levelset.faces==v1,2),and(any(levelset.faces==v2,2),any(levelset.faces==v3,2))),1));
   if not(test)%length(test>0)
       weights(i)=boundary_weight;
       %weights=[weights,boundary_weight];
   else
       weights(i)=1;
       %weights=[weights,1];
   end
end

%% testing ends




G=graph(dual(:,1),dual(:,2),weights);
Tree=minspantree(G);
max_ind=find(draw==max(draw),1);
TestVertex=[X(max_ind) Y(max_ind) Z(max_ind)];



kikka=zeros(size(delT,1),1);
index=0;
for i= 1:size(delT,1)
    %tetra=levelset.vertices(delT(i,:),:);
    tetra=supplement(delT(i,:),:);
    %if inpolygon(TestVertex(1),TestVertex(2),tetra(:,1),tetra(:,2))...
    %        & inpolygon(TestVertex(1),TestVertex(3),tetra(:,1),tetra(:,3))...
    %        & inpolygon(TestVertex(2),TestVertex(3),tetra(:,2),tetra(:,3))
    test=tetrahedron_membership_test(tetra,TestVertex);
    kikka(i)=test;
    if (test==1)
        index=i;
    end
end


keepers=zeros(max(max(G.Edges.EndNodes)),1);
for i=1:max(max(G.Edges.EndNodes))
   [~,d]=shortestpath(Tree,i,index);
   test=ceil(d/boundary_weight);
   keepers(i)=mod(test,2);
end



vertices=levelset.vertices;
tetrahedra=delT(find(keepers),:);

triangles=unique([tetrahedra(:,[1 2 3]);tetrahedra(:,[1 2 4]);...
    tetrahedra(:,[1 3 4]);tetrahedra(:,[2 3 4])],'rows');

edges=unique([triangles(:,[1 2]);triangles(:,[2 3]);triangles(:,[1 3])],'rows');
[~,idx]=unique(sort(edges')','rows','stable');
Ed=edges(idx,:);

[~,idx]=unique(sort(triangles')','rows','stable');
Tr=triangles(idx,:);

[~,idx]=unique(sort(tetrahedra')','rows','stable');
Tet=tetrahedra(idx,:);
end
