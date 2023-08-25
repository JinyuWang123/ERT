function ECT=compute_ECT_3d(vertices,edges,triangles,tetrahedra,directions,heights)
p=size(heights,1);
V=vertices';
n=size(directions,1);
ECT=zeros(n,p);
for i=1:n
   direction=directions(i,:);
   evals=(direction*V);
   for j=1:p
      test=evals-heights(j)<=0;
      inds=find(test);
      v=sum(test);
      e=sum(and(ismember(edges(:,1),inds),ismember(edges(:,2),inds)));
      t=sum(and(ismember(triangles(:,1),inds),...
          and(ismember(triangles(:,2),inds),...
          ismember(triangles(:,3),inds))));
      f= sum(and(and(ismember(tetrahedra(:,1),inds),...
       ismember(tetrahedra(:,2),inds)),...
       and(ismember(tetrahedra(:,3),inds),...
       ismember(tetrahedra(:,4),inds))));
      ECT(i,j)=v-e+t-f;
   end
   
end
end