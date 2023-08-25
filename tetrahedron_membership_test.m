function is_member=tetrahedron_membership_test(tetrahedron,point)
% Test is a point is contained in a closed tetrahedron
% This is achieved by a hyperplane containment argument
% Right side found by the omitted point
is_member=0;
v0=tetrahedron(1,:);v1=tetrahedron(2,:);v2=tetrahedron(3,:);v3=tetrahedron(4,:);
s01=v1-v0;s02=v2-v0;s03=v3-v0;
s12=v2-v1;s13=v3-v1;s23=v3-v2;
n012=cross(s01,s02);
n013=cross(s01,s03);
n023=cross(s02,s03);
n123=cross(s12,s13);
if dot(n012,v3-v0)*dot(n012,point-v0)>=0 && dot(n013,v2-v0)*dot(n013,point-v0)>=0 && ...
        dot(n023,v1-v0)*dot(n023,point-v0)>=0 && dot(n123,v0-v1)*dot(n123,point-v1)>=0
    is_member=1;
end

end