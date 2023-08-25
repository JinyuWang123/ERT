% Orchestration script

% Simulation setup
%
rng('default')
tic
parpool(8)
n=10; % grid size
samples=1; % number of samples per category
categories=4;
draws=zeros(n,n,n,samples*categories);
draws_format=zeros(n^3,samples*categories);
x=linspace(-1,1,n)'; y=x; z=x;
[X,Y,Z]=meshgrid(x,y,z);
N=size(y,1);
%scaling: 0: No scaling (bad idea)
%1: The whole shape collection is scaled to [0,1] (standard)
%2: Each group is scaled to [0,1] (more overlap between groups, but messes up interpretation)
scaling=1
add_noise=0
noise_level=0.1;
%draws=[]

% Creating the grid

tmp1=repmat(1:N,1,N^2);
tmp2=repmat(repelem(1:N,N),1,N);
tmp3=repelem(1:N,N^2);
M=[tmp3' tmp2' tmp1'];
m=0;%n/2;

% Sampling stuff
tic
coefficients=zeros(categories*samples,3);
for kategory=1:categories
    % Would be better to replace this loop
    for sample=1:samples
        index=(samples)*(kategory-1)+sample;
        coefficients(index,:)=0.5+0.5*rand(1,3);coefs=coefficients(index,:);
    draw=zeros(size(M,1),1);
    if kategory==1
       signs=[1 1 1]; 
    end
    if kategory==2
        signs=[1 1 -1];
    end
    if kategory==3
       signs=[1 -1 -1];
    end
    for i=1:size(M,1)
        row=[X(i) Y(i) Z(i)];
        %row=M(i,:);
        draw(i)=coefs(1)*((row(1)-m))^2+coefs(2)*signs(2)*((row(2)-m))^2+coefs(3)*signs(3)*((row(3)-m))^2;
    end
    
    if kategory==4
        for i=1:size(M,1)
            row=M(i,:);
            row=[X(i) Y(i) Z(i)];
            tmp_r=sqrt(coefs(1)*((row(1)-m))^2+coefs(2)*((row(2)-m))^2);
            draw(i)=exp(-(tmp_r-0.4-0.2*rand(1))^2+coefs(3)*((row(3)-m))^2);
        end
    end

        Q=reshape(draw,n,n,n);
        draws(:,:,:,index)=Q;
        draws_format(:,index)=draw';
    end
    if add_noise==1
        noise=normrnd(0,noise_level,numel(draws),1);
        a_noise=reshape(noise,N,N,N,categories*samples);
        draws=draws+a_noise;
        draws_format=reshape(draws,N^3,categories*samples);%,draws_format+reshape(noise,N^3,categories*samples);
    end
end

toc

tic

if scaling==2
    for k=1:categories
        ind_s=1+(k-1)*samples;
        ind_e=5+(k-1)*samples;
        maxx=max(draws(:,:,:,ind_s:ind_e),[],'all');
        minn=min(draws(:,:,:,ind_s:ind_e),[],'all');
        scale=maxx-minn;
        draws(:,:,:,ind_s:ind_e)=(draws(:,:,:,ind_s:ind_e)-minn)/scale;
        draws_format(:,ind_s:ind_e)=(draws_format(:,ind_s:ind_e)-minn)/scale;
    end
end


if scaling==1
    maxx=max(draws,[],'all');
    minn=min(draws,[],'all');
    scale=maxx-minn;
    draws=(draws-minn)/scale;
    draws_format=(draws_format-minn)/scale;
end


toc



thresholds=(1:30)/30;
thresholds(end)=thresholds(end)-0.001;

Vertex_cell=cell(samples*categories,size(thresholds,1));
edge_cell=cell(samples*categories,size(thresholds,1));
triangle_cell=cell(samples*categories,size(thresholds,1));
tetra_cell=cell(samples*categories,size(thresholds,1));
flags=zeros(samples*categories,size(thresholds,1));


tmp_n=length(thresholds);
tic
parfor s=1:size(draws_format,2)
    scaled_draw=draws_format(:,s);
    Q=draws(:,:,:,s);
    for t=1:tmp_n
        threshold=thresholds(t);
        % Suggest: Q=reshape(scaled_draw,20,20,20);
        [vertices,edges,triangles,tetrahedra,flag]=Compute_superlevelset(threshold,Q,X,Y,Z,scaled_draw);
        Vertex_cell{s,t}=vertices;
        edge_cell{s,t}=edges;
        triangle_cell{s,t}=triangles;
        tetra_cell{s,t}=tetrahedra;
        flags(s,t)=flag;
    end
end
sound(sin(1:3000));
toc

tic
directions=load("directions326.csv");
n_h=100;
radius=sqrt(3);%sqrt(8/3)%*n;
heights=linspace(m-radius,m+radius,n_h)';
ECTs=zeros(size(directions,1),n_h,length(thresholds),categories*samples);

parfor s=1:size(draws_format,2)
    ECT=zeros(size(directions,1),n_h,length(thresholds));
    for t=1:tmp_n
        if flags(s,t)==-1
            continue
        else
        ECT(:,:,t)=compute_ECT_3d(Vertex_cell{s,t},edge_cell{s,t},triangle_cell{s,t},tetra_cell{s,t},directions,heights);
        end
    end
    ECTs(:,:,:,s)=ECT;
    disp(s)
end
sound(0.5*sin(1:3000));
toc
