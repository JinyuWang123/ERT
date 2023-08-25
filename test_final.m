

% Orchestration script

% Simulation setup
%
rng(0) % Used seed: 16, 0, 1, 3, 7, 11, 13, 17, 19, 23, 29
tic
parpool(8)
n=10; % grid size
samples=1; % number of samples per category
categories=2;
alpha_type_one_er=0.05;
draws=zeros(n,n,n,samples*categories);
draws_format=zeros(n^3,samples*categories);
x=linspace(-1,1,n)'; y=x; z=x;
[X,Y,Z]=meshgrid(x,y,z);
N=size(y,1);
%scaling: 0: No scaling (bad idea)
%1: The whole shape collection is scaled to [0,1] (standard)
%2: Each group is scaled to [0,1] (more overlap between groups, but messes up interpretation)
scaling=1;
add_noise=0;
delta=1; %%its epsilon
noise_level=0.1;
number_studies=30;
%draws=[]

% Creating the grid

tmp1=repmat(1:N,1,N^2);   
tmp2=repmat(repelem(1:N,N),1,N);
tmp3=repelem(1:N,N^2);
M=[tmp3' tmp2' tmp1'];
m=0;%n/2;


for w=1:number_studies
    coefficients=zeros(samples,4);
    for i=1:samples
        coefficients(i,1)=0.5+0.5*rand(1);
        coefficients(i,2)=0.5+0.5*rand(1);
        coefficients(i,3)=0.5+0.5*rand(1);
        coefficients(i,4)=0.4+0.2*rand(1);
    end
    
    
    
    coefficients_2=zeros(samples,4);
    for i=1:samples
        coefficients_2(i,1)=0.5+0.5*rand(1);
        coefficients_2(i,2)=0.5+0.5*rand(1);
        coefficients_2(i,3)=0.5+0.5*rand(1);
        coefficients_2(i,4)=0.4+0.2*rand(1);
    end
    
    % Sampling stuff
    tic
    
    for kategory=1:categories
        % Would be better to replace this loop
        for sample=1:samples
            index=(samples)*(kategory-1)+sample;
            draw=zeros(size(M,1),1);
        
        if kategory==1
            for i=1:size(M,1)
                row=M(i,:);
                row=[X(i) Y(i) Z(i)];
                tmp_r=sqrt(coefficients(sample,1)*((row(1)-m))^2+coefficients(sample,2)*((row(2)-m))^2);
                draw(i)=((tmp_r-coefficients(sample,4))^2+coefficients(sample,3)*((row(3)-m))^2);
            end
        end
    
        if kategory==2
             for i=1:size(M,1)
                row=M(i,:);
                row=[X(i) Y(i) Z(i)];
                tmp_r_2=sqrt(coefficients_2(sample,1)/delta*((row(1)-m))^2+coefficients_2(sample,2)*delta*((row(2)-m))^2);
                draw(i)=(((tmp_r_2-coefficients_2(sample,4))^2+coefficients_2(sample,3)*((row(3)-m))^2));
            end
        end
    
            Q=reshape(draw,n,n,n);
            draws(:,:,:,index)=Q;
            draws_format(:,index)=draw';
        end
       
    end
    
    toc
    
    
    %data_sert_dir1=SERTs(:,:,1) writematrix(data_sert_dir1)
  
    
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
    
    
    
    tic
    LECTs=zeros(size(directions,1),n_h,length(thresholds),categories*samples);
    
    parfor s=1:(samples*categories)
        donut=ECTs(:,:,:,s)
        lect=zeros(size(directions,1),n_h,length(thresholds))
            for t=1:30
                if t<30
                    lect(:,:,t)=donut(:,:,t)-donut(:,:,t+1);
                else
                    lect(:,:,t)=donut(:,:,t);
                end
            end
        LECTs(:,:,:,s)=lect;
    end
    toc
    
    tic
    MECs=sum(ECTs,3);
    ERTs=zeros(size(directions,1),n_h,(samples*categories));
    marginals_ECT=sum(ECTs,3);
    marginals_LECT=sum(LECTs,3);
    for i=1:(samples*categories)
        ERTs(:,:,i)=marginals_ECT(:,:,:,i)-0.5*(marginals_LECT(:,:,:,i));
    end
    
    T=4;
    t_SERT=linspace(0,4,n_h);
    SERTs=zeros(n_h,size(directions,1),(samples*categories));
    for i=1:(samples*categories)
        ERT_2D_results=ERTs(:,:,i).';
        PERT=zeros(n_h,size(directions,1));
        for s=1:size(directions,1)
            PERT(:,s)=(T/n_h)*cumsum(ERT_2D_results(:,s));
        end
        for t=1:size(directions,1)
            SERTs(:,t,i)=PERT(:,t)-((PERT(n_h,t)/T)*t_SERT).';
        end
    end
    toc
    
    %%Algorithm 3

    F_statistic_original_ert=0;
    F_statistic_original_sert=0;
    F_statistic_original_select=0;
    F_statistic_original_mec=0;

    %%ERT
    DistanceM_ert=zeros(categories*samples,categories*samples);
    for i=1:categories*samples
       ERTi=ERTs(:,:,i);
       for j=1:categories*samples
           ERTj=ERTs(:,:,j);
           DistanceM_ert(i,j)=sqrt(sum((ERTi-ERTj).^2,'all'));
       end
    end

    for i=1:samples
        for j=1:samples
            F_statistic_original_ert=F_statistic_original_ert+DistanceM_ert(i,j)+DistanceM_ert(i+samples,j+samples);
        end
    end

    %%SERT
    DistanceM_sert=zeros(categories*samples,categories*samples);
    for i=1:categories*samples
       SERTi=SERTs(:,:,i);
       for j=1:categories*samples
           SERTj=SERTs(:,:,j);
           DistanceM_sert(i,j)=sqrt(sum((SERTi-SERTj).^2,'all'));
       end
    end

    for i=1:samples
        for j=1:samples
            F_statistic_original_sert=F_statistic_original_sert+DistanceM_sert(i,j)+DistanceM_sert(i+samples,j+samples);
        end
    end

    %%SECELT
    DistanceM_select=zeros(categories*samples,categories*samples);
    for i=1:categories*samples
        ECTi=ECTs(:,:,:,i);
        for j=1:categories*samples
            ECTj=ECTs(:,:,:,j);
            DistanceM_select(i,j)=sqrt(sum((ECTi-ECTj).^2,'all'));
        end
    end

    for i=1:samples
        for j=1:samples
            F_statistic_original_select=F_statistic_original_select+DistanceM_select(i,j)+DistanceM_select(i+samples,j+samples);
        end
    end

    %%MEC
    DistanceM_mec=zeros(categories*samples,categories*samples);
    for i=1:categories*samples
        MECTi=MECs(:,:,:,i);
        for j=1:categories*samples
            MECTj=MECs(:,:,:,j);
            DistanceM_mec(i,j)=sqrt(sum((MECTi-MECTj).^2,'all'));
        end
    end

    for i=1:samples
        for j=1:samples
            F_statistic_original_mec=F_statistic_original_mec+DistanceM_mec(i,j)+DistanceM_mec(i+samples,j+samples);
        end
    end

    %%Compute the mean functions of the two SECT collections
    
    tic
    M_1=zeros(n_h,size(directions,1));
    M_2=zeros(n_h,size(directions,1));
    for i=1:samples
        M_1=M_1+SERTs(:,:,i);
        M_2=M_2+SERTs(:,:,i+samples);
    end
    M_1=(1/samples)*M_1;
    M_2=(1/samples)*M_2;
    toc
    
    %%Select the distinguishing direction
    
    tic 
    supnorm=zeros(1,size(directions,1));
    for i=1:size(directions,1)
        supnorm(1,i)=max(abs(M_1(:,i)-M_2(:,i)));
    end
    [max_supnorm,distinguishing_direction_index]=max(supnorm);
    toc
    
    %%Karhunen-Loeve decomposition
    
    tic
    SERT_1_distinguishing_direction=zeros(n_h,samples);
    SERT_2_distinguishing_direction=zeros(n_h,samples);
    for i=1:samples
        SERT_1_distinguishing_direction(:,i)=SERTs(:,distinguishing_direction_index,i);
    end
    
    for i=1:samples
        SERT_2_distinguishing_direction(:,i)=SERTs(:,distinguishing_direction_index,i+samples);
    end
    
    demean_all=zeros(n_h,(samples*categories));
    for i=1:samples
        demean_all(:,i)=SERT_1_distinguishing_direction(:,i)-M_1(:,distinguishing_direction_index);
    end
    
    for i=1:samples
        demean_all(:,i+samples)=SERT_2_distinguishing_direction(:,i)-M_2(:,distinguishing_direction_index);
    end
    CovKer=cov((demean_all).');
    eigen_results_val=flip(eig(CovKer));
    [eigen_results_vec,ww]=eig(CovKer);
    eigen_results_vec=flip(eigen_results_vec,2);
    toc
    
    %%We only care about the components having 95% cumulative variance.
    
    tic
    L_coef=(cumsum(eigen_results_val)/sum(eigen_results_val));
    if L_coef(1)>0.97
        L=1;
    else
        L=max(find((cumsum(eigen_results_val)/sum(eigen_results_val))<0.97));
    end
    Eigen_values=(T/n_h)*eigen_results_val(1:L);
    Eigen_vectors=sqrt(n_h/T)*eigen_results_vec(:,1:L);
    Xis=zeros(samples,L);
    
    for i=1:samples
        if L==1
            Xi=(1/sqrt(2*Eigen_values(1)))*(T/n_h)*sum((SERT_1_distinguishing_direction(:,i)-SERT_2_distinguishing_direction(:,i)).*Eigen_vectors);
        else
            for l=1:L
                Xi(l)=(1/sqrt(2*Eigen_values(l)))*(T/n_h)*sum((SERT_1_distinguishing_direction(:,i)-SERT_2_distinguishing_direction(:,i)).*Eigen_vectors(:,l));
            end
        end
        Xis(i,:)=Xi;
    end
    
    toc
    
    
    
    tic
    chisq_statistic_original=sum((sqrt(samples)*mean(Xis,1)).^2);
    df_seq=L;
    p_value_seq_al1=1-chi2cdf(sum((sqrt(samples)*mean(Xis,1)).^2), L);
    
    toc
    
    p_value_seq_all_al1(w)=p_value_seq_al1;

    


    %%Permutation test
    number_permutations=1000;
    for perm=1:number_permutations
        all_indices=1:(2*samples);
        sampled_indices=randsample(all_indices, samples);
        all_indices(sampled_indices)=[];
        unsampled_indices=all_indices;

        F_statistic_original_ert_perm=0;
        F_statistic_original_sert_perm=0;
        F_statistic_original_select_perm=0;
        F_statistic_original_mec_perm=0;

        for i=1:samples
            for j=1:samples
                F_statistic_original_ert_perm=F_statistic_original_ert_perm+DistanceM_ert(sampled_indices(i),sampled_indices(j))+DistanceM_ert(unsampled_indices(i),unsampled_indices(j));
            end
        end        
        F_statistic_seq_ert(perm)=F_statistic_original_ert_perm;


        for i=1:samples
            for j=1:samples
                F_statistic_original_sert_perm=F_statistic_original_sert_perm+DistanceM_sert(sampled_indices(i),sampled_indices(j))+DistanceM_sert(unsampled_indices(i),unsampled_indices(j));
            end
        end        
        F_statistic_seq_sert(perm)=F_statistic_original_sert_perm;


        for i=1:samples
            for j=1:samples
                F_statistic_original_select_perm=F_statistic_original_select_perm+DistanceM_select(sampled_indices(i),sampled_indices(j))+DistanceM_select(unsampled_indices(i),unsampled_indices(j));
            end
        end        
        F_statistic_seq_select(perm)=F_statistic_original_select_perm;



        for i=1:samples
            for j=1:samples
                F_statistic_original_mec_perm=F_statistic_original_mec_perm+DistanceM_mec(sampled_indices(i),sampled_indices(j))+DistanceM_mec(unsampled_indices(i),unsampled_indices(j));
            end
        end        
        F_statistic_seq_mec(perm)=F_statistic_original_mec_perm;



        M_1_perm=zeros(n_h,size(directions,1));
        M_2_perm=zeros(n_h,size(directions,1));
        for i=1:samples
            M_1_perm=M_1_perm+SERTs(:,:,sampled_indices(i));
            M_2_perm=M_2_perm+SERTs(:,:,unsampled_indices(i));
        end
        M_1_perm=(1/samples)*M_1_perm;
        M_2_perm=(1/samples)*M_2_perm;

        supnorm_perm=zeros(1,size(directions,1));
        for i=1:size(directions,1)
            supnorm_perm(1,i)=max(abs(M_1_perm(:,i)-M_2_perm(:,i)));
        end
        [max_supnorm_perm,distinguishing_direction_index_perm]=max(supnorm_perm);

        SERT_1_distinguishing_direction_perm=zeros(n_h,samples);
        SERT_2_distinguishing_direction_perm=zeros(n_h,samples);
        for i=1:samples
            SERT_1_distinguishing_direction_perm(:,i)=SERTs(:,distinguishing_direction_index_perm,sampled_indices(i));
        end
        
        for i=1:samples
            SERT_2_distinguishing_direction_perm(:,i)=SERTs(:,distinguishing_direction_index_perm,unsampled_indices(i));
        end
        
        demean_all_perm=zeros(n_h,(samples*categories));
        for i=1:samples
            demean_all_perm(:,i)=SERT_1_distinguishing_direction_perm(:,i)-M_1_perm(:,distinguishing_direction_index_perm);
        end
        
        for i=1:samples
            demean_all_perm(:,i+samples)=SERT_2_distinguishing_direction_perm(:,i)-M_2_perm(:,distinguishing_direction_index_perm);
        end
        CovKer_perm=cov((demean_all_perm).');
        eigen_results_val_perm=flip(eig(CovKer_perm));
        [eigen_results_vec_perm,ww]=eig(CovKer_perm);
        eigen_results_vec_perm=flip(eigen_results_vec_perm,2);

        L_coef_perm=(cumsum(eigen_results_val_perm)/sum(eigen_results_val_perm));
        if L_coef_perm(1)>0.97
            L_perm=1;
        else
            L_perm=max(find((cumsum(eigen_results_val_perm)/sum(eigen_results_val_perm))<0.97));
        end
        Eigen_values_perm=(T/n_h)*eigen_results_val_perm(1:L_perm);
        Eigen_vectors_perm=sqrt(n_h/T)*eigen_results_vec_perm(:,1:L_perm);
        Xis_perm=zeros(samples,L_perm);
        
        for i=1:samples
            if L_perm==1
                Xi=(1/sqrt(2*Eigen_values_perm(1)))*(T/n_h)*sum((SERT_1_distinguishing_direction_perm(:,i)-SERT_2_distinguishing_direction_perm(:,i)).*Eigen_vectors_perm);
            else
                for l=1:L_perm
                    Xi(l)=(1/sqrt(2*Eigen_values_perm(l)))*(T/n_h)*sum((SERT_1_distinguishing_direction_perm(:,i)-SERT_2_distinguishing_direction_perm(:,i)).*Eigen_vectors_perm(:,l));
                end
            end
            Xis_perm(i,:)=Xi;
        end
        chisq_statistic_seq(perm)=sum((sqrt(samples)*mean(Xis_perm,1)).^2);
    end
    increasing_chisq_stat=sort(chisq_statistic_seq);

    if chisq_statistic_original>max(increasing_chisq_stat)
        p_value_seq_all_al2(w)=0;
    else
        p_value_seq_all_al2(w)=1-max(find(chisq_statistic_original>increasing_chisq_stat))/number_permutations;
    end
    
    decreasing_F_stat_ert=sort(F_statistic_seq_ert,'descend');
    if F_statistic_original_ert<min(decreasing_F_stat_ert)
        p_value_seq_all_al3_ert(w)=0;
    else
        p_value_seq_all_al3_ert(w)=1-min(find(decreasing_F_stat_ert<F_statistic_original_ert))/number_permutations;
    end

    
    decreasing_F_stat_sert=sort(F_statistic_seq_sert,'descend');
    if F_statistic_original_sert<min(decreasing_F_stat_sert)
        p_value_seq_all_al3_sert(w)=0;
    else
        p_value_seq_all_al3_sert(w)=1-min(find(decreasing_F_stat_sert<F_statistic_original_sert))/number_permutations;
    end

    decreasing_F_stat_select=sort(F_statistic_seq_select,'descend');
    if F_statistic_original_select<min(decreasing_F_stat_select)
        p_value_seq_all_al3_select(w)=0;
    else
        p_value_seq_all_al3_select(w)=1-min(find(decreasing_F_stat_select<F_statistic_original_select))/number_permutations;
    end

    decreasing_F_stat_mec=sort(F_statistic_seq_mec,'descend');
    if F_statistic_original_mec<min(decreasing_F_stat_mec)
        p_value_seq_all_al3_mec(w)=0
    else
        p_value_seq_all_al3_mec(w)=1-min(find(decreasing_F_stat_mec<F_statistic_original_mec))/number_permutations;
    end
end


%%p_value_seq_all_al1 p_value_seq_all_al2  p_value_seq_all_al3_ert
%%p_value_seq_all_al3_sert  p_value_seq_all_al3_select  p_value_seq_all_al3_mec











