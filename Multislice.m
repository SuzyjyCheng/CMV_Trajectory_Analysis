filenameM1 = "cell1k_mat_d0.csv";
filenameM2 = "cell1k_mat_d1.csv";
filenameM3 = "cell1k_mat_d4.csv";
filenameM4 = "cell1k_mat_d7.csv";
filenameM5 = "cell1k_mat_d14.csv";


A1 = readmatrix(filenameM1,"NumHeaderLines", 1);
A2 = readmatrix(filenameM2, "NumHeaderLines", 1);
A3 = readmatrix(filenameM3, "NumHeaderLines", 1);
A4 = readmatrix(filenameM4, "NumHeaderLines", 1);
A5 = readmatrix(filenameM5, "NumHeaderLines", 1);

A = {A1, A2, A3, A4, A5};


omega = linspace(0.001,40,100);
n = 100 ;

gamma = linspace(0.001,10,100);
com=[];
para_mat = [];
for j = 1:n
    for i = 1:length(gamma)
        N=length(A{1});
        T=length(A);
        B=spalloc(N*T,N*T,N*N*T+2*N*T);
        twomu=0;
        for s=1:T
            k=sum(A{s});
            twom=sum(k);
            twomu=twomu+twom;
            indx=[1:N]+(s-1)*N;
            B(indx,indx)=A{s}-gamma(i)*k'*k/twom;
        end
        twomu=twomu+2*omega*N*(T-1);
        B = B + omega(j)*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
        [S,Q] = genlouvain(B);
        com=[com; S];
    end
    para_mat = [para_mat; S];
end
      
community = reshape(com, 5000, 10000);
max_mat = max(community);
pmat = reshape(max_mat, 100, 100);
writematrix(pmat,'pmat_1k.csv')


