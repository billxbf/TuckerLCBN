TimeListExp = [];
TimeListImp = [];
SpeedUpList = [];
FlopListExp = [];
FlopListImp = []; 
AList = [100:100:2500]';

for a = [100:100:2500]

    W = 60;

    % Generate normal distributed tensor
    A = normrnd(2,3,[a,147,6]);
    
    
    %profile on
    %tic
    [U1,~,~,etime1] = fast_HOSVD(A,W);
    %etime1 = toc;
    %profileStruct = profile('info');
    %[flopTotal,~] = FLOPS('fast_HOSVD','fast_HOSVDMat',profileStruct);
    TimeListImp = [TimeListImp; etime1];
    %FlopListImp = [FlopListImp; flopTotal];
    

    %profile on
    %tic 
    [U2,~,~,etime2] = corr_HOSVD(A,W);
    %etime2 = toc;
    %profileStruct = profile('info');
    %[flopTotal2,~] = FLOPS('corr_HOSVD','corr_HOSVDMat',profileStruct);
    TimeListExp = [TimeListExp; etime2];
    %FlopListExp = [FlopListExp; flopTotal2];
    
    speedup = etime2/etime1;
    SpeedUpList = [SpeedUpList; speedup];
    
    
    disp("time for implicit a = " + a + ": " + etime1);
    disp("time for explicit a = " + a + ": " + etime2);
    disp("time speedup for a = " + a + ": " + speedup);
    disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    
end

T = table(AList, TimeListExp, TimeListImp, SpeedUpList)
writetable(T,'ExperimentData1.csv')