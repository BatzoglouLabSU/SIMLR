function [cluster_result, classification_result] = GB_get_metric(F,ytrue)

%% clustering results
F = F./repmat(sqrt(mean(F.^2)),size(F,1),1);
if isempty(F)
    cluster_result=zeros(1,3);
    classification_result=zeros(1,3);
else
    
    c = max(ytrue);
    cluster_result = zeros(1,3);
    for t = 1:10
        predY = litekmeans(F,c,'replicates',100);
        cluster_result =  cluster_result + ClusteringMeasure(ytrue, predY);
    end
    cluster_result = cluster_result/t;
    
    %% Classification results
%     KNNMdl = fitcknn(F,ytrue,...
%         'NumNeighbors',1,'Standardize',1);
%     %CVKNNMdl = crossval(KNNMdl,'Leaveout','on');
%     CVKNNMdl = crossval(KNNMdl);
%     classError = kfoldLoss(CVKNNMdl);
%     classification_result(1) = classError;
%   
D = dist2(F);
for t = 1:20
    nfold = 10;
    
    for i = 1:max(ytrue)
        ind = find(ytrue==i);
        randomind{i} = ind(randperm(length(ind)));
    end
    for i = 1:nfold
        testind = [];
        for k = 1:max(ytrue)
            step = round(length(randomind{k})/nfold);
            testind=[testind;randomind{k}((i-1)*step+1:min(i*step,length(randomind{k})))];
        end
        trainind = setdiff([1:length(ytrue)],testind);
        trainData = F(trainind,:);trainY = ytrue(trainind);valData=F(testind,:);valY = ytrue(testind);
        [a,b] = sort(D(testind,trainind),2);
        predY = trainY(b(:,1));
        err(t,i) = 1-sum(predY==valY)/length(valY);
    end
end
    classification_result(1) = mean(err(:));
    clear err
    
    % % KNNMdl = fitcknn(F,ytrue,...
    % %     'NumNeighbors',5,'Standardize',1);
    % % CVKNNMdl = crossval(KNNMdl,'Leaveout','on');
    % % %CVKNNMdl = crossval(KNNMdl);
    % % classError = kfoldLoss(CVKNNMdl);
    % % classification_result(2) = classError;
    
    nfold = 10;
    
    for i = 1:max(ytrue)
        ind = find(ytrue==i);
        randomind{i} = ind(randperm(length(ind)));
    end
    
    
    for i = 1:nfold
        testind = [];
        for k = 1:max(ytrue)
            step = round(length(randomind{k})/nfold);
            testind=[testind;randomind{k}((i-1)*step+1:min(i*step,length(randomind{k})))];
        end
        trainind = setdiff([1:length(ytrue)],testind);
        trainData = F(trainind,:);trainY = ytrue(trainind);valData=F(testind,:);valY = ytrue(testind);
        for t = 1:20
            [model.w,model.b] = vl_svmtrain(trainData', trainY, 20,'epsilon',1e-3,'Solver','SDCA');
            testscore = valData * model.w + model.b;
            %[~,~,info] = vl_pr(valY,);
            auc = colAUC(testscore,valY);
            err(i,t) = 1-mean(auc);
        end
    end
    classification_result(2) = mean(err(:));
    
    
    
    % linclass = fitcdiscr(F,ytrue,'Prior','empirical','discrimType', 'pseudoLinear');
    % %cvmodel = crossval(linclass,'Leaveout','on');
    %
    % cvmodel = crossval(linclass);
    % cverror = kfoldLoss(cvmodel);
    % classification_result(2) = cverror;
    
    
    linclass = fitcdiscr(F,ytrue,'Prior','empirical','discrimType', 'pseudoQuadratic');
    %cvmodel = crossval(linclass,'Leaveout','on');
    
    cvmodel = crossval(linclass);
    cverror = kfoldLoss(cvmodel);
    classification_result(3) = cverror;
    
    
end
end




