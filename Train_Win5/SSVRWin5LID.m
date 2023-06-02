clear
clc


ReGabor1 = load ('Feature-Win5-LID\3D_LogGaborRe.mat');
ReGabor1 = ReGabor1.featuress1(:,1:60);
ReGabor2 = load ('Feature-Win5-LID\3D_LogGaborIm.mat');
ReGabor2 = ReGabor2.featuress3(:,1:60);
ReGabor5 = (ReGabor1 + ReGabor2)/2;

ReGabor3 = load ('Feature-Win5-LID\3D_LogGaborRe1.mat');
ReGabor3 = ReGabor3.featuress1(:,1:60);
ReGabor4 = load ('Feature-Win5-LID\3D_LogGaborIm1.mat');
ReGabor4 = ReGabor4.featuress3(:,1:60);
ReGabor6 = (ReGabor3 + ReGabor4)/2;
ReGabor = [ReGabor5; ReGabor6];





EPI_Gabor_width = load ('Feature-Win5-LID\EPI_Gabor_width.mat');
EPI_Gabor_width = EPI_Gabor_width.EPI_Gabor_width;
EPI_Gabor_width = permute(EPI_Gabor_width,[1,3,2]);
EPI_Gabor_height = load ('Feature-Win5-LID\EPI_Gabor_height.mat');
EPI_Gabor_height = EPI_Gabor_height.EPI_Gabor_height;
EPI_Gabor_height = permute(EPI_Gabor_height,[1,3,2]);
EPI_Gabor = cat(3, EPI_Gabor_width, EPI_Gabor_height);
EPI_Gabor = mean(EPI_Gabor,3);

EPI_Gabor_width = load ('Feature-Win5-LID\EPI_Gabor_width1.mat');
EPI_Gabor_width = EPI_Gabor_width.EPI_Gabor_width;
EPI_Gabor_width = permute(EPI_Gabor_width,[1,3,2]);
EPI_Gabor_height = load ('Feature-Win5-LID\EPI_Gabor_height1.mat');
EPI_Gabor_height = EPI_Gabor_height.EPI_Gabor_height;
EPI_Gabor_height = permute(EPI_Gabor_height,[1,3,2]);
EPI_Gabor1 = cat(3, EPI_Gabor_width, EPI_Gabor_height);
EPI_Gabor1 = mean(EPI_Gabor1,3);

EPI_Gabor = [EPI_Gabor; EPI_Gabor1];




svd_proportion = load( 'Feature-Win5-LID\svd_proportion.mat' );
svd_proportion = svd_proportion.svd_proportion;
Block_Size = 8;
Block_Size_Stride = 8;
M = 625;
N = 434;
M1 = size([Block_Size : Block_Size_Stride : M],2);
N1 = size([Block_Size : Block_Size_Stride : N],2);
fg1 = 6;
fg2 = 6;
svd_proportion = reshape (svd_proportion, [132, M1, N1]);
locfg1 = floor(linspace(1,M1,fg1+1));
locfg2 = floor(linspace(1,N1,fg2+1));
for i = 1 : 132
    for j = 1 : fg1
        for z = 1 : fg2
        Feature_SVD(i,(j-1)*fg2+z) = mean2(svd_proportion( i, locfg1(j):locfg1(j+1), locfg2(z):locfg2(z+1) ));
        end
    end
end


svd_proportion = load( 'Feature-Win5-LID\svd_proportion1.mat' );
svd_proportion = svd_proportion.svd_proportion;
Block_Size = 8;
Block_Size_Stride = 8;
M = 512;
N = 512;
M1 = size([Block_Size : Block_Size_Stride : M],2);
N1 = size([Block_Size : Block_Size_Stride : N],2);
fg1 = 6;
fg2 = 6;
svd_proportion = reshape (svd_proportion, [88, M1, N1]);
locfg1 = floor(linspace(1,M1,fg1+1));
locfg2 = floor(linspace(1,N1,fg2+1));
for i = 1 : 88
    for j = 1 : fg1
        for z = 1 : fg2
        Feature_SVD1(i,(j-1)*fg2+z) = mean2(svd_proportion( i, locfg1(j):locfg1(j+1), locfg2(z):locfg2(z+1) ));
        end
    end
end
Feature_SVD = [Feature_SVD; Feature_SVD1];




LBP_features_width = load ('Feature-Win5-LID\LBP_features_width.mat');
LBP_features_width = LBP_features_width.LBP_features_width;
LBP_features_width = reshape(LBP_features_width, 132, []);
LBP_features_height = load ('Feature-Win5-LID\LBP_features_height.mat');
LBP_features_height = LBP_features_height.LBP_features_height;
LBP_features_height = reshape(LBP_features_height, 132, []);
LBP_features1 = [LBP_features_width, LBP_features_height];

LBP_features_width = load ('Feature-Win5-LID\LBP_features_width1.mat');
LBP_features_width = LBP_features_width.LBP_features_width;
LBP_features_width = reshape(LBP_features_width, 88, []);
LBP_features_height = load ('Feature-Win5-LID\LBP_features_height1.mat');
LBP_features_height = LBP_features_height.LBP_features_height;
LBP_features_height = reshape(LBP_features_height, 88, []);
LBP_features2 = [LBP_features_width, LBP_features_height];

LBP_features = [LBP_features1; LBP_features2];


real_mos = load('Feature-Win5-LID\real_mos.mat');
real_mos = real_mos.real_mos;

synthetis_mos = load('Feature-Win5-LID\synthetis_mos.mat');
synthetis_mos = synthetis_mos.synthetis_mos;
jodmos = [real_mos(:,2); synthetis_mos(:,2)];


warning off


for z = 1 :1000
    z
    ind = randperm(length(jodmos));
    trainind = ind(1:floor(0.8*length(jodmos)));
    testind = ind(floor(0.8*length(jodmos))+1:length(jodmos));

    [coeff,score,latent] = pca(LBP_features(trainind,:),'Centered', false);
    Contribution_rate=cumsum(latent)./sum(latent);
    LBP_features_train = LBP_features(trainind,:)*coeff(:,1:40);
    LBP_features_teat = LBP_features(testind,:)*coeff(:,1:40); 
    
    features_train = [  Feature_SVD(trainind,:), ReGabor(trainind,:), EPI_Gabor(trainind,:), LBP_features_train   ];
    features_test = [  Feature_SVD(testind,:), ReGabor(testind,:), EPI_Gabor(testind,:), LBP_features_teat   ];
    


    [trainX, testX] = featNormalize(features_train, features_test);

    
    c_str = sprintf('%.0f',2^9);
    g_str = sprintf('%.8f',2^-8);  
    libsvm_options = ['-s 3 -t 2 -g ',g_str,' -c ',c_str];
    svm_model = svmtrain(jodmos(trainind), trainX, libsvm_options);
        
    

    [score, ~,~] = svmpredict(jodmos(testind), testX, svm_model);
    
    [ pearson_cc, spearmn_srocc, rmse, mae, dmosp ] = RegressionIQA( score, jodmos(testind) );
    PLCC(z)=pearson_cc;  SRCC(z)=spearmn_srocc;  RMSE(z)=rmse;
    
end
median(PLCC)
median(SRCC)
median(RMSE)

warning on


