function [featuress1,featuress3] = ThreeLogGabor(pos, lf_num, rows, cols, volume)

angtheta_set = [0, pi/4, pi/2, pi*3/4, pi];  % -pi -- pi
angphi_set = [0, pi/4, pi/3, pi/2];  % -pi/2 -- pi/2  
norient = size(angtheta_set,2);
nscale = size(angphi_set, 2);
minWaveLength = 25;
sigmaOnf = 14.3/180*pi;
thetaSigma = 0.25;

filter = Newthree1_gaborconvolve(rows, cols, volume, nscale, minWaveLength, ...
                                sigmaOnf, thetaSigma, angtheta_set, angphi_set, norient);

j_idx_temp = 9: -1: 1;
j_idx = [];
for i = 1 : 9
    j_idx_temp = flip(j_idx_temp);
    j_idx = [j_idx; j_idx_temp];
end
                           
                            
for k = 1:lf_num
    k
    j_count = 0;
    for i = 1:9
        for j = j_idx(i,:)
            
            j_count = j_count+1;
            im = imread(strcat(pos, '\',num2str(k),'\',num2str(i),num2str(j),'.bmp'));%SAI
            ALL_LF(:,:,j_count) = double(rgb2gray(im));
            
        end
    end
    
    Y = double(ALL_LF);
    imagefft = fftn(Y);
    
    featuress1 = []; 
    featuress3 = []; 
    
    for numi  = 1 : nscale
        for numj = 1 : norient
            EO =ifftn(imagefft.* filter{numi,numj});            
            
            [structdis,~,~,~]= divisiveNormalization3D(real(EO));
            [alpha leftstd rightstd] = estimateaggdparam( structdis(:) );         
            featuress1 = [featuress1 alpha leftstd rightstd];
            
            [structdis,~,~,~]= divisiveNormalization3D(imag(EO));
            [alpha leftstd rightstd] = estimateaggdparam( structdis(:)  );
            featuress3 = [featuress3 alpha leftstd rightstd];
            
            
        end
    end

    featuress1(k,:) = [featuress1];
    featuress3(k,:) = [featuress3];

      
end

save 3D_LogGaborRe.mat featuress1
save 3D_LogGaborIm.mat featuress3
