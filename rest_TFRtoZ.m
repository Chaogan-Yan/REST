function [Z P Header] = rest_TFRtoZ(ImgFile,OutputName,Flag,Df1,Df2,Header)
% FORMAT [Z P Header] = y_TFRtoZ(ImgFile,OutputName,Flag,Df1,Df2,Header)
%   Input:
%     ImgFile    - T, F or R statistical image which wanted to be converted to Z statistical value
%     OutputName - The output name
%     Flag       - 'T', 'F' or 'R'. Indicate the type of the input statsical image
%     Df1        - the degree of freedom of the statistical image. For F 
%                statistical image, there is also Df2
%     Df2        - the second degree of freedom of F statistical image
%     Header     - If ImgFile is a MATRIX, the Nifti Header should be specified 
%   Output:
%     Z          - Z statistical image. Also output as .img/.hdr.
%     P          - The corresponding P value
%     Header     - the output Nifti Header
%___________________________________________________________________________
% Written by YAN Chao-Gan 100424.
% State Key Laboratory of Cognitive Neuroscience and Learning, Beijing Normal University, China, 100875
% ycg.yan@gmail.com
% Last revised by YAN Chao-Gan 100814. Changed to call spm_t2z for T and R images to use approximation in case of big T values.


if ischar(ImgFile)
    [Data VoxelSize Header]=rest_readfile(ImgFile);
else
    Data = ImgFile;
    VoxelSize=sqrt(sum(Header.mat(1:3,1:3).^2)); 
end

[nDim1,nDim2,nDim3]=size(Data);


if strcmpi(Flag,'F')
    fprintf('Convert F to Z...\n');

    Z = norminv(fcdf(Data,Df1,Df2)); %YAN Chao-Gan 100814. Use one-tail because F value is positive and one-tail.
    Z(Data==0) = 0;
    P = 1-fcdf(Data,Df1,Df2);
    
else % T image or R image: YAN Chao-Gan 100814. Changed to call spm_t2z to use approximation in case of big T values.
    
    if strcmpi(Flag,'R')
        fprintf('Convert R to T...\n');
        Data = Data .* sqrt(Df1 / (1 - Data.*Data));
    end
    
    fprintf('Convert T to Z...\n');
    
    P = 2*(1-tcdf(abs(Data),Df1)); %Two-tailed P value
    
    Tol = 1E-16; %minimum tail probability for direct computation
    Z = spm_t2z(Data,Df1,Tol);
    Z = reshape(Z,[nDim1,nDim2,nDim3]);
end


%YAN Chao-Gan 100814. Changed to call spm_t2z to use approximation in case of big T values.
% Z=zeros(nDim1,nDim2,nDim3);
% P=ones(nDim1,nDim2,nDim3);
% fprintf('\n\tT/F/R to Z Calculating...\n');
% for i=1:nDim1
%     fprintf('.');
%     for j=1:nDim2
%         for k=1:nDim3
%             if Data(i,j,k)~=0;
%                 if strcmp(Flag,'F')
%                     F=Data(i,j,k);
%                     PTemp =1-fcdf(abs(F),Df1,Df2);
%                     ZTemp=norminv(1 - PTemp/2)*sign(F);
%                     P(i,j,k)=PTemp;
%                     Z(i,j,k)=ZTemp;
%                 elseif strcmp(Flag,'T')
%                     T=Data(i,j,k);
%                     PTemp=2*(1-tcdf(abs(T),Df1));
%                     ZTemp=norminv(1 - PTemp/2)*sign(T);
%                     P(i,j,k)=PTemp;
%                     Z(i,j,k)=ZTemp;
%                 elseif strcmp(Flag,'R')
%                     R=Data(i,j,k);
%                     PTemp=2*(1-tcdf(abs(R).*sqrt((Df1)./(1-R.*R)),Df1));
%                     ZTemp=norminv(1 - PTemp/2)*sign(R);
%                     P(i,j,k)=PTemp;
%                     Z(i,j,k)=ZTemp;
%                 end
%                 
%             end
% 
%         end
%     end
% end

Z(isnan(Z))=0;
P(isnan(P))=1;


Header.descrip=sprintf('{Z_[%.1f]}',1);
if ~strcmpi(OutputName,'DO NOT OUTPUT IMAGE')%Added by Sandy to make it compatible with Image matrix
    rest_writefile(Z,OutputName,[nDim1,nDim2,nDim3],VoxelSize, Header,'double');
end

fprintf('\n\tT/F/R to Z Calculation finished.\n');