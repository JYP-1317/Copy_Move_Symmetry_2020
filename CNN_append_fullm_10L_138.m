%% Min Variation + New Descriptor(138 with 10LBP) 
%% SIFT normalization 확인
%% SIFT, Min variation Delete
clear all
close all

clc;
clf;

%tic

list = dir('D:\Dataset\CNN_Dataset\Coverage\Ori\'); 
list = list(3:size(list,1),:);
numfiles = size(list,1);

%%
for o = 1:numfiles %numfiles %Grip10번파일 40번부터 엑셀값다시 - g해결햇음
    try

        fn = strcat('D:\Dataset\CNN_Dataset\Coverage\Ori','\',list(o).name);
        [image, descriptors, locs] = sift(fn);

        fdes = descriptors(:,1:128);

        [cols, rows]=size(image);

        loc5 = round(locs(:,1:2));
        locs = locs(:,1:2);

        descriptors_n = descriptors;

        for i = 1:size(descriptors,1)
            d_n = descriptors(i,:) / sqrt(sum(descriptors(i,:).^2));
            descriptors_n(i, :) = d_n;
        end

        for m = 1:size(loc5,1)
           if (loc5(m,1) < 9) || (loc5(m,1) > cols-9) || (loc5(m,2) < 9) || (loc5(m,2) > rows-9) 
             loc5(m,:) = 9;
           end
        end    

    %% LBP

        LBP = LocalBinary(image);

        for k=1:size(loc5,1) %디스크립터 112로 바뀌면 수치 변경

            f = loc5(k,1)-7;  
            f1 = loc5(k,1)+7;
            g = loc5(k,2)-7;
            g1 = loc5(k,2)+7;

            h0=zeros(1,59);
            h1=zeros(1,10);

            for i = f:f1
                for j = g:g1
                    image_key = LBP(i,j);      
                    if image_key == 0
                        h0(1) = h0(1) + 1;
                    elseif image_key == 1
                        h0(2) = h0(2) + 1;
                    elseif image_key == 2
                        h0(3) = h0(3) + 1;
                    elseif image_key == 3
                        h0(4) = h0(4) + 1;
                    elseif image_key == 4
                        h0(5) = h0(5) + 1;
                    elseif image_key == 6
                        h0(6) = h0(6) + 1;
                    elseif image_key == 7
                        h0(7) = h0(7) + 1;
                    elseif image_key == 8
                        h0(8) = h0(8) + 1;
                    elseif image_key == 12
                        h0(9) = h0(9) + 1;
                    elseif image_key == 14
                        h0(10) = h0(10) + 1;
                    elseif image_key == 15
                        h0(11) = h0(11) + 1;
                    elseif image_key == 16
                        h0(12) = h0(12) + 1;    
                    elseif image_key == 24
                        h0(13) = h0(13) + 1;   
                    elseif image_key == 28
                        h0(14) = h0(14) + 1;  
                    elseif image_key == 30
                        h0(15) = h0(15) + 1;    
                    elseif image_key == 31
                        h0(16) = h0(16) + 1; 
                    elseif image_key == 32
                        h0(17) = h0(17) + 1; 
                    elseif image_key == 48
                        h0(18) = h0(18) + 1; 
                    elseif image_key == 56
                        h0(19) = h0(19) + 1;    
                    elseif image_key == 60
                        h0(20) = h0(20) + 1;   
                    elseif image_key == 62
                        h0(21) = h0(21) + 1; 
                    elseif image_key == 63
                        h0(22) = h0(22) + 1;  
                    elseif image_key == 64
                        h0(23) = h0(23) + 1;    
                    elseif image_key == 96
                        h0(24) = h0(24) + 1;   
                    elseif image_key == 112
                        h0(25) = h0(25) + 1; 
                    elseif image_key == 120
                        h0(26) = h0(26) + 1; 
                    elseif image_key == 124
                        h0(27) = h0(27) + 1; 
                    elseif image_key == 126
                        h0(28) = h0(28) + 1; 
                    elseif image_key == 127
                        h0(29) = h0(29) + 1; 
                    elseif image_key == 128
                        h0(30) = h0(30) + 1; 
                    elseif image_key == 129
                        h0(31) = h0(31) + 1;    
                    elseif image_key == 131
                        h0(32) = h0(32) + 1;  
                    elseif image_key == 135
                        h0(33) = h0(33) + 1;
                    elseif image_key == 143
                        h0(34) = h0(34) + 1;
                    elseif image_key == 13
                        h0(35) = h0(35) + 1;  
                    elseif image_key == 191
                        h0(36) = h0(36) + 1;
                    elseif image_key == 192
                        h0(37) = h0(37) + 1;   
                    elseif image_key == 193
                        h0(38) = h0(38) + 1;    
                    elseif image_key == 195
                        h0(39) = h0(39) + 1;    
                    elseif image_key == 199
                        h0(40) = h0(40) + 1;    
                    elseif image_key == 207
                        h0(41) = h0(41) + 1;       
                    elseif image_key == 223
                        h0(42) = h0(42) + 1;      
                    elseif image_key == 224
                        h0(43) = h0(43) + 1;        
                    elseif image_key == 225
                        h0(44) = h0(44) + 1;
                    elseif image_key == 227
                        h0(45) = h0(45) + 1;       
                    elseif image_key == 231
                        h0(46) = h0(46) + 1;    
                    elseif image_key == 239
                        h0(47) = h0(47) + 1;        
                    elseif image_key == 240
                        h0(48) = h0(48) + 1;
                    elseif image_key == 241
                        h0(49) = h0(49) + 1;     
                    elseif image_key == 243
                        h0(50) = h0(50) + 1; 
                    elseif image_key == 247
                        h0(51) = h0(51) + 1;        
                    elseif image_key == 248
                        h0(52) = h0(52) + 1;
                    elseif image_key == 249
                        h0(53) = h0(53) + 1;
                    elseif image_key == 251
                        h0(54) = h0(54) + 1;
                    elseif image_key == 252
                        h0(55) = h0(55) + 1;
                    elseif image_key == 253
                        h0(56) = h0(56) + 1;
                    elseif image_key == 254
                        h0(57) = h0(57) + 1;
                    elseif image_key == 255
                        h0(58) = h0(58) + 1;
                    else
                        h0(59) = h0(59) + 1;
                    end
                end
            end
            h1(1) = h0(1);
            h1(2) = h0(2)+h0(3)+h0(5)+h0(8)+h0(12)+h0(17)+h0(23)+h0(30);
            h1(3) = h0(4)+h0(6)+h0(9)+h0(13)+h0(18)+h0(24)+h0(31)+h0(37);
            h1(4) = h0(7)+h0(10)+h0(14)+h0(19)+h0(25)+h0(32)+h0(38)+h0(43);
            h1(5) = h0(11)+h0(15)+h0(20)+h0(26)+h0(33)+h0(39)+h0(44)+h0(48);
            h1(6) = h0(16)+h0(21)+h0(27)+h0(34)+h0(40)+h0(45)+h0(49)+h0(52);
            h1(7) = h0(22)+h0(28)+h0(35)+h0(41)+h0(46)+h0(50)+h0(53)+h0(55);
            h1(8) = h0(29)+h0(36)+h0(42)+h0(47)+h0(51)+h0(54)+h0(56)+h0(57);
            h1(9) = h0(58);
            h1(10) = h0(59);

            descriptors_n(k,129:138) = h1;
        end

        ddescriptors_n = descriptors_n(:,129:138);

        for i = 1:size(ddescriptors_n,1)
            d_n = ddescriptors_n(i,:) / sqrt(sum(ddescriptors_n(i,:).^2));
            descriptors_n(i, 129:138) = d_n;
        end


        if (image==0)
            p1=[];
            p2=[];
            tp=[];

        else
            p1=[];
            p2=[];
            num=0;

         % load data
            loc1 = locs(:,1:2);
            des1 = descriptors_n; %descriptors;
            dr2=0.65; %%나머지이미지 0.65 
            % descriptor are normalized with norm-2
            if (size(des1,1)<30000)
                des1 = des1./repmat(sqrt(diag(des1*des1')),1,size(des1,2));
            else
                des1_norm = des1; 
                for j= 1 : size(des1,2)
                    des1_j = des1_norm(j,:);
                    des1_norm(j,:) = des1_j/norm(des1_j,2); 
                end
                des1 = des1_norm;
            end


            % sift matching
            des2t = des1';   % precompute matrix transpose
            if size(des1,1) > 1 % start the matching procedure iff there are at least 2 points
                for i = 1 : size(des1,1)
                %논문사진 체크i=358,359,528;
                    dotprods = des1(i,:) * des2t;        % Computes vector of dot products
                    [vals,indx] = sort(acos(dotprods));  % Take inverse cosine and sort results


                    j=2;
                    while vals(j)<dr2* vals(j+1) 
                        j=j+1;
                    end
                    for k = 2 : j-1
                        match(i) = indx(k); 
                        if pdist([loc1(i,1) loc1(i,2); loc1(match(i),1) loc1(match(i),2)],'euclidean') > 30  
                            p1 = [p1 [loc1(i,2); loc1(i,1); 1]];
                            p2 = [p2 [loc1(match(i),2); loc1(match(i),1); 1]];
    %                           line([loc1(i,2) loc1(match(i),2)], ...
    %                                 [loc1(i,1) loc1(match(i),1)],'Color','r');

                            num=num+1;
                        end
                    end
                end
            end
       end

        assignin('base','Bp1',p1);
        assignin('base','Bp2',p2);

    %     cen_mat = [];
    %     
    %     for i = 1:size(Bp1,2)
    %         cen_mat(:,i) = (Bp1(:,i) + Bp2(:,i))/2;
    %     end
    % 
    %     avg_cen = sum(cen_mat,2)/size(cen_mat,2);
    % 
    %     T = 150; %%사이즈커지면 좀 더 큰값 CMH-300 MICC_150
    %     nBp1_temp = [];
    %     nBp2_temp = [];
    % 
    %     for i = 1:size(cen_mat,2)
    %         if  abs(cen_mat(:,i) - avg_cen(:,1)) < T
    %             cen_mat(:,i) = cen_mat(:,i);
    %         else
    %             cen_mat(:,i) = 0;
    %         end 
    %     end
    % 
    %     for j = 1:size(cen_mat,2)
    %         if  cen_mat(:,j) ~= 0
    %             nBp1_temp = [nBp1_temp; Bp1(:,j)'];
    %             nBp2_temp = [nBp2_temp; Bp2(:,j)'];
    %         end
    %     end 
    % 
    %     nBp1 = nBp1_temp';
    %     nBp2 = nBp2_temp';
    % 
    %     if size(nBp1,1) == 0 || size(nBp2,1) == 0
    %         nBp1 = Bp1;
    %         Bp1 = nBp1;
    %         nBp2 = Bp2;
    %         Bp2 = nBp2;
    %     else
    %         Bp1 = nBp1;
    %         Bp2 = nBp2;
    %     end

        p=[Bp1(1:2,:) Bp2(1:2,:)]';
          distance_p=pdist(p);
            Z = linkage(distance_p,'ward');
            %c = cluster(Z,'cutoff',2.2,'depth',4);
            c = cluster(Z,'maxclust',4);
            % show an image depicting clusters and matches
               %figure;
               %imshow(fn);
               hold on
                for i = 1: size(Bp1,2)
                    line([Bp1(1,i)' Bp2(1,i)'], [Bp1(2,i)' Bp2(2,i)'], 'Color', 'r');
                end
                gscatter(p(:,1),p(:,2),c)

        inliers1 = [];
        inliers2 = [];
        c_max = max(c);  %%c cluster 수 
            if(c_max > 1)
                n_combination_cluster = nchoosek(1:c_max,2);

                for i=1:1:size(n_combination_cluster,1)
                    k=n_combination_cluster(i,1);
                    j=n_combination_cluster(i,2);

                    z1=[];
                    z2=[];
                    %k=1; j=2;
                    for r=1:1:size(Bp1,2)
                        if c(r)==k && c(r+size(Bp1,2))==j
                            z1 = [z1; [p(r,:) 1]];
                            z2 = [z2; [p(r+size(Bp1,2),:) 1]];
                        end
                        if c(r)==j && c(r+size(Bp1,2))==k
                            z1 = [z1; [p(r+size(Bp1,2),:) 1]];
                            z2 = [z2; [p(r,:) 1]];
                        end
                    end 
                    %z1 are coordinates of points in the first cluster 
                    %z2 are coordinates of points in the second cluster            
                    if (size(z1,1) > 3 && size(z2,1) > 3)
                       % run ransacfithomography for affine homography
                       [H, inliers, dx, dy, xc, yc] = ransacfithomography2(z1', z2', 0.05); %%inliers number of cluster (클러스터번호)
                       %[H, inliers] = ransacfithomography(z1', z2', 0.05);

                        if size(H,1)==0
                            %num_gt = num_gt;
                        else
                            H = H / H(3,3);
                            %num_gt = num_gt+1;
                            inliers1 = [inliers1; [z1(inliers,1) z1(inliers,2)]];
                            inliers2 = [inliers2; [z2(inliers,1) z2(inliers,2)]]; %%clu 두개여서 노상관
                            %show_inliers(fn,z1',z2',inliers);
                        end
                    end
                end  
            end
    %end
    %% correlation 완판
    im = imread(fn);
    rim = rgb2gray(im);
    intensity = double(rim);

    %% x - Hx correlation 7by7
     for e = 4:1:size(intensity,2)-4 % e(열)
         for f = 4:1:size(intensity,1)-4 % f(행)

               %exy = [exy; [e,f]];
               %e=3; f=3;
                Hx = round(H(1,1)*e + H(1,2)*f+H(1,3)); % x좌표 이동 열의 값 
                Hy = round(H(2,1)*e + H(2,2)*f+H(2,3)); % y좌표 이동 행의 값
                %Hxy = [Hxy; [Hx, Hy]];

                if (Hy <= cols-4) && (Hy >= 4) && (Hx <= rows-4) && (Hx >= 4) %회전시 이미지 범위에 안맞으면 ㅠㅠ먼가이상해
                %597,797
                   is1(f,e) = sum(sum(intensity(f-3:f+3,e-3:e+3)))/49;
                   is(f,e) = sum(sum((intensity(f-3:f+3,e-3:e+3)-is1(f,e)).^2));

                   it1(Hy,Hx) = sum(sum(intensity(Hy-3:Hy+3,Hx-3:Hx+3))) / 49;
                   it(Hy,Hx) = sum(sum((intensity(Hy-3:Hy+3,Hx-3:Hx+3)-it1(Hy,Hx)).^2));

                   its(f,e) = sum(sum((intensity(f-3:f+3,e-3:e+3)-is1(f,e)).*(intensity(Hy-3:Hy+3,Hx-3:Hx+3)-it1(Hy,Hx))));

                   cor_cof(f,e) = its(f,e)/sqrt(it(Hy,Hx)*is(f,e));
               else
                   cor_cof(f,e) = 0;
               end      
         end
     end

    %% x - Hx correlation map postporcessing
    G = fspecial('gaussian',[7 7],0.5);
    im2 = imfilter(cor_cof,G,'same');
    bw = im2bw(im2,0.25); %0.4~0.6, 0.4 best
    %figure, imshow(bw);
    img_out = imfill(bw,'holes');
    img_out = bwmorph(img_out,'erode');
    img_out = bwmorph(img_out,'dilate');

    bound = bwboundaries(bw);
        inModel_x = [ Bp1(1,inliers) Bp2(1,inliers)]';
        inModel_y = [ Bp1(2,inliers) Bp2(2,inliers)]';

        img_out = false(size(im,1),size(im,2));    
        for k=1:size(bound,1)
          b= bound{k};
            in = inpolygon(inModel_x,inModel_y,b(:,2),b(:,1));
            if ~isempty(find(in,1))
                bw_b =false(size(im,1),size(im,2)); 
                bw_b = roipoly(bw_b,b(:,2),b(:,1));
                img_out = img_out | bw_b;
            end
        end
    img_out = bwmorph(img_out,'close');    
    img_out = bwmorph(img_out,'open');
    area=round(0.005*size(intensity,1)*size(intensity,2));
    img_out=bwareaopen(img_out,area,8);
    img_out = imfill(img_out,'holes');
    %figure, imshow(img_out);

    %% x - H^(-1)x correlation
    for g = 4:1:size(intensity,2)-4
            for h = 4:1:size(intensity,1)-4 
                A = H(1:2,1:2);
                inH = inv(A);
    %            H1 = inv(H);
    %            Hx1 = round(H1(1,1)*g + H1(1,2)*h+H1(1,3)); % x좌표 이동 열의 값 
    %            Hy1 = round(H1(2,1)*g + H1(2,2)*h+H1(2,3)); % y좌표 이동 행의 값
                Hx1 = round(inH(1,1)*(g-H(1,3)) + inH(1,2)*(h-H(2,3)));
                Hy1 = round(inH(2,1)*(g-H(1,3)) + inH(2,2)*(h-H(2,3)));

               if (Hy1 <= cols-4) && (Hy1 >= 4) && (Hx1 <= rows-4) && (Hx1 >= 4) 
               %497/797
               is1(h,g) = sum(sum(intensity(h-3:h+3,g-3:g+3))) / 49;
               is(h,g) = sum(sum((intensity(h-3:h+3,g-3:g+3)-is1(h,g)).^2));

               it1(Hy1,Hx1) = sum(sum(intensity(Hy1-3:Hy1+3,Hx1-3:Hx1+3))) / 49;
               it(Hy1,Hx1) = sum(sum((intensity(Hy1-3:Hy1+3,Hx1-3:Hx1+3)-it1(Hy1,Hx1)).^2));

               its(h,g) = sum(sum((intensity(h-3:h+3,g-3:g+3)-is1(h,g)).*(intensity(Hy1-3:Hy1+3,Hx1-3:Hx1+3)-it1(Hy1,Hx1))));

                  hcor_cof(h,g) = its(h,g)/sqrt(it(Hy1,Hx1)*is(h,g));
               else
                  hcor_cof(h,g) = 0;
               end 
            end
    end

    %% x - H^(-1)x correlation map postporcessing

    im3 = imfilter(hcor_cof,G,'same');
    bw1 = im2bw(hcor_cof,0.25);
    %figure, imshow(bw1); 
    himg_out = imfill(bw1,'holes');
    himg_out = bwmorph(himg_out,'erode');
    himg_out = bwmorph(himg_out,'dilate');
    bound = bwboundaries(bw1);
        inModel_x = [ Bp1(1,inliers) Bp2(1,inliers)]';
        inModel_y = [ Bp1(2,inliers) Bp2(2,inliers)]';

        himg_out = false(size(im,1),size(im,2));    
        for k=1:size(bound,1)
          b= bound{k};
            in = inpolygon(inModel_x,inModel_y,b(:,2),b(:,1));
            if ~isempty(find(in,1))
                bw_b1 =false(size(im,1),size(im,2)); 
                bw_b1 = roipoly(bw_b1,b(:,2),b(:,1));
                himg_out = himg_out | bw_b1;
            end
        end
    himg_out = bwmorph(himg_out,'close');    
    himg_out = bwmorph(himg_out,'open');     
    area1=round(0.005*size(intensity,1)*size(intensity,2));
    himg_out=bwareaopen(himg_out,area1,8);
    himg_out = imfill(himg_out,'holes');
    %figure, imshow(himg_out);

    final_image = zeros(size(image));
    for n = 3:1:size(image,2)-3
        for m= 3:1:size(image,1)-3
            if  (img_out(m,n) == 1) || (himg_out(m,n) == 1)
                final_image(m,n) = 1;
            end
        end
    end

    %figure, imshow(final_image);
    %imagename = strcat(int2str(o),'.png');
    [path, name, ext] = fileparts(list(o).name);
    name1 = strcat('D:\CNNpaper\IEEE Access\cov_cslbp','\',name,'.png');
    imwrite(final_image,name1)
    %name1 = strcat('C:\Users\JYP\Desktop\jpggt\','\',list(o).name);
    %imwrite(final_image,name1);
    catch 
        warning('not H, imageneame');
        H=0;
    end
end
