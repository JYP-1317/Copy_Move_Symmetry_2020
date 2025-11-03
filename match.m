function num = match(image1)

% Find SIFT keypoints for each image
[im1, descs, locs] = sift(image1);


if (im1==0)
    p1=[];
    p2=[];
    tp=[];
    
    
else
    p1=[];
    p2=[];
    num=0;

 % load data
    loc1 = locs(:,1:2);
    %scale1 = locs(:,3);
    %ori1 = locs(:,4);
    des1 = descs;

    % descriptor are normalized with norm-2
    if (size(des1,1)<15000)
        des1 = des1./repmat(sqrt(diag(des1*des1')),1,size(des1,2));
    else
        des1_norm = des1; 
        for j= 1 : size(des1,2)
            des1_j = des1_norm(j,:);
            des1_norm(j,:) = des1_j/norm(des1_j,2); 
        end
        des1 = des1_norm;
    end
    dr2=0.4; 
    
    % sift matching
    des2t = des1';   % precompute matrix transpose
    if size(des1,1) > 1 % start the matching procedure iff there are at least 2 points
        for i = 1 : size(des1,1)
            dotprods = des1(i,:) * des2t;        % Computes vector of dot products
            [vals,indx] = sort(acos(dotprods));  % Take inverse cosine and sort results

            j=2;
            while vals(j)<dr2* vals(j+1) 
                j=j+1;
            end
            for k = 2 : j-1
                match(i) = indx(k); 
                if pdist([loc1(i,1) loc1(i,2); loc1(match(i),1) loc1(match(i),2)],'euclidean') > 20  
                    p1 = [p1 [loc1(i,2); loc1(i,1); 1]];
                    p2 = [p2 [loc1(match(i),2); loc1(match(i),1); 1]];
%                     line([loc1(i,2) loc1(match(i),2)], ...
%                           [loc1(i,1) loc1(match(i),1)],'Color','r');
                   
                    num=num+1;
                end
            end
        end
    end
    
    assignin('base','Bp1',p1);
    assignin('base','Bp2',p2);  
    
    if size(p1,1)==0
        fprintf('Found %d matches.\n', num);
    else
        p=[p1(1:2,:)' p2(1:2,:)'];
        p=unique(p,'rows');
        p1=[p(:,1:2)'; ones(1,size(p,1))];
        p2=[p(:,3:4)'; ones(1,size(p,1))];
        num=size(p1,2);
        fprintf('Found %d matches.\n', num);
    end
    
    figure(1);
    imshow(image1)
    %imagesc(im1);
    %colormap('gray');
      hold on;
      for i = 1: size(des1,1)
         if (match(i) > 0)
             line([loc1(i,2) loc1(match(i),2)], ...
                  [loc1(i,1) loc1(match(i),1)],'Color','r');
         end
      end
  
      hold off;
end  

