function output = LBP(Image)
    if size(Image,3)==3
        GrayImage = rgb2gray(Image);
    else
        GrayImage = Image;
    end
    
    GrayImage = cat(1, GrayImage(1, :), GrayImage, GrayImage(end, :));
    GrayImage = cat(2, GrayImage(:, 1), GrayImage, GrayImage(:, end));

    ColImg = im2col(GrayImage, [3 3] ,'sliding');
    ColImg5 = ColImg(5,:);
    ColImg5 = repmat(ColImg5,[8,1]);
    ColImg(5,:)=[];

    temp1 = zeros(size(ColImg));
    temp2 = ones(size(ColImg));
    temp1(ColImg>=ColImg5) = temp2(ColImg>=ColImg5);

    Mult = zeros(8,1,8);
    Mult(:,1,1) = [1;	2;   4;   128; 8;   64;  32;  16 ];
    Mult(:,1,2) = [128; 1;	 2;   64;  4;   32;  16;  8  ];
    Mult(:,1,3) = [64;	128; 1;   32;  2;   16;  8;   4  ];
    Mult(:,1,4) = [32;  64;  128; 16;  1;   8;   4;   2  ];
    Mult(:,1,5) = [16;  32;  64;  8;   128; 4;   2;   1  ];
    Mult(:,1,6) = [8;   16;  32;  4;   64;  2;   1;   128];
    Mult(:,1,7) = [4;   8;   16;  2;   32;  1;   128; 64 ];
    Mult(:,1,8) = [2;   4;   8;   1;   16;  128; 64;  32 ];
    
    Mult = repmat(Mult, [1, size(ColImg,2), 1]);
    
    tempLBP = zeros(1, size(temp1,2), 8);
    for i=1:8
        tempLBP(:,:,i) = sum(temp1.*Mult(:,:,i));
    end
    
    MaxLBP = max(tempLBP,[],3);    
    LBPImage = col2im(MaxLBP, [1 1], [size(GrayImage,1)-2, size(GrayImage,2)-2], 'distinct');

    output = LBPImage;
end