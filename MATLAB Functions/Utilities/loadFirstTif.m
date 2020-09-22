function img = loadFirstTif(fname,bits, frame)

    if nargin < 3
        frame = 1;
    end
    
    switch bits
            case 8
                bitstr = 'uint8';
            case 16
                bitstr = 'uint16';
            case 32
                bitstr = 'single';
            otherwise
                bitstr = 'uint8';
                disp('May not be loading tiff properly');
    end
    
    infoImage=imfinfo(fname);
    mImage=infoImage(1).Width;
    nImage=infoImage(1).Height;
    img=zeros(nImage,mImage,1,bitstr);
    
    
    img(:,:,1) = imread(fname,frame);

end
