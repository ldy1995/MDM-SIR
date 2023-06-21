function tv_x = diff_x(x,sizeD)

tenX     = reshape(x, sizeD);
dfx1     = diff(tenX, 1, 1); % diff along the x
dfx      = zeros(sizeD);
dfx(1:end-1,:,:) = dfx1;
dfx(end,:,:)     =  tenX(1,:,:) - tenX(end,:,:); % the first diff of real output is the diff between first and end
tv_x=dfx(:);

