% problem 1b solution
tx = 3;
ty = 4;
x = linspace(0,2*tx,100);
% x = linespace(0,2*tx,50); size_x
y = linspace(0,2*ty,100);
% y = linespace(0,2*ty,50); size_y
[xx,yy] = ndgrid(x,y);

% Define rect and mod function
rect = @(t) single(abs(t)<1/2);
mod_function = @(x,T) mod(x+T/2,T)-T/2;
% Define g(x,y)
g_xy = 5*rect(sqrt(mod_function(xx,tx).^2+mod_function(yy,ty).^2) / 2);


clf, pl = @(i)subplot(221+i);
pl(0), imagesc(x,y,g_xy'), axis image, colorbar
xlabel x, ylabel y
title "g(x,y)"
xt = (0:1:2)*tx;
yt = (0:1:2)*ty;
xticks(xt), yticks(yt)

% Define c_kl
ckl_function = @(p)20/tx/ty * jinc(2*p);
for i=1:3
    level=4+(i-1)*10; % for levels of truncation K=0,5,15 and 25
    g=0;
    for k=-level:level
        for l=-level:level
            ckl = ckl_function(sqrt((k/tx)^2+(l/ty)^2));
            g = g+ckl * cos(2*pi*(k*xx/tx+l*yy/ty));
        end
    end

    pl(i), imagesc(x,y,g'), colorbar
    xlabel x, ylabel y
    xticks(xt), yticks(yt)
    title(sprintf("g(x,y) with level K=%d", level))
end

colormap gray

