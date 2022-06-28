
% inputs: pinned points (X1, y1) and (x2, y2); and length of MTs (l)
% At critical load the beam equation read for pinned support
% w(x) = B*sin(x*pi./l)
% w(x) is the hight of the beam at position x
% l is the length of the beam

% given l and the distance between pinned point (d), 
% B has a unique value
% l is the arc length of the curve B*sin(x*pi./l) rewrite
% thus l = int( sqrt( 1 + B*pi*cos(x) ) , 0, pi) which is not integrable 
% Thus we guess B via minimizing {discrete sums of integral -l}

% x1= 0.2 ; y1=0.8 ;
% x2= 6.8 ; y2=-0.4;
% l=8.5;

function [GuessB] = FUNBuckle_ampli_approx(x1, y1, x2, y2, l)

d = pdist2( ([x1,y1]),  ([x2,y2]) ) ;

if l>d
    N=10;
    GuessB= sqrt((l./2).^2 - (d./2).^2) ;
    
    for k=1:20
        I=0:N;
        GAr =sqrt((ones(1,N).*(d./N)).^2  + (GuessB.*sin(  I(2:end).*(pi/N)) - GuessB.*sin( I(1:end-1).*(pi/N))).^2)  ;
        SumG = sum(GAr) ;
        DelErr = l - SumG ;
        
        GuessB = GuessB + DelErr./(2+rand) ;
    end
    
else
    GuessB=0;
end

end

