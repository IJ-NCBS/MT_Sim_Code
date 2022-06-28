

function [B_predicted_max, L_MT] = FUNBuckle_Bcric(Xa, Ya, Xb, Yb , xc, yc , xd, yd )

% Inputs
XaI = [ Xa; Ya ] ; % MTOC MT
XbI = [ Xb; Yb ] ; % EDGE MT
XcI = [ xc; yc ] ;  % Center 
XdI = [ xd; yd ] ;  % EDGE

%--- Rotate and translation ---
Theta_M = atand(((XcI(2,1)-XdI(2,1))./(XcI(1,1)-XdI(1,1))))   ; % rotation so that xc xd are in a line 
Rot_Mat = [cosd(-Theta_M),  -sind(-Theta_M)  ; sind(-Theta_M) , cosd(-Theta_M) ] ;
XaR = Rot_Mat*XaI  ;
XbR = Rot_Mat*XbI  ;
XcR = Rot_Mat*XcI  ;
XdR = Rot_Mat*XdI  ;

% translation so that line (XaR, XbR) and (XcR, XdR) intersect on origin
% intersection point
m1 = ((XbR(2,1)-XaR(2,1))./(XbR(1,1)-XaR(1,1)))  ;
m2 = ((XdR(2,1)-XcR(2,1))./(XdR(1,1)-XcR(1,1)))  ;

A = [-m1, 1 ; -m2, 1 ] ;
c1 = XaR(2,1)-m1*XaR(1,1) ;
c2 = XcR(2,1) -m2*XcR(1,1);
B = [ c1 ; c2 ] ;
Xo = linsolve(A,B) ;


% translation
XaT = XaR-Xo  ;
XbT = XbR-Xo  ;
XcT = XcR-Xo  ;
XdT = XdR-Xo  ;


m1_line = ((XbT(2,1)-XaT(2,1))./(XbT(1,1)-XaT(1,1)))  ; % atand(m1 ) 
%m2 = ((XdT(2,1)-XcT(2,1))./(XdT(1,1)-XcT(1,1)))  ;

A_c = norm(XcT(1,1)) ;
B_c = norm(XdT(1,1)) ;

B_predicted_max =  (abs((m1_line ./4).*(A_c + B_c + sqrt( (A_c + B_c).^2 -((( A_c - B_c ).^2)*8./pi.^2)  )   )) ) ;

N=20; I=0:N; GuessB = B_predicted_max  ; d = pdist2( (XaI'),  (XbI') ) ;
GAr =sqrt((ones(1,N).*(d./N)).^2  + (GuessB.*sin(  I(2:end).*(pi/N)) - GuessB.*sin( I(1:end-1).*(pi/N))).^2)  ;
L_MT =  (sum(GAr)) ;



