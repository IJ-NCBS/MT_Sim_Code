
% Cell
Cell_Width  = 3.2 ; %
Cell_length_in = 14;% CELLSIZEVAR ; WT Cell_length = 14;
Vary_Cell_size = 1; %<-----------------------------------

% MT growth params
Vg_o_ini =  0.074 ;% mum/sec
Vs_o_ini =  -0.26 ;
Vg_f_ini =  0.048 ;
% Gamma dis para
% Cdc25
N_gamma_ini = 3.63  ; %  Ng=3.63, kg=0.029 from Mal3 data
k_gamma_ini = 0.029  ; %
%
Cat_time_stall_ini = 25;% CATTIMEVAR ; % in sec.; -1 == no force dependent cat. % 25 sec. from jason 2003

Num_MT_ini=10; %NUMARVAR % Num MT ;
vary_SPB_num=1; % < ---------------------------------
Binomial = -1 ; % if -1, MT are distributed equally; otherwise (1) from a random binomial

Ori_SD_ini= deg2rad(9.5) ; % WT Ori_STD 9.5 degree

% MT mech
f_mt=4 ; %
EI = 1.5 ;
restrictive_buckle = 1 ; % 1 on/else off



eta_cell = 0.9 ; % 0.9 pN s/µm2
%--- Nucleus constants ---
Radius_Nuc = 1.3 ;
Drag_nuc_tras=6*pi*eta_cell*Radius_Nuc   ;% --
Drag_nuc_rot =8*pi*eta_cell*Radius_Nuc.^3;

% distribution of MT on SPB
Even_MT_at_SPB=1 ; %=0 random distribution ; 1 even MTs at SPB and directed in opposite direction
Initial_state_offset = "Tip" ; % "Tip" or "Center"

Sim_Time_ini_tot = 7200 ; % sec
Sim_Time_ini = Sim_Time_ini_tot ; % sec
tot_iter = 100 ;
%file_name = sprintf('%s_Sim_Grand_Ar_Vary_NucRate_CellSz_fsR_XXXX_NumMT%d_CellLen%d_f_mt%1.1f_EI%1.1f.mat', Mean_MT_Array{N_dis,1}, Num_MT_ini , Cell_length_in , f_mt, EI)


%----------------
Cell_Ar = cell(tot_iter,1) ;
MT_Ar   = cell(tot_iter,1) ;
Length_Number   = cell(tot_iter,1) ;
for iter=1:tot_iter  %iter=1

    %Num_SPB = Num_SPB_ini ;
    N_MT = Num_MT_ini ;  %

    Delta_T = 1 ;
    Sim_Time = Sim_Time_ini ;
    Sim_TimeSer = 1:Delta_T:Sim_Time ;

    Vg_o =  Vg_o_ini;% mum/sec -- growth rate
    Vs_o = Vs_o_ini ; % mum/sec -- shrinkage rate
    Vg_f = Vg_f_ini ;
    Vs = Drag_nuc_tras; %Fs_p  ;%1.67 ; % pN force senstivity factor
    % kBT
    %KbT = 4.25.*10^-3 ; % pN mum    % 4.114 at 298

    Cat_time_stall = Cat_time_stall_ini ;
    N_gamma = N_gamma_ini ;
    k_gamma = k_gamma_ini ;
    Ori_SD = Ori_SD_ini ; % Ori_SD_ini = 0


    %--------- direction of MT
    if Binomial==-1
        Direction_State_AR = zeros(N_MT , 1) ;
        Direction_State_AR(1:2:N_MT,1) = 1 ;
        Direction_State_AR(2:2:N_MT,1) = 0 ;
    else
        Direction_State_AR = zeros(N_MT , 1) ;
        Direction_State_AR = randi([0,1] , N_MT,1)  ; % sum(Direction_State_AR)
    end

    %------------------------------
    % [time, Len_MT' , MT_Age',  Growth_State' , theta_dis', Direction, Cat_time_force, Dwell, Bend, Asso_SPB, SPB_Angle, L_E] 1 or 0
    Cell_Time_ser_MTs = cell(N_MT , 1) ;
    for j=1:N_MT
        Cell_Time_ser_MTs{j,1} = zeros(Sim_Time,12 ) ;
        Cell_Time_ser_MTs{j,1}(:,1) = 1:Sim_Time ;
    end



    %------------------- -------------------------------------------
    X_nuc = []; Y_nuc = [] ;
    if Vary_Cell_size==1
        Cell_length = Cell_length_in + 1.5*randn ;  % for WT
    else
        Cell_length = Cell_length_in ;
    end
    Bound_right = Cell_length./2 ;
    Bound_left = -Cell_length./2 ;
    Bound_up   = Cell_Width./2 ;
    Bound_down = -Cell_Width./2 ;
    if Initial_state_offset=="Tip"
        X_nuc = Bound_right-Radius_Nuc;
        Y_nuc = 0;
    elseif Initial_state_offset=="Center" % Initial_state_offset="Center"
        X_nuc = 0;
        Y_nuc = 0;
    end
    X_nuc_Ar = zeros(size(Sim_TimeSer,2),2) ;
    ki=1;
    X_nuc_Ar(1,:) = [X_nuc, Y_nuc ] ; % position and

    %-------- SPB angle and associated MT
    if vary_SPB_num==1
        Num_SPB = ceil(3.7 + 0.5*randn); % for WT
    else
        Num_SPB = floor(N_MT./4) ;
    end

    Temp_theta = zeros(Num_SPB,2) ;
    for ti=1:Num_SPB
        if mod(ti,2) == 0
            Rpi=0;
        else
            Rpi=pi;
        end
        Temp_theta(ti,1) = Rpi ;% + (rand*pi/4 )*sign(rand*2-1) ;  %<----- distribution
    end
    Temp_theta(:,2)=2;
    if Even_MT_at_SPB==0 % random MT distribution on SPB
        for ti=1:N_MT-2*Num_SPB
            t_ran = randi(Num_SPB);
            Temp_theta(t_ran,2) = Temp_theta(t_ran,2)+1 ;
        end
        Theta_AR = [] ;
        for ti=1:Num_SPB
            Theta_AR = [Theta_AR ;  [ti.*ones(Temp_theta(ti,2),1) , Temp_theta(ti,1).*ones(Temp_theta(ti,2),1)] ];
        end
    elseif Even_MT_at_SPB==1 % even MTs at SPB and directed in opposite direction
        for ti=1:(N_MT-2*Num_SPB)./2
            t_ran = randi(Num_SPB);
            Temp_theta(t_ran,2) = Temp_theta(t_ran,2)+2 ;
        end
        Theta_AR = [] ;
        for ti=1:Num_SPB
            Theta_AR = [Theta_AR ;  [ti.*ones(Temp_theta(ti,2),1) , Temp_theta(ti,1).*ones(Temp_theta(ti,2),1)] ];
        end
    else

    end

    %--- Nuc Angle
    T_nuc_Ar = zeros(size(Sim_TimeSer,2),1) ;
    T_Nuc = pi./2 ; % nucleus angle
    T_nuc_Ar(1,:) = T_Nuc ;



    %------ sim sttart
    tic
    for J=1:Sim_Time % J=J+1  % J=50

        for xn=1:length(Cell_Time_ser_MTs) % xn=19
            % [J.*Delta_T , Len_MT , MT_Age,  Growth_State , theta_dis, direction_state , force_cat_time, touch ] ;

            if J==1
                Cell_Time_ser_MTs{xn,1}(J,2) = Vg_o ; % length
                Cell_Time_ser_MTs{xn,1}(J,3) = 1 ;    % Age
                Cell_Time_ser_MTs{xn,1}(J,4) = 1 ;    % Growth state
                Cell_Time_ser_MTs{xn,1}(J,5) = Ori_SD*randn(1,1) ; % theta
                Cell_Time_ser_MTs{xn,1}(1:Sim_Time,6) = Direction_State_AR(xn, 1) ; % Direction state
                Cell_Time_ser_MTs{xn,1}(J,7) = 0 ; % force_on
                %Cell_Time_ser_MTs{xn,1}(J,8) = 0 ; % touching
                %Cell_Time_ser_MTs{xn,1}(J,9) = 0 ; % bending
                Cell_Time_ser_MTs{xn,1}(1:Sim_Time,10) = Theta_AR(xn,1) ; % Associated SPB
                Cell_Time_ser_MTs{xn,1}(1:Sim_Time,11) = Theta_AR(xn,2) ; % SPB Angle


            else % if J!=1
                if Cell_Time_ser_MTs{xn,1}(J-1,2)<=0
                    Cell_Time_ser_MTs{xn,1}(J,2) = Vg_o ; % length
                    Cell_Time_ser_MTs{xn,1}(J,3) = 1 ;    % Age
                    Cell_Time_ser_MTs{xn,1}(J,4) = 1 ;    % Growth state
                    Cell_Time_ser_MTs{xn,1}(J,5) = Ori_SD*randn(1,1) ; % theta
                    %Cell_Time_ser_MTs{xn,1}(J,6) = Direction_State_AR(xn, 1) ; % Direction state
                    Cell_Time_ser_MTs{xn,1}(J,7) = 0 ; % force_on
                else
                    if Cell_Time_ser_MTs{xn,1}(J-1,4)==1

                        Ft = Cell_Time_ser_MTs{xn,1}(J-1,7);
                        Vg_t = Vg_o.*(1-Ft./f_mt) ;
                        if Ft>=f_mt
                            Vg_t=0;
                        end

                        if Cell_Time_ser_MTs{xn,1}(J-1,8)==1
                            Cell_Time_ser_MTs{xn,1}(J,2) = Cell_Time_ser_MTs{xn,1}(J-1,2) + Vg_f ; % length with force
                        else
                            Cell_Time_ser_MTs{xn,1}(J,2) = Cell_Time_ser_MTs{xn,1}(J-1,2) + Vg_o ; % length
                        end
                        Cell_Time_ser_MTs{xn,1}(J,3) = Cell_Time_ser_MTs{xn,1}(J-1,3) + 1 ;    % Age

                        % catastrophe event %--- changes growth state
                        Age = Cell_Time_ser_MTs{xn,1}(J-1,3) ; % Age = 250
                        Cat_time_mean_o= 1./(( pdf('gam',Age, N_gamma, 1/k_gamma)  )./(1-cdf('gam',Age, N_gamma, 1/k_gamma )) );
                        if Cat_time_stall < 0
                            Cat_time_t=Cat_time_mean_o;
                        else
                            b_cat = (Cat_time_mean_o-Cat_time_stall)./Vg_o;
                            Cat_time_t = Cat_time_stall + b_cat.*Vg_t ;
                            if Cat_time_t>Cat_time_mean_o
                                Cat_time_t=Cat_time_mean_o ;
                            end
                        end

                        R_cat = 1/Cat_time_t ;
                        if rand<= 1- exp(-R_cat) % 1- exp(-R_cat.*DelT) DelT=1
                            Cell_Time_ser_MTs{xn,1}(J,4) = -1 ;  % Shrinkage state
                        else
                            Cell_Time_ser_MTs{xn,1}(J,4) = 1 ;  %  Growth state
                        end
                        Cell_Time_ser_MTs{xn,1}(J,5) = Cell_Time_ser_MTs{xn,1}(J-1,5) ; % theta
                        %Cell_Time_ser_MTs{xn,1}(J,6) = Direction_State_AR(xn, 1) ; % Direction state
                        %Cell_Time_ser_MTs{xn,1}(J,7) = 0 ; % force_on

                    elseif Cell_Time_ser_MTs{xn,1}(J-1,4)==-1
                        Cell_Time_ser_MTs{xn,1}(J,2) = Cell_Time_ser_MTs{xn,1}(J-1,2) + Vs_o ; % length
                        Cell_Time_ser_MTs{xn,1}(J,3) = Cell_Time_ser_MTs{xn,1}(J-1,3) + 1 ;    % Age
                        Cell_Time_ser_MTs{xn,1}(J,4) = Cell_Time_ser_MTs{xn,1}(J-1,4) ;  % Growth state

                        Cell_Time_ser_MTs{xn,1}(J,5) = Cell_Time_ser_MTs{xn,1}(J-1,5) ; % theta
                        %Cell_Time_ser_MTs{xn,1}(J,6) = Direction_State_AR(xn, 1) ; % Direction state
                        Cell_Time_ser_MTs{xn,1}(J,7) = 0 ; % force_on

                    else
                        "error"
                    end % is if growing or shrinking

                end  % if length >0 or <0
            end  % if J==1

        end % for each MT

        %

        

        %--------- force calculation % J=122
        D_s = Vs ; % 
        F_r_x=0 ; F_l_x=0;
        F_r_y=0 ; F_l_y=0;  % random force =zero
        T_r = [0 , 0, 0]; T_l = [0 , 0, 0];

        % [time, Len_MT' , MT_Age',  Growth_State' , theta_dis', Direction, Cat_time_force, Dwell, Bend, Asso_SPB, SPB_Angle] 1 or 0
        F_MTs = zeros(length(Cell_Time_ser_MTs),2) ;

        for I=1:length(Cell_Time_ser_MTs) % I=I+1  %I=40
            if Cell_Time_ser_MTs{I,1}(J,4)==1 % growth state
                if Cell_Time_ser_MTs{I,1}(J,6)==1  % direction

                    SPB_X = X_nuc + Radius_Nuc.*cos(T_Nuc+ Cell_Time_ser_MTs{I,1}(J,11) )  ;
                    SPB_Y = Y_nuc + Radius_Nuc.*sin( T_Nuc+ Cell_Time_ser_MTs{I,1}(J,11) )   ;

                    %theta wedged
                    Theta_range = [ atan((Bound_up-SPB_Y)./(Bound_right-(SPB_X))) , atan((Bound_down+SPB_Y)./(Bound_right-(SPB_X )))] ;
                    Theta_MT = Cell_Time_ser_MTs{I,1}(J,5) ;
                    if Theta_MT <= min(Theta_range)
                        Theta_MT = min(Theta_range)  ;
                    elseif Theta_MT >= max(Theta_range)
                        Theta_MT = min(Theta_range)  ;
                    end

                    if Cell_Time_ser_MTs{I,1}(J,2)>abs(Bound_right-(SPB_X))./cos(Theta_MT)

                        Fo = (1- Vg_f./Vg_o)*(f_mt)  ;
                        if Fo>f_mt
                            Fo=f_mt ;
                        end
                        %Fo = Cell_Time_ser_MTs{I,1}(J-1,7);
                        %------------>
                        if restrictive_buckle==1
                            Xa=SPB_X; Ya=SPB_Y ; L_MT = Cell_Time_ser_MTs{I,1}(J,2) ;
                            Xb=Bound_right ; Yb = (Bound_right-SPB_X)*tan(Theta_MT) + SPB_Y ;
                            if Theta_MT>=0
                                xd = Bound_right ; yd= Bound_up ;
                                xc = 0 ;           yc= Bound_up ;
                            else
                                xd = Bound_right ; yd= Bound_down ;
                                xc = 0 ;           yc= Bound_down ;
                            end
                            B_predicted      = FUNBuckle_ampli_approx(Xa, Ya, Xb, Yb, L_MT)   ; % w(x) = B*sin(x*pi/l)
                            [B_max , L_MT_co]= FUNBuckle_Bcric(Xa, Ya, Xb, Yb , xc, yc , xd, yd )  ;
                            if B_max>B_predicted
                                Fe = EI*pi^2./Cell_Time_ser_MTs{I,1}(J,2)^2 ;
                            else
                                Fe = EI*pi^2./L_MT_co^2 ;
                            end

                        else
                            Fe = EI*pi^2./Cell_Time_ser_MTs{I,1}(J,2)^2 ; %  Fe = EI*pi^.2./7 ;
                        end
                        if Fe<Fo
                            Vo=Fe ;
                            Cell_Time_ser_MTs{I,1}(J,9) = 1;

                            if Cell_Time_ser_MTs{I,1}(J,2)>L_MT_co
                                Cell_Time_ser_MTs{I,1}(J,12) = L_MT_co ;
                            else
                                Cell_Time_ser_MTs{I,1}(J,12) = Cell_Time_ser_MTs{I,1}(J,2) ;
                            end
                        else
                            Vo=Fo ;
                        end
                        F_r_x = F_r_x + Vo*cos(pi-Theta_MT);
                        F_r_y = F_r_y + Vo*sin(pi-Theta_MT);
                        F_MTs(I,1) = Vo; % Vo*cos(pi-Theta_MT);
                        T_r =  T_r + cross( [ SPB_X, SPB_Y, 0   ] , [Vo*cos(pi-Theta_MT), Vo*sin(pi-Theta_MT), 0 ] ) ;
                        Cell_Time_ser_MTs{I,1}(J,7) = Vo ;
                        Cell_Time_ser_MTs{I,1}(J,8)=1;

                    end
                elseif Cell_Time_ser_MTs{I,1}(J,6)==0 % direction

                    SPB_X = X_nuc + Radius_Nuc.*cos(T_Nuc+ Cell_Time_ser_MTs{I,1}(J,11) )  ;
                    SPB_Y = Y_nuc + Radius_Nuc.*sin( T_Nuc+ Cell_Time_ser_MTs{I,1}(J,11) )   ;

                    %theta wedged
                    Theta_range = [ atan((Bound_up-SPB_Y)./(abs(Bound_left)+(SPB_X))) , atan((Bound_down+SPB_Y)./(abs(Bound_left)+(SPB_X)))] ;
                    Theta_MT = Cell_Time_ser_MTs{I,1}(J,5) ;
                    if Theta_MT <= min(Theta_range)
                        Theta_MT = min(Theta_range)  ;
                    elseif Theta_MT >= max(Theta_range)
                        Theta_MT = min(Theta_range)  ;
                    end

                    if Cell_Time_ser_MTs{I,1}(J,2)>abs(Bound_left-(SPB_X))./cos(Theta_MT)

                        Fo = (1- Vg_f./Vg_o)*(f_mt)  ;
                        if Fo>f_mt
                            Fo=f_mt ;
                        end
                        %Fo = Cell_Time_ser_MTs{I,1}(J-1,7);
                        %------------>
                        if restrictive_buckle==1
                            Xa=SPB_X; Ya=SPB_Y ; L_MT = Cell_Time_ser_MTs{I,1}(J,2) ;
                            Xb=Bound_left ; Yb = (Bound_left-SPB_X)*tan(Theta_MT) + SPB_Y ;
                            if Theta_MT>=0
                                xd = Bound_left ; yd= Bound_up ;
                                xc = 0 ;           yc= Bound_up ;
                            else
                                xd = Bound_left ; yd= Bound_down ;
                                xc = 0 ;           yc= Bound_down ;
                            end
                            B_predicted      = FUNBuckle_ampli_approx(Xa, Ya, Xb, Yb, L_MT)   ; % w(x) = B*sin(x*pi/l)
                            [B_max , L_MT_co]= FUNBuckle_Bcric(Xa, Ya, Xb, Yb , xc, yc , xd, yd )  ;
                            if B_max > B_predicted
                                Fe = EI*pi^2./Cell_Time_ser_MTs{I,1}(J,2)^2 ;
                            else
                                Fe = EI*pi^2./L_MT_co^2 ;
                            end

                        else
                            Fe = EI*pi^2./Cell_Time_ser_MTs{I,1}(J,2)^2 ; %  Fe = EI*pi^.2./7 ;
                        end
                        if Fe<Fo
                            Vo=Fe ;
                            Cell_Time_ser_MTs{I,1}(J,9) = 1;
                            if Cell_Time_ser_MTs{I,1}(J,2)>L_MT_co
                                Cell_Time_ser_MTs{I,1}(J,12) = L_MT_co ;
                            else
                                Cell_Time_ser_MTs{I,1}(J,12) = Cell_Time_ser_MTs{I,1}(J,2) ;
                            end
                        else
                            Vo=Fo ;
                        end
                        F_l_x = F_l_x + Vo*cos(Theta_MT);
                        F_l_y = F_l_y + Vo*sin(Theta_MT);
                        F_MTs(I,2) = Vo; % Vo*cos(Theta_MT)
                        T_l =  T_l + cross( [ SPB_X, SPB_Y, 0   ] , [Vo*cos(Theta_MT), Vo*sin(Theta_MT), 0 ] ) ;
                        Cell_Time_ser_MTs{I,1}(J,7) = Vo ;
                        Cell_Time_ser_MTs{I,1}(J,8)=1;
                    end
                end
            end
        end
        % [time, Len_MT' , MT_Age',  Growth_State' , theta_dis', Direction, Cat_time_force, Dwell, Bend, Asso_SPB, SPB_Angle, ForceApp] 1 or 0

        % Translational
        % MT Drag
        Arr_MT = cell2mat(cellfun(@(x) x(J,[2,5,6,10]),Cell_Time_ser_MTs,'UniformOutput',false)); % length of MT, angle, associated SPB

        Max_len_Angle = zeros(max(Arr_MT(:,4))*2,2);
        Ij1=1;
        for Ij=1:max(Arr_MT(:,4))
            Max_len_Angle(Ij1,1) = max(Arr_MT( find((1:length(Arr_MT))'.*((Arr_MT(:,4)==Ij).*(Arr_MT(:,3)==1))), 1 ));
            Max_len_Angle(Ij1,2) = mean(Arr_MT( find((1:length(Arr_MT))'.*((Arr_MT(:,4)==Ij).*(Arr_MT(:,3)==1))), 2 ));
            Ij1=Ij1+1;
            Max_len_Angle(Ij1,1) = max(Arr_MT( find((1:length(Arr_MT))'.*((Arr_MT(:,4)==Ij).*(Arr_MT(:,3)==0))), 1 ));
            Max_len_Angle(Ij1,2) = mean(Arr_MT( find((1:length(Arr_MT))'.*((Arr_MT(:,4)==Ij).*(Arr_MT(:,3)==0))), 2 ));
            Ij1=Ij1+1;
        end

        R_MT = 0.025  ; % Radius of MT ~ R_bundle
        %eta_cell = 0.9 ;  % 0.9 pN s/µm2
        drag_MT_p=(2*pi*eta_cell*sum(Max_len_Angle(:,1).*cos(Max_len_Angle(:,2) )) )./(log(sum(Max_len_Angle(:,1).*cos(Max_len_Angle(:,2) )) ./(2*R_MT)) -0.2 ); % from Howard 2001
        drag_MT_l=(4*pi*eta_cell*sum(abs(Max_len_Angle(:,1).*sin(Max_len_Angle(:,2) ))) )./(log(sum(abs(Max_len_Angle(:,1).*sin(Max_len_Angle(:,2) ))) ./(2*R_MT)) +0.84 ); % from Howard 2001
        drag_MT_x = drag_MT_p+drag_MT_l ;

        drag_MT_p=(2*pi*eta_cell*sum(abs(Max_len_Angle(:,1).*sin(Max_len_Angle(:,2) ))) )./(log(sum(abs(Max_len_Angle(:,1).*sin(Max_len_Angle(:,2) ))) ./(2*R_MT)) -0.2 ); % from Howard 2001
        drag_MT_l=(4*pi*eta_cell*sum(abs(Max_len_Angle(:,1).*cos(Max_len_Angle(:,2) ))) )./(log(sum(abs(Max_len_Angle(:,1).*cos(Max_len_Angle(:,2) ))) ./(2*R_MT)) +0.84 ); % from Howard 2001
        drag_MT_y = drag_MT_p+drag_MT_l ;

        drag_Nuc=D_s;
        TotDrag_x = drag_Nuc*drag_MT_x./( drag_Nuc+drag_MT_x ) ;
        dx = (F_l_x+F_r_x)./(TotDrag_x ) ;
        X_nuc = X_nuc + dx  ;
        if X_nuc > Bound_right-Radius_Nuc
            X_nuc = Bound_right-Radius_Nuc ;
        end
        if X_nuc < Bound_left+Radius_Nuc
            X_nuc = Bound_left+Radius_Nuc ;
        end
        TotDrag_y = drag_Nuc*drag_MT_y./( drag_Nuc+drag_MT_y ) ;
        dy = (F_l_y+F_r_y)./(TotDrag_y) ;
        Y_nuc = Y_nuc + dy ;
        if Y_nuc > Bound_up-Radius_Nuc
            Y_nuc = Bound_up-Radius_Nuc ;
        end
        if Y_nuc < Bound_down+Radius_Nuc
            Y_nuc = Bound_down+Radius_Nuc ;
        end

        %Rotational
        drag_Rot=Drag_nuc_rot;
        TotDrag_Rot = drag_Rot;
        T_Nuc = T_Nuc + ( T_r(1,3) + T_l(1,3)  )./(TotDrag_Rot) ;
        if T_Nuc<0
            T_Nuc=0;
        elseif T_Nuc>pi
            T_Nuc=pi;
        end



        ki=ki+1;
        
        X_nuc_Ar(ki,:) = [X_nuc, Y_nuc ] ;
        T_nuc_Ar(ki,:) = T_Nuc ;

        %
    end % time array
    
    toc
    %figure, plot(X_nuc_Ar(:,2))

    X_nuc_Ar_T = [ [1:length(X_nuc_Ar)]'   , X_nuc_Ar(:,:) , T_nuc_Ar(:,1) ] ;
    X_nuc_Ar_T2 = X_nuc_Ar_T(1:5:length(X_nuc_Ar_T) , :) ;



    Cell_Ar{iter,1} = X_nuc_Ar_T2  ;
    clear X_nuc_Ar_T2 ;
    clear X_nuc_Ar_T ;
    clear X_nuc_Ar ;

    Cell_Time_ser_MTs_T = cell(0 , 1) ;
    for j=1:N_MT % j=50
        if sum(Cell_Time_ser_MTs{j,1}(:,2))>0
            Cell_Time_ser_MTs_T{j,1} = Cell_Time_ser_MTs{j,1}(1:5:length(Cell_Time_ser_MTs{j,1}) , :) ;
        end
    end



    MT_Ar{iter,1} = Cell_Time_ser_MTs_T ;

    clear Cell_Time_ser_MTs_T  ;
    clear Cell_Time_ser_MTs ;

    Length_Number{iter,1} = [Cell_length , Cell_Width, Num_SPB ] ;

end




