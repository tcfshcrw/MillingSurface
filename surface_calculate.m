classdef surface_calculate < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        MillingcalculationUIFigure      matlab.ui.Figure
        KtEditFieldLabel                matlab.ui.control.Label
        KtEditField                     matlab.ui.control.NumericEditField
        CalculateButton                 matlab.ui.control.Button
        KrEditFieldLabel                matlab.ui.control.Label
        KrEditField                     matlab.ui.control.NumericEditField
        AxialDepthmmEditFieldLabel      matlab.ui.control.Label
        daEditField                     matlab.ui.control.NumericEditField
        RadialDepthmmEditFieldLabel     matlab.ui.control.Label
        drEditField                     matlab.ui.control.NumericEditField
        ToolFluteDiametermmLabel        matlab.ui.control.Label
        ToolREditField                  matlab.ui.control.NumericEditField
        ToolaxialimersionangleLabel     matlab.ui.control.Label
        ToolBetaEditField               matlab.ui.control.NumericEditField
        FeedmmminLabel                  matlab.ui.control.Label
        FeedEditField                   matlab.ui.control.NumericEditField
        RotationspeedRPMEditFieldLabel  matlab.ui.control.Label
        RPMEditField                    matlab.ui.control.NumericEditField
        UIAxesCuttingForce              matlab.ui.control.UIAxes
        UIAxesSGW                       matlab.ui.control.UIAxes
        UIAxesError                     matlab.ui.control.UIAxes
        ToolLengthmmEditFieldLabel      matlab.ui.control.Label
        ToolLengthEditField             matlab.ui.control.NumericEditField
        YoungsMoldulesofToolGPaEditFieldLabel  matlab.ui.control.Label
        YoungEditField                  matlab.ui.control.NumericEditField
        MillingTypeButtonGroup          matlab.ui.container.ButtonGroup
        DownMillingButton               matlab.ui.control.RadioButton
        UpMillingButton                 matlab.ui.control.RadioButton
        ToolFluteEditFieldLabel         matlab.ui.control.Label
        ToolFluteEditField              matlab.ui.control.NumericEditField
        ToolShankdiamatermmLabel        matlab.ui.control.Label
        ToolShankDiameter               matlab.ui.control.NumericEditField
        ToolFluteLengthmmEditFieldLabel  matlab.ui.control.Label
        ToolFluteLength                 matlab.ui.control.NumericEditField
        ProcessstandarddeviationmmEditFieldLabel  matlab.ui.control.Label
        ProcessStd                      matlab.ui.control.NumericEditField
        MateialKtKr                     matlab.ui.control.Label
        TextShow                        matlab.ui.control.Label
    end

    methods (Access = private)

        % Button pushed function: CalculateButton
        function CalculateButtonPushed(app, event)

da=app.daEditField.Value;
dr=app.drEditField.Value;
alpha_theta=app.ToolBetaEditField.Value/360*2*pi;
tool_r=app.ToolREditField.Value/2;
tool_l=app.ToolLengthEditField.Value;
tool_m=tool_r^4*pi/4;

tool_sr=app.ToolShankDiameter.Value/2;
tool_lf=app.ToolFluteLength.Value;
tool_ls=tool_l-tool_lf;
tool_ms=tool_sr^4*pi/4;

E=app.YoungEditField.Value*10^3;
feed=app.FeedEditField.Value;
N=app.ToolFluteEditField.Value;
speed_RPM=app.RPMEditField.Value;

tx=feed/speed_RPM/N;
kt=app.KtEditField.Value;
kr=app.KrEditField.Value;
process_std=app.ProcessStd.Value;
process_zone=3*process_std;

if app.DownMillingButton.Value==1
milling_geo=1;%1=down, 2=up
end
if app.UpMillingButton.Value==1
milling_geo=2;%1=down, 2=up
end
f_p=3600;%一周期點
   
        theta=linspace(0,2*pi,f_p);
        theta_space=linspace(0,540,round(f_p+f_p/2));
        beta_a=round(da*tan(alpha_theta)/tool_r/2/pi*f_p);
        beta_p=360/N;
        cut_window=zeros(1,f_p);
        if (milling_geo == 1)
            cut_window(round((f_p/2-round(acos((tool_r-dr)/tool_r)/2/pi*f_p)+1)):f_p/2)=1;
            theta_sgw_s=180;
        end
        if (milling_geo == 2)
            cut_window(1:round((acos((tool_r-dr)/tool_r)/pi*f_p/2+1)))=1;
            theta_sgw_s=0;
        end
        

        for ii=1:f_p
        fx_base(ii)=kt*tx*(sin(2*theta(ii))/2+kr*(1-cos(2*theta(ii)))/2);
        fy_base(ii)=kt*tx*(-1*kr*sin(2*theta(ii))/2+(1-cos(2*theta(ii)))/2);
        end
        fx_base_w(:)=cut_window.*fx_base;
        fy_base_w(:)=cut_window.*fy_base;
        
        cwd=zeros(1,f_p);
        cwd(1:beta_a)=tool_r/tan(alpha_theta);
        ts=zeros(1,f_p);
        ts(1)=1;
        if(N>1)
            for jj=1:N-1
            ts(round(1+f_p/N*jj))=1;
            end
        end
        %%%%%%milling force conv開始
        cwdc=conv(cwd,ts);
        fx_sim(:)=conv(cwdc,fx_base_w(:))*2*pi/f_p;
        fy_sim(:)=conv(cwdc,fy_base_w(:))*2*pi/f_p;
    %%%calculate the surface generation window    
    sgw=zeros(1,round(f_p+f_p/2));
    for jj=0:N-1
        sgw(round((theta_sgw_s+jj*beta_p)/180*f_p/2)+1:round((theta_sgw_s+jj*beta_p)/180*f_p/2)+beta_a)=1;
    end
    surface_force(:)=sgw(1:round(f_p*1.5)).*fy_sim(1:round(f_p+f_p/2));
    beta_theta_s=zeros(1,round(1.5*f_p));
    for jj=1:f_p
        beta_theta_s(round(theta_sgw_s/360*f_p)+jj)=theta(jj);
        for nn=1:N
        if(beta_theta_s(round(theta_sgw_s/360*f_p)+jj)>=2*pi/N)
            beta_theta_s(round(theta_sgw_s/360*f_p)+jj)=beta_theta_s(round(theta_sgw_s/360*f_p)+jj)-2*pi/N;
        end
        end
    end
    h_theta=zeros(1,round(1.5*f_p));
    if (milling_geo==2)
        h_theta(1:f_p)=beta_theta_s(1:f_p).*cwdc(1:f_p);
    end
    if(milling_geo==1)
        h_theta(round(theta_sgw_s/360*f_p+1):round(1.5*f_p))=beta_theta_s(round(theta_sgw_s/360*f_p+1):round(1.5*f_p)).*cwdc(1:f_p);
    end
    %K=-6*E*tool_m./(-1*h_theta.^3+3*tool_l^2.*h_theta-2*tool_l^3);
    
    z=h_theta;
    Ks=6*E*tool_ms./((z-tool_lf).^3-tool_ls^3+3*tool_ls^2.*(tool_l-z));
    Ksf=6*E*tool_m./(-1*(tool_lf-z).^3+3*(tool_lf-z).^2.*(tool_l-z));
    Kf=2*E*tool_ms./(2.*(tool_l-z)*tool_ls-tool_ls^2)./(tool_lf-z);
    %Dis=surface_force./K;
    Dis=surface_force./Ks+surface_force./Ksf+surface_force./Kf;
    
    cla(app.UIAxesCuttingForce);
    cla(app.UIAxesSGW);
    cla(app.UIAxesError);
    
    plot(app.UIAxesCuttingForce,theta_space,fx_sim(1:round(f_p+f_p/2)));
    hold(app.UIAxesCuttingForce);
    plot(app.UIAxesCuttingForce,theta_space,fy_sim(1:round(f_p+f_p/2)),'r');
    legend(app.UIAxesCuttingForce,'Fx','Fy');
    xlim(app.UIAxesCuttingForce,[0 540]);
    ylim(app.UIAxesCuttingForce,[1.1*min([min(fx_sim) min(fy_sim)]) 1.1*max([max(fx_sim) max(fy_sim)])]);
    grid(app.UIAxesCuttingForce);
    hold(app.UIAxesCuttingForce);
    %figure(2)
    plot(app.UIAxesSGW,theta_space,surface_force);
    hold(app.UIAxesSGW);
    xlim(app.UIAxesSGW,[0 540]);
    ylim(app.UIAxesSGW,[1.1*min(surface_force) 1.1*max(surface_force)]);
    grid(app.UIAxesSGW);
    %figure(3)
    hold(app.UIAxesError);
    plot(app.UIAxesError,h_theta(round((theta_sgw_s+beta_p)/180*f_p/2)+1:round((theta_sgw_s+beta_p)/180*f_p/2)+beta_a),Dis(round((theta_sgw_s+1*beta_p)/180*f_p/2)+1:round((theta_sgw_s+1*beta_p)/180*f_p/2)+beta_a));
    hold(app.UIAxesError);
    hold(app.UIAxesError);
    plot(app.UIAxesError,h_theta(round((theta_sgw_s+beta_p)/180*f_p/2)+1:round((theta_sgw_s+beta_p)/180*f_p/2)+beta_a),Dis(round((theta_sgw_s+1*beta_p)/180*f_p/2)+1:round((theta_sgw_s+1*beta_p)/180*f_p/2)+beta_a)+process_zone,'r--');
    hold(app.UIAxesError);
    hold(app.UIAxesError);
    plot(app.UIAxesError,h_theta(round((theta_sgw_s+beta_p)/180*f_p/2)+1:round((theta_sgw_s+beta_p)/180*f_p/2)+beta_a),Dis(round((theta_sgw_s+1*beta_p)/180*f_p/2)+1:round((theta_sgw_s+1*beta_p)/180*f_p/2)+beta_a)-process_zone,'r--');
    hold(app.UIAxesError);
    hold(app.UIAxesError);
    xlim(app.UIAxesError,[0 da]);
    hold(app.UIAxesError);
    grid(app.UIAxesError);
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create MillingcalculationUIFigure
            app.MillingcalculationUIFigure = uifigure;
            app.MillingcalculationUIFigure.Position = [100 100 1203 665];
            app.MillingcalculationUIFigure.Name = 'Milling calculation';

            % Create KtEditFieldLabel
            app.KtEditFieldLabel = uilabel(app.MillingcalculationUIFigure);
            app.KtEditFieldLabel.HorizontalAlignment = 'right';
            app.KtEditFieldLabel.Position = [153 613 25 22];
            app.KtEditFieldLabel.Text = 'Kt';

            % Create KtEditField
            app.KtEditField = uieditfield(app.MillingcalculationUIFigure, 'numeric');
            app.KtEditField.Position = [193 613 100 22];
            app.KtEditField.Value = 1659;

            % Create CalculateButton
            app.CalculateButton = uibutton(app.MillingcalculationUIFigure, 'push');
            app.CalculateButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateButtonPushed, true);
            app.CalculateButton.Position = [193 11 100 22];
            app.CalculateButton.Text = 'Calculate';

            % Create KrEditFieldLabel
            app.KrEditFieldLabel = uilabel(app.MillingcalculationUIFigure);
            app.KrEditFieldLabel.HorizontalAlignment = 'right';
            app.KrEditFieldLabel.Position = [153 576 25 22];
            app.KrEditFieldLabel.Text = 'Kr';

            % Create KrEditField
            app.KrEditField = uieditfield(app.MillingcalculationUIFigure, 'numeric');
            app.KrEditField.Position = [193 576 100 22];
            app.KrEditField.Value = 0.48;

            % Create AxialDepthmmEditFieldLabel
            app.AxialDepthmmEditFieldLabel = uilabel(app.MillingcalculationUIFigure);
            app.AxialDepthmmEditFieldLabel.HorizontalAlignment = 'right';
            app.AxialDepthmmEditFieldLabel.Position = [84 539 94 22];
            app.AxialDepthmmEditFieldLabel.Text = 'Axial Depth(mm)';

            % Create daEditField
            app.daEditField = uieditfield(app.MillingcalculationUIFigure, 'numeric');
            app.daEditField.Position = [193 539 100 22];
            app.daEditField.Value = 2.1;

            % Create RadialDepthmmEditFieldLabel
            app.RadialDepthmmEditFieldLabel = uilabel(app.MillingcalculationUIFigure);
            app.RadialDepthmmEditFieldLabel.HorizontalAlignment = 'right';
            app.RadialDepthmmEditFieldLabel.Position = [76 502 102 22];
            app.RadialDepthmmEditFieldLabel.Text = 'Radial Depth(mm)';

            % Create drEditField
            app.drEditField = uieditfield(app.MillingcalculationUIFigure, 'numeric');
            app.drEditField.Position = [193 502 100 22];
            app.drEditField.Value = 0.5;

            % Create ToolFluteDiametermmLabel
            app.ToolFluteDiametermmLabel = uilabel(app.MillingcalculationUIFigure);
            app.ToolFluteDiametermmLabel.HorizontalAlignment = 'right';
            app.ToolFluteDiametermmLabel.Position = [39 465 139 22];
            app.ToolFluteDiametermmLabel.Text = 'Tool Flute Diameter (mm)';

            % Create ToolREditField
            app.ToolREditField = uieditfield(app.MillingcalculationUIFigure, 'numeric');
            app.ToolREditField.Position = [193 465 100 22];
            app.ToolREditField.Value = 4;

            % Create ToolaxialimersionangleLabel
            app.ToolaxialimersionangleLabel = uilabel(app.MillingcalculationUIFigure);
            app.ToolaxialimersionangleLabel.HorizontalAlignment = 'right';
            app.ToolaxialimersionangleLabel.Position = [35 317 143 22];
            app.ToolaxialimersionangleLabel.Text = 'Tool  axial imersion angle ';

            % Create ToolBetaEditField
            app.ToolBetaEditField = uieditfield(app.MillingcalculationUIFigure, 'numeric');
            app.ToolBetaEditField.Position = [193 317 100 22];
            app.ToolBetaEditField.Value = 45;

            % Create FeedmmminLabel
            app.FeedmmminLabel = uilabel(app.MillingcalculationUIFigure);
            app.FeedmmminLabel.HorizontalAlignment = 'right';
            app.FeedmmminLabel.Position = [92 207 86 22];
            app.FeedmmminLabel.Text = 'Feed (mm/min)';

            % Create FeedEditField
            app.FeedEditField = uieditfield(app.MillingcalculationUIFigure, 'numeric');
            app.FeedEditField.Position = [193 207 100 22];
            app.FeedEditField.Value = 360;

            % Create RotationspeedRPMEditFieldLabel
            app.RotationspeedRPMEditFieldLabel = uilabel(app.MillingcalculationUIFigure);
            app.RotationspeedRPMEditFieldLabel.HorizontalAlignment = 'right';
            app.RotationspeedRPMEditFieldLabel.Position = [55 171 123 22];
            app.RotationspeedRPMEditFieldLabel.Text = 'Rotation speed (RPM)';

            % Create RPMEditField
            app.RPMEditField = uieditfield(app.MillingcalculationUIFigure, 'numeric');
            app.RPMEditField.Position = [193 171 100 22];
            app.RPMEditField.Value = 2600;

            % Create UIAxesCuttingForce
            app.UIAxesCuttingForce = uiaxes(app.MillingcalculationUIFigure);
            title(app.UIAxesCuttingForce, 'Cutting Force')
            xlabel(app.UIAxesCuttingForce, 'Tool Rotation Angle')
            ylabel(app.UIAxesCuttingForce, 'Cutting Force (N)')
            app.UIAxesCuttingForce.Position = [315 372 412 263];

            % Create UIAxesSGW
            app.UIAxesSGW = uiaxes(app.MillingcalculationUIFigure);
            title(app.UIAxesSGW, 'Surface Generation Force')
            xlabel(app.UIAxesSGW, 'Tool Rotation Angle')
            ylabel(app.UIAxesSGW, 'Force (N)')
            app.UIAxesSGW.Position = [770 372 412 263];

            % Create UIAxesError
            app.UIAxesError = uiaxes(app.MillingcalculationUIFigure);
            title(app.UIAxesError, 'Surface Dimension Error (mm)')
            xlabel(app.UIAxesError, 'Distance from bottom side (mm)')
            ylabel(app.UIAxesError, 'Error (mm)')
            app.UIAxesError.Position = [315 92 867 263];

            % Create ToolLengthmmEditFieldLabel
            app.ToolLengthmmEditFieldLabel = uilabel(app.MillingcalculationUIFigure);
            app.ToolLengthmmEditFieldLabel.HorizontalAlignment = 'right';
            app.ToolLengthmmEditFieldLabel.Position = [80 391 98 22];
            app.ToolLengthmmEditFieldLabel.Text = 'Tool Length (mm)';

            % Create ToolLengthEditField
            app.ToolLengthEditField = uieditfield(app.MillingcalculationUIFigure, 'numeric');
            app.ToolLengthEditField.Position = [193 391 100 22];
            app.ToolLengthEditField.Value = 32;

            % Create YoungsMoldulesofToolGPaEditFieldLabel
            app.YoungsMoldulesofToolGPaEditFieldLabel = uilabel(app.MillingcalculationUIFigure);
            app.YoungsMoldulesofToolGPaEditFieldLabel.HorizontalAlignment = 'right';
            app.YoungsMoldulesofToolGPaEditFieldLabel.Position = [5 243 173 22];
            app.YoungsMoldulesofToolGPaEditFieldLabel.Text = 'Young''s Moldules of Tool (GPa)';

            % Create YoungEditField
            app.YoungEditField = uieditfield(app.MillingcalculationUIFigure, 'numeric');
            app.YoungEditField.Position = [193 243 100 22];
            app.YoungEditField.Value = 212;

            % Create MillingTypeButtonGroup
            app.MillingTypeButtonGroup = uibuttongroup(app.MillingcalculationUIFigure);
            app.MillingTypeButtonGroup.TitlePosition = 'centertop';
            app.MillingTypeButtonGroup.Title = 'Milling Type';
            app.MillingTypeButtonGroup.Position = [170 47 123 74];

            % Create DownMillingButton
            app.DownMillingButton = uiradiobutton(app.MillingTypeButtonGroup);
            app.DownMillingButton.Text = 'Down Milling';
            app.DownMillingButton.Position = [11 28 91 22];
            app.DownMillingButton.Value = true;

            % Create UpMillingButton
            app.UpMillingButton = uiradiobutton(app.MillingTypeButtonGroup);
            app.UpMillingButton.Text = 'Up Milling';
            app.UpMillingButton.Position = [11 6 76 22];

            % Create ToolFluteEditFieldLabel
            app.ToolFluteEditFieldLabel = uilabel(app.MillingcalculationUIFigure);
            app.ToolFluteEditFieldLabel.HorizontalAlignment = 'right';
            app.ToolFluteEditFieldLabel.Position = [121 280 57 22];
            app.ToolFluteEditFieldLabel.Text = 'Tool Flute';

            % Create ToolFluteEditField
            app.ToolFluteEditField = uieditfield(app.MillingcalculationUIFigure, 'numeric');
            app.ToolFluteEditField.Position = [193 280 100 22];
            app.ToolFluteEditField.Value = 3;

            % Create ToolShankdiamatermmLabel
            app.ToolShankdiamatermmLabel = uilabel(app.MillingcalculationUIFigure);
            app.ToolShankdiamatermmLabel.HorizontalAlignment = 'right';
            app.ToolShankdiamatermmLabel.Position = [33 428 145 22];
            app.ToolShankdiamatermmLabel.Text = 'Tool Shank diamater (mm)';

            % Create ToolShankDiameter
            app.ToolShankDiameter = uieditfield(app.MillingcalculationUIFigure, 'numeric');
            app.ToolShankDiameter.Position = [193 428 100 22];
            app.ToolShankDiameter.Value = 6;

            % Create ToolFluteLengthmmEditFieldLabel
            app.ToolFluteLengthmmEditFieldLabel = uilabel(app.MillingcalculationUIFigure);
            app.ToolFluteLengthmmEditFieldLabel.HorizontalAlignment = 'right';
            app.ToolFluteLengthmmEditFieldLabel.Position = [50 354 128 22];
            app.ToolFluteLengthmmEditFieldLabel.Text = 'Tool Flute Length (mm)';

            % Create ToolFluteLength
            app.ToolFluteLength = uieditfield(app.MillingcalculationUIFigure, 'numeric');
            app.ToolFluteLength.Position = [193 354 100 22];
            app.ToolFluteLength.Value = 12;

            % Create ProcessstandarddeviationmmEditFieldLabel
            app.ProcessstandarddeviationmmEditFieldLabel = uilabel(app.MillingcalculationUIFigure);
            app.ProcessstandarddeviationmmEditFieldLabel.HorizontalAlignment = 'right';
            app.ProcessstandarddeviationmmEditFieldLabel.Position = [-6 135 182 22];
            app.ProcessstandarddeviationmmEditFieldLabel.Text = 'Process standard deviation (mm)';

            % Create ProcessStd
            app.ProcessStd = uieditfield(app.MillingcalculationUIFigure, 'numeric');
            app.ProcessStd.Position = [191 135 100 22];
            app.ProcessStd.Value = 0.001;

            % Create MateialKtKr
            app.MateialKtKr = uilabel(app.MillingcalculationUIFigure);
            app.MateialKtKr.Position = [13 560 99 98];
            app.MateialKtKr.Text = {'Al 6061 T6'; 'Kt=1659'; 'Kr=0.48'; '============='; 'Al  7075 T6'; 'Kt~=5000'; 'Kr~=0.4'; ''};

            % Create TextShow
            app.TextShow = uilabel(app.MillingcalculationUIFigure);
            app.TextShow.Position = [359 11 756 64];
            app.TextShow.Text = ' ';
        end
    end

    methods (Access = public)

        % Construct app
        function app = surface_calculate

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.MillingcalculationUIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.MillingcalculationUIFigure)
        end
    end
end