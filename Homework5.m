clear all
clc

mass=288778.1   % mass          = 飞机质量, kg
s=510.95        % S             = 机翼面积, m^2
c_ref=8.32      % c_ref         = 参考气动弦长, m
b_ref=59.74     % b_ref         = 参考机翼展长, m
Ix=24679200     % Ix            = 滚转惯矩, kg.m^2
Iy=44883600     % Iy            = 俯仰惯矩, kg.m^2
Iz=67393200     % Iz            = 偏航惯矩, kg.m^2;
Ixz=1315320     % Ixz           = 交叉惯矩, kg.m^2

% x_ref         = 气动参考位置, 从参考气动弦前缘向后为正，按c_ref无量纲化
x_cg=0.25       % x_cg          = 实际中心位置, 从参考气动弦前缘向后为正，按c_ref无量纲化
% z_p           = 发动机推力线偏离质心距离, 在质心下方为正, 按c_ref无量纲化

% Ma            = 配平飞行时马赫数
% H             = 配平飞行高度，m
% rho           = 配平飞行时大气密度, kg/m^3
% c             = 配平飞行时音速, m/s
% rho_H         = 大气密度随高度变化率，kg/m^4
% c_H           = 音速随高度变化率，1/s
% alpha_star    = 配平迎角, rad
% gamma_star    = 配平爬升角, rad
% theta_star    = 配平俯仰角=alpha_star + gamma_star, rad

% %% 纵向气动导数
% CLstar,CDstar,Cmstar       = 配平升力、阻力、俯仰力矩系数
% CLalpha,CDalpha,Cmalpha    = 迎角导数, 1/rad
% CLq,Cmq                    = 俯仰角速度导数
% CLde,CDde,Cmde             = 升降舵操纵导数, 1/rad %升力/阻力，纵向力矩系数对elevator
% CLalphadot,Cmalphadot      = 迎角变化率导数
% CLMa,CDMa,CmMa             = 马赫数影响，是用在对高度大导数的计算公式中，是升力、阻力、力矩系数对马赫数的导数
% T_dp                       = 推力随油门变化(油门变化范围0-100), N/%
% T_V                        = 推力随速度变化, N/(m/s)
% T_H                        = 推力随高度变化, N/(m/s)
% %% 横航向气动导数
% Clbeta,Cnbeta,Ccbeta       = 侧滑导数, 1/rad
% Clp,Cnp,Ccp                = 滚转角速度导数
% Clr,Cnr,Ccr                = 偏航角速度导数
% Clda,Cnda,Ccda             = 副翼操纵导数, 1/rad
% Cldr,Cndr,Ccdr             = 方向舵操纵导数, 1/rad
D_V              =        rho * V_star   * S * (0.5 * CDMa * Ma + CDstar);




X_V           = (T_V     * cos(alpha_star) - D_V    ) / mass; 



A1            = [ X_V                          X_alpha+g                          0                      -g        X_H;
                  -Z_V                         -Z_alpha                           1                      0         -Z_H;
                  Mbar_V-Mbar_alphadot*Z_V     Mbar_alpha-Mbar_alphadot*Z_alpha   Mbar_q+Mbar_alphadot   0         Mbar_H-Mbar_alphadot*Z_H;
                  0                            0                                  1                      0         0;
                  0                            -V_star                             0                     V_star   0];
B1            = [ X_de/R2D                           X_dp*20;
                  -Z_de/R2D                          -Z_dp*20;
                  (Mbar_de-Mbar_alphadot*Z_de)/R2D   (Mbar_dp-Mbar_alphadot*Z_dp)*20;
                  0                                  0;
                  0                                  0];                   % 升降舵输入转换成“每度”,油门输入按每20%


%% 解释大气计算函数语法Compute atmospheric parameters and their rates with altitude.Mbar_dp
function [At,At2h] = Atmo2h(alt,j)
%   density(J=1),pressure(J=2),temprature(J=3),speed of sound(J=4)
%   Altitude range: 0 to 32000 meters
%   Units         : (Kg/m**3),(Paska),(degrees Kelvin) or (m/sec)

% rho0             = 1.225;                                                  % 海平面标准密度
% T0               = 288.15;                                                 % 海平面标准温度
% if H <= 11000                                                              % 计算状态点密度
%     T = T0 - 0.0065 * H;
%     rho = rho0 * (T / T0)^4.25588;
% elseif H > 11000 && H <= 20000
%     T = 216.65;
%     rho = 0.36392 * exp((-H + 11000) / 6341.62);
% else
%     T = 216.65 + 0.001 * (H - 20000);
%     rho = 0.088035 * (T / 216.65)^-35.1632;
% end
g = 9.81;
if alt < 0, h = 0;
elseif alt > 32000, h = 32000;
else h = alt;
end

if(h <= 11000.)
    At = 288.15 -.0065 * h; At2h = -.0065;
    switch j
        case 1
            At   = (1. - 2.25577 * 1.e-5 * h)^4.25588 * g * 0.12492;
            At2h = 4.25588 * (1. - 2.25577 * 1.e-5 * h)^3.25588 * g * 0.12492 * (-2.25577 * 1.e-5);
        case 2
            At   = 10332.3 * (1. - 2.25577 * 1.e-5*h)^5.2588*g;
            At2h = 5.2588 * 10332.3 * (1. - 2.25577 * 1.e-5 * h)^4.2588 * g * (-2.25577 * 1.e-5);
        case 4
            At   = sqrt(1.4 * 287. * At);
            At2h = 0.5 / At * 1.4 * 287. * (-.0065);
    end
elseif(h <= 20000.)
    At = 216.68; At2h = 0;
    switch j
        case 1
            At   = 0.037109 * exp((11000. - h) / 6341.62) * g;
            At2h = At / (-6341.62);
        case 2
            At   = 2307.8 * exp((11000. - h) / 6341.62) * g;
            At2h = At / (-6341.62);
        case 4
            At   = 295.069; At2h = 0;
    end
else
    At = 216.65 + 0.001 * (h - 20000.); At2h = 0.001;
    switch j
        case 1
            At   = 0.008977 * g / (1. + 4.61574 * 1.e-6 * (h - 20000.))^35.1632;
            At2h = (-35.1632) * 0.008977 * g /(1. + 4.61574 * 1.e-6 * (h - 20000.))^36.1632 * 4.61574 * 1.e-6;
        case 2
            At   = 558.28 / (1. + 4.61574 * 1.e-6 * (h - 20000.))^35.1632;
            At2h = (-35.1632) * 558.28 / (1. + 4.61574 * 1.e-6 * (h - 20000.))^36.1632 * 4.61574 * 1.e-6;
        case 4
            At   = sqrt(1.4 * 287. * At);
            At2h = 0.5 / At * 1.4 * 287. * 0.001;
    end
end
end
