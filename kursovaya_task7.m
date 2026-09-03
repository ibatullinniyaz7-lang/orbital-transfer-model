function kursovaya_task7()
% Курсовая работа. Задача 7. Межорбитальный перелет КЛА.
% Расчет траектории управляемого перелета между двумя круговыми орбитами.
% В файле используются методы, изученные в семестре:
% 1) ode45 - встроенный метод Рунге-Кутта 4-5 порядка;
% 2) rk4_method - собственная реализация метода Рунге-Кутта 4 порядка;
% 3) euler_method - собственная реализация метода Эйлера.

clc;
clear;
close all;

% ----------------------- Исходные данные -----------------------
% Так как в задании конкретные численные значения не заданы,
% используются безразмерные расчетные данные.
p.mu = 1.0;          % гравитационный параметр центрального тела
p.P  = 0.08;         % постоянная тяга ЭРД
p.q  = 0.0005;       % секундный расход массы
p.m0 = 1.0;          % начальная масса КЛА

p.r0 = 1.0;          % радиус начальной круговой орбиты
p.rk = 1.5;          % радиус конечной круговой орбиты
p.v0 = sqrt(p.mu/p.r0);  % скорость на начальной круговой орбите
p.vk = sqrt(p.mu/p.rk);  % скорость на конечной круговой орбите
p.theta0 = 0.0;      % начальный угол наклона скорости к горизонту
p.thetak = 0.0;      % конечный угол наклона скорости к горизонту
p.phi0 = 0.0;        % начальный полярный угол

x0 = [p.v0; p.theta0; p.r0; p.phi0];

% ----------------------- Управление ----------------------------
% Угол тяги аппроксимируется кубическим полиномом:
% alpha(t) = a0 + a1*tau + a2*tau^2 + a3*tau^3, tau = t/tk.
% Вектор u = [a0 a1 a2 a3 tk]. Начальное приближение уже близко к решению,
% а fminsearch уточняет его методом пристрелки.
u0 = [-0.62853419, -1.33301805, -0.14829783, 0.52414630, 6.09257065];

options = optimset('Display', 'iter', ...
                   'MaxIter', 500, ...
                   'MaxFunEvals', 2000, ...
                   'TolX', 1e-10, ...
                   'TolFun', 1e-12);

u = fminsearch(@(u) objective_function(u, p, x0), u0, options);

fprintf('\nНайденные параметры управления:\n');
fprintf('a0 = %.10f\n', u(1));
fprintf('a1 = %.10f\n', u(2));
fprintf('a2 = %.10f\n', u(3));
fprintf('a3 = %.10f\n', u(4));
fprintf('tk = %.10f\n', u(5));

% ----------------------- Основной расчет ode45 ------------------
tk = u(5);
ode_options = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
[t_ode, x_ode] = ode45(@(t, x) rhs_orbit(t, x, u, p), [0 tk], x0, ode_options);

% ----------------------- Ручные методы --------------------------
N = 300;
[t_euler, x_euler] = euler_method(@(t, x) rhs_orbit(t, x, u, p), [0 tk], x0, N);
[t_rk4, x_rk4] = rk4_method(@(t, x) rhs_orbit(t, x, u, p), [0 tk], x0, N);

% ----------------------- Итоговые значения ----------------------
xk_ode = x_ode(end, :);
xk_rk4 = x_rk4(end, :);
xk_euler = x_euler(end, :);

fprintf('\nЦелевые значения:\n');
fprintf('V_k = %.10f, theta_k = %.10f, r_k = %.10f\n', p.vk, p.thetak, p.rk);

fprintf('\nИтог ode45:\n');
fprintf('V = %.10f, theta = %.10f, r = %.10f, phi = %.10f\n', xk_ode(1), xk_ode(2), xk_ode(3), xk_ode(4));

fprintf('\nИтог RK4:\n');
fprintf('V = %.10f, theta = %.10f, r = %.10f, phi = %.10f\n', xk_rk4(1), xk_rk4(2), xk_rk4(3), xk_rk4(4));

fprintf('\nИтог Эйлера:\n');
fprintf('V = %.10f, theta = %.10f, r = %.10f, phi = %.10f\n', xk_euler(1), xk_euler(2), xk_euler(3), xk_euler(4));

% ----------------------- Таблица результатов --------------------
alpha_ode = alpha_control(t_ode, u);
mass_ode = p.m0 - p.q*t_ode;
result = [t_ode, x_ode(:,1), x_ode(:,2), x_ode(:,3), x_ode(:,4), alpha_ode, mass_ode];

fid = fopen('kursovaya_task7_results.csv', 'w');
fprintf(fid, 't,V,theta,r,phi,alpha,m\n');
fprintf(fid, '%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n', result');
fclose(fid);

% ----------------------- Графики --------------------------------
% 1. Траектория в декартовых координатах.
figure('Name', 'Траектория межорбитального перелета');
phi = x_ode(:,4);
r = x_ode(:,3);
X = r.*cos(phi);
Y = r.*sin(phi);
plot(X, Y, 'LineWidth', 1.5);
hold on;
ang = linspace(0, 2*pi, 400);
plot(p.r0*cos(ang), p.r0*sin(ang), '--');
plot(p.rk*cos(ang), p.rk*sin(ang), '--');
plot(X(1), Y(1), 'o', 'MarkerSize', 7, 'LineWidth', 1.5);
plot(X(end), Y(end), 's', 'MarkerSize', 7, 'LineWidth', 1.5);
axis equal;
grid on;
xlabel('x');
ylabel('y');
title('Траектория КЛА между круговыми орбитами');
legend('траектория КЛА', 'начальная орбита', 'конечная орбита', 'старт', 'финиш', 'Location', 'best');

% 2. Радиус от времени.
figure('Name', 'Изменение радиуса');
plot(t_ode, x_ode(:,3), 'LineWidth', 1.5);
hold on;
plot(t_rk4, x_rk4(:,3), '--', 'LineWidth', 1.0);
plot(t_euler, x_euler(:,3), ':', 'LineWidth', 1.0);
grid on;
xlabel('t');
ylabel('r(t)');
title('Изменение радиуса орбиты');
legend('ode45', 'Рунге-Кутта 4', 'Эйлер', 'Location', 'best');

% 3. Скорость от времени.
figure('Name', 'Изменение скорости');
plot(t_ode, x_ode(:,1), 'LineWidth', 1.5);
hold on;
plot([0 tk], [p.vk p.vk], '--');
grid on;
xlabel('t');
ylabel('V(t)');
title('Изменение скорости КЛА');
legend('V(t)', 'целевая скорость', 'Location', 'best');

% 4. Угол theta от времени.
figure('Name', 'Угол наклона скорости');
plot(t_ode, x_ode(:,2), 'LineWidth', 1.5);
hold on;
plot([0 tk], [p.thetak p.thetak], '--');
grid on;
xlabel('t');
ylabel('\theta(t), рад');
title('Изменение угла наклона скорости к местному горизонту');
legend('\theta(t)', 'целевое значение', 'Location', 'best');

% 5. Закон управления alpha(t).
figure('Name', 'Закон управления');
plot(t_ode, alpha_ode, 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('\alpha(t), рад');
title('Закон изменения направления тяги ЭРД');

% 6. Изменение массы.
figure('Name', 'Изменение массы');
plot(t_ode, mass_ode, 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('m(t)');
title('Изменение массы КЛА');

end

% =================================================================
function dx = rhs_orbit(t, x, u, p)
% Правая часть системы дифференциальных уравнений межорбитального перелета.
V = x(1);
theta = x(2);
r = x(3);
phi = x(4);

alpha = alpha_control(t, u);
m = p.m0 - p.q*t;

if m <= 0 || V <= 0 || r <= 0
    dx = [0; 0; 0; 0];
    return;
end

dV = p.P*cos(alpha)/m - p.mu*sin(theta)/(r^2);
dtheta = p.P*sin(alpha)/(m*V) + (V/r - p.mu/(r^2*V))*cos(theta);
dr = V*sin(theta);
dphi = V*cos(theta)/r;

dx = [dV; dtheta; dr; dphi];
end

% =================================================================
function alpha = alpha_control(t, u)
% Кубический закон управления направлением тяги.
tk = u(5);
tau = t./tk;
alpha = u(1) + u(2).*tau + u(3).*tau.^2 + u(4).*tau.^3;
end

% =================================================================
function J = objective_function(u, p, x0)
% Целевая функция метода пристрелки.
% Минимизируются отклонения V(tk), theta(tk), r(tk) от заданных значений.
tk = u(5);

if tk <= 0 || (p.m0 - p.q*tk) <= 0.1
    J = 1e20;
    return;
end

try
    ode_options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);
    [~, x] = ode45(@(t, x) rhs_orbit(t, x, u, p), [0 tk], x0, ode_options);
    xf = x(end, :);
catch
    J = 1e20;
    return;
end

if any(~isfinite(xf)) || xf(1) <= 0 || xf(3) <= 0
    J = 1e20;
    return;
end

eV = (xf(1) - p.vk)/p.vk;
etheta = xf(2) - p.thetak;
er = (xf(3) - p.rk)/p.rk;

% Большие коэффициенты заставляют метод точно выполнить граничные условия.
J = 1e6*eV^2 + 1e6*etheta^2 + 1e6*er^2 + 1e-4*sum(u(1:4).^2);
end

% =================================================================
function [t, x] = euler_method(f, interval, x0, N)
% Собственная реализация метода Эйлера для системы ДУ.
t0 = interval(1);
tk = interval(2);
h = (tk - t0)/N;
t = linspace(t0, tk, N + 1)';
x = zeros(N + 1, length(x0));
x(1, :) = x0';

for i = 1:N
    dx = f(t(i), x(i, :)');
    x(i + 1, :) = x(i, :) + h*dx';
end
end

% =================================================================
function [t, x] = rk4_method(f, interval, x0, N)
% Собственная реализация классического метода Рунге-Кутта 4 порядка.
t0 = interval(1);
tk = interval(2);
h = (tk - t0)/N;
t = linspace(t0, tk, N + 1)';
x = zeros(N + 1, length(x0));
x(1, :) = x0';

for i = 1:N
    xi = x(i, :)';
    ti = t(i);
    k1 = f(ti, xi);
    k2 = f(ti + h/2, xi + h*k1/2);
    k3 = f(ti + h/2, xi + h*k2/2);
    k4 = f(ti + h, xi + h*k3);
    x(i + 1, :) = (xi + h*(k1 + 2*k2 + 2*k3 + k4)/6)';
end
end
