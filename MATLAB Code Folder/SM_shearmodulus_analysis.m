%% linear regression analysis of current shear modulus stuff
p0c = 3.813;
%x = p0c*ones(size(p0sm_r1_r2(:,1))) - p0sm_r1_r2(:,1);
x = p0sm(:,1);%p0sm_r1_r2(:,1);
y = 0.05 * p0sm(:,2);

[p,~] = polyfit(x,y,1);
best_fit_line = polyval(p,x);
figure
plot(x,y,'o');
hold on
plot(x, best_fit_line);
xlabel('$P_0 / \sqrt{A_0}$', 'Interpreter', 'latex');
ylabel('$r \langle G_{xy} \rangle$', 'Interpreter', 'latex');
title('Shear Modulus vs Shape Index')
hold off;
%% find projected critical shape index
projected_critical_shape_index = p(2) / p(1)*(-1);