%% debugging combined probability
data = load('comb_debug_list.txt');
test1 = data(:,1);
test2 = data(:,2);
test3 = data(:,3);
X = linspace(0,pi,100);
histogram(test1, 'facecolor', 'green', 'BinEdges', X);
hold on;
histogram(test2, 'facecolor', 'red', 'BinEdges', X);
histogram(test3, 'facecolor', 'blue', 'BinEdges', X);
hold off;

%% second attempt at debugging anisotropic strength
data = load('ani_debug_3row.txt');
test10 = data(:,1);
test100 = data(:,2);
test10000 = data(:,3);
X = linspace(0,pi,100);
const = 3.14/2.0 * ones(1000);
histogram(test10,'facecolor','green', 'BinEdges', X);
hold on;
histogram(test100, 'facecolor', 'yellow', 'BinEdges', X);
histogram(test10000,'facecolor', 'blue', 'BinEdges', X);
hold off;


%% first attempt at debugging anisotropic strength
test10 = load('ani_debug_strength_10.txt');
test100 = load('ani_debug_strength_100.txt');
test10000 = load('ani_debug_strength_10000.txt');
X = linspace(0,pi,1000);
% plot(X, 1000*exp(-(cos(X) .* cos(X) * 10))/(exp(-5) * pi * besseli(0,5)));
figure();
histogram(test10, 'facecolor', 'blue');
figure();
histogram(test100, 'facecolor', 'green');
figure();
histogram(test10000, 'facecolor', 'red');
