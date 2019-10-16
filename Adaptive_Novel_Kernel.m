%% Cite: "A Novel Adaptive Kernel for the RBF Neural Networks"
% Khan, S., Naseem, I., Togneri, R. et al. 
% Circuits Syst Signal Process (2017) 36: 1639. doi:10.1007/s00034-016-0375-7 
% https://link.springer.com/article/10.1007/s00034-016-0375-7

% Function approximation using "A Novel Adaptive Kernel for the RBF Neural Networks",
% with fixed centers, and fixed spread

clc;
clear all;
close all;

% Function  to be approximated 
% $f=e^{-(x^2+y)}$

%% Training Phase

% Generate training data
x=-1:.2:1; 
y=-1:.2:1;
f=[];
P=[];
for i=1:length(y)
    f=[f exp(-(x.^2+y(i)))];
    P=[P; x' repmat(y(i),length(x),1)];
end
[m n] = size(P);

%% Simulation parameters (RBF)
% Initialize mixing parameters
alpha1=0.5;
alpha2=0.5;

% Initialize training parameters
epoch = 100;% number of training epochs
eta = 1e-1; % Global learning rate
eta1 = eta; % Learning rate for Alpha 1
eta2 = eta; % Learning rate for Alpha 2
NC = 10;    % Number of centers (neurons)

%% Fixed centers of RBFNN using Kmean clustering
[~,c,beeta]=kmeans(P,NC);

%% Initialize weights and bias
n1 = length(c);
w = randn(1,n1);
b = randn();

%% Start
for k=1:epoch
    I(k)=0;
    ind=randperm(m); % Shuffle
    
    for i1=1:m
        for i2=1:n1
            % Cosine Kernel
            CD(i2)=(P(ind(i1),:)*c(i2,:)')/(norm(P(ind(i1),:))*norm(c(i2,:))+1e-50);
            % Euclidean Kernel
            ED(i2)=exp((-(norm(P(ind(i1),:)-c(i2,:))^2))/beeta(i2)^2);
            % Kernel mixing
            phi(i1,i2)=((abs(alpha1)*CD(i2))+(abs(alpha2)*ED(i2)))/sum([abs(alpha1) abs(alpha2)]);
        end
        
        % Calculate output and error
        y(i1)=w*phi(i1,:)' + b; % Output of RBF
        d(i1)=f(ind(i1));       % Desired output
        e=d(i1)-y(i1);          % Estimation error
        I(k)=I(k)+e*e';         % Objective Function

        
        % Gradient Descent learning algorithm
        % Weigth and Bias update 
        w=w+eta*e*phi(i1,:);
        b=b+eta*e;
        
        % Update mizing parameters
        R=(abs(alpha1)+abs(alpha2))^2;
        alpha1=alpha1+((eta1*e*w*(CD'-ED'))*(((alpha1)*(alpha2))/(R*alpha1)));
        alpha2=alpha2+((eta2*e*w*(ED'-CD'))*(((alpha1)*(alpha2))/(R*alpha2)));
    
    end
    
end

%% Performance measures
figure
semilogy(I) % MSE (curve)
grid minor
I(end)  % MSE
title('Cost Fuction');
xlabel('Epochs');
ylabel('MSE (dB)')

figure
plot(f(ind)); % desired value
hold on;
plot(y,'*-r'); % approximated value
legend('desired value','approximated value');


%% Test phase

% Generate Test samples
x=-.9:.2:.9;
y=-.9:.2:.9;
f=[];
P=[];

for i=1:length(y)
    f=[f exp(-(x.^2+y(i)))];
    P=[P; x' repmat(y(i),1,length(x))'];
end
[m n] = size(P);

%% Testing
for i1=1:m
    for i2=1:n1
        CD(i2)=(P(i1,:)*c(i2,:)')/(norm(P(i1,:))*norm(c(i2,:))+1e-50);
        ED(i2)=exp((-(norm(P(i1,:)-c(i2,:))^2))/beeta(i2)^2);
        phi(i1,i2)=((abs(alpha1(end))*CD(i2))+(abs(alpha2(end))*ED(i2)))/sum([abs(alpha1(end)) abs(alpha2(end))]);
    end
    y(i1)=w*phi(i1,:)' + b;
end

%% Performance measures
figure
plot(f);
hold on;
plot(y,'*-r');

% Surface
X=[-.9 -.7 -.5 -.3 -.1 .1 .3 .5 .7 .9]';
Y=[-.9 -.7 -.5 -.3 -.1 .1 .3 .5 .7 .9]';
for i=1:10
   Z(i,:)=f(10*i-9:10*i);
   Zout(i,:)=y(10*i-9:10*i);
end

figure
subplot(2,1,1)
surfc(X,Y,Z)
title('Desired Values')
subplot(2,1,2)
surfc(X,Y,Zout)
title('Approximated Values')

