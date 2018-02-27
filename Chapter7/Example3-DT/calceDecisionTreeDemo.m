% Deomstrate decision tree classification
dim=2;
X1=[randn(100,dim)*0.5+2*ones(100,dim);randn(100,dim)-2*ones(100,dim)];
y1=zeros(200,1);
X2=randn(100,dim)*0.75;
y2=ones(100,1);

% Training data
X=[X1;X2];  % Predictor
y=[y1;y2];  % Response (label)

% Test data
Xtest=randn(10,dim);

% Plot original data using dots with different colors
dim_x=1;
dim_y=2;
plot(X(y==0,dim_x),X(y==0,dim_y),'b.')
hold on
plot(X(y==1,dim_x),X(y==1,dim_y),'r.')

% Fit the model and predict
model=DecisionTreeFit(X,y);
% Predict the label of test data using the fit model
yhat=DecisionTreePredict(Xtest,model);

% Plot the estimated class using circles with different colors
plot(Xtest(yhat==0,dim_x),Xtest(yhat==0,dim_y),'bo')
plot(Xtest(yhat==1,dim_x),Xtest(yhat==1,dim_y),'ro')
xlabel('x')
ylabel('y')
