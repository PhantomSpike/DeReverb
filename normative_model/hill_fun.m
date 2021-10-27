function out_data = hill_fun(data,n,c,theta)
%%out_data = hill_fun(data,n,c,theta)
%This function applies the hill transformation to the data of the form
%h(x) = c*x^n/(theta + c*x^n)
%data is a vector of the data to be processed
%n - the power of x
%c - normalizing constant to scale the data in the appropraite range
%theta - additive constant

out_data = (c*data.^n)./(theta + c*data.^n);
