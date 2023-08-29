function [resp_est,X1,LN1,linear_out1,X2,LN2,linear_out2] = compute_est(a,stim)
%a: the list of the parameters, either 27 for one pathway or 54 for two
%stim: stimulus
L_mode = length(a); %number of parameters in a
%compute basis for filter
t=0.001:0.001:1;
t2 = 2*t-t.^2;
m1 = t2'*ones(1,15);
m2 = ones(1000,1)*(1:15);
A = sin(m1.*m2*pi);
basis = orth(A);


if(L_mode==27)
    %convolve with linear filter
    h = basis*a(1:15)';
    linear_out1 = conv(stim,h);
    linear_out1 = linear_out1(1:length(stim));
    %apply nonlinearity
    LN1 = (a(16).^(erf(linear_out1+a(17))+1))+a(18);%LN1 is output of L and N
    LN_out2 = a(19)*(LN1) +a(20); %scaled nonlinearity output
    %compute output of kinetics block
    X1 = simulate_4state(a(21:26),[0;0;0;100],LN1,LN_out2); %X1 is set of four states 
    resp_est = a(27)*(X1(2,:)); % X1(2,:) is active (second) state, the model output
    resp_est = resp_est-mean(resp_est(50001:end));
elseif(L_mode==54) %same for two pathways
    h1 = basis*a(1:15)';
    h2 = basis*a(28:42)';
    linear_out1 = conv(stim,h1);
    linear_out1 = linear_out1(1:length(stim));
    linear_out2 = conv(stim,h2);
    linear_out2 = linear_out2(1:length(stim));
    LN1 = (a(16).^(erf(linear_out1+a(17))+1))+a(18);
    LN2 = (a(43).^(erf(linear_out2+a(44))+1))+a(45);
    LN1_2 = a(19)*(LN1) +a(20);
    LN2_2 = a(46)*(LN2) +a(47);
    X1 = simulate_4state(a(21:26),[0;0;0;100],LN1,LN1_2);
    X2 = simulate_4state(a(48:53),[0;0;0;100],LN2,LN2_2);
    resp_est = a(27)*(X1(2,:))+ a(54)*(X2(2,:));
    resp_est = resp_est-mean(resp_est(50001:end));
    
end


end