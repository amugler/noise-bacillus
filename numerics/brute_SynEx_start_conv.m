function [pnm,P,djs,C] = brute_SynEx_start_conv(gn,qn,rm,sbar,N,M,I,Pstart)

nvec = (0:N)';

% C = -[diagonal part of L]^{-1}*[nondiagonal part of L]
Cp1 = []; Cm1 = []; CpN = []; CmN = [];
for m = 0:M
  D = gn + rm(m+1)*nvec + qn + sbar*m;
  Cp1 = [Cp1; [rm(m+1)*(nvec(1:end-1)+1);0]./D];
  Cm1 = [Cm1; [0;gn(1:end-1)]./D];
  CpN = [CpN; sbar*(m+1)./D];
  CmN = [CmN; qn./D];
end
Cp1 = Cp1(1:end-1);
Cm1 = Cm1(1+1:end);
CpN = CpN(1:end-(N+1));
CmN = CmN(1+(N+1):end);

C = diag(Cp1,1)+diag(Cm1,-1)+diag(CpN,N+1)+diag(CmN,-(N+1));

P = Pstart;
P = P/sum(P);
for i = 1:I
%  i
  Pnew = C*P;
  Pnew = Pnew/sum(Pnew);
  djs(i) = JS_div(P,Pnew);
  P = Pnew;
end

pnm = reshape(P,[N+1 M+1]);