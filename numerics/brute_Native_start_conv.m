function [pnm,P,djs,C] = brute_Native_start_conv(gn,qn,rnm,snm,N,M,I,Pstart)

nvec = (0:N)';

% C = -[diagonal part of L]^{-1}*[nondiagonal part of L]
Cp1 = []; Cm1 = []; CpN = []; CmN = [];
for m = 0:M
  D = gn + rnm(:,m+1).*nvec + qn + snm(:,m+1)*m;
  Cp1 = [Cp1; [rnm(2:end,m+1).*nvec(2:end);0]./D];
  Cm1 = [Cm1; [0;gn(1:end-1)]./D];
  if m ~= M; CpN = [CpN; snm(:,m+1+1)*(m+1)./D]; end
  CmN = [CmN; qn./D];
end
Cp1 = Cp1(1:end-1);
Cm1 = Cm1(1+1:end);
%CpN = CpN(1:end-(N+1));
CmN = CmN(1+(N+1):end);

C = diag(Cp1,1)+diag(Cm1,-1)+diag(CpN,N+1)+diag(CmN,-(N+1));

P = Pstart;
P = P/sum(P);
for i = 1:I
  %i
  Pnew = C*P;
  Pnew = Pnew/sum(Pnew);
  djs(i) = JS_div(P,Pnew);
  P = Pnew;
end

pnm = reshape(P,[N+1 M+1]);