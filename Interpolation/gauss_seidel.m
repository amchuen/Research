function [m_old] = gauss_seidel(Abar, b, m_guess, error)

k = 0;
nvals = length(b);
m_old = m_guess;
resid = max(abs(Abar * m_guess - b));

while resid >= error
   for i = 1:nvals
      sig_new = 0;
      sig_old = 0;
      
      for j = 1:i
         sig_new = sig_new + A(i,j) * m_guess(j);
      end
       
      for j = (i+1):n
         sig_old =  sig_old + A(i,j) * m_guess(j);
      end
      
      m_guess(i) = (b(i) - sig_new - sig_old)/A(i,j);
   end
   
   % compute residual and correction for each iteration
   resid = max(abs((Abar * m_guess) - b));
   corr = max(abs(m_guess - m_old));
   
   % assign "new" old guess
   m_old = m_guess;
   
   % increase iteration count
   k = k+1;
end

end