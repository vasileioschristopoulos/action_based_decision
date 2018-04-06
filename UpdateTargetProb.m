function biasStim1Coeff = UpdateTargetProb(biasStim1Coeff,biasStim2Coeff,learn_rate,Seq)

Pcued          = biasStim1Coeff/(biasStim1Coeff + biasStim2Coeff);

for jj = 1:Seq
    Pcued = Pcued + learn_rate*(1-Pcued);
end

biasStim1Coeff    =  Pcued*biasStim2Coeff/(1-Pcued);