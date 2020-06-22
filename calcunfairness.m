function unfairness = calcunfairness(wg, pi_g, alpha0g, alpha1g, a0, a1)


nonzero0 = pi_g ~= 1;
nonzero1 = pi_g ~= 0;
unfairness = max(1-alpha0g(nonzero0)/a0, 1-(1-alpha0g(nonzero0))/(1-a0)) * ((1-pi_g(nonzero0)) .* wg(nonzero0)) +...
                    max(1-alpha1g(nonzero1)/a1, 1-(1-alpha1g(nonzero1))/(1-a1)) * (pi_g(nonzero1) .* wg(nonzero1));

end

