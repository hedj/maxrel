module Constants
        const c = 2.997e8  # m/s
        const epsilon_0 = 8.85418782e-12 # m^3 kg^-1 s^4 A^2
        const mu_0 = 1.25663706e-6 # m kg s^-2 A^-2
        const invc = 1.0 / c
        const invcsquared = 1.0 / c^2
        const q_e = -1.60217657e-19 # C
        const m_e = 9.10938291e-31 # kg
        const zero_vec = [0.0, 0.0, 0.0]

        function DebyeLength(eV, n_per_cc, charge)
          return sqrt( 2.0 * epsilon_0 * eV * abs(q_e) / (3.0 * n_per_cc * 10^6 * charge^2) ) 
        end

        function ParticlesInADebyeSphere(eV, n_per_cc, charge)
           lambda = DebyeLength(eV, n_per_cc, charge) 
           return 4*pi*10^6*n_per_cc * lambda^3 / 3.0
        end

        export c,epsilon_0,mu_0,invc,invcsquared,q_e,m_e,zero_vec
end
