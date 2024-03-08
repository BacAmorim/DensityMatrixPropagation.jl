function propagate(ph::PeierlsHamiltonian, rho0, tstep, tmax, A, alg = Tsit5())
    
    function eom!(drho, rho, p, t)
        
        ht = ham(ph, A(t))
        
        applyham!(drho, ht, rho, -im, 0)
        hermitianpart!(drho)
        
        return nothing
        
    end
    
    saved_values = SavedValues(Float64, Tuple{ComplexF64, ComplexF64})
    save_callback = SavingCallback((u, t, integrator) -> (tr(u), tr(u*current(ph, A(t))[1])), saved_values)
    
    
    prob = ODEProblem(eom!, rho0, 
        (0.0, tmax), dt = tstep, dtmax = tstep, tstops = [tstep*n for n in 0:(Int(fld(tmax, tstep))-1)],
        callback = save_callback, save_on = false, save_start = false, save_end = false)
    
    sol = solve(prob, alg)
    
    return sol, saved_values
    
end 
