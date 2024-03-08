#############################
# Auxiliar functions
#############################
function sin2switchON(t, tf, ti)
    
    if t < ti
        return 0.0
    elseif t < tf
        return sin(0.5*pi*(t-ti)/(tf-ti))^2
    else 
        return 1.0
    end
    
end

function sin2switchOFF(t, tf, ti)
        
    if t < ti
        return 1.0
    elseif t < tf
        return 1.0 - sin(0.5*pi*(t-ti)/(tf-ti))^2
    else
        return 0.0
    end
            
end

function sin2pulse(t, tf, ti)
    
    tcenter = (tf + ti)/2
    
    if t < ti
        return 0.0
    elseif t < tcenter
        return sin(0.5*pi*(t-ti)/(tcenter-ti))^2
    elseif t < tf
        return 1.0 - sin(0.5*pi*(t-tcenter)/(tf-tcenter))^2
    else
        return 0.0
    end

end

function sin2window(t, t4, t3, t2, t1)
        
    if t < t1
        return 0.0
    elseif t < t2
        return sin(0.5*pi*(t-t1)/(t2-t1))^2
    elseif t < t3
        return 1.0
    elseif t < t4
        return 1.0 - sin(0.5*pi*(t-t3)/(t4-t3))^2
    else
        return 0.0
    end

end

function lorentzian(y, eta)
    return eta/(y^2 + eta^2)/pi
end

function gaussian(y, sigma)
    return exp(-(y/sigma)^2/2)/(sqrt(2*pi)*sigma)
end

#############################

abstract type AbstractLaserField end

# Monochromatic field
struct MonocromaticSinField <: AbstractLaserField
    amplitude::Float64
    omega::Float64
    tswitch::Float64
end

function (MF::MonocromaticSinField)(t)
    
    A = MF.amplitude
    tswitch = MF.tswitch
    omega = MF.omega
    
    return A*sin(omega*t)*sin2switchON(t, tswitch, 0.0)
    
end

# sin-square pulse
struct SinSquaredPulseField <: AbstractLaserField
    amplitude::Float64
    tcenter::Float64
    fwhm::Float64
end

function (SP::SinSquaredPulseField)(t)
    
    dt = pi*SP.fwhm/2/asin(1/sqrt(2))
    
    t1 = SP.tcenter - dt/2
    t2 = SP.tcenter + dt/2
    A = SP.amplitude
    
    return A*sin2pulse(t, t2, t1)*2/dt
    
end

# Lorentzian pulse
struct LorentzianPulseField <: AbstractLaserField
    amplitude::Float64
    tcenter::Float64
    fwhm::Float64
    duration::Float64
    tswitch::Float64
end

function (LP::LorentzianPulseField)(t)
    
    t1 = LP.tcenter - LP.duration/2
    t2 = LP.tcenter - LP.duration/2 + LP.tswitch
    t3 = LP.tcenter + LP.duration/2 - LP.tswitch
    t4 = LP.tcenter + LP.duration/2
    
    tcenter = LP.tcenter
    eta = LP.fwhm/2
    A = LP.amplitude
    
    return A*sin2window(t, t4, t3, t2, t1)*lorentzian(t - tcenter, eta)
    
end

# Gaussian pulse
struct GaussianPulseField <: AbstractLaserField
    amplitude::Float64
    tcenter::Float64
    fwhm::Float64
    duration::Float64
    tswitch::Float64
end

function (GP::GaussianPulseField)(t)
    
    t1 = GP.tcenter - GP.duration/2
    t2 = GP.tcenter - GP.duration/2 + GP.tswitch
    t3 = GP.tcenter + GP.duration/2 - GP.tswitch
    t4 = GP.tcenter + GP.duration/2
    
    tcenter = GP.tcenter
    sigma = GP.fwhm/2/sqrt(2*log(2))
    A = GP.amplitude
    
    return A*sin2window(t, t4, t3, t2, t1)*gaussian(t - tcenter, sigma)
    
end

# Step function
struct StepFunctionField <: AbstractLaserField
    amplitude::Float64
    tcenter::Float64
    tswitch::Float64
end

function (SF::StepFunctionField)(t)
    
    t1 = SF.tcenter - SF.tswitch/2
    t2 = SF.tcenter + SF.tswitch/2
    A = SF.amplitude
    
    return A*sin2switchON(t, t2, t1)
    
end