"""
`module  JAC.PhotoIonization`
... a submodel of JAC that contains all methods for computing photoionization properties between some initial
    and final-state multiplets.
"""
module PhotoIonization
using Plots

using Printf, ..AngularMomentum, ..Basics, ..Continuum, ..Defaults, ..Radial, ..Nuclear, ..ManyElectron, ..PhotoEmission,
              ..TableStrings

"""
`struct  PhotoIonization.Settings  <:  AbstractProcessSettings`  ... defines a type for the details and parameters of computing photoionization lines.

    + multipoles                    ::Array{EmMultipole}  ... Specifies the multipoles of the radiation field that are to be included.
    + gauges                        ::Array{UseGauge}     ... Specifies the gauges to be included into the computations.
    + photonEnergies                ::Array{Float64,1}    ... List of photon energies [in user-selected units].
    + electronEnergies              ::Array{Float64,1}    ... List of electron energies; usually only one of these lists are utilized.
    + thetas                        ::Array{Float64,1}    ... List of theta-values if angle-differential CS are calculated explicitly.
    + phis                          ::Array{Float64,1}    ... List of phi]-values if angle-differential CS are calculated explicitly.
    + calcAnisotropy                ::Bool                ... True, if the beta anisotropy parameters are to be calculated and false otherwise (o/w).
    + calcPartialCs                 ::Bool                ... True, if partial cross sections are to be calculated and false otherwise.
    + calcTimeDelay                 ::Bool                ... True, if time-delays are to be calculated and false otherwise.
    + calcNonE1AngleDifferentialCS  ::Bool                ... True, if non-E1 angle-differential CS are be calculated and false otherwise.
    + calcTensors                   ::Bool                ... True, if statistical tensors of the excited atom are to be calculated and false o/w.
    + printBefore                   ::Bool                ... True, if all energies and lines are printed before their evaluation.
    + lineSelection                 ::LineSelection       ... Specifies the selected levels, if any.
    + stokes                        ::ExpStokes           ... Stokes parameters of the incident radiation.
    + freeElectronShift             ::Float64             ... An overall energy shift of all free-electron energies [user-specified units].
    + lValues                       ::Array{Int64,1}      ... Orbital angular momentum of free-electrons, for which partial waves are considered.
"""
struct Settings  <:  AbstractProcessSettings
    multipoles                      ::Array{EmMultipole}
    gauges                          ::Array{UseGauge}
    photonEnergies                  ::Array{Float64,1}
    electronEnergies                ::Array{Float64,1}
    thetas                          ::Array{Float64,1}
    phis                            ::Array{Float64,1}
    calcAnisotropy                  ::Bool
    calcPartialCs                   ::Bool
    calcTimeDelay                   ::Bool
    calcNonE1AngleDifferentialCS    ::Bool
    calcTensors                     ::Bool
    printBefore                     ::Bool
    lineSelection                   ::LineSelection
    stokes                          ::ExpStokes
    freeElectronShift               ::Float64
    lValues                         ::Array{Int64,1}
end


"""
`PhotoIonization.Settings()`  ... constructor for the default values of photoionization line computations
"""
function Settings()
    Settings(Basics.EmMultipole[E1], Basics.UseGauge[Basics.UseCoulomb, Basics.UseBabushkin], Float64[], Float64[], Float64[], Float64[],
                false, false, false, false, false, false, LineSelection(), Basics.ExpStokes(), 0., [0,1,2,3,4,5])
end


"""
`PhotoIonization.Settings(set::PhotoIonization.Settings;`

        multipoles=..,                      gauges=..,                  photonEnergies=..,          electronEnergies=..,
        thetas=..,                          calcAnisotropy=..,          calcPartialCs..,            calcTimeDelay=..,
        calcNonE1AngleDifferentialCS=..,    calcTensors=..,             printBefore=..,             lineSelection=..,
        stokes=..,                          freeElectronShift=..,       lValues=.. )

    ... constructor for modifying the given PhotoIonization.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::PhotoIonization.Settings;
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,                gauges::Union{Nothing,Array{UseGauge,1}}=nothing,
    photonEnergies::Union{Nothing,Array{Float64,1}}=nothing,                electronEnergies::Union{Nothing,Array{Float64,1}}=nothing,
    thetas::Union{Nothing,Array{Float64,1}}=nothing,                        phis::Union{Nothing,Array{Float64,1}}=nothing,
    calcAnisotropy::Union{Nothing,Bool}=nothing,                            calcPartialCs::Union{Nothing,Bool}=nothing,
    calcTimeDelay::Union{Nothing,Bool}=nothing,                             calcNonE1AngleDifferentialCS::Union{Nothing,Bool}=nothing,
    calcTensors::Union{Nothing,Bool}=nothing,                               printBefore::Union{Nothing,Bool}=nothing,
    lineSelection::Union{Nothing,LineSelection}=nothing,                    stokes::Union{Nothing,ExpStokes}=nothing,
    freeElectronShift::Union{Nothing,Float64}=nothing,                      lValues::Union{Nothing,Array{Int64,1}}=nothing)


    if  multipoles        == nothing   multipolesx        = set.multipoles        else  multipolesx        = multipoles         end
    if  gauges            == nothing   gaugesx            = set.gauges            else  gaugesx            = gauges             end
    if  photonEnergies    == nothing   photonEnergiesx    = set.photonEnergies    else  photonEnergiesx    = photonEnergies     end
    if  electronEnergies  == nothing   electronEnergiesx  = set.electronEnergies  else  electronEnergiesx  = electronEnergies   end
    if  thetas            == nothing   thetasx            = set.thetas            else  thetasx            = thetas             end
    if  phis              == nothing   phisx              = set.phis              else  phisx              = phis               end
    if  calcAnisotropy    == nothing   calcAnisotropyx    = set.calcAnisotropy    else  calcAnisotropyx    = calcAnisotropy     end
    if  calcPartialCs     == nothing   calcPartialCsx     = set.calcPartialCs     else  calcPartialCsx     = calcPartialCs      end
    if  calcTimeDelay     == nothing   calcTimeDelayx     = set.calcTimeDelay     else  calcTimeDelayx     = calcTimeDelay      end
    if  calcNonE1AngleDifferentialCS   == nothing   calcNonE1AngleDifferentialCSx = set.calcNonE1AngleDifferentialCS        else
        calcNonE1AngleDifferentialCSx  = calcNonE1AngleDifferentialCS                                                           end
    if  calcTensors       == nothing   calcTensorsx       = set.calcTensors       else  calcTensorsx       = calcTensors        end
    if  printBefore       == nothing   printBeforex       = set.printBefore       else  printBeforex       = printBefore        end
    if  lineSelection     == nothing   lineSelectionx     = set.lineSelection     else  lineSelectionx     = lineSelection      end
    if  stokes            == nothing   stokesx            = set.stokes            else  stokesx            = stokes             end
    if  freeElectronShift == nothing   freeElectronShiftx = set.freeElectronShift else  freeElectronShiftx = freeElectronShift  end
    if  lValues           == nothing   lValuesx           = set.lValues           else  lValuesx           = lValues            end

    Settings( multipolesx, gaugesx, photonEnergiesx, electronEnergiesx, thetasx, phisx, calcAnisotropyx, calcPartialCsx, calcTimeDelayx,
                calcNonE1AngleDifferentialCSx, calcTensorsx, printBeforex, lineSelectionx, stokesx, freeElectronShiftx, lValuesx)
end


# `Base.show(io::IO, settings::PhotoIonization.Settings)`  ... prepares a proper printout of the variable settings::PhotoIonization.Settings.
function Base.show(io::IO, settings::PhotoIonization.Settings)
    println(io, "multipoles:                    $(settings.multipoles)  ")
    println(io, "gauges:                        $(settings.gauges)  ")
    println(io, "photonEnergies:                $(settings.photonEnergies)  ")
    println(io, "electronEnergies:              $(settings.electronEnergies)  ")
    println(io, "thetas:                        $(settings.thetas)  ")
    println(io, "phis:                          $(settings.phis)  ")
    println(io, "calcAnisotropy:                $(settings.calcAnisotropy)  ")
    println(io, "calcPartialCs:                 $(settings.calcPartialCs)  ")
    println(io, "calcTimeDelay:                 $(settings.calcTimeDelay)  ")
    println(io, "calcNonE1AngleDifferentialCS:  $(settings.calcNonE1AngleDifferentialCS)  ")
    println(io, "calcTensors:                   $(settings.calcTensors)  ")
    println(io, "printBefore:                   $(settings.printBefore)  ")
    println(io, "lineSelection:                 $(settings.lineSelection)  ")
    println(io, "stokes:                        $(settings.stokes)  ")
    println(io, "freeElectronShift:             $(settings.freeElectronShift)  ")
    println(io, "lValues:                       $(settings.lValues)  ")
end


"""
`struct  PhotoIonization.PlasmaSettings  <:  Basics.AbstractLineShiftSettings`
    ... defines a type for the details and parameters of computing photoionization rates with plasma interactions.

    + multipoles             ::Array{Basics.EmMultipole}     ... Specifies the multipoles of the radiation field that are to be included.
    + gauges                 ::Array{Basics.UseGauge}        ... Specifies the gauges to be included into the computations.
    + photonEnergies         ::Array{Float64,1}              ... List of photon energies.
    + printBefore            ::Bool                          ... True, if all energies and lines are printed before their evaluation.
    + lineSelection          ::LineSelection                 ... Specifies the selected levels, if any.
"""
struct PlasmaSettings  <:  Basics.AbstractLineShiftSettings
    multipoles               ::Array{Basics.EmMultipole}
    gauges                   ::Array{Basics.UseGauge}
    photonEnergies           ::Array{Float64,1}
    printBefore              ::Bool
    lineSelection            ::LineSelection
end


"""
`PhotoIonization.PlasmaSettings()`  ... constructor for a standard instance of PhotoIonization.PlasmaSettings.
"""
function PlasmaSettings()
    PlasmaSettings([E1], [Basics.UseCoulomb], Float64[], true, LineSelection() )
end


# `Base.show(io::IO, settings::PhotoIonization.PlasmaSettings)`  ... prepares a proper printout of the settings::PhotoIonization.PlasmaSettings.
function Base.show(io::IO, settings::PhotoIonization.PlasmaSettings)
    println(io, "multipoles:              $(settings.multipoles)  ")
    println(io, "gauges:                  $(settings.gauges)  ")
    println(io, "photonEnergies:          $(settings.photonEnergies)  ")
    println(io, "printBefore:             $(settings.printBefore)  ")
    println(io, "lineSelection:           $(settings.lineSelection)  ")
end


"""
`struct  PhotoIonization.Channel`
    ... defines a type for a photoionization channel to help characterize a single multipole and scattering (continuum) state
        of many electron-states with a single free electron.

    + multipole      ::EmMultipole          ... Multipole of the photon absorption.
    + gauge          ::EmGauge              ... Gauge for dealing with the (coupled) radiation field.
    + kappa          ::Int64                ... partial-wave of the free electron
    + symmetry       ::LevelSymmetry        ... total angular momentum and parity of the scattering state
    + phase          ::Float64              ... phase of the partial wave
    + amplitude      ::Complex{Float64}     ... Photoionization amplitude associated with the given channel.
"""
struct  Channel
    multipole        ::EmMultipole
    gauge            ::EmGauge
    kappa            ::Int64
    symmetry         ::LevelSymmetry
    phase            ::Float64
    amplitude        ::Complex{Float64}
end


"""
`struct  Line`  ... defines a type for a photoionization line that may include the definition of channels.

    + initialLevel   ::Level                  ... initial-(state) level
    + finalLevel     ::Level                  ... final-(state) level
    + electronEnergy ::Float64                ... Energy of the (outgoing free) electron.
    + photonEnergy   ::Float64                ... Energy of the absorbed photon.
    + crossSection   ::EmProperty             ... Cross section for this photoionization.
    + angularBeta    ::EmProperty             ... beta -parameter for unpolarized targets with J=0, 1/2, 1
    + coherentDelay  ::EmProperty             ... coherent time-delay due to the selected averaging of phases.
    + incoherentDelay::EmProperty             ... incoherent time-delay due to the selected averaging of phases.
    + channels       ::Array{PhotoIonization.Channel,1}  ... List of PhotoIonization.Channels of this line.
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    electronEnergy   ::Float64
    photonEnergy     ::Float64
    crossSection     ::EmProperty
    angularBeta      ::EmProperty
    coherentDelay    ::EmProperty
    incoherentDelay  ::EmProperty
    channels         ::Array{PhotoIonization.Channel,1}
end


"""
`PhotoIonization.Line(initialLevel::Level, finalLevel::Level, crossSection::Float64)`
    ... constructor for an photoionization line between a specified initial and final level.
"""
function Line(initialLevel::Level, finalLevel::Level, crossSection::EmProperty)
    Line(initialLevel, finalLevel, totalRate, 0., 0., crossSection, EmProperty(0.), EmProperty(0.), EmProperty(0.),
            EmProperty(0.), PhotoChannel[] )
end


# `Base.show(io::IO, line::PhotoIonization.Line)`  ... prepares a proper printout of the variable line::PhotoIonization.Line.
function Base.show(io::IO, line::PhotoIonization.Line)
    println(io, "initialLevel:      $(line.initialLevel)  ")
    println(io, "finalLevel:        $(line.finalLevel)  ")
    println(io, "electronEnergy:    $(line.electronEnergy)  ")
    println(io, "photonEnergy:      $(line.photonEnergy)  ")
    println(io, "crossSection:      $(line.crossSection)  ")
    println(io, "angularBeta:       $(line.angularBeta)  ")
    println(io, "coherentDelay:     $(line.coherentDelay)  ")
    println(io, "incoherentDelay:   $(line.incoherentDelay)  ")
end


"""
`PhotoIonization.amplitude(kind::String, channel::PhotoIonization.Channel, omega::Float64, continuumLevel::Level,
                                initialLevel::Level, grid::Radial.Grid)`
    ... to compute the kind = (photoionization) amplitude  <(alpha_f J_f, epsilon kappa) J_t || O^(photoionization) || alpha_i J_i>
        due to the electron-photon interaction for the given final and initial level, the partial wave of the outgoing
        electron as well as the given multipole and gauge. A value::ComplexF64 is returned.
"""
function amplitude(kind::String, channel::PhotoIonization.Channel, omega::Float64, continuumLevel::Level, initialLevel::Level, grid::Radial.Grid)
    if      kind in [ "photoionization"]
    #-----------------------------------
        amp = PhotoEmission.amplitude("absorption", channel.multipole, channel.gauge, omega, continuumLevel, initialLevel, grid,
                                        display=false, printout=false)
        l         = Basics.subshell_l(Subshell(101, channel.kappa))
        amplitude = (1.0im)^(-l) * exp( -im*channel.phase ) * amp

    else    error("stop b")
    end

    return( amplitude )
end


"""
`PhotoIonization.angularFunctionK(L1::Int64, L2::Int64, X::Int64, Ji::AngularJ64, Jf::AngularJ64,
                                  kappa1::Int64, J1::AngularJ64, kappa2::Int64, J2::AngularJ64)`
    ... to compute angular function K(...) as defined for the non-E1 angle-differential cross sections by
        Nishita Hosea (2025). No tests are made that the triangular conditions of the quantum numbers
        are fulfilled. A wa::Float64 is returned.
"""
function angularFunctionK(L1::Int64, L2::Int64, X::Int64, Ji::AngularJ64, Jf::AngularJ64,
                          kappa1::Int64, J1::AngularJ64, kappa2::Int64, J2::AngularJ64)
    s1 = Subshell(20,kappa1);   j1 = Basics.subshell_j(s1)
    s2 = Subshell(20,kappa2);   j2 = Basics.subshell_j(s2)
    wb = (2L1+1) * (Basics.twice(j1)+1) * (Basics.twice(J1)+1) * (2L2+1) * (Basics.twice(j2)+1) * (Basics.twice(J2)+1)
    #wa = (2X+1) / (Basics.twice(Ji)+1)  *
    wa = (2X+1) * sqrt(wb) * AngularMomentum.phaseFactor([Ji, -1, Jf, 1, AngularJ64(1//2)])
    wa = wa * AngularMomentum.Wigner_6j(J2, J1, X, j1, j2, Jf) * AngularMomentum.Wigner_6j(J2, J1, X, L1, L2, Ji)

    return( wa )
end


"""
`PhotoIonization.angularFunctionW(theta::Float64, L1::Int64, L2::Int64, X::Int64, lambda1::Int64, lambda2::Int64,
                                  kappa1::Int64, mu1::Rational{Int64}, kappa2::Int64, mu2::Rational{Int64})`
    ... to compute angular function W(theta; ...) as defined for the non-E1 angle-differential cross sections by
        Nishita Hosea (2025). No tests are made that the triangular conditions of the quantum numbers
        are fulfilled. A  wa::Float64 is returned.
"""
function angularFunctionW(theta::Float64, L1::Int64, L2::Int64, X::Int64, lambda1::Int64, lambda2::Int64,
                          kappa1::Int64, mu1::Rational{Int64}, kappa2::Int64, mu2::Rational{Int64})
    if  abs(lambda1) > L1  ||   abs(lambda2) > L2  ||   abs(lambda2-lambda1) > X   return( 0.)    end
    s1 = Subshell(20,kappa1);   j1 = Basics.subshell_j(s1)
    s2 = Subshell(20,kappa2);   j2 = Basics.subshell_j(s2)
    wa = AngularMomentum.Wigner_dmatrix(X, lambda2-lambda1, Int64(mu2-mu1), theta) *
         AngularMomentum.Wigner_3j(j2, j1, X, abs(mu2), -abs(mu1), mu2-mu1) *
         AngularMomentum.Wigner_3j(L2, L1, X, -lambda2, lambda1, lambda2-lambda1)

    return( wa )
end

"""
`PhotoIonization.computeDisplayNonE1AngleDifferentialCS(stream::IO, lines::Array{PhotoIonization.Line,1},
                                                        settings::PhotoIonization.Settings)`
    ... to compute & display the non-E1 angle-differential photoionization cross sections for all PhotoIonization.Line's
        and at all angles theta as defined in the settings. The general formula by Nishita Hosea (2025) is applied here.
        A neat table is printed for each line but nothing is returned otherwise.
"""
function computeDisplayNonE1AngleDifferentialCS(stream::IO, lines::Array{PhotoIonization.Line,1}, settings::PhotoIonization.Settings)
    function spinDensityMatrix(lambda1::Int64, lambda2::Int64, Stokes::ExpStokes)
        # Convert the Stokes parameters of the incoming light into a spin-density matrix on the indices lambda = +-1
        # For linearly polarised light sqrt{S_1^2 + S_2^2 + S_3^2} = 1
        # We take gamma = 0

        #RCPL
        # p = 1.
        # alpha = pi/4
        # gamma = 0.

        #LCPL
        # p = 1.
        # alpha = -pi/4.
        # gamma = 0.

        #LPL
        # p = 1.
        # alpha = 0.
        # gamma = 0.

        #UPL
        # p = 0.
        # alpha = -pi/4.
        # gamma = 0.

        # S1 = - p * cos( 2 * alpha ) * cos( 2 * gamma )
        # S2 = - p * cos( 2 * alpha ) * sin( 2 * gamma )
        # S3 =   p * sin( 2 * alpha )

        # S1 = - cos( 2 * phi )
        # S2 = - sin( 2 * phi )
        # S3 =   0

        # Linear
        S1 = -1.
        S2 = 0.
        S3 = 0.

        # RCPL
        # S1 = 0.
        # S2 = 0.
        # S3 = 1.

        if      lambda1 == lambda2  == 1               return( (1.0 + Stokes.P3)/2. )
        elseif  lambda1 ==  1   &&   lambda2  == -1    return( (Stokes.P1 - Stokes.P2*im)/2. )
        elseif  lambda1 == -1   &&   lambda2  ==  1    return( (Stokes.P1 + Stokes.P2*im)/2. )
        elseif  lambda1 == lambda2  == -1              return( (1.0 - Stokes.P3)/2. )
        else    error("stop a")
        end
    end

    f = open("Beta-test.dat", "w")
    g = open("DCS-test.dat", "w")

    # Beta, Gamma1, Gamma3, Lambda1, Lambda3, Delta1
    angularParams = Dict{Float64, Vector{Basics.EmProperty}}()

    nx = 50
    # Define the 2x2 spins
    # Loop about all lines; a table is printed independently for each line
    lineCount = 0;
    for  line in lines
        lineCount += 1
        @printf( "Working on energy, E = %0.5F H; ", line.photonEnergy )
        println( "$(line.initialLevel.J) -- $(line.finalLevel.J)" )
        # readline()

        angCS = Tuple{AngularJ64, Float64, Float64, ComplexF64, ComplexF64}[]  # finalLevel.J, theta, phi, angCs.Coulomb, angCs.Babushkin
        sigmaBarC = 0.
        sigmaBarB = 0.

        println( "Line $(lineCount)" )
        chCount = 0;

        for  cha in line.channels
            # println( "Amplitude ", abs( cha.amplitude ) ^ 2, "gauge ", cha.gauge, "kappa ", cha.kappa )
            if  cha.gauge == Basics.Coulomb
                sigmaBarC += abs( cha.amplitude ) ^ 2
            elseif  cha.gauge == Basics.Babushkin
                sigmaBarB += abs( cha.amplitude ) ^ 2
            elseif cha.gauge == Basics.Magnetic
                sigmaBarC += abs( cha.amplitude ) ^ 2
                sigmaBarB += abs( cha.amplitude ) ^ 2
            end
        end

        # @printf( "sigmaBar (Coulomb)   = %25.10e\n", sigmaBarC )
        # @printf( "sigmaBar (Babushkin) = %25.10e\n", sigmaBarB )
        # println( "==========")

        angCsCoeff =  2. * pi^3 * 137.03599 / ( line.photonEnergy * ( Basics.twice(line.initialLevel.J) + 1 ) )

        # Loop over all angles theta, phi
        for theta in settings.thetas, phi in settings.phis
            csCoulomb = 1.;   csBabushkin = 1.
            angCsParams = Dict{Tuple, Vector{ComplexF64}}()
            for  X = 1:20  # Test for triangular conditions for X and continue otherwise
                # Loop twice about all channels but distinguish the two gauges
                for  cha in line.channels
                    for  chb in line.channels
                        if  cha.gauge == Basics.Coulomb  &&   chb.gauge == Basics.Babushkin   continue    end
                        if  chb.gauge == Basics.Coulomb  &&   cha.gauge == Basics.Babushkin   continue    end

                        s1 = Subshell(20,cha.kappa);   j1 = Basics.subshell_j(s1)
                        s2 = Subshell(20,chb.kappa);   j2 = Basics.subshell_j(s2)
                        s1_l   = Subshell( 101, cha.kappa )
                        ell1 = Basics.subshell_l( s1_l )
                        s1_2   = Subshell( 101, chb.kappa )
                        ell2 = Basics.subshell_l( s1_2 )

                        if  !AngularMomentum.isTriangle(cha.multipole.L, chb.multipole.L, X)                continue
                        elseif !AngularMomentum.isTriangle(cha.symmetry.J,  chb.symmetry.J, AngularJ64(X) ) continue
                        elseif !AngularMomentum.isTriangle(j1,  j2, AngularJ64(X) )                         continue
                        elseif mod( ell1 + ell2 + X, 2 ) != 0                                               continue
                        end

                        # Always use key1
                        key1 = (cha.multipole.L, chb.multipole.L, X, cha.multipole.electric, chb.multipole.electric)
                        key2 = (chb.multipole.L, cha.multipole.L, X, chb.multipole.electric, cha.multipole.electric)

                        # If both keys are not in the dict, then create Key1.
                        if !haskey( angCsParams, key1 ) && !haskey( angCsParams, key2 )
                            angCsParams[ key1 ] = [0.0 + 0.0im, 0.0 + 0.0im]
                        end

                        # If key2 is in dict, convert key1 to key2
                        if haskey( angCsParams, key2 )
                            key1 = key2
                        end

                        # Multiply by 2, because mu is fixed as 1//2
                        K = 2. * PhotoIonization.angularFunctionK(
                            cha.multipole.L, chb.multipole.L, X,
                            line.initialLevel.J, line.finalLevel.J,
                            cha.kappa, cha.symmetry.J, chb.kappa, chb.symmetry.J)

                        # Compute the summation over lambda's and mu's
                        mu = 1//2
                        for  lambda1 = -1:2:1,   lambda2 = -1:2:1
                            W = PhotoIonization.angularFunctionW(
                                theta,
                                cha.multipole.L, chb.multipole.L, X,
                                lambda1, lambda2,
                                cha.kappa, mu, chb.kappa, mu)

                            # PartW = W

                            W = W * spinDensityMatrix(lambda1, lambda2, settings.stokes) * exp( 1.0im * (lambda1 - lambda2) * phi )
                            W = W * (1.0im)^(chb.multipole.L - cha.multipole.L) * (lambda1*lambda2) / 2.
                            W = W * AngularMomentum.phaseMultipole(1.0im*lambda1, cha.multipole)
                            W = W * AngularMomentum.phaseMultipole(-1.0im*lambda2, chb.multipole)

                            betaPart  = K * W * cha.amplitude * conj( chb.amplitude )

                            if cha.gauge == Basics.Coulomb || chb.gauge == Basics.Coulomb
                                csCoulomb += ( betaPart / sigmaBarC )

                                if lambda1 * lambda2 == 1 && theta == phi == 0.
                                    if ( cha.gauge == Basics.Magnetic || chb.gauge == Basics.Magnetic && cha.gauge != chb.gauge )
                                        # chCount += 1
                                        # println("$(chCount) =>    $(cha.kappa)  --  $(chb.kappa)    $(j1)    $(j2)    $(cha.multipole.L)    $(chb.multipole.L)    $(X)    $(cha.multipole.electric)    $(chb.multipole.electric)")
                                        # println( key1 )
                                        # println("K = ", K)
                                        # println("W = ", W)
                                        # println("KW = ", K * W )
                                        # println("DD*  = ", cha.amplitude * conj( chb.amplitude ) )
                                        # println("KWDD*  = ", betaPart )
                                        # println("==============")
                                    end
                                    angCsParams[ key1 ][ 1 ] += ( betaPart / sigmaBarC )
                                end

                            elseif cha.gauge == Basics.Babushkin || chb.gauge == Basics.Babushkin
                                csBabushkin += ( betaPart / sigmaBarB )

                                if lambda1 * lambda2 == 1 && theta == phi == 0.
                                    angCsParams[ key1 ][ 2 ] += ( betaPart / sigmaBarB )
                                end
                            end
                        end
                    end
                end
            end

            if theta == phi == 0.0
                angularParams[line.photonEnergy * 27.2114079527] = [
                    EmProperty(0.),   # Beta_1
                    EmProperty(0.),   # Gamma1
                    EmProperty(0.),   # Gamma3
                    EmProperty(0.),   # Pi2
                    EmProperty(0.),   # Pi4
                    EmProperty(0.),   # Delta1
                    EmProperty(0.),   # Lambda2
                    EmProperty(0.),    # Lambda4
                    EmProperty(0.)    # Upsilon
                ]
                for (key, value) in angCsParams
                    print( f, "$(key[1])        $(key[2])        $(key[3])        $(key[4] ? 1 : 0)        $(key[5] ? 1 : 0)        " )
                    @printf( f, "%15.5e        %15.5e\n", value[1].re, value[2].re )

                    if key[1] == 1 && key[2] == 1 && key[3] == 2 && key[4] && key[5]
                        angularParams[line.photonEnergy * 27.2114079527][1] = EmProperty( value[1].re * -2., value[2].re * -2. )

                    elseif key[1] == 1 && key[2] == 2 && key[3] == 1 && key[4] && key[5]
                        angularParams[line.photonEnergy * 27.2114079527][2] = EmProperty( value[1].re, value[2].re )

                    elseif key[1] == 1 && key[2] == 2 && key[3] == 3 && key[4] && key[5]
                        angularParams[line.photonEnergy * 27.2114079527][3] = EmProperty( value[1].re, value[2].re )

                    elseif key[1] == 2 && key[2] == 2 && key[3] == 2 && key[4] && key[5]
                        angularParams[line.photonEnergy * 27.2114079527][4] = EmProperty( value[1].re, value[2].re )

                    elseif key[1] == 2 && key[2] == 2 && key[3] == 4 && key[4] && key[5]
                        angularParams[line.photonEnergy * 27.2114079527][5] = EmProperty( value[1].re, value[2].re )

                    elseif key[1] == 1 && key[2] == 1 && key[3] == 1 && key[4] && !key[5]
                        angularParams[line.photonEnergy * 27.2114079527][6] = EmProperty( value[1].re, value[2].re )

                    elseif key[1] == 1 && key[2] == 3 && key[3] == 2 && key[4] && key[5]
                        angularParams[line.photonEnergy * 27.2114079527][7] = EmProperty( value[1].re, value[2].re )

                    elseif key[1] == 1 && key[2] == 3 && key[3] == 4 && key[4] && key[5]
                        angularParams[line.photonEnergy * 27.2114079527][8] = EmProperty( value[1].re, value[2].re )

                    elseif key[1] == 1 && key[2] == 2 && key[3] == 2 && key[4] && !key[5]
                        angularParams[line.photonEnergy * 27.2114079527][9] = EmProperty( value[1].re, value[2].re )

                    end
                end
            end

            # push!(angCS, (line.finalLevel.J, theta, phi, csCoulomb * angCsCoeff, csBabushkin * angCsCoeff))
            print( g, "$(line.initialLevel.J) --  $(line.finalLevel.J)        " )
            @printf( g, "%0.5f\t\t",  theta )
            @printf( g, "%0.5f\t\t",  phi )
            @printf( g, "%15.5e\t\t", csCoulomb.re * angCsCoeff )
            @printf( g, "%15.5e\t\t", csBabushkin.re * angCsCoeff )
            @printf( g, "%15.8e\n",   line.photonEnergy * 27.2114079527 )
        end
        # println()

        println( g, "\n" )

        # angCS_axis = [[],[],[],[],[],[]]

        # for elm in angCS
        #     X_axis_C = elm[4] * sin( elm[2] ) * cos( elm[3] )
        #     X_axis_B = elm[5] * sin( elm[2] ) * cos( elm[3] )
        #     Y_axis_C = elm[4] * sin( elm[2] ) * sin( elm[3] )
        #     Y_axis_B = elm[5] * sin( elm[2] ) * sin( elm[3] )
        #     Z_axis_C = elm[4] * cos( elm[2] )
        #     Z_axis_B = elm[5] * cos( elm[2] )
        #     push!( angCS_axis[1], X_axis_C.re)
        #     push!( angCS_axis[2], Y_axis_C.re)
        #     push!( angCS_axis[3], Z_axis_C.re)
        #     push!( angCS_axis[4], X_axis_B.re)
        #     push!( angCS_axis[5], Y_axis_B.re)
        #     push!( angCS_axis[6], Z_axis_B.re)
        # end
        # Plot3DGraphDifferentialCrossSection( line, angCS_axis )
    end

    ## We're printing only Beta, Gamma1, Gamma3, Pi2 , Pi4
    parameterNames = ["Beta_1", "Gamma1", "Gamma3", "Pi2", "Pi4", "Delta1", "Lambda2", "Lambda4", "Upsilon2" ]
    for p in 1:9
        nx = 120
        println(stream, " ")
        println(stream, "  Angular " * parameterNames[p] * "-parameters for unpolarized target atoms with Ji = 0, 1/2, 1:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "omega"     ; na=4)
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "Energy e_p"; na=3)
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(30, "Cou -- angular " * parameterNames[p] * " -- Bab"; na=3);       sb = sb * TableStrings.hBlank(33)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
        #
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                            fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
            en = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "   "
            sa = sa * @sprintf("% .6e", angularParams[line.photonEnergy * 27.2114079527][p].Coulomb)     * "   "
            sa = sa * @sprintf("% .6e", angularParams[line.photonEnergy * 27.2114079527][p].Babushkin)   * "   "
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end

    close(f)
    close(g)

    return( nothing )
end

function ManualGamma1(lines::Array{PhotoIonization.Line,1})
    G1 = open("Gamma1.dat", "w")

    mu1 = mu2 = 1//2

    lcount = 0
    for  line in lines
        lcount += 1
        @printf( "Working on energy, E = %0.5F H; ", line.photonEnergy )
        println( "$(line.initialLevel.J) -- $(line.finalLevel.J)" )
        # readline()

        angCS = Tuple{AngularJ64, Float64, Float64, ComplexF64, ComplexF64}[]  # finalLevel.J, theta, phi, angCs.Coulomb, angCs.Babushkin
        sigmaBarC = 0.
        sigmaBarB = 0.

        for  cha in line.channels
            # println( "Amplitude ", abs( cha.amplitude ) ^ 2, "gauge ", cha.gauge, "kappa ", cha.kappa )
            if  cha.gauge == Basics.Coulomb
                sigmaBarC += abs( cha.amplitude ) ^ 2
            elseif  cha.gauge == Basics.Babushkin
                sigmaBarB += abs( cha.amplitude ) ^ 2
            end
        end

        Gamma1_C = 0.0 + 0.0im
        Gamma1_B = 0.0 + 0.0im
        println( "Line $(lcount)" )
        chCount = 0
        for  cha in line.channels, chb in line.channels
            if  cha.gauge == Basics.Coulomb  &&   chb.gauge == Basics.Babushkin   continue    end
            if  chb.gauge == Basics.Coulomb  &&   cha.gauge == Basics.Babushkin   continue    end
            s1 = Subshell(20,cha.kappa);   j1 = Basics.subshell_j(s1)
            s2 = Subshell(20,chb.kappa);   j2 = Basics.subshell_j(s2)
            s1_l   = Subshell( 101, cha.kappa )
            ell1 = Basics.subshell_l( s1_l )
            s1_2   = Subshell( 101, chb.kappa )
            ell2 = Basics.subshell_l( s1_2 )

            if  !AngularMomentum.isTriangle(cha.multipole.L, chb.multipole.L, 1)                continue
            elseif !AngularMomentum.isTriangle(cha.symmetry.J,  chb.symmetry.J, AngularJ64(1) ) continue
            elseif !AngularMomentum.isTriangle(j1,  j2, AngularJ64(1) )                         continue
            elseif mod( ell1 + ell2 + 1, 2 ) != 0                                               continue
            elseif !cha.multipole.electric || !chb.multipole.electric                           continue
            end

            Gamma1a = (Basics.twice(j1)+1) * (Basics.twice(cha.symmetry.J)+1) * (Basics.twice(j2)+1) * (Basics.twice(chb.symmetry.J)+1 )
            Gamma1h = 3.0 * sqrt(15) * sqrt(Gamma1a)
            Gamma1b = AngularMomentum.Wigner_3j(chb.multipole.L, cha.multipole.L, 1.0, -1.0, 1.0, 0) * AngularMomentum.Wigner_3j(j2, j1, 1.0, abs(mu2), -abs(mu1), mu2-mu1)
            Gamma1c = AngularMomentum.Wigner_6j(chb.symmetry.J, cha.symmetry.J, 1, j1, j2, line.finalLevel.J)
            Gamma1d = AngularMomentum.Wigner_6j(chb.symmetry.J, cha.symmetry.J, 1, cha.multipole.L, chb.multipole.L, line.initialLevel.J)
            Gamma1e = (1.0im)^(chb.multipole.L - cha.multipole.L )
            Gamma1f = cha.amplitude * conj( chb.amplitude )
            Gamma1g = AngularMomentum.phaseFactor([line.initialLevel.J, -1, line.finalLevel.J, 1, AngularJ64(1//2)])


            # if ( cha.gauge == chb.gauge == Basics.Coulomb )
            #     chCount += 1
            #     println("$(chCount) =>    $(cha.kappa)  --  $(chb.kappa)    $(j2)    $(j1)")
            #     println(Gamma1a)
            #     println(Gamma1b)
            #     println(Gamma1c)
            #     println(Gamma1d)
            #     println(Gamma1e)
            #     println(Gamma1f)
            #     println(Gamma1g)
            #     println(Gamma1a * Gamma1b * Gamma1c * Gamma1d * Gamma1e * Gamma1f * Gamma1g / sigmaBarC)
            #     println("----------------------")
            # end

            if ( cha.gauge == chb.gauge == Basics.Coulomb )
                Gamma1_C += ( Gamma1h * Gamma1b * Gamma1c * Gamma1d * Gamma1e * Gamma1f * Gamma1g )

            elseif ( cha.gauge == chb.gauge == Basics.Babushkin )
                Gamma1_B += ( Gamma1h * Gamma1b * Gamma1c * Gamma1d * Gamma1e * Gamma1f * Gamma1g )
            end
        end

        print( G1, "$(line.initialLevel.J) --  $(line.finalLevel.J)        " )
        @printf( G1, "%15.5e\t\t", Gamma1_C.re/sigmaBarC )
        @printf( G1, "%15.5e\t\t", Gamma1_B.re/sigmaBarB )
        @printf( G1, "%15.8e\n",   line.photonEnergy * 27.2114079527 )

    end
    close(G1)
end

function ManualGamma3(lines::Array{PhotoIonization.Line,1})
    G3 = open("Gamma3.dat", "w")

    mu1 = mu2 = 1//2

    X = 3
    lcount = 0
    for  line in lines
        lcount += 1
        @printf( "Working on energy, E = %0.5F H; ", line.photonEnergy )
        println( "$(line.initialLevel.J) -- $(line.finalLevel.J)" )
        # readline()

        angCS = Tuple{AngularJ64, Float64, Float64, ComplexF64, ComplexF64}[]  # finalLevel.J, theta, phi, angCs.Coulomb, angCs.Babushkin
        sigmaBarC = 0.
        sigmaBarB = 0.

        for  cha in line.channels
            # println( "Amplitude ", abs( cha.amplitude ) ^ 2, "gauge ", cha.gauge, "kappa ", cha.kappa )
            if  cha.gauge == Basics.Coulomb
                sigmaBarC += abs( cha.amplitude ) ^ 2
            elseif  cha.gauge == Basics.Babushkin
                sigmaBarB += abs( cha.amplitude ) ^ 2
            end
        end

        Gamma3_C = 0.0 + 0.0im
        Gamma3_B = 0.0 + 0.0im
        println( "Line $(lcount)" )
        chCount = 0
        for  cha in line.channels, chb in line.channels
            if  cha.gauge == Basics.Coulomb  &&   chb.gauge == Basics.Babushkin   continue    end
            if  chb.gauge == Basics.Coulomb  &&   cha.gauge == Basics.Babushkin   continue    end
            s1 = Subshell(20,cha.kappa);   j1 = Basics.subshell_j(s1)
            s2 = Subshell(20,chb.kappa);   j2 = Basics.subshell_j(s2)
            s1_l   = Subshell( 101, cha.kappa )
            ell1 = Basics.subshell_l( s1_l )
            s1_2   = Subshell( 101, chb.kappa )
            ell2 = Basics.subshell_l( s1_2 )

            if  !AngularMomentum.isTriangle(cha.multipole.L, chb.multipole.L, X)                continue
            elseif !AngularMomentum.isTriangle(cha.symmetry.J,  chb.symmetry.J, AngularJ64(X) ) continue
            elseif !AngularMomentum.isTriangle(j1,  j2, AngularJ64(X) )                         continue
            elseif mod( ell1 + ell2 + X, 2 ) != 0                                               continue
            elseif !cha.multipole.electric || !chb.multipole.electric                           continue
            end

            Gamma3a = (Basics.twice(j1)+1) * (Basics.twice(cha.symmetry.J)+1) * (Basics.twice(j2)+1) * (Basics.twice(chb.symmetry.J)+1 )
            Gamma3h = (2X+1) * sqrt(2cha.multipole.L+1) * sqrt(2chb.multipole.L+1) * sqrt(Gamma3a)
            Gamma3b = AngularMomentum.Wigner_3j(chb.multipole.L, cha.multipole.L, X, -1.0, 1.0, 0) * AngularMomentum.Wigner_3j(j2, j1, X, abs(mu2), -abs(mu1), mu2-mu1)
            Gamma3c = AngularMomentum.Wigner_6j(chb.symmetry.J, cha.symmetry.J, X, j1, j2, line.finalLevel.J)
            Gamma3d = AngularMomentum.Wigner_6j(chb.symmetry.J, cha.symmetry.J, X, cha.multipole.L, chb.multipole.L, line.initialLevel.J)
            Gamma3e = (1.0im)^(chb.multipole.L - cha.multipole.L )
            Gamma3f = cha.amplitude * conj( chb.amplitude )
            Gamma3g = AngularMomentum.phaseFactor([line.initialLevel.J, -1, line.finalLevel.J, 1, AngularJ64(1//2)])

            # if ( cha.gauge == chb.gauge == Basics.Coulomb )
            #     chCount += 1
            #     println("$(chCount) =>    $(cha.kappa)  --  $(chb.kappa)    $(j2)    $(j1)")
            #     println(Basics.twice(j1)+1)
            #     println(Basics.twice(j2)+1)
            #     println((Basics.twice(cha.symmetry.J)+1))
            #     println((Basics.twice(chb.symmetry.J)+1))
            #     println(3.0 * sqrt(15) * sqrt(Gamma3a))
            #     println(Gamma3c)
            #     println(Gamma3d)
            #     println("K = ", 3.0 * sqrt(15) * sqrt(Gamma3a) * Gamma3c * Gamma3d)
            #     println("W = ", Gamma3b * Gamma3e)
            #     println(" D * Ddash")
            #     println(Gamma3f)
            #     println(Gamma3g)
            #     println("----------------------")
            # end

            if ( cha.gauge == chb.gauge == Basics.Coulomb )
                Gamma3_C += ( Gamma3h * Gamma3b * Gamma3c * Gamma3d * Gamma3e * Gamma3f * Gamma3g )

            elseif ( cha.gauge == chb.gauge == Basics.Babushkin )
                Gamma3_B += ( Gamma3h * Gamma3b * Gamma3c * Gamma3d * Gamma3e * Gamma3f * Gamma3g )
            end
        end

        print( G3, "$(line.initialLevel.J) --  $(line.finalLevel.J)        " )
        @printf( G3, "%15.5e\t\t", Gamma3_C.re / sigmaBarC )
        @printf( G3, "%15.5e\t\t", Gamma3_B.re / sigmaBarB )
        @printf( G3, "%15.8e\n",   line.photonEnergy * 27.2114079527 )

    end
    close(G3)
end

function ManualPi2(lines::Array{PhotoIonization.Line,1})
    P2 = open("Pi2.dat", "w")

    mu1 = mu2 = 1//2

    X = 2
    lcount = 0
    for  line in lines
        lcount += 1
        @printf( "Working on energy, E = %0.5F H; ", line.photonEnergy )
        println( "$(line.initialLevel.J) -- $(line.finalLevel.J)" )
        # readline()

        angCS = Tuple{AngularJ64, Float64, Float64, ComplexF64, ComplexF64}[]  # finalLevel.J, theta, phi, angCs.Coulomb, angCs.Babushkin
        sigmaBarC = 0.
        sigmaBarB = 0.

        for  cha in line.channels
            # println( "Amplitude ", abs( cha.amplitude ) ^ 2, "gauge ", cha.gauge, "kappa ", cha.kappa )
            if  cha.gauge == Basics.Coulomb
                sigmaBarC += abs( cha.amplitude ) ^ 2
            elseif  cha.gauge == Basics.Babushkin
                sigmaBarB += abs( cha.amplitude ) ^ 2
            end
        end

        Pi2_C = 0.0 + 0.0im
        Pi2_B = 0.0 + 0.0im
        println( "Line $(lcount)" )
        chCount = 0
        for  cha in line.channels, chb in line.channels
            if  cha.gauge == Basics.Coulomb  &&   chb.gauge == Basics.Babushkin   continue    end
            if  chb.gauge == Basics.Coulomb  &&   cha.gauge == Basics.Babushkin   continue    end
            s1 = Subshell(20,cha.kappa);   j1 = Basics.subshell_j(s1)
            s2 = Subshell(20,chb.kappa);   j2 = Basics.subshell_j(s2)
            s1_l   = Subshell( 101, cha.kappa )
            ell1 = Basics.subshell_l( s1_l )
            s1_2   = Subshell( 101, chb.kappa )
            ell2 = Basics.subshell_l( s1_2 )

            if  !AngularMomentum.isTriangle(cha.multipole.L, chb.multipole.L, X)                continue
            elseif !AngularMomentum.isTriangle(cha.symmetry.J,  chb.symmetry.J, AngularJ64(X) ) continue
            elseif !AngularMomentum.isTriangle(j1,  j2, AngularJ64(X) )                         continue
            elseif mod( ell1 + ell2 + X, 2 ) != 0                                               continue
            elseif cha.multipole.L == chb.multipole.L != 2                                      continue
            elseif !cha.multipole.electric || !chb.multipole.electric                           continue
            end

            Pi2a = (Basics.twice(j1)+1) * (Basics.twice(cha.symmetry.J)+1) * (Basics.twice(j2)+1) * (Basics.twice(chb.symmetry.J)+1 )
            Pi2h = (2X+1) * sqrt(2cha.multipole.L+1) * sqrt(2chb.multipole.L+1) * sqrt(Pi2a)
            Pi2b = AngularMomentum.Wigner_3j(chb.multipole.L, cha.multipole.L, X, -1.0, 1.0, 0) * AngularMomentum.Wigner_3j(j2, j1, X, abs(mu2), -abs(mu1), mu2-mu1)
            Pi2c = AngularMomentum.Wigner_6j(chb.symmetry.J, cha.symmetry.J, X, j1, j2, line.finalLevel.J)
            Pi2d = AngularMomentum.Wigner_6j(chb.symmetry.J, cha.symmetry.J, X, cha.multipole.L, chb.multipole.L, line.initialLevel.J)
            Pi2e = (1.0im)^(chb.multipole.L - cha.multipole.L )
            Pi2f = cha.amplitude * conj( chb.amplitude )
            Pi2g = AngularMomentum.phaseFactor([line.initialLevel.J, -1, line.finalLevel.J, 1, AngularJ64(1//2)])

            # if ( cha.gauge == chb.gauge == Basics.Coulomb )
            #     chCount += 1
            #     println("$(chCount) =>    $(cha.kappa)  --  $(chb.kappa)    $(j2)    $(j1)")
            #     println(2X+1)
            #     println("L1 = ", cha.multipole.L)
            #     println("L2 = ", chb.multipole.L )
            #     println(sqrt(2(cha.multipole.L)+1) )
            #     println(sqrt(2chb.multipole.L+1))
            #     println(Basics.twice(j1)+1)
            #     println(Basics.twice(j2)+1)
            #     println((Basics.twice(cha.symmetry.J)+1))
            #     println((Basics.twice(chb.symmetry.J)+1))
            #     println( (2X + 1) * sqrt(2cha.multipole.L+1) * sqrt(2chb.multipole.L+1) * sqrt(Pi2a))
            #     println(Pi2c)
            #     println(Pi2d)
            #     println("K = ", (2X+1) * sqrt(2cha.multipole.L+1) * sqrt(2chb.multipole.L+1) * sqrt(Pi2a) * Pi2c * Pi2d)
            #     println("W = ", Pi2b * Pi2e)
            #     println(" D * Ddash")
            #     println(Pi2f)
            #     println(Pi2g)
            #     println("----------------------")
            # end

            if ( cha.gauge == chb.gauge == Basics.Coulomb )
                Pi2_C += ( Pi2h * Pi2b * Pi2c * Pi2d * Pi2e * Pi2f * Pi2g )

            elseif ( cha.gauge == chb.gauge == Basics.Babushkin )
                Pi2_B += ( Pi2h * Pi2b * Pi2c * Pi2d * Pi2e * Pi2f * Pi2g )
            end
        end

        print( P2, "$(line.initialLevel.J) --  $(line.finalLevel.J)        " )
        @printf( P2, "%15.5e\t\t", Pi2_C.re / sigmaBarC )
        @printf( P2, "%15.5e\t\t", Pi2_B.re / sigmaBarB )
        @printf( P2, "%15.8e\n",   line.photonEnergy * 27.2114079527 )

    end
    close(P2)
end

function ManualPi4(lines::Array{PhotoIonization.Line,1})
    P4 = open("Pi4.dat", "w")

    mu1 = mu2 = 1//2

    X = 4
    lcount = 0
    for  line in lines
        lcount += 1
        @printf( "Working on energy, E = %0.5F H; ", line.photonEnergy )
        println( "$(line.initialLevel.J) -- $(line.finalLevel.J)" )
        # readline()

        angCS = Tuple{AngularJ64, Float64, Float64, ComplexF64, ComplexF64}[]  # finalLevel.J, theta, phi, angCs.Coulomb, angCs.Babushkin
        sigmaBarC = 0.
        sigmaBarB = 0.

        for  cha in line.channels
            # println( "Amplitude ", abs( cha.amplitude ) ^ 2, "gauge ", cha.gauge, "kappa ", cha.kappa )
            if  cha.gauge == Basics.Coulomb
                sigmaBarC += abs( cha.amplitude ) ^ 2
            elseif  cha.gauge == Basics.Babushkin
                sigmaBarB += abs( cha.amplitude ) ^ 2
            end
        end

        Pi4_C = 0.0 + 0.0im
        Pi4_B = 0.0 + 0.0im
        println( "Line $(lcount)" )
        chCount = 0
        for  cha in line.channels, chb in line.channels
            if  cha.gauge == Basics.Coulomb  &&   chb.gauge == Basics.Babushkin   continue    end
            if  chb.gauge == Basics.Coulomb  &&   cha.gauge == Basics.Babushkin   continue    end
            s1 = Subshell(20,cha.kappa);   j1 = Basics.subshell_j(s1)
            s2 = Subshell(20,chb.kappa);   j2 = Basics.subshell_j(s2)
            s1_l   = Subshell( 101, cha.kappa )
            ell1 = Basics.subshell_l( s1_l )
            s1_2   = Subshell( 101, chb.kappa )
            ell2 = Basics.subshell_l( s1_2 )

            if  !AngularMomentum.isTriangle(cha.multipole.L, chb.multipole.L, X)                continue
            elseif !AngularMomentum.isTriangle(cha.symmetry.J,  chb.symmetry.J, AngularJ64(X) ) continue
            elseif !AngularMomentum.isTriangle(j1,  j2, AngularJ64(X) )                         continue
            elseif mod( ell1 + ell2 + X, 2 ) != 0                                               continue
            elseif cha.multipole.L == chb.multipole.L != 2                                      continue
            elseif !cha.multipole.electric || !chb.multipole.electric                           continue
            end

            Pi4a = (Basics.twice(j1)+1) * (Basics.twice(cha.symmetry.J)+1) * (Basics.twice(j2)+1) * (Basics.twice(chb.symmetry.J)+1 )
            Pi4h = (2X+1) * sqrt(2cha.multipole.L+1) * sqrt(2chb.multipole.L+1) * sqrt(Pi4a)
            Pi4b = AngularMomentum.Wigner_3j(chb.multipole.L, cha.multipole.L, X, -1.0, 1.0, 0) * AngularMomentum.Wigner_3j(j2, j1, X, abs(mu2), -abs(mu1), mu2-mu1)
            Pi4c = AngularMomentum.Wigner_6j(chb.symmetry.J, cha.symmetry.J, X, j1, j2, line.finalLevel.J)
            Pi4d = AngularMomentum.Wigner_6j(chb.symmetry.J, cha.symmetry.J, X, cha.multipole.L, chb.multipole.L, line.initialLevel.J)
            Pi4e = (1.0im)^(chb.multipole.L - cha.multipole.L )
            Pi4f = cha.amplitude * conj( chb.amplitude )
            Pi4g = AngularMomentum.phaseFactor([line.initialLevel.J, -1, line.finalLevel.J, 1, AngularJ64(1//2)])

            # if ( cha.gauge == chb.gauge == Basics.Coulomb )
            #     chCount += 1
            #     println("$(chCount) =>    $(cha.kappa)  --  $(chb.kappa)    $(j2)    $(j1)")
            #     println(Basics.twice(j1)+1)
            #     println(Basics.twice(j2)+1)
            #     println((Basics.twice(cha.symmetry.J)+1))
            #     println((Basics.twice(chb.symmetry.J)+1))
            #     println(3.0 * sqrt(15) * sqrt(Gamma3a))
            #     println(Gamma3c)
            #     println(Gamma3d)
            #     println("K = ", 3.0 * sqrt(15) * sqrt(Gamma3a) * Gamma3c * Gamma3d)
            #     println("W = ", Gamma3b * Gamma3e)
            #     println(" D * Ddash")
            #     println(Gamma3f)
            #     println(Gamma3g)
            #     println("----------------------")
            # end

            if ( cha.gauge == chb.gauge == Basics.Coulomb )
                Pi4_C += ( Pi4h * Pi4b * Pi4c * Pi4d * Pi4e * Pi4f * Pi4g )

            elseif ( cha.gauge == chb.gauge == Basics.Babushkin )
                Pi4_B += ( Pi4h * Pi4b * Pi4c * Pi4d * Pi4e * Pi4f * Pi4g )
            end
        end

        print( P4, "$(line.initialLevel.J) --  $(line.finalLevel.J)        " )
        @printf( P4, "%15.5e\t\t", Pi4_C.re / sigmaBarC )
        @printf( P4, "%15.5e\t\t", Pi4_B.re / sigmaBarB )
        @printf( P4, "%15.8e\n",   line.photonEnergy * 27.2114079527 )

    end
    close(P4)
end

function ManualDelta1(lines::Array{PhotoIonization.Line,1})
    D1 = open("Delta1.dat", "w")

    mu1 = mu2 = 1//2

    X = 1
    lcount = 0
    for  line in lines
        lcount += 1
        @printf( "Working on energy, E = %0.5F H; ", line.photonEnergy )
        println( "$(line.initialLevel.J) -- $(line.finalLevel.J)" )
        # readline()

        angCS = Tuple{AngularJ64, Float64, Float64, ComplexF64, ComplexF64}[]  # finalLevel.J, theta, phi, angCs.Coulomb, angCs.Babushkin
        sigmaBarC = 0.
        sigmaBarB = 0.

        for  cha in line.channels
            if  cha.gauge == Basics.Coulomb
                sigmaBarC += abs( cha.amplitude ) ^ 2
            elseif  cha.gauge == Basics.Babushkin
                sigmaBarB += abs( cha.amplitude ) ^ 2
            end
        end

        Delta1_C = 0.0 + 0.0im
        Delta1_B = 0.0 + 0.0im
        println( "Line $(lcount)" )
        chCount = 0
        for  cha in line.channels, chb in line.channels
            if  cha.gauge == Basics.Coulomb  &&   chb.gauge == Basics.Babushkin   continue    end
            if  chb.gauge == Basics.Coulomb  &&   cha.gauge == Basics.Babushkin   continue    end
            s1 = Subshell(20,cha.kappa);   j1 = Basics.subshell_j(s1)
            s2 = Subshell(20,chb.kappa);   j2 = Basics.subshell_j(s2)
            s1_l = Subshell( 101, cha.kappa )
            ell1 = Basics.subshell_l( s1_l )
            s1_2 = Subshell( 101, chb.kappa )
            ell2 = Basics.subshell_l( s1_2 )

            if  !AngularMomentum.isTriangle(cha.multipole.L, chb.multipole.L, X)                     continue
            elseif !AngularMomentum.isTriangle(cha.symmetry.J,  chb.symmetry.J, AngularJ64(X) )      continue
            elseif !AngularMomentum.isTriangle(j1,  j2, AngularJ64(X) )                              continue
            elseif mod( ell1 + ell2 + X, 2 ) != 0                                                    continue
            elseif cha.multipole.L == chb.multipole.L != 1                                           continue
            elseif cha.multipole.electric == chb.multipole.electric                                  continue
            end

            Delta1a = (Basics.twice(j1)+1) * (Basics.twice(cha.symmetry.J)+1) * (Basics.twice(j2)+1) * (Basics.twice(chb.symmetry.J)+1 )
            Delta1h = (2X+1) * sqrt(2cha.multipole.L+1) * sqrt(2chb.multipole.L+1) * sqrt(Delta1a)
            Delta1b = AngularMomentum.Wigner_3j(chb.multipole.L, cha.multipole.L, X, -1.0, 1.0, 0) * AngularMomentum.Wigner_3j(j2, j1, X, abs(mu2), -abs(mu1), mu2-mu1)
            Delta1c = AngularMomentum.Wigner_6j(chb.symmetry.J, cha.symmetry.J, X, j1, j2, line.finalLevel.J)
            Delta1d = AngularMomentum.Wigner_6j(chb.symmetry.J, cha.symmetry.J, X, cha.multipole.L, chb.multipole.L, line.initialLevel.J)
            Delta1e = (1.0im)^(chb.multipole.L - cha.multipole.L)
            Delta1f = cha.amplitude * conj( chb.amplitude )
            Delta1g = AngularMomentum.phaseFactor([line.initialLevel.J, -1, line.finalLevel.J, 1, AngularJ64(1//2)])
            Delta1i = (-1) ^ ( cha.multipole.electric ? 1 : 0 )

            if ( cha.gauge == Basics.Coulomb || chb.gauge == Basics.Coulomb )
                chCount += 1
                # println("$(chCount) =>    $(cha.kappa)  --  $(chb.kappa)    $(j1)    $(j2)    $(cha.multipole.L)    $(chb.multipole.L)    $(cha.multipole.electric)    $(chb.multipole.electric)")
                # println(2X+1)
                # println(sqrt(2cha.multipole.L+1))
                # println(sqrt(2chb.multipole.L+1))
                # println(Basics.twice(j1)+1)
                # println(Basics.twice(j2)+1)
                # println((Basics.twice(cha.symmetry.J)+1))
                # println((Basics.twice(chb.symmetry.J)+1))
                # println((2X+1) * sqrt(2cha.multipole.L+1) * sqrt(2chb.multipole.L+1) * sqrt(Delta1a))
                # println(Delta1b)
                # println(Delta1c)
                # println(Delta1d)
                # println("K = ", (2X+1) * sqrt(2cha.multipole.L+1) * sqrt(2chb.multipole.L+1) * sqrt(Delta1a) * Delta1c * Delta1d)
                # println("W = ", Delta1b * Delta1e)
                # println("DD* ", Delta1f)
                # println("KWDD* = ", Delta1h * Delta1b * Delta1c * Delta1d * Delta1e * Delta1f * Delta1g)
                # println("----------------------")
                Delta1_C += ( Delta1h * Delta1b * Delta1c * Delta1d * Delta1e * Delta1f * Delta1g * Delta1i )

            elseif ( cha.gauge == Basics.Babushkin || chb.gauge == Basics.Babushkin )
                Delta1_B += ( Delta1h * Delta1b * Delta1c * Delta1d * Delta1e * Delta1f * Delta1g * Delta1i )

            end
        end

        print( D1, "$(line.initialLevel.J) --  $(line.finalLevel.J)        " )
        @printf( D1, "%15.5e\t\t", Delta1_C.im / sigmaBarC )
        @printf( D1, "%15.5e\t\t", Delta1_B.im / sigmaBarB )
        @printf( D1, "%15.8e\n",   line.photonEnergy * 27.2114079527 )

    end
    close(D1)
end

function ManualLambda2(lines::Array{PhotoIonization.Line,1})
    L2 = open("Lambda2.dat", "w")

    mu1 = mu2 = 1//2

    X = 2
    lcount = 0
    for  line in lines
        lcount += 1
        @printf( "Working on energy, E = %0.5F H; ", line.photonEnergy )
        println( "$(line.initialLevel.J) -- $(line.finalLevel.J)" )
        # readline()

        angCS = Tuple{AngularJ64, Float64, Float64, ComplexF64, ComplexF64}[]  # finalLevel.J, theta, phi, angCs.Coulomb, angCs.Babushkin
        sigmaBarC = 0.
        sigmaBarB = 0.

        for  cha in line.channels
            # println( "Amplitude ", abs( cha.amplitude ) ^ 2, "gauge ", cha.gauge, "kappa ", cha.kappa )
            if  cha.gauge == Basics.Coulomb
                sigmaBarC += abs( cha.amplitude ) ^ 2
            elseif  cha.gauge == Basics.Babushkin
                sigmaBarB += abs( cha.amplitude ) ^ 2
            end
        end

        Lambda2_C = 0.0 + 0.0im
        Lambda2_B = 0.0 + 0.0im
        println( "Line $(lcount)" )
        chCount = 0
        for  cha in line.channels, chb in line.channels
            if  cha.gauge == Basics.Coulomb  &&   chb.gauge == Basics.Babushkin   continue    end
            if  chb.gauge == Basics.Coulomb  &&   cha.gauge == Basics.Babushkin   continue    end
            s1 = Subshell(20,cha.kappa);   j1 = Basics.subshell_j(s1)
            s2 = Subshell(20,chb.kappa);   j2 = Basics.subshell_j(s2)
            s1_l   = Subshell( 101, cha.kappa )
            ell1 = Basics.subshell_l( s1_l )
            s1_2   = Subshell( 101, chb.kappa )
            ell2 = Basics.subshell_l( s1_2 )

            if  !AngularMomentum.isTriangle(cha.multipole.L, chb.multipole.L, X)                continue
            elseif !AngularMomentum.isTriangle(cha.symmetry.J,  chb.symmetry.J, AngularJ64(X) ) continue
            elseif !AngularMomentum.isTriangle(j1,  j2, AngularJ64(X) )                         continue
            elseif mod( ell1 + ell2 + X, 2 ) != 0                                               continue
            elseif cha.multipole.L * chb.multipole.L != 3                                       continue
            elseif !cha.multipole.electric || !chb.multipole.electric                           continue
            end

            Lambda2a = (Basics.twice(j1)+1) * (Basics.twice(cha.symmetry.J)+1) * (Basics.twice(j2)+1) * (Basics.twice(chb.symmetry.J)+1 )
            Lambda2h = (2X+1) * sqrt(2cha.multipole.L+1) * sqrt(2chb.multipole.L+1) * sqrt(Lambda2a)
            Lambda2b = AngularMomentum.Wigner_3j(chb.multipole.L, cha.multipole.L, X, -1.0, 1.0, 0) * AngularMomentum.Wigner_3j(j2, j1, X, abs(mu2), -abs(mu1), mu2-mu1)
            Lambda2c = AngularMomentum.Wigner_6j(chb.symmetry.J, cha.symmetry.J, X, j1, j2, line.finalLevel.J)
            Lambda2d = AngularMomentum.Wigner_6j(chb.symmetry.J, cha.symmetry.J, X, cha.multipole.L, chb.multipole.L, line.initialLevel.J)
            Lambda2e = (1.0im)^(chb.multipole.L - cha.multipole.L )
            Lambda2f = cha.amplitude * conj( chb.amplitude )
            Lambda2g = AngularMomentum.phaseFactor([line.initialLevel.J, -1, line.finalLevel.J, 1, AngularJ64(1//2)])

            # if ( cha.gauge == chb.gauge == Basics.Coulomb )
            #     chCount += 1
            #     println("$(chCount) =>    $(cha.kappa)  --  $(chb.kappa)    $(j2)    $(j1)")
            #     println(2X+1)
            #     println("L1 = ", cha.multipole.L)
            #     println("L2 = ", chb.multipole.L )
            #     println(sqrt(2(cha.multipole.L)+1) )
            #     println(sqrt(2chb.multipole.L+1))
            #     println(Basics.twice(j1)+1)
            #     println(Basics.twice(j2)+1)
            #     println((Basics.twice(cha.symmetry.J)+1))
            #     println((Basics.twice(chb.symmetry.J)+1))
            #     println( (2X + 1) * sqrt(2cha.multipole.L+1) * sqrt(2chb.multipole.L+1) * sqrt(Pi2a))
            #     println(Pi2c)
            #     println(Pi2d)
            #     println("K = ", (2X+1) * sqrt(2cha.multipole.L+1) * sqrt(2chb.multipole.L+1) * sqrt(Pi2a) * Pi2c * Pi2d)
            #     println("W = ", Pi2b * Pi2e)
            #     println(" D * Ddash")
            #     println(Pi2f)
            #     println(Pi2g)
            #     println("----------------------")
            # end

            if ( cha.gauge == chb.gauge == Basics.Coulomb )
                Lambda2_C += ( Lambda2h * Lambda2b * Lambda2c * Lambda2d * Lambda2e * Lambda2f * Lambda2g )

            elseif ( cha.gauge == chb.gauge == Basics.Babushkin )
                Lambda2_B += ( Lambda2h * Lambda2b * Lambda2c * Lambda2d * Lambda2e * Lambda2f * Lambda2g )
            end
        end

        print( L2, "$(line.initialLevel.J) --  $(line.finalLevel.J)        " )
        @printf( L2, "%15.5e\t\t", Lambda2_C.re / sigmaBarC )
        @printf( L2, "%15.5e\t\t", Lambda2_B.re / sigmaBarB )
        @printf( L2, "%15.8e\n",   line.photonEnergy * 27.2114079527 )

    end
    close(L2)
end

function ManualLambda4(lines::Array{PhotoIonization.Line,1})
    L4 = open("Lambda4.dat", "w")

    mu1 = mu2 = 1//2

    X = 4
    lcount = 0
    for  line in lines
        lcount += 1
        @printf( "Working on energy, E = %0.5F H; ", line.photonEnergy )
        println( "$(line.initialLevel.J) -- $(line.finalLevel.J)" )
        # readline()

        angCS = Tuple{AngularJ64, Float64, Float64, ComplexF64, ComplexF64}[]  # finalLevel.J, theta, phi, angCs.Coulomb, angCs.Babushkin
        sigmaBarC = 0.
        sigmaBarB = 0.

        for  cha in line.channels
            # println( "Amplitude ", abs( cha.amplitude ) ^ 2, "gauge ", cha.gauge, "kappa ", cha.kappa )
            if  cha.gauge == Basics.Coulomb
                sigmaBarC += abs( cha.amplitude ) ^ 2
            elseif  cha.gauge == Basics.Babushkin
                sigmaBarB += abs( cha.amplitude ) ^ 2
            end
        end

        Lambda4_C = 0.0 + 0.0im
        Lambda4_B = 0.0 + 0.0im
        println( "Line $(lcount)" )
        chCount = 0
        for  cha in line.channels, chb in line.channels
            if  cha.gauge == Basics.Coulomb  &&   chb.gauge == Basics.Babushkin   continue    end
            if  chb.gauge == Basics.Coulomb  &&   cha.gauge == Basics.Babushkin   continue    end
            s1 = Subshell(20,cha.kappa);   j1 = Basics.subshell_j(s1)
            s2 = Subshell(20,chb.kappa);   j2 = Basics.subshell_j(s2)
            s1_l   = Subshell( 101, cha.kappa )
            ell1 = Basics.subshell_l( s1_l )
            s1_2   = Subshell( 101, chb.kappa )
            ell2 = Basics.subshell_l( s1_2 )

            if  !AngularMomentum.isTriangle(cha.multipole.L, chb.multipole.L, X)                continue
            elseif !AngularMomentum.isTriangle(cha.symmetry.J,  chb.symmetry.J, AngularJ64(X) ) continue
            elseif !AngularMomentum.isTriangle(j1,  j2, AngularJ64(X) )                         continue
            elseif mod( ell1 + ell2 + X, 2 ) != 0                                               continue
            elseif cha.multipole.L * chb.multipole.L != 3                                     continue
            elseif !cha.multipole.electric || !chb.multipole.electric                           continue
            end

            Lambda4a = (Basics.twice(j1)+1) * (Basics.twice(cha.symmetry.J)+1) * (Basics.twice(j2)+1) * (Basics.twice(chb.symmetry.J)+1 )
            Lambda4h = (2X+1) * sqrt(2cha.multipole.L+1) * sqrt(2chb.multipole.L+1) * sqrt(Lambda4a)
            Lambda4b = AngularMomentum.Wigner_3j(chb.multipole.L, cha.multipole.L, X, -1.0, 1.0, 0) * AngularMomentum.Wigner_3j(j2, j1, X, abs(mu2), -abs(mu1), mu2-mu1)
            Lambda4c = AngularMomentum.Wigner_6j(chb.symmetry.J, cha.symmetry.J, X, j1, j2, line.finalLevel.J)
            Lambda4d = AngularMomentum.Wigner_6j(chb.symmetry.J, cha.symmetry.J, X, cha.multipole.L, chb.multipole.L, line.initialLevel.J)
            Lambda4e = (1.0im)^(chb.multipole.L - cha.multipole.L )
            Lambda4f = cha.amplitude * conj( chb.amplitude )
            Lambda4g = AngularMomentum.phaseFactor([line.initialLevel.J, -1, line.finalLevel.J, 1, AngularJ64(1//2)])

            # if ( cha.gauge == chb.gauge == Basics.Coulomb )
            #     chCount += 1
            #     println("$(chCount) =>    $(cha.kappa)  --  $(chb.kappa)    $(j2)    $(j1)")
            #     println(Basics.twice(j1)+1)
            #     println(Basics.twice(j2)+1)
            #     println((Basics.twice(cha.symmetry.J)+1))
            #     println((Basics.twice(chb.symmetry.J)+1))
            #     println(3.0 * sqrt(15) * sqrt(Gamma3a))
            #     println(Gamma3c)
            #     println(Gamma3d)
            #     println("K = ", 3.0 * sqrt(15) * sqrt(Gamma3a) * Gamma3c * Gamma3d)
            #     println("W = ", Gamma3b * Gamma3e)
            #     println(" D * Ddash")
            #     println(Gamma3f)
            #     println(Gamma3g)
            #     println("----------------------")
            # end

            if ( cha.gauge == chb.gauge == Basics.Coulomb )
                Lambda4_C += ( Lambda4h * Lambda4b * Lambda4c * Lambda4d * Lambda4e * Lambda4f * Lambda4g )

            elseif ( cha.gauge == chb.gauge == Basics.Babushkin )
                Lambda4_B += ( Lambda4h * Lambda4b * Lambda4c * Lambda4d * Lambda4e * Lambda4f * Lambda4g )
            end
        end

        print( L4, "$(line.initialLevel.J) --  $(line.finalLevel.J)        " )
        @printf( L4, "%15.5e\t\t", Lambda4_C.re / sigmaBarC )
        @printf( L4, "%15.5e\t\t", Lambda4_B.re / sigmaBarB )
        @printf( L4, "%15.8e\n",   line.photonEnergy * 27.2114079527 )

    end
    close(L4)
end

function ManualUpsilon2(lines::Array{PhotoIonization.Line,1})
    U2 = open("Upsilon2.dat", "w")

    mu1 = mu2 = 1//2

    X = 2
    lcount = 0
    for  line in lines
        lcount += 1
        # @printf( "Working on energy, E = %0.5F H; ", line.photonEnergy )
        # println( "$(line.initialLevel.J) -- $(line.finalLevel.J)" )
        # readline()

        angCS = Tuple{AngularJ64, Float64, Float64, ComplexF64, ComplexF64}[]  # finalLevel.J, theta, phi, angCs.Coulomb, angCs.Babushkin
        sigmaBarC = 0.
        sigmaBarB = 0.

        for  cha in line.channels
            println( "Amplitude ", abs( cha.amplitude ) ^ 2, "gauge ", cha.gauge, "kappa ", cha.kappa )
            if  cha.gauge == Basics.Coulomb
                sigmaBarC += abs( cha.amplitude ) ^ 2
            elseif  cha.gauge == Basics.Babushkin
                sigmaBarB += abs( cha.amplitude ) ^ 2
            elseif  cha.gauge == Basics.Babushkin
                sigmaBarB += abs( cha.amplitude ) ^ 2
            end
        end

        Upsilon2_C = 0.0 + 0.0im
        Upsilon2_B = 0.0 + 0.0im
        println( "Line $(lcount)" )
        chCount = 0
        for  cha in line.channels, chb in line.channels
            if  cha.gauge == Basics.Coulomb  &&   chb.gauge == Basics.Babushkin   continue    end
            if  chb.gauge == Basics.Coulomb  &&   cha.gauge == Basics.Babushkin   continue    end
            s1 = Subshell(20,cha.kappa);   j1 = Basics.subshell_j(s1)
            s2 = Subshell(20,chb.kappa);   j2 = Basics.subshell_j(s2)
            s1_l   = Subshell( 101, cha.kappa )
            ell1 = Basics.subshell_l( s1_l )
            s1_2   = Subshell( 101, chb.kappa )
            ell2 = Basics.subshell_l( s1_2 )

            if  !AngularMomentum.isTriangle(cha.multipole.L, chb.multipole.L, X)                continue
            elseif !AngularMomentum.isTriangle(cha.symmetry.J,  chb.symmetry.J, AngularJ64(X) ) continue
            elseif !AngularMomentum.isTriangle(j1,  j2, AngularJ64(X) )                         continue
            elseif mod( ell1 + ell2 + X, 2 ) != 0                                               continue
            elseif (cha.multipole.L != 1 || !cha.multipole.electric )                           continue
            elseif (chb.multipole.L != 2 || chb.multipole.electric )                            continue
            end

            Upsilon2a = (Basics.twice(j1)+1) * (Basics.twice(cha.symmetry.J)+1) * (Basics.twice(j2)+1) * (Basics.twice(chb.symmetry.J)+1 )
            Upsilon2h = (2X+1) * sqrt(2cha.multipole.L+1) * sqrt(2chb.multipole.L+1) * sqrt(Upsilon2a)
            Upsilon2b = AngularMomentum.Wigner_3j(chb.multipole.L, cha.multipole.L, X, -1.0, 1.0, 0) * AngularMomentum.Wigner_3j(j2, j1, X, abs(mu2), -abs(mu1), mu2-mu1)
            Upsilon2c = AngularMomentum.Wigner_6j(chb.symmetry.J, cha.symmetry.J, X, j1, j2, line.finalLevel.J)
            Upsilon2d = AngularMomentum.Wigner_6j(chb.symmetry.J, cha.symmetry.J, X, cha.multipole.L, chb.multipole.L, line.initialLevel.J)
            Upsilon2e = (1.0im)^(chb.multipole.L - cha.multipole.L )
            Upsilon2f = cha.amplitude * conj( chb.amplitude )
            Upsilon2g = AngularMomentum.phaseFactor([line.initialLevel.J, -1, line.finalLevel.J, 1, AngularJ64(1//2)])
            Upsilon2i = (-1im) ^ ( chb.multipole.electric ? 1 : 0 ) * (1im) ^ ( cha.multipole.electric ? 1 : 0 )

            # if ( cha.gauge == chb.gauge == Basics.Coulomb )
            #     chCount += 1
            #     println("$(chCount) =>    $(cha.kappa)  --  $(chb.kappa)    $(j2)    $(j1)")
            #     println(2X+1)
            #     println("L1 = ", cha.multipole.L)
            #     println("L2 = ", chb.multipole.L )
            #     println(sqrt(2(cha.multipole.L)+1) )
            #     println(sqrt(2chb.multipole.L+1))
            #     println(Basics.twice(j1)+1)
            #     println(Basics.twice(j2)+1)
            #     println((Basics.twice(cha.symmetry.J)+1))
            #     println((Basics.twice(chb.symmetry.J)+1))
            #     println( (2X + 1) * sqrt(2cha.multipole.L+1) * sqrt(2chb.multipole.L+1) * sqrt(Pi2a))
            #     println(Pi2c)
            #     println(Pi2d)
            #     println("K = ", (2X+1) * sqrt(2cha.multipole.L+1) * sqrt(2chb.multipole.L+1) * sqrt(Pi2a) * Pi2c * Pi2d)
            #     println("W = ", Pi2b * Pi2e)
            #     println(" D * Ddash")
            #     println(Pi2f)
            #     println(Pi2g)
            #     println("----------------------")
            # end

            if ( cha.gauge == Basics.Coulomb || chb.gauge == Basics.Coulomb )
                Upsilon2_C += ( Upsilon2h * Upsilon2b * Upsilon2c * Upsilon2d * Upsilon2e * Upsilon2f * Upsilon2g * Upsilon2i )

            elseif ( cha.gauge == Basics.Babushkin || chb.gauge == Basics.Babushkin )
                Upsilon2_B += ( Upsilon2h * Upsilon2b * Upsilon2c * Upsilon2d * Upsilon2e * Upsilon2f * Upsilon2g * Upsilon2i )
            end
        end
        println( "-------------------" )

        # We're multiplying by 2 because we're considering only cha->electric and chb->magnetic
        print( U2, "$(line.initialLevel.J) --  $(line.finalLevel.J)        " )
        @printf( U2, "%15.5e\t\t", 2 * Upsilon2_C.re / sigmaBarC )
        @printf( U2, "%15.5e\t\t", 2 * Upsilon2_B.re / sigmaBarB )
        @printf( U2, "%15.8e\n",   line.photonEnergy * 27.2114079527 )

    end
    close(U2)
end

"""
`PhotoIonization.computeAmplitudesProperties(line::PhotoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64,
                                                    settings::PhotoIonization.Settings; printout::Bool=false)`
    ... to compute all amplitudes and properties of the given line; a line::PhotoIonization.Line is returned for which the amplitudes and
        properties are now evaluated.
"""
function  computeAmplitudesProperties(line::PhotoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64,
                                        settings::PhotoIonization.Settings; printout::Bool=false)
    nChannels    = PhotoIonization.Channel[];   nxChannels    = PhotoIonization.Channel[];
    contSettings = Continuum.Settings(false, nrContinuum);      csC = csB = dtC = dtB = 0.
    for channel in line.channels
        newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, line.initialLevel.basis.subshells)
        newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, newiLevel.basis.subshells)
        newiLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, channel.kappa), newiLevel)
        cOrbital, phase  = Continuum.generateOrbitalForLevel(line.electronEnergy, Subshell(101, channel.kappa), newfLevel, nm, grid, contSettings)
        #
        newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, channel.symmetry, newfLevel)
        nChannel   = PhotoIonization.Channel(channel.multipole, channel.gauge, channel.kappa, channel.symmetry, phase, 0.)
        amplitude  = PhotoIonization.amplitude("photoionization", nChannel, line.photonEnergy, newcLevel, newiLevel, grid)
        push!( nChannels, PhotoIonization.Channel(nChannel.multipole, nChannel.gauge, nChannel.kappa, nChannel.symmetry,
                                                    nChannel.phase, amplitude) )
        if       channel.gauge == Basics.Coulomb     csC = csC + abs(amplitude)^2
        elseif   channel.gauge == Basics.Babushkin   csB = csB + abs(amplitude)^2
        elseif   channel.gauge == Basics.Magnetic    csB = csB + abs(amplitude)^2;   csC = csC + abs(amplitude)^2
        end
        #
        if  settings.calcTimeDelay
            cOrbitalx, phasex  = Continuum.generateOrbitalForLevel(line.electronEnergy+0.01, Subshell(101, channel.kappa),
                                                                    newfLevel, nm, grid, contSettings)
            newcLevelx = Basics.generateLevelWithExtraElectron(cOrbitalx, channel.symmetry, newfLevel)
            nxChannel  = PhotoIonization.Channel(channel.multipole, channel.gauge, channel.kappa, channel.symmetry, phasex, 0.)
            amplitude  = PhotoIonization.amplitude("photoionization", nxChannel, line.photonEnergy+0.01, newcLevelx, newiLevel, grid)
            push!( nxChannels, PhotoIonization.Channel(nxChannel.multipole, nxChannel.gauge, nxChannel.kappa, nxChannel.symmetry,
                                                       nxChannel.phase, amplitude) )
        end
    end
    Ji2 = Basics.twice(line.initialLevel.J)
    ##x csFactor     = 4 * pi^2 * Defaults.getDefaults("alpha") * line.photonEnergy / (2*(Ji2 + 1))
    ##x csFactor     = 4 * pi^2 * Defaults.getDefaults("alpha") / line.photonEnergy / (Ji2 + 1)
    ##x csFactor     = 4 * pi^2 / Defaults.getDefaults("alpha") / line.photonEnergy / (Ji2 + 1)
    csFactor     = 8 * pi^3 / Defaults.getDefaults("alpha") / line.photonEnergy
    ##  csFactor     = csFactor / 2.   # Not fully clear, arises likely from the Rydberg normalization
    ##  Correct for energy normalization
    ##  if  line.electronEnergy < 2.0   csFactor = csFactor * (line.electronEnergy/2.0)^1.5     end
    crossSection = EmProperty(csFactor * csC, csFactor * csB)
    if    settings.calcAnisotropy
            angularBeta  = PhotoIonization.computeAngularBeta(line.initialLevel, line.finalLevel, nChannels)
    else  angularBeta  = EmProperty(0.)
    end
    if    settings.calcTimeDelay
          coherentDelay, incoherentDelay = PhotoIonization.computeTimeDelays(nChannels, nxChannels, 0.01, line.finalLevel.J)
    else  coherentDelay = EmProperty(0.);     incoherentDelay = EmProperty(0.)
    end
    #
    nLine = PhotoIonization.Line( line.initialLevel, line.finalLevel, line.electronEnergy, line.photonEnergy,
                                    crossSection, angularBeta, coherentDelay, incoherentDelay, nChannels)

    return( nLine )
end


"""
`PhotoIonization.computeAmplitudesPropertiesPlasma(line::PhotoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                                   settings::PhotoIonization.PlasmaSettings)`
    ... to compute all amplitudes and properties of the given line but for the given plasma model;
        a line::PhotoIonization.Line is returned for which the amplitudes and properties are now evaluated.
"""
function  computeAmplitudesPropertiesPlasma(line::PhotoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, settings::PhotoIonization.PlasmaSettings)
    newChannels = PhotoIonization.Channel[];;   contSettings = Continuum.Settings(false, grid.NoPoints-50);    csC = csB = 0.
    for channel in line.channels
        newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel)
        newiLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, channel.kappa), newiLevel)
        newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel)
        @warn "Adapt a proper continuum orbital for the plasma potential"
        cOrbital, phase  = Continuum.generateOrbitalForLevel(line.electronEnergy, Subshell(101, channel.kappa), newfLevel, nm, grid, contSettings)
        newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, channel.symmetry, newfLevel)
        newChannel = PhotoIonization.Channel(channel.multipole, channel.gauge, channel.kappa, channel.symmetry, phase, 0., 0.)
        @warn "Adapt a proper Auger amplitude for the plasma e-e interaction"
        amplitude = 1.0
        # amplitude  = PhotoIonization.amplitude("photoionization", channel, line.photonEnergy, newcLevel, newiLevel, grid)
        push!( newChannels, PhotoIonization.Channel(newChannel.multipole, newChannel.gauge, newChannel.kappa, newChannel.symmetry,
                                                    newChannel.phase, amplitude) )
        if       channel.gauge == Basics.Coulomb     csC = csC + abs(amplitude)^2
        elseif   channel.gauge == Basics.Babushkin   csB = csB + abs(amplitude)^2
        elseif   channel.gauge == Basics.Magnetic    csB = csB + abs(amplitude)^2;   csC = csC + abs(amplitude)^2
        end
    end
    Ji2 = Basics.twice(line.initialLevel.J)
    csFactor     = 4 * pi^2 * Defaults.getDefaults("alpha") * line.photonEnergy / (2*(Ji2 + 1))
    crossSection = EmProperty(csFactor * csC, csFactor * csB)
    println("plasma-photo cs = $crossSection")
    newline = PhotoIonization.Line( line.initialLevel, line.finalLevel, line.electronEnergy, line.photonEnergy,
                                    crossSection, EmProperty(0.), EmProperty(0.), EmProperty(0.), newChannels)
    return( newline )
end


"""
`PhotoIonization.computeAngularBeta(iLevel::Level, fLevel::Level, channels::Array{PhotoIonization.Channel,1})`
    ... to compute the beta anisotropy parameter for the photoionization transition i -> f with the given channels;
        here, the formula from Balashov (1994, Eq. 2.135) has been utilized. A beta::EmProperty parameter is returned.
        These (gauge-dependent) beta parameters are set to -9., if no amplitudes are calculated for the given gauge.
"""
function  computeAngularBeta(iLevel::Level, fLevel::Level, channels::Array{PhotoIonization.Channel,1})
    wnC = wnB = waC = waB = 0.
    for  ch in channels
        if  ch.multipole != E1   continue    end      # These beta parameters are valid only in E1 approximation
        if       ch.gauge == Basics.Coulomb     wnC = wnC + conj(ch.amplitude) * ch.amplitude
        elseif   ch.gauge == Basics.Babushkin   wnB = wnB + conj(ch.amplitude) * ch.amplitude
        else
        end
    end
    Ji = iLevel.J;    Jf = fLevel.J;
    for  ch  in channels
        if  ch.multipole != E1   continue    end  # These beta parameters are valid only in E1 approximation
        j = AngularMomentum.kappa_j(ch.kappa);    l = AngularMomentum.kappa_l(ch.kappa);   Jt = ch.symmetry.J
        for  chp  in channels
            if  ch.gauge !=  chp.gauge    continue     end
            jp = AngularMomentum.kappa_j(chp.kappa);    lp = AngularMomentum.kappa_l(chp.kappa);   Jtp = chp.symmetry.J
            wa = AngularMomentum.phaseFactor([Jf, -1, Ji, -1, AngularJ64(1//2)]) *
                    sqrt( AngularMomentum.bracket([Jt, Jtp, j, jp, l, lp]) ) *
                        AngularMomentum.ClebschGordan(l, AngularM64(0), lp, AngularM64(0), AngularJ64(2), AngularM64(0)) *
                        AngularMomentum.Wigner_6j(j, l, AngularJ64(1//2), lp, jp, AngularJ64(2)) *
                        AngularMomentum.Wigner_6j(j, Jt, Jf, Jtp, jp, AngularJ64(2)) *
                        AngularMomentum.Wigner_6j(AngularJ64(1), Jt, Ji, Jtp, AngularJ64(1), AngularJ64(2)) *
                        ch.amplitude * conj(chp.amplitude)
            if      chp.gauge == Basics.Coulomb      waC = waC + wa
            elseif  chp.gauge == Basics.Babushkin    waB = waB + wa
            else    error("stop a")
            end
        end
    end

    if  wnC == 0.   waC = ComplexF64(-9.0)    else    waC = sqrt(6.0) * waC / wnC      end
    if  wnB == 0.   waB = ComplexF64(-9.0)    else    waB = sqrt(6.0) * waB / wnB      end
    @show waC, waB

    return( EmProperty(waC.re, waB.re) )
end

"""
`PhotoIonization.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                    settings::PhotoIonization.Settings; output::Bool=true)`
    ... to compute the photoIonization transition amplitudes and all properties as requested by the given settings.
        A list of lines::Array{PhotoIonization.Lines} is returned.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                        settings::PhotoIonization.Settings; output::Bool=true)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")

    println("")
    printstyled("PhotoIonization.computeLines(): The computation of photo-ionization and properties starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    lines = PhotoIonization.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    PhotoIonization.displayLines(stdout, lines)    end
    # Determine maximum energy and check for consistency of the grid
    maxEnergy = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    # Calculate all amplitudes and requested properties
    newLines = PhotoIonization.Line[]
    for  line in lines
        println("\n>> Calculate photoionization amplitudes and properties for line: $(line.initialLevel.index) - $(line.finalLevel.index) " *
                "for the photon energy $(Defaults.convertUnits("energy: from atomic", line.photonEnergy)) " * Defaults.GBL_ENERGY_UNIT)
        newLine = PhotoIonization.computeAmplitudesProperties(line, nm, grid, nrContinuum, settings)
        push!( newLines, newLine)
    end
    # Print all results to screen
    PhotoIonization.displayPhases(newLines)
    PhotoIonization.displayResults(stdout, newLines, settings)
    ## PhotoIonization.displayTimeDelay(stdout, newLines, settings)
    if  printSummary   PhotoIonization.displayResults(iostream, newLines, settings)
                        ## PhotoIonization.displayTimeDelay(iostream, newLines, settings)
    end
    #
    # Add printout about non-E1 angle-differential cross sections, if required
    if  settings.calcNonE1AngleDifferentialCS
        PhotoIonization.computeDisplayNonE1AngleDifferentialCS(stdout, newLines, settings)
        PhotoIonization.ManualGamma1(newLines)
        PhotoIonization.ManualGamma3(newLines)
        PhotoIonization.ManualPi2(newLines)
        PhotoIonization.ManualPi4(newLines)
        PhotoIonization.ManualDelta1(newLines)
        PhotoIonization.ManualLambda2(newLines)
        PhotoIonization.ManualLambda4(newLines)
        PhotoIonization.ManualUpsilon2(newLines)
    end
    #
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end



"""
`PhotoIonization.computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                     settings::PhotoIonization.Settings, initialLevelSelection::LevelSelection;
                                     output=true, printout::Bool=true)`
    ... to compute the photoionization transition amplitudes and all properties as requested by the given settings. The computations
        and printout is adapted for large cascade computations by including only lines with at least one channel and by sending
        all printout to a summary file only. A list of lines::Array{PhotoIonization.Lines} is returned.
"""
function  computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                              settings::PhotoIonization.Settings, initialLevelSelection::LevelSelection; output=true, printout::Bool=true)

    lines = PhotoIonization.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    # if  settings.printBefore    PhotoIonization.displayLines(stdout, lines)    end
    # Determine maximum energy and check for consistency of the grid
    maxEnergy = 0.;   for  en in settings.photonEnergies     maxEnergy = max(maxEnergy, Defaults.convertUnits("energy: to atomic", en))   end
                      for  en in settings.electronEnergies   maxEnergy = max(maxEnergy, Defaults.convertUnits("energy: to atomic", en))   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    # Calculate all amplitudes and requested properties
    newLines = PhotoIonization.Line[]
    for  (i,line)  in  enumerate(lines)
        if  rem(i,10) == 0    println("> Photo line $i:  ... calculated ")    end
        # Do not compute line if initial level is not in initialLevelSelection()
        ##x @show Basics.selectLevel(line.initialLevel, initialLevelSelection)
        if  !Basics.selectLevel(line.initialLevel, initialLevelSelection)     continue    end
        #
        newLine = PhotoIonization.computeAmplitudesProperties(line, nm, grid, nrContinuum, settings, printout=printout)
        push!( newLines, newLine)
    end
    # Print all results to a summary file, if requested
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotoIonization.displayResults(iostream, newLines, settings)   end
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`PhotoIonization.computeLinesPlasma(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                    settings::PhotoIonization.PlasmaSettings; output::Bool=true)`
    ... to compute the photoIonization transition amplitudes and all properties as requested by the given settings.
        A list of lines::Array{PhotoIonization.Lines} is returned.
"""
function  computeLinesPlasma(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                settings::PhotoIonization.PlasmaSettings; output::Bool=true)
    println("")
    printstyled("PhotoIonization.computeLinesPlasma(): The computation of photo-ionization cross sections starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    photoSettings = PhotoIonization.Settings(settings.multipoles, settings.gauges, settings.photonEnergies, settings.electronEnergies,
                        false, false, false, settings.printBefore, settings.selectLines, settings.selectedLines, Basics.ExpStokes() )

    lines = PhotoIonization.determineLines(finalMultiplet, initialMultiplet, photoSettings)
    # Display all selected lines before the computations start
    if  settings.printBefore    PhotoIonization.displayLines(stdout, lines)    end
    # Calculate all amplitudes and requested properties
    newLines = PhotoIonization.Line[]
    for  line in lines
        newLine = PhotoIonization.computeAmplitudesPropertiesPlasma(line, nm, grid, settings)
        push!( newLines, newLine)
    end
    # Print all results to screen
    PhotoIonization.displayResults(stdout, lines, photoSettings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotoIonization.displayResults(iostream, lines, photoSettings)     end
    #
    if    output    return( lines )
    else            return( nothing )
    end
end



"""
`PhotoIonization.computePartialCrossSectionUnpolarized(gauge::EmGauge, Mf::AngularM64, line::PhotoIonization.Line)`
    ... to compute the partial photoionization cross section for initially unpolarized atoms by unpolarized plane-wave photons.
        A value::Float64 is returned.
"""
function  computePartialCrossSectionUnpolarized(gauge::EmGauge, Mf::AngularM64, line::PhotoIonization.Line)
    # Define an internal Racah expression for the summation over t, lambda
    function Racahexpr(kappa::Int64, Ji::AngularJ64, Jf::AngularJ64, Mf::AngularM64, J::AngularJ64, Jp::AngularJ64,
                        L::Int64, Lp::Int64, p::Int64, pp::Int64)
        # Determine the allowed values of t
        t1 = Basics.oplus( AngularJ64(Lp), Jf);    t2 = Basics.oplus( AngularJ64(L), Jf);    tList = intersect(t1, t2)
        wb = 0.
        for  t  in tList
            for  lambda = -1:2:1
                j = AngularMomentum.kappa_j(kappa);   Mf_lambda = Basics.add(AngularM64(lambda), Mf)
                wb = wb + (1.0im * lambda)^p * (-1.0im * lambda)^pp *
                        AngularMomentum.ClebschGordan( AngularJ64(Lp), AngularM64(lambda), Jf, Mf, t, Mf_lambda) *
                        AngularMomentum.ClebschGordan( AngularJ64(L),  AngularM64(lambda), Jf, Mf, t, Mf_lambda) *
                        AngularMomentum.Wigner_9j(j, Jp, Jf, J, Ji, AngularJ64(L), Jf, AngularJ64(Lp), t)
            end
        end

        return( wb )
    end

        wa = 0.0im;    Ji = line.initialLevel.J;    Jf = line.finalLevel.J;
    kappaList = PhotoIonization.getLineKappas(line)
    for  kappa in kappaList
        for  cha  in line.channels
            if  kappa != cha.kappa  ||  (gauge != cha.gauge  &&  gauge != Basics.Magnetic)   continue    end
            J = cha.symmetry.J;    L = cha.multipole.L;    if  cha.multipole.electric   p = 1   else    p = 0   end
            for  chp  in line.channels
                if  kappa != chp.kappa  ||  (gauge != cha.gauge  &&  gauge != Basics.Magnetic)    continue    end
                Jp = chp.symmetry.J;    Lp = chp.multipole.L;    if  chp.multipole.electric   pp = 1   else    pp = 0   end
                wa = wa + 1.0im^(L - Lp) * (-1)^(L + Lp) * AngularMomentum.bracket([AngularJ64(L), AngularJ64(Lp), J, Jp]) *
                        Racahexpr(kappa, Ji, Jf, Mf, J, Jp, L, Lp, p, pp) * cha.amplitude * conj(chp.amplitude)
            end
        end
    end
    csFactor = 8 * pi^3 * Defaults.getDefaults("alpha") / (2*line.photonEnergy * (Basics.twice(Ji) + 1))
    wa       = csFactor * wa

    return( wa )
end



"""
`PhotoIonization.computeStatisticalTensorUnpolarized(k::Int64, q::Int64, gauge::EmGauge, line::PhotoIonization.Line,
                                                            settings::PhotoIonization.Settings)`
    ... to compute the statistical tensor of the photoion in its final level after the photoionization of initially unpolarized atoms
        by plane-wave photons with given Stokes parameters (density matrix). A value::ComplexF64 is returned.
"""
function  computeStatisticalTensorUnpolarized(k::Int64, q::Int64, gauge::EmGauge, line::PhotoIonization.Line,
                                                settings::PhotoIonization.Settings)
    wa = 0.0im;    Ji = line.initialLevel.J;    Jf = line.finalLevel.J
    kappaList = PhotoIonization.getLineKappas(line);    P1 = settings.stokes.P1;   P2 = settings.stokes.P2;   P3 = settings.stokes.P3
    for  kappa in kappaList
        j = AngularMomentum.kappa_j(kappa)
        for  cha  in line.channels
            if  kappa != cha.kappa  ||  (gauge != cha.gauge  &&  gauge != Basics.Magnetic)   continue    end
            J = cha.symmetry.J;    L = cha.multipole.L;    if  cha.multipole.electric   p = 1   else    p = 0   end
            for  chp  in line.channels
                if  kappa != chp.kappa  ||  (gauge != cha.gauge  &&  gauge != Basics.Magnetic)    continue    end
                Jp = chp.symmetry.J;    Lp = chp.multipole.L;    if  chp.multipole.electric   pp = 1   else    pp = 0   end
                #
                for  lambda = -1:2:1
                    for  lambdap = -1:2:1
                        if  lambda == lambdap   wb = (1.0 + 0.0im +lambda*P3)    else    wb = P1 - lambda * P2 * im    end
                        wa = wa + wb * 1.0im^(L - Lp + p - pp) * lambda^p * lambdap^pp *
                                sqrt( AngularMomentum.bracket([AngularJ64(L), AngularJ64(Lp), J, Jp]) ) *
                                AngularMomentum.phaseFactor([J, +1, Jp, +1, Jf, +1, Ji, +1, j, +1, AngularJ64(1)]) *
                                AngularMomentum.ClebschGordan( AngularJ64(L),  AngularM64(lambda), AngularJ64(Lp),  AngularM64(-lambda),
                                                                AngularJ64(k),  AngularM64(q)) *
                                AngularMomentum.Wigner_6j(Jf, j, Jp, J, AngularJ64(k), Jf) *
                                AngularMomentum.Wigner_6j(Jp, Ji, AngularJ64(Lp), AngularJ64(L), AngularJ64(k), J) *
                                cha.amplitude * conj(chp.amplitude)
                    end
                end
            end
        end
    end

    wa = pi / (Basics.twice(Ji) + 1) * wa
    return( wa )
end


"""
`PhotoIonization.computeTimeDelays(channels::Array{PhotoIonization.Channel,1}, xchannels::Array{PhotoIonization.Channel,1},
                                    deltaE::Float64, Jf::AngularJ64)`
    ... to compute the -- coherent and incoherent -- time delay from the channels as calculated for two neighboured photon
        energies (deltaE = xE - E). Two tuple of two time delays (coherentDelay::EmProperty, incoherentDelay::EmProperty)
        is returned.
"""
function  computeTimeDelays(channels::Array{PhotoIonization.Channel,1}, xchannels::Array{PhotoIonization.Channel,1},
                            deltaE::Float64, Jf::AngularJ64)
    #== Calculate the coherent and incorent time delays separately; start with the coherent delays.
    # First attempt
    MeffC = MeffCx = MeffB = MeffBx = ComplexF64(0.)
    for  ch in channels
        if       ch.gauge == Basics.Coulomb     MeffC  = MeffC  + ch.amplitude
        elseif   ch.gauge == Basics.Babushkin   MeffB  = MeffB  + ch.amplitude
        end
    end
    for  ch in xchannels
        if       ch.gauge == Basics.Coulomb     MeffCx = MeffCx + ch.amplitude
        elseif   ch.gauge == Basics.Babushkin   MeffBx = MeffBx + ch.amplitude
        end
    end
    DeffC = angle(MeffC);    DeffCx = angle(MeffCx);    DeffB = angle(MeffB);    DeffBx = angle(MeffBx)
    coherentDelay = EmProperty( (DeffCx - DeffC) / deltaE,  - (DeffBx - DeffB) / deltaE )
    #
    # Calculate the coherent time delays
    # Second attempt
    nomTauC = denTauC = nomTauB = denTauB= ComplexF64(0.)
    for  (ic, ch) in  enumerate(channels)
        @show  "***", ch.gauge, ch.amplitude, xchannels[ic].amplitude, (xchannels[ic].amplitude - ch.amplitude) / deltaE
        if       ch.gauge == Basics.Coulomb
            nomTauC = nomTauC + conj(ch.amplitude) * (xchannels[ic].amplitude - ch.amplitude) / deltaE
            denTauC = denTauC + conj(ch.amplitude) * ch.amplitude
        elseif   ch.gauge == Basics.Babushkin
            nomTauB = nomTauB + conj(ch.amplitude) * (xchannels[ic].amplitude - ch.amplitude) / deltaE
            denTauB = denTauB + conj(ch.amplitude) * ch.amplitude
        end
    end
    @warn "Multiply coherentTauC by 20. for mean energy calibration !!!"
    coherentTauC = -im * nomTauC / denTauC / deltaE * 20.;   coherentTauB = im * nomTauB / denTauB
    @show coherentTauC, coherentTauB
    coherentDelay = EmProperty( coherentTauC.im, coherentTauB.im )    ==#
    #
    # Calculate the coherent time delays
    # Third attempt, explicit derivation by Nikolay, April 2024
    ## @warn "l0 = 1 ... for p_1/2, 3/2 splitting"
    @warn "l0 = 2 ... for d_3/2, 5/2 splitting"
    ampC = ampCx = ampB = ampBx = ComplexF64(0.)
    for  (ic, ch) in  enumerate(channels)
        j  = AngularMomentum.kappa_j(ch.kappa);   l  = AngularMomentum.kappa_l(ch.kappa);   l0 = AngularJ64(2);   Jc = Jf
        @show  "***", j, l, l0, Jc
        @show  "***", ch.gauge, ch.amplitude, xchannels[ic].amplitude, (xchannels[ic].amplitude - ch.amplitude) / deltaE
        ## factor = (1.0im)^( Basics.twice(l)/2 ) * sqrt(3/(4pi))  * AngularMomentum.phaseFactor([j, +1, l, +1, AngularJ64(1//2)]) *
        factor = sqrt(3/(4pi))  * AngularMomentum.phaseFactor([j, +1, l, +1, AngularJ64(1//2)]) *
                 AngularMomentum.ClebschGordan( l0, AngularM64(0), AngularJ64(1), AngularM64(0), l,  AngularM64(0)) *
                 sqrt(Basics.twice(j)+1) * AngularMomentum.Wigner_6j(Jc, AngularJ64(1//2), l0, l, AngularJ64(1), j)
        if       ch.gauge == Basics.Coulomb
            ampC  = ampC  + factor * ch.amplitude
            ampCx = ampCx + factor * xchannels[ic].amplitude
        elseif   ch.gauge == Basics.Babushkin
            ampB  = ampB  + factor * ch.amplitude
            ampBx = ampBx + factor * xchannels[ic].amplitude
        end
    end
    coherentTauC = (ampCx - ampC) / deltaE / ampC
    coherentTauB = (ampBx - ampB) / deltaE / ampB
    ## phiEffB      = log(ampB);       phiEffBx = log(ampBx);    coherentTaulnB = (phiEffBx.re - phiEffB.re) / deltaE
    @show coherentTauC, coherentTauB
    coherentDelay = EmProperty( coherentTauC.im, coherentTauB.im )
    #
    println("\n\nChannel amplitudes M_lj for photon energy: \n")
    for channel in channels
        println("   $(Subshell(11, channel.kappa))   $(channel.gauge)     phase=$(channel.phase)    M_lj=$(channel.amplitude) ")
    end
    # Incoherent time delays
    nomDeffC = nomDeffCx = nomDeffB = nomDeffBx = 0.
    denDeffC = denDeffCx = denDeffB = denDeffBx = 0.
    for  ch in channels
        if       ch.gauge == Basics.Coulomb     nomDeffC  = nomDeffC  + abs(ch.amplitude)^2 * ch.phase
                                                denDeffC  = denDeffC  + abs(ch.amplitude)^2
        elseif   ch.gauge == Basics.Babushkin   nomDeffB  = nomDeffB  + abs(ch.amplitude)^2 * ch.phase
                                                denDeffB  = denDeffB  + abs(ch.amplitude)^2
        end
    end
    for  ch in xchannels
        if       ch.gauge == Basics.Coulomb     nomDeffCx = nomDeffCx + abs(ch.amplitude)^2 * ch.phase
                                                denDeffCx = denDeffCx + abs(ch.amplitude)^2
        elseif   ch.gauge == Basics.Babushkin   nomDeffBx = nomDeffBx + abs(ch.amplitude)^2 * ch.phase
                                                denDeffBx = denDeffBx + abs(ch.amplitude)^2
        end
    end
    DeffC  = nomDeffC  / denDeffC;     DeffB  = nomDeffB  / denDeffB
    DeffCx = nomDeffCx / denDeffCx;    DeffBx = nomDeffBx / denDeffBx
    incoherentDelay = EmProperty( (DeffCx - DeffC) / deltaE,  (DeffBx - DeffB) / deltaE )

    return( coherentDelay, incoherentDelay)
end


"""
`PhotoIonization.determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoIonization.Settings)`
    ... to determine a list of photoionization Channel for a transitions from the initial to final level and by taking into account
        the particular settings of for this computation; an Array{PhotoIonization.Channel,1} is returned.
"""
function determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoIonization.Settings)
    channels = PhotoIonization.Channel[];
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity)
    if  Basics.UseCoulomb  in  settings.gauges   gaugeM = Basics.UseCoulomb    else   gaugeM = Basics.UseBabushkin    end
    for  mp in settings.multipoles
        symList = AngularMomentum.allowedMultipoleSymmetries(symi, mp)
        ##x println("mp = $mp   symi = $symi   symList = $symList")
        for  symt in symList
            kappaList = AngularMomentum.allowedKappaSymmetries(symt, symf)
            for  kappa in kappaList
                if  !(Basics.subshell_l(Subshell(10,kappa)) in  settings.lValues)    continue     end
                for  gauge in settings.gauges
                    # Include further restrictions if appropriate
                    if     string(mp)[1] == 'E'  &&   gauge == Basics.UseCoulomb
                        push!(channels, PhotoIonization.Channel(mp, Basics.Coulomb,   kappa, symt, 0., Complex(0.)) )
                    elseif string(mp)[1] == 'E'  &&   gauge == Basics.UseBabushkin
                        push!(channels, PhotoIonization.Channel(mp, Basics.Babushkin, kappa, symt, 0., Complex(0.)) )
                    elseif string(mp)[1] == 'M'  &&   gauge == gaugeM
                        push!(channels, PhotoIonization.Channel(mp, Basics.Magnetic,  kappa, symt, 0., Complex(0.)) )
                    end
                end
            end
        end
    end
    return( channels )
end


"""
`PhotoIonization.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoIonization.Settings)`
    ... to determine a list of PhotoIonization.Line's for transitions between levels from the initial- and final-state multiplets,
        and  by taking into account the particular selections and settings for this computation; an Array{PhotoIonization.Line,1}
        is returned. Apart from the level specification, all physical properties are set to zero during the initialization process.
"""
function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoIonization.Settings)
    lines    = PhotoIonization.Line[]
    shift_au = Defaults.convertUnits("energy: to atomic", settings.freeElectronShift)
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                # Add lines for all photon energies
                for  omega in settings.photonEnergies
                    # Photon energies are still in 'pre-defined' units; convert to Hartree
                    omega_au = Defaults.convertUnits("energy: to atomic", omega)
                    energy   = omega_au - (fLevel.energy - iLevel.energy) + shift_au
                    if  energy < 0.    continue   end
                    channels = PhotoIonization.determineChannels(fLevel, iLevel, settings)
                    push!( lines, PhotoIonization.Line(iLevel, fLevel, energy, omega_au, EmProperty(0.), EmProperty(0.),
                                                        EmProperty(0.), EmProperty(0.), channels) )
                end
                # Add lines for all electron energies
                for  en in settings.electronEnergies
                    # Electron energies are still in 'pre-defined' units; convert to Hartree
                    energy_au = Defaults.convertUnits("energy: to atomic", en) + shift_au
                    omega     = energy_au + (fLevel.energy - iLevel.energy)
                    if  energy_au < 0.    continue   end
                    channels = PhotoIonization.determineChannels(fLevel, iLevel, settings)
                    push!( lines, PhotoIonization.Line(iLevel, fLevel, energy_au, omega, EmProperty(0.), EmProperty(0.),
                                                        EmProperty(0.), EmProperty(0.), channels) )
                end
            end
        end
    end
    return( lines )
end


"""
`PhotoIonization.displayLineData(stream::IO, lines::Array{PhotoIonization.Line,1})`
    ... to display the calculated data, ordered by the initial levels and the photon energies involved.
        Neat tables of all initial levels and photon energies as well as all associated cross sections are printed
        but nothing is returned otherwise.
"""
function  displayLineData(stream::IO, lines::Array{PhotoIonization.Line,1})
    # Extract and display all initial levels by their total energy and symmetry
    energies = Float64[]
    for  line  in  lines    push!(energies, line.initialLevel.energy)  end
    energies = unique(energies);    energies = sort(energies)
    println(stream, "\n  Initial levels, available in the given photoionization line data:")
    println(stream, "\n  ", TableStrings.hLine(42))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                              sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                              sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(12, "Level energy"   ; na=3);
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(42))
    for  energy in energies
        isNew = true
        for  line in lines
            if  line.initialLevel.energy == energy  &&  isNew
                sa  = "  ";    sym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity )
                sa = sa * TableStrings.center(10, TableStrings.level(line.initialLevel.index); na=2)
                sa = sa * TableStrings.center(10, string(sym); na=4)
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy)) * "    "
                println(stream, sa);   isNew = false
            end
        end
    end
    #
    # Extract and display photon energies for which line data are available
    omegas = Float64[]
    for  line  in  lines    push!(omegas, line.photonEnergy)   end;
    omegas = unique(omegas);    omegas = sort(omegas)
    sa  = "    ";
    for   omega in omegas  sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", omega)) * "   "      end
    #
    println(stream, "\n  Photon energies $(TableStrings.inUnits("energy")) for which given photoionization line data are given:")
    println(stream, "\n" * sa)
    #
    # Extract and display total cross sections, ordered by the initial levels and photon energies, for which line data are available
    println(stream, "\n  Initial level, photon energies and total cross sections for given photoionization line data: \n")
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                              sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                              sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(12, "Omega"   ; na=3);
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(26, "Cou--Total CS--Bab"   ; na=5);
    sb = sb * TableStrings.center(26,TableStrings.inUnits("cross section"); na=5)
    sa = sa * TableStrings.center(18, "Electron energies"   ; na=3);
    sb = sb * TableStrings.center(18,TableStrings.inUnits("energy"); na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(107))
    #
    if  length(lines) == 0    println("  >>> No photoionization data in lines here ... return ");    return(nothing)    end
    #
    wLine = lines[1]
    for  energy in energies
        for  omega in omegas
            cs = Basics.EmProperty(0.);   electronEnergies = Float64[]
            for  line in lines
                if  line.initialLevel.energy == energy  &&   line.photonEnergy == omega     cs = cs + line.crossSection
                    push!(electronEnergies, line.electronEnergy); wLine = line   end
            end
            sa  = "  ";    sym = LevelSymmetry( wLine.initialLevel.J, wLine.initialLevel.parity )
            sa = sa * TableStrings.center(10, TableStrings.level(wLine.initialLevel.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", omega)) * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", cs.Coulomb))   * "  " *
                        @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", cs.Babushkin)) * "     "
            for en in electronEnergies  sa = sa * @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", en)) * "  " end
            println(stream, sa)
        end
    end


    return( nothing )
end


"""
`PhotoIonization.displayLines(stream::IO, lines::Array{PhotoIonization.Line,1})`
    ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all selected
        transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines(stream::IO, lines::Array{PhotoIonization.Line,1})
    nx = 175
    println(stream, " ")
    println(stream, "  Selected photoionization lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(10, "Energy_fi"; na=3);
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(10, "omega"; na=3);
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(12, "Energy e_p"; na=3);
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.flushleft(57, "List of multipoles, gauges, kappas and total symmetries"; na=4)
    sb = sb * TableStrings.flushleft(57, "partial (multipole, gauge, total J^P)                  "; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    #
    nchannels = 0
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
        energy = line.finalLevel.energy - line.initialLevel.energy
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", energy))              * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
        kappaMultipoleSymmetryList = Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry}[]
        for  i in 1:length(line.channels)
            push!( kappaMultipoleSymmetryList, (line.channels[i].kappa, line.channels[i].multipole, line.channels[i].gauge,
                                                line.channels[i].symmetry) )
            nchannels = nchannels + 1
        end
        ##x println("PhotoIonization-diplayLines-ad: kappaMultipoleSymmetryList = ", kappaMultipoleSymmetryList)
        wa = TableStrings.kappaMultipoleSymmetryTupels(85, kappaMultipoleSymmetryList)
        sb = sa * wa[1];    println(stream,  sb )
        for  i = 2:length(wa)
            sb = TableStrings.hBlank( length(sa) ) * wa[i];    println(stream,  sb )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n")
    println(stream, "  A total of $nchannels channels need to be calculated. \n")
    #
    return( nothing )
end


"""
`PhotoIonization.displayPhases(lines::Array{PhotoIonization.Line,1})`
    ... to display a list of lines, channels and phases of the continuum wave that have been selected due to the prior settings.
        A neat table of all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayPhases(lines::Array{PhotoIonization.Line,1})
    nx = 185
    println(" ")
    println("  Selected photoionization lines and phases:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(10, "Energy_fi"; na=3);
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(10, "omega"; na=3);
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(12, "Energy e_p"; na=3);
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.flushleft(57, "List of multipoles, gauges, kappas and total symmetries"; na=4)
    sb = sb * TableStrings.flushleft(57, "partial (multipole, gauge, total J^P, phase)           "; na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx))
    #
    nchannels = 0
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
        energy = line.finalLevel.energy - line.initialLevel.energy
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", energy))              * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
        kappaMultipoleSymmetryPhaseList = Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry,Float64}[]
        for  i in 1:length(line.channels)
            push!( kappaMultipoleSymmetryPhaseList, (line.channels[i].kappa, line.channels[i].multipole, line.channels[i].gauge,
                                                        line.channels[i].symmetry, line.channels[i].phase) )
            nchannels = nchannels + 1
        end
        wa = TableStrings.kappaMultipoleSymmetryPhaseTupels(85, kappaMultipoleSymmetryPhaseList)
        sb = sa * wa[1];    println( sb )
        for  i = 2:length(wa)
            sb = TableStrings.hBlank( length(sa) ) * wa[i];    println( sb )
        end
    end
    println("  ", TableStrings.hLine(nx), "\n")
    println("  A total of $nchannels channels has been calculated. \n")
    #
    return( nothing )
end



"""
`PhotoIonization.displayResults(stream::IO, lines::Array{PhotoIonization.Line,1}, settings::PhotoIonization.Settings)`
    ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but nothing
        is returned otherwise.
"""
function  displayResults(stream::IO, lines::Array{PhotoIonization.Line,1}, settings::PhotoIonization.Settings)
    nx = 130
    println(stream, " ")
    println(stream, "  Total photoionization cross sections for initially unpolarized atoms by unpolarized plane-wave photons:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(12, "omega"     ; na=4)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(12, "Energy e_p"; na=3)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(10, "Multipoles"; na=1);                              sb = sb * TableStrings.hBlank(13)
    sa = sa * TableStrings.center(30, "Cou -- Cross section -- Bab"; na=3)
    sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section") * "          " *
                                            TableStrings.inUnits("cross section"); na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    #
    wx = 1.0 ## wx = 1.0  (Schippers, August'23; wx = 2.0; wx = pi/2)
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
        en = line.finalLevel.energy - line.initialLevel.energy
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
        multipoles = EmMultipole[]
        for  ch in line.channels
            multipoles = push!( multipoles, ch.multipole)
        end
        multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
        sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=2)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wx * line.crossSection.Coulomb))     * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wx * line.crossSection.Babushkin))   * "                 "
        ##x sa = sa * @sprintf("%.6e", line.crossSection.Coulomb)     * "    "
        ##x sa = sa * @sprintf("%.6e", line.crossSection.Babushkin)   * "    "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    #==
    # Print, if useful, the sum of all cross sections (Schippers, August'23)
    tcs = Basics.EmProperty(0.);    for  line in lines   tcs = tcs + line.crossSection   end
    sa = repeat(" ", 102)
    sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wx * tcs.Coulomb))     * "    "
    sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wx * tcs.Babushkin))   * "                 "
    println(stream, sa)  ==#
    #
    #
    if  settings.calcAnisotropy
        nx = 120
        println(stream, " ")
        println(stream, "  Angular beta-parameters in E1 approximation for unpolarized target atoms with Ji = 0, 1/2, 1:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "omega"     ; na=4)
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "Energy e_p"; na=3)
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(30, "Cou -- angular beta_2 -- Bab"; na=3);       sb = sb * TableStrings.hBlank(33)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
        #
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                            fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
            en = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "   "
            sa = sa * @sprintf("% .6e", line.angularBeta.Coulomb)     * "   "
            sa = sa * @sprintf("% .6e", line.angularBeta.Babushkin)   * "   "
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    #
    if  settings.calcTimeDelay
        nx = 150
        println(stream, " ")
        println(stream, "  (Averaged) Time-delays of individual lines:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "omega"     ; na=4)
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "Energy e_p"; na=3)
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(30, "Cou -- coherent delay -- Bab"; na=3)
        sb = sb * TableStrings.center(30, TableStrings.inUnits("time"); na=4)
        sa = sa * TableStrings.center(30, "Cou -- incoherent delay -- Bab"; na=3);
        sb = sb * TableStrings.center(30, TableStrings.inUnits("time"); na=4)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
        #
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                            fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
            en = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "   "
            sa = sa * @sprintf("% .6e", Defaults.convertUnits("time: from atomic", line.coherentDelay.Coulomb))     * "   "
            sa = sa * @sprintf("% .6e", Defaults.convertUnits("time: from atomic", line.coherentDelay.Babushkin))   * "     "
            sa = sa * @sprintf("% .6e", Defaults.convertUnits("time: from atomic", line.incoherentDelay.Coulomb))   * "   "
            sa = sa * @sprintf("% .6e", Defaults.convertUnits("time: from atomic", line.incoherentDelay.Babushkin)) * "     "
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    #
    if  settings.calcPartialCs
        nx = 144
        println(stream, " ")
        println(stream, "  Partial cross sections for initially unpolarized atoms by unpolarized plane-wave photons:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
            sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "omega"     ; na=4)
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "Energy e_p"; na=3)
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(10, "Multipoles"; na=1);                              sb = sb * TableStrings.hBlank(13)
        sa = sa * TableStrings.center( 7, "M_f"; na=1);                                     sb = sb * TableStrings.hBlank(11)
        sa = sa * TableStrings.center(30, "Cou -- Partial cross section -- Bab"; na=3)
        sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section") * "          " *
                                            TableStrings.inUnits("cross section"); na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
        #
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                            fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
            en = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
            multipoles = EmMultipole[]
            for  ch in line.channels
                multipoles = push!( multipoles, ch.multipole)
            end
            multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
            sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=2)
            println(stream, sa)
            MfList = Basics.projections(line.finalLevel.J)
            for  Mf in MfList
                sb  = TableStrings.hBlank(97)
                wac = PhotoIonization.computePartialCrossSectionUnpolarized(Basics.Coulomb, Mf, line)
                wab = PhotoIonization.computePartialCrossSectionUnpolarized(Basics.Babushkin, Mf, line)
                sb  = sb * TableStrings.flushright( 8, string(Mf))                             * "       "
                sb  = sb * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wac.re))     * "    "
                sb  = sb * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wab.re))     * "    "
                println(stream, sb)
            end
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    #
    if  settings.calcTensors
        nx = 144
        println(stream, " ")
        println(stream, "  Reduced statistical tensors of the photoion in its final level after the photoionization ")
        println(stream, "  of initially unpolarized atoms by plane-wave photons with given Stokes parameters (density matrix):")
        println(stream, "\n     + tensors are printed for k = 0, 1, 2 and if non-zero only.")
        println(stream,   "     + Stokes parameters are:  P1 = $(settings.stokes.P1),  P2 = $(settings.stokes.P2),  P3 = $(settings.stokes.P3) ")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
            sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "omega"     ; na=4)
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "Energy e_p"; na=3)
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(10, "Multipoles"; na=1);                              sb = sb * TableStrings.hBlank(13)
        sa = sa * TableStrings.center(10, "k    q"; na=4);                                  sb = sb * TableStrings.hBlank(11)
        sa = sa * TableStrings.center(30, "Cou --  rho_kq (J_f)  -- Bab"; na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
        #
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                            fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
            en = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
            multipoles = EmMultipole[]
            for  ch in line.channels
                multipoles = push!( multipoles, ch.multipole)
            end
            multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
            sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=2)
            println(stream, sa)
            for  k = 0:2
                for q = -k:k
                    sb   = TableStrings.hBlank(102)
                    rhoc = PhotoIonization.computeStatisticalTensorUnpolarized(k, q, Basics.Coulomb,   line, settings)
                    rhob = PhotoIonization.computeStatisticalTensorUnpolarized(k, q, Basics.Babushkin, line, settings)
                    sb   = sb * string(k) * " " * TableStrings.flushright( 4, string(q))             * "       "
                    sb   = sb * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", rhoc.re))     * "    "
                    sb   = sb * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", rhob.re))     * "    "
                    println(stream, sb)
                end
            end
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    return( nothing )
end

#==
"""
`PhotoIonization.displayResultsDetailed(stream::IO, line::PhotoIonization.Line, settings::PhotoIonization.Settings)`
    ... to list the detailed results, energies, etc. for the given line. A neat table is printed but nothing
        is returned otherwise.
"""
function  displayResultsDetailed(stream::IO, line::PhotoIonization.Line, settings::PhotoIonization.Settings)
    symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)
    symf = LevelSymmetry(line.finalLevel.J, line.finalLevel.parity)
    aeffC = ComplexF64(0.);     aeffB = ComplexF64(0.);
    beffC = Float64(0.);        beffB = Float64(0.);       deffC = Float64(0.);     deffB = Float64(0.)
    for  channel in line.channels
        if      channel.gauge == Basics.Coulomb     aeffC = aeffC + channel.amplitude
                                                    beffC = beffC + abs(channel.amplitude)^2
                                                    deffC = deffC + abs(channel.amplitude)^2 * atan(channel.amplitude.im, channel.amplitude.re)
                                                    ## deffC = deffC + abs(channel.amplitude)^2 * angle(channel.amplitude)
        elseif  channel.gauge == Basics.Babushkin   aeffB = aeffB + channel.amplitude
                                                    beffB = beffB + abs(channel.amplitude)^2
                                                    deffB = deffB + abs(channel.amplitude)^2 * atan(channel.amplitude.im, channel.amplitude.re)
                                                    ## deffB = deffB + abs(channel.amplitude)^2 * angle(channel.amplitude)
        else
        end
    end
    deffC = deffC / beffC;    deffB = deffB / beffB
    #
    sa   = "\n  Results for PI line from the transition  $(line.initialLevel.index) -  $(line.finalLevel.index):  " *
            "  $symi  - $symf "

    println(stream, sa, "\n  ", TableStrings.hLine(length(sa)-3))
    println(stream, "\n  Photon energy               = " * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.photonEnergy)) *
                        "  " * TableStrings.inUnits("energy") )
    println(stream,   "  Electron energy             = " * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) *
                        "  " * TableStrings.inUnits("energy") )
    println(stream,   "  Total cross section         = " * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", line.crossSection.Coulomb)) *
                        "  " * TableStrings.inUnits("cross section") * " (Coulomb gauge)" )
    println(stream,   "                                " * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", line.crossSection.Babushkin)) *
                        "  " * TableStrings.inUnits("cross section") * " (Babushkin gauge)" )
    println(stream,   "  A^eff (Re, Im, abs, phi)    = " * @sprintf("%.4e", aeffC.re)   * "   " * @sprintf("%.4e", aeffC.im) * "   " *
                                                            @sprintf("%.4e", abs(aeffC)) * "   " * @sprintf("%.4e", angle(aeffC)) * "  (Coulomb gauge)" )
    println(stream,   "                              = " * @sprintf("%.4e", aeffB.re)   * "   " * @sprintf("%.4e", aeffB.im) * "   " *
                                                            @sprintf("%.4e", abs(aeffB)) * "   " * @sprintf("%.4e", angle(aeffB)) * "  (Babushkin gauge)" )
    println(stream,   "  A^r_eff, del_eff (non-coh)  = " * @sprintf("%.4e", sqrt(beffC))* "   " * @sprintf("%.4e", deffC)        * "  (Coulomb gauge)" )
    println(stream,   "                              = " * @sprintf("%.4e", sqrt(beffB))* "   " * @sprintf("%.4e", deffB)        * "  (Babushkin gauge)" )
    println(stream,   "  A^r_eff, del_eff (coherent) = " * @sprintf("%.4e", abs(aeffC)) * "   " * @sprintf("%.4e", atan(aeffC.im, aeffC.re)) * "  (Coulomb gauge)" )
    println(stream,   "                              = " * @sprintf("%.4e", abs(aeffB)) * "   " * @sprintf("%.4e", atan(aeffB.im, aeffB.re)) * "  (Babushkin gauge)" )
    println(stream, "\n  Kappa    Total J^P  Mp     Gauge               Amplitude          Real-Amplitude  Cross section (b)   Phase    atan()" *
                    "\n  " * TableStrings.hLine(112) )

    for  channel in line.channels
        Ji2 = Basics.twice(line.initialLevel.J)
        csFactor = 4 * pi^3 / Defaults.getDefaults("alpha") / line.photonEnergy
        # csFactor     = 4 * pi^2 * Defaults.getDefaults("alpha") * line.photonEnergy / (2*(Ji2 + 1))
        cs = csFactor * abs(channel.amplitude)^2
        sg = string(channel.gauge) * "      "
        sk = "    " * string(channel.kappa)
        sb = "    " * @sprintf("%.4e", channel.amplitude.re)
        sc = "    " * @sprintf("%.4e", channel.amplitude.im)
        sp = "    " * @sprintf("%.4e", channel.phase)
        sx = "    " * @sprintf("%.4e", atan(channel.amplitude.im, channel.amplitude.re) )
        sa = " " * sk[end-3:end] * "          " * string(channel.symmetry) * "    " * string(channel.multipole) * "     " * sg[1:11] *
                sb[end-12:end] * sc[end-12:end] * "    " * @sprintf("%.4e", abs(channel.amplitude)) *
                "       " * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", cs)) * "   " * sp[end-12:end] *
                "       " * sx[end-12:end]
        println(stream, sa)
    end

    return( nothing )
end   ==#



"""
`PhotoIonization.extractPhotonEnergies(lines::Array{PhotoIonization.Line,1})`
    ... to extract all photon energies for which photoionization data and cross sections are provided by lines;
        an list of energies::Array{Float64,1} is returned.
"""
function  extractPhotonEnergies(lines::Array{PhotoIonization.Line,1})
    pEnergies = Float64[]
    for  line  in  lines    push!(pEnergies, line.photonEnergy)  end
    pEnergies = unique(pEnergies)
    pEnergies = sort(pEnergies)

    return( pEnergies )
end


"""
`PhotoIonization.extractLines(lines::Array{PhotoIonization.Line,1}, omega::Float64)`
    ... to extract from lines all those that refer to the given omega;
        a reduced list rLines::Array{PhotoIonization.Line,1} is returned.
"""
function  extractLines(lines::Array{PhotoIonization.Line,1}, omega::Float64)
    rLines = PhotoIonization.Line[]
    for  line  in  lines
        if  line.photonEnergy == omega   @show  "aa", omega;    push!(rLines, line)  end
    end

    return( rLines )
end


"""
`PhotoIonization.extractCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, initialLevel)`
    ... to extract from lines the total PI cross section that refer to the given omega and initial level;
        a cross section cs::EmProperty is returned.
"""
function  extractCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, initialLevel)
    cs = Basics.EmProperty(0.)
    for  line  in  lines
        if  line.initialLevel.index == initialLevel.index  &&  line.initialLevel.energy == initialLevel.energy  &&
            line.photonEnergy       == omega
            cs = cs + line.crossSection
        end
    end

    return( cs )
end


"""
`PhotoIonization.extractCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, shell::Shell, initialLevel)`
    ... to extract from lines the total PI cross section that refer to the given omega and initial level and to
        to the ionization of an electron from shell; a cross section cs::EmProperty is returned.
"""
function  extractCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, shell::Shell, initialLevel)
    cs = Basics.EmProperty(0.)
    for  line  in  lines
        if  line.initialLevel.index == initialLevel.index  &&  line.initialLevel.energy == initialLevel.energy  &&
            line.photonEnergy       == omega
            # Now determined of whether the photoionization refers to the given shell
            confi     = Basics.extractLeadingConfiguration(line.initialLevel)
            conff     = Basics.extractLeadingConfiguration(line.finalLevel)
            shellOccs = Basics.extractShellOccupationDifference(confi::Configuration, conff::Configuration)
            if  length(shellOccs) > 1  ||  shellOccs[1][2] < 0      error("stop a")   end
            if  shellOccs[1][1] == shell   cs = cs + line.crossSection      end
        end
    end

    return( cs )
end


"""
`PhotoIonization.interpolateCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, initialLevel)`
    ... to interpolate (or extrapolate) from lines the total PI cross section for any given omega and
        initial level. The procedure applies a linear interpolation/extrapolation by just using the
        cross sections from the two nearest (given) omega points; a cross section cs::EmProperty is returned.
"""
function  interpolateCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, initialLevel)
    # First determine for which omegas cross sections are available and which associated cross sections
    # are to be applied in the interpolation
    omegas     = PhotoIonization.extractPhotonEnergies(lines)
    oms, diffs = Basics.determineNearestPoints(omega, 2, omegas)
    if  oms[1] < oms[2]     om1 = oms[1];   om2 = oms[2];   diff = - diffs[1]
    else                    om1 = oms[2];   om2 = oms[1];   diff = - diffs[2]
    end
    cs1        = PhotoIonization.extractCrossSection(lines, om1, initialLevel)
    cs2        = PhotoIonization.extractCrossSection(lines, om2, initialLevel)
    if  omega < om1  &&  omega < om2   ||  omega > om1  &&  omega > om2
        @warn("No extrapolation of cross sections; cs = 0.");
        return( Basics.EmProperty(0.) )
    end
    # Interpolate/extrapolate linearly with the cross section data for the two omegas
    wm         = 1 / (om2 - om1) * (cs2 - cs1)
    cs         = cs1  +  diff * wm
    ## @show om1, om2, diff, cs1.Coulomb, cs2.Coulomb, cs.Coulomb

    return( cs )
end


"""
`PhotoIonization.getLineKappas(line::PhotoIonization.Line)`
    ... returns a list of kappa-values (partial waves) which contribute to the given line, to which one or several channels are
        assigned. An kappaList::Array{Int64,1} is returned.
"""
function getLineKappas(line::PhotoIonization.Line)
    kappaList = Int64[]
    for  ch  in line.channels
        kappaList = union(kappaList, ch.kappa)
    end

    return( kappaList )
end

end # module
