using Distances: euclidean
#using Kalman
# using GaussianDistributions
# using GaussianDistributions: ⊕ # independent sum of Gaussian r.v.
import ProgressMeter
# using TensorCore

function calc_dists(X, q, w)

    n_samples, n_features = size(X);
    n_neighbors = ceil(Int, (1-q)*n_samples);

    dists = zeros(eltype(X), n_samples, n_neighbors)
    ind_neighbors = zeros(Int, n_samples, n_neighbors)

    δs = [copy(dists[:,1]) for _ in 1:Threads.nthreads()]
    Threads.@threads for j in eachindex(X)
        δ = δs[Threads.threadid()]
        @inbounds map!(x -> euclidean(x, X[j]), δ, vec(X))
        logdist = _calc_logdist(δ, w, n_samples, n_neighbors+1);
        # Here `logdist` is already the -log(euclidean) distance of one point
        # to all other points in the set.
        # Extract the threshold corresponding to the quantile defined
        thresh = quantile(logdist, q)
        # Filter to obtain Peaks Over Threshold (PoTs)
        # PoTs = logdist[findall(≥(thresh), logdist)]
        PoTs = filter(≥(thresh), logdist)
        ind_neighbors[j,:] = findall(≥(thresh), logdist)
        # We need to filter, because one entry will have infinite value,
        # because one entry has 0 Euclidean distance in the set.
        filter!(isfinite, PoTs)
        #exceedances = PoTs .- thresh
        dists[j,:] = exp.(PoTs)
    end
    return dists, ind_neighbors
end

function _calc_logdist(δ, w, n_samples, nneighbors)
    ind = []
    for i in 1:n_samples
        if i<=w
            istart = 1;
        else
            istart = i-w;
        end
        if i+w>n_samples
            iend = n_samples;
        else
            iend = i+w;
        end
        if (sum(δ[istart:i-1] .> δ[i]) + sum(δ[i+1:iend] .> δ[i])) == 0
            append!(ind,i);
        end
    end
    if size(ind)[1] < nneighbors
        nneighbors = size(ind)[1]
    end
    logdist = -ones(size(δ))*Inf;
    logdist[ind[sortperm(δ[ind])[1:nneighbors]]] = -log.(δ[ind[sortperm(δ[ind])[1:nneighbors]]])
    return logdist
end


function estimate_delay_multivariate(ss, τₘₐₓ, method::String)
    n_samples, n_features = size(ss)
    τ = zeros(Int, n_features)
    # Construct a ParforProgressbar object
    progress = ProgressMeter.Progress(
        n_features; desc = "Calculating τ", enabled = true
        )
    Threads.@threads for j=1:n_features
        τ[j] = estimate_delay(ss[:,j],method, 0:τₘₐₓ)
        ProgressMeter.next!(progress)
    end
    return τ
end

function filter_X(X, dists, ind_neighbors)
    n_samples, n_features = size(X)
    Xweighted = zeros(n_samples, n_features)
    for t=1:n_samples
        ind_t = [t; ind_neighbors[t,:]]
        # weights of analogues
        wa = [1; exp.(-dists[t,:])]
        # analogues
        a = X[ind_t,:]
        # weighted data
        Xweighted[t,:] = (transpose(wa) * a) ./ sum(wa)
    end
    return Xweighted
end

function _finite_diff(X, dt)
    n_samples, n_features = size(X)
    dX = zeros(n_samples-2, n_features)
    for i=1:n_features
        dX[:,i] = (X[3:end,i] .- 2*X[2:end-1,i] .+ X[1:end-2,i]) / (2*dt)
    end
    return dX
end




function filt_kalm(X, ind_neighbors, L, cov_L, decomp, ind_comps)
    ########################
    ### USEFUL VARIABLES ###
    ########################
    S = decomp["S"][ind_comps,ind_comps]
    V = decomp["V"][:,ind_comps]
    var_V = decomp["var_V"][:,ind_comps]
    # Extract the number of patches P
    P = Int(size(L,1)/2);
    # Extract the number of components inverted N and the number of epochs T
    N,T = size(transpose(V));
    # Extract the spatial slip component L relative to the strike and dip
    L_strike = L[1:P,:];
    L_dip    = L[P+1:2*P,:];
    # Initialize the covariances and variances variables to speed up the loop
    cov_L_strike_L_strike = zeros(P,P,N)
    cov_L_strike_L_dip    = zeros(P,P,N)
    cov_L_dip_L_dip       = zeros(P,P,N)
    L_strike_L_strike     = zeros(P,P,N)
    L_strike_L_dip        = zeros(P,P,N)
    L_dip_L_dip           = zeros(P,P,N)
    var_L_strike = zeros(P,N);
    var_L_dip    = zeros(P,N);
    # For every inverted component...
    for nn=1:N
        # Calculate the covariance matrices between the strike and strike
        # directions, between the strike and dip directions, and between the
        # dip and dip directions
        cov_L_strike_L_strike[:,:,nn] = cov_L[nn][1:P,1:P];
        cov_L_strike_L_dip[:,:,nn]    = cov_L[nn][P+1:2*P,1:P];
        cov_L_dip_L_dip[:,:,nn]       = cov_L[nn][P+1:2*P,P+1:2*P];

        L_strike_L_strike[:,:,nn] = L_strike[:,nn] * L_strike[:,nn]';
        L_strike_L_dip[:,:,nn]    = L_strike[:,nn] * L_dip[:,nn]';
        L_dip_L_dip[:,:,nn]       = L_dip[:,nn] * L_dip[:,nn]';
        # To speed up future calculation, store the variance in a simple matrix
        # instead of being the diagonal of matrices in different cells
        var_L_strike[:,nn]  = diag(cov_L_strike_L_strike[:,:,nn]);
        var_L_dip[:,nn]     = diag(cov_L_dip_L_dip[:,:,nn]);
    end
    
    P0 = Matrix(1.0I, 2*P, 2*P)
    p = Gaussian(X[1,:], P0)
    ps = [p] # vector of filtered Gaussians
    H = Matrix(1.0I, 2*P, 2*P)
    # Construct a ParforProgressbar object
    progress = ProgressMeter.Progress(
        T; desc = "test", enabled = true
        )
    #Threads.@threads for t=1:10
    for t=1:T
        # current state
        xnow = X[t,:]
        # next state
        xnext = X[t+1,:]
        # current neighbors
        xnow_neigh = X[ind_neighbors[t,:],:]
        # neighbors of next state
        xnext_neigh = X[ind_neighbors[t+1,:],:]
        # evolution of current neighbors
        xnowneigh_next = X[ind_neighbors[t,:].+1,:]
        # current distances of neighbors from state
        δnow = euclidean.(transpose(xnow_neigh), xnow)
        # distances of neighbors from state at next time step
        δnext = euclidean.(transpose(xnext_neigh), xnext)
        # distances of evolved neighbors from evolved state (i.e., state at next
        # time step)
        δnowneigh_next = euclidean.(transpose(xnowneigh_next), xnext)
        F = δnowneigh_next*pinv(δnow)
        
        cov_slip_strike_slip_strike = zeros(P,P)
        cov_slip_strike_slip_dip = zeros(P,P)
        cov_slip_dip_slip_dip = zeros(P,P)
        
        cov_slip_strike_slip_strike = 
            boxdot(cov_L_strike_L_strike, diag(S.^2 .* V[t,:].^2)) .+
            boxdot(L_strike_L_strike, diag(S.^2 .* var_V[t,:]))
        cov_slip_strike_slip_dip = 
            boxdot(cov_L_strike_L_dip, diag(S.^2 .* V[t,:].^2)) .+
            boxdot(L_strike_L_dip, diag(S.^2 .* var_V[t,:]))
        cov_slip_dip_slip_dip = 
            boxdot(cov_L_dip_L_dip, diag(S.^2 .* V[t,:].^2)) .+
            boxdot(L_dip_L_dip, diag(S.^2 .* var_V[t,:]))

        R = [cov_slip_strike_slip_strike   cov_slip_strike_slip_dip;
             cov_slip_strike_slip_dip      cov_slip_dip_slip_dip];
        cov_slip_strike_slip_strike = nothing
        cov_slip_strike_slip_dip    = nothing
        cov_slip_dip_slip_dip       = nothing

        #global p
        # predict
        p = F*p ⊕ Gaussian(zero(xnow), F * P0) #same as Gaussian(Φ*p.μ, Φ*p.Σ*Φ' + Q)
        # correct
        p, xnowres, _ = Kalman.correct(Kalman.JosephForm(), p,
                                        (Gaussian(xnow, R), H))
        push!(ps, p) # save filtered density
        
        ProgressMeter.next!(progress)
    end

    # n_samples, n_features = size(X)
    # for t = 1:n_samples-1
    #     # current state
    #     xnow = X[t,:]
    #     # next state
    #     xnext = X[t+1,:]
    #     # current neighbors
    #     xnow_neigh = X[ind_neighbors[t,:],:]
    #     # neighbors of next state
    #     xnext_neigh = X[ind_neighbors[t+1,:],:]
    #     # evolution of current neighbors
    #     xnowneigh_next = X[ind_neighbors[t,:].+1,:]
    #     # current distances of neighbors from state
    #     println(size(xnow_neigh))
    #     println(size(xnow))
    #     δnow = euclidean.(transpose(xnow_neigh), xnow)
    #     # distances of neighbors from state at next time step
    #     δnext = euclidean.(transpose(xnext_neigh), xnext)
    #     # distances of evolved neighbors from evolved state (i.e., state at next
    #     # time step)
    #     δnowneigh_next = euclidean.(transpose(xnowneigh_next), xnext)


    #     # observational noise characterized by how spread are the neighbors at a
    #     # given time
    #     var(δnow)
    #     Qnow = zeros(n_features, n_features)
        
    #     var(δnext)
    #     # dynamical noise characterized by distance of evolved neighbors
    #     var(δnowneigh_next)

    #     wnow = Gaussian(zero(xnow), Qnow)
    #     wnext = Gaussian(zero(xnext), Qnext)
    # end
    return ps
end



# function filt_kalm(X, ind_neighbors, L, cov_L, decomp, ind_comps)
#     ########################
#     ### USEFUL VARIABLES ###
#     ########################
#     S = decomp["S"][ind_comps,ind_comps]
#     V = decomp["V"][:,ind_comps]
#     var_V = decomp["var_V"][:,ind_comps]
#     # Extract the number of patches P
#     P = Int(size(L,1)/2);
#     # Extract the number of components inverted N and the number of epochs T
#     N,T = size(transpose(V));
#     # Extract the spatial slip component L relative to the strike and dip
#     L_strike = L[1:P,:];
#     L_dip    = L[P+1:2*P,:];
#     # Initialize the covariances and variances variables to speed up the loop
#     cov_L_strike_L_strike = zeros(N,P,P)
#     cov_L_strike_L_dip = zeros(N,P,P)
#     cov_L_dip_L_dip = zeros(N,P,P)
#     var_L_strike = zeros(P,N);
#     var_L_dip    = zeros(P,N);
#     # For every inverted component...
#     for nn=1:N
#         # Calculate the covariance matrices between the strike and strike
#         # directions, between the strike and dip directions, and between the
#         # dip and dip directions
#         cov_L_strike_L_strike[nn,:,:] = cov_L[nn][1:P,1:P];
#         cov_L_strike_L_dip[nn,:,:]    = cov_L[nn][P+1:2*P,1:P];
#         cov_L_dip_L_dip[nn,:,:]       = cov_L[nn][P+1:2*P,P+1:2*P];
#         # To speed up future calculation, store the variance in a simple matrix
#         # instead of being the diagonal of matrices in different cells
#         var_L_strike[:,nn]  = diag(cov_L_strike_L_strike[nn,:,:]);
#         var_L_dip[:,nn]     = diag(cov_L_dip_L_dip[nn,:,:]);
#     end
#     vnow = nothing
    
#     # Construct a ParforProgressbar object
#     progress = ProgressMeter.Progress(
#         T; desc = "test", enabled = true
#         )
    
#     Threads.@threads for tt=1:T
#     #for tt=1:T
#         cov_slip_strike_slip_strike = zeros(P,P)
#         cov_slip_strike_slip_dip = zeros(P,P)
#         cov_slip_dip_slip_dip = zeros(P,P)
#         for nn=1:N
#             cov_slip_strike_slip_strike_comp = 
#                 S[nn,nn]^2 * repeat([V[tt,nn]^2], P, P) .* 
#                 cov_L_strike_L_strike[nn,:,:] .+
#                 S[nn,nn]^2 * repeat([var_V[tt,nn]], P, P) *
#                 transpose(L_strike[:,nn])*L_strike[:,nn]
#             cov_slip_strike_slip_strike = cov_slip_strike_slip_strike .+
#                 cov_slip_strike_slip_strike_comp
#             cov_slip_strike_slip_dip_comp = nothing

#             cov_slip_strike_slip_dip_comp = 
#                 S[nn,nn]^2 * repeat([V[tt,nn]^2], P, P) .* 
#                 cov_L_strike_L_dip[nn,:,:] .+
#                 S[nn,nn]^2 * repeat([var_V[tt,nn]], P, P) *
#                 transpose(L_strike[:,nn])*L_dip[:,nn]
#             cov_slip_strike_slip_dip = cov_slip_strike_slip_dip .+
#                 cov_slip_strike_slip_dip_comp
#             cov_slip_strike_slip_dip_comp = nothing

#             cov_slip_dip_slip_dip_comp = 
#                 S[nn,nn]^2 * repeat([V[tt,nn]^2], P, P) .* 
#                 cov_L_dip_L_dip[nn,:,:] .+
#                 S[nn,nn]^2 * repeat([var_V[tt,nn]], P, P) *
#                 transpose(L_dip[:,nn])*L_dip[:,nn]
#             cov_slip_dip_slip_dip = cov_slip_dip_slip_dip .+
#                 cov_slip_dip_slip_dip_comp
#             cov_slip_dip_slip_dip_comp = nothing
#         end
#         R = [cov_slip_strike_slip_strike   cov_slip_strike_slip_dip;
#              cov_slip_strike_slip_dip      cov_slip_dip_slip_dip];
#         cov_slip_strike_slip_strike = nothing
#         cov_slip_strike_slip_dip    = nothing
#         cov_slip_dip_slip_dip       = nothing

#         vnow = Gaussian(zeros(2*P,1), R)
#         ProgressMeter.next!(progress)
#     end

#     # n_samples, n_features = size(X)
#     # for t = 1:n_samples-1
#     #     # current state
#     #     xnow = X[t,:]
#     #     # next state
#     #     xnext = X[t+1,:]
#     #     # current neighbors
#     #     xnow_neigh = X[ind_neighbors[t,:],:]
#     #     # neighbors of next state
#     #     xnext_neigh = X[ind_neighbors[t+1,:],:]
#     #     # evolution of current neighbors
#     #     xnowneigh_next = X[ind_neighbors[t,:].+1,:]
#     #     # current distances of neighbors from state
#     #     println(size(xnow_neigh))
#     #     println(size(xnow))
#     #     δnow = euclidean.(transpose(xnow_neigh), xnow)
#     #     # distances of neighbors from state at next time step
#     #     δnext = euclidean.(transpose(xnext_neigh), xnext)
#     #     # distances of evolved neighbors from evolved state (i.e., state at next
#     #     # time step)
#     #     δnowneigh_next = euclidean.(transpose(xnowneigh_next), xnext)


#     #     # observational noise characterized by how spread are the neighbors at a
#     #     # given time
#     #     var(δnow)
#     #     Qnow = zeros(n_features, n_features)
        
#     #     var(δnext)
#     #     # dynamical noise characterized by distance of evolved neighbors
#     #     var(δnowneigh_next)

#     #     wnow = Gaussian(zero(xnow), Qnow)
#     #     wnext = Gaussian(zero(xnext), Qnext)
#     # end
#     return vnow
# end


    # println("Initializing variables for covariance...")
    # # var_slip_strike = zeros(P,T)
    # # var_slip_dip = zeros(P,T)
    # cov_slip_strike_slip_strike = zeros(T,P,P)
    # cov_slip_strike_slip_dip = zeros(T,P,P)
    # cov_slip_dip_slip_dip = zeros(T,P,P)
    # println("Done")
    # # Construct a ParforProgressbar object
    # progress = ProgressMeter.Progress(
    #     N; desc = "Calculating slip strike and slip dip  covariance",
    #     enabled = true)
    # for nn=1:N
    #     # der_slip_strike_L_strike = 
    #     #     S[nn,nn].*repeat(transpose(V[:,nn]),P,1);
    #     # der_slip_strike_V =  repeat(L_strike[:,nn],1,T) .* S[nn,nn];
        
    #     # der_slip_dip_L_strike = 
    #     #     S[nn,nn].*repeat(transpose(V[:,nn]),P,1);
    #     # der_slip_dip_V =  repeat(L_dip[:,nn],1,T) .* S[nn,nn];
    #     # # Calculate the variance related to the component nn for the strike
    #     # # slip. The calculation is performed simultaneously for all the patches
    #     # # and epochs. This is thus a matrix PxT.
    #     # var_slip_strike_comp = 
    #     #     (der_slip_strike_L_strike.^2).*repeat(var_L_strike[:,nn],1,T) .+ 
    #     #     (der_slip_strike_V.^2).*repeat(transpose(var_V[:,nn]),P,1)
    #     # var_slip_dip_comp = 
    #     #     (der_slip_dip_L_strike.^2).*repeat(var_L_dip[:,nn],1,T) .+ 
    #     #     (der_slip_dip_V.^2).*repeat(transpose(var_V[:,nn]),P,1)
    #     # var_slip_strike = var_slip_strike + var_slip_strike_comp;
    #     # var_slip_dip = var_slip_dip + var_slip_dip_comp;
    #     # var_slip_strike_comp = nothing
    #     # var_slip_dip_comp = nothing

    #     cov_slip_strike_slip_strike_comp = 
    #         S[nn,nn]^2 * repeat(V[:,nn]^2, 1, P, P) .*
    #         permutedims(repeat(cov_L_strike_L_strike[nn,:,:], 1, 1, T),
    #             (3,1,2)) .+
    #         S[nn,nn]^2 * repeat(var_V[:,nn], 1, P, P) *
    #         permutedims(repeat(transpose(L_strike[:,nn])*L_strike[:,nn],1,1,T),
    #             (3,1,2))
    #     cov_slip_strike_slip_strike = cov_slip_strike_slip_strike .+
    #             cov_slip_strike_slip_strike_comp
    #     cov_slip_strike_slip_strike_comp = nothing

    #     cov_slip_strike_slip_dip_comp = 
    #         S[nn,nn]^2 * repeat(V[:,nn]^2, 1, P, P) .*
    #         permutedims(repeat(cov_L_strike_L_dip[nn,:,:], 1, 1, T),(3,1,2)) .+
    #         S[nn,nn]^2 * repeat(var_V[:,nn], 1, P, P) *
    #         permutedims(repeat(transpose(L_strike[:,nn])*L_dip[:,nn], 1, 1, T),
    #             (3,1,2))
    #     cov_slip_strike_slip_dip = cov_slip_strike_slip_dip .+
    #         cov_slip_strike_slip_dip_comp
    #     cov_slip_strike_slip_dip_comp = nothing

    #     cov_slip_dip_slip_dip_comp = 
    #         S[nn,nn]^2 * repeat(V[:,nn]^2, 1, P, P) .*
    #         permutedims(repeat(cov_L_dip_L_dip[nn,:,:], 1, 1, T),(3,1,2)) .+
    #         S[nn,nn]^2 * repeat(var_V[:,nn], 1, P, P) *
    #         permutedims(repeat(transpose(L_dip[:,nn])*L_dip[:,nn], 1, 1, T),
    #             (3,1,2))
    #     cov_slip_dip_slip_dip = cov_slip_dip_slip_dip .+
    #         cov_slip_dip_slip_dip_comp
    #     cov_slip_dip_slip_dip_comp = nothing

    #     ProgressMeter.next!(progress)
    # end