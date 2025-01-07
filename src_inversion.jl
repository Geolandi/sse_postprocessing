using Distances
using LinearAlgebra
import ProgressMeter
using GLM
include("src_coordinates.jl")

function create_priors(n_patches, options, ind_sigma, d_patches)
    lambda_strike   = options["inversion"]["lambda_strike"];
    lambda_dip      = options["inversion"]["lambda_dip"];
    lambda0         = options["inversion"]["lambda0"];
    etm0s           = options["inversion"]["sigma"];
    which_smoothing = options["inversion"]["which_smoothing"];

    m0 = options["inversion"]["m0"]*ones(2*n_patches,1)
    
    etm0 = etm0s[ind_sigma]
    
    if options["inversion"]["verbose"] == true
        println("sigma0 = ", etm0);
    end

    Em0 = ones(1,n_patches) .* etm0;

    # normalize standard deviation on model depending on the smoothing
    # (see Valette, 05/11/2009)
    Em_strike = Em0 .*lambda0 ./ lambda_strike;
    Em_dip    = Em0 .*lambda0 ./ lambda_dip;

    if which_smoothing == "gaussian"
        Cm0_strike = (Em_strike.^2) .* exp.(-0.5 * d_patches.^2 / lambda_strike)
        Cm0_dip = (Em_dip.^2) .* exp.(-0.5 * d_patches.^2 / lambda_dip)
    elseif which_smoothing == "exponential"
        Cm0_strike = (Em_strike.^2) .* exp.(-d_patches / lambda_strike);
        Cm0_dip = (Em_dip.^2) .* exp.(-d_patches / lambda_dip);
    else
        error(which_smoothing*"is not a valid option for "*
            "smoothing choice. Choices are "*
            "`expotential` or `gaussian`'");
    end

    Cm0 = [Cm0_strike                   zeros(n_patches,n_patches);
           zeros(n_patches,n_patches)   Cm0_dip];

    return m0, Cm0
end

function invert_comp(G,d,Cd,m0,Cm0)
    G_t = transpose(G);
    B = inv(G * Cm0 * G_t + Cd);
    m = m0 + (Cm0 * G_t / (G * Cm0 * G_t + Cd) * (d - G * m0));
    Cm = Cm0 - Cm0 * G_t * B * G * Cm0;
    return m, Cm
end

function invert_comps(ICA, ind_comps, fault, G, ind_sigma0_comps, options)
    n_ICs2invert = length(ind_comps)

    n_patches = length(fault["lon"][:,1])

    dobs = ICA["U"][:,ind_comps]
    Cdobs = Matrix{Float64}[]
    for i = 1:n_ICs2invert
        push!(Cdobs, Diagonal(ICA["var_U"][:,ind_comps[i]]))
    end

    # compute distances between patches in the local coordinates
    print("Calculating patch pairwise distances... ")
    d_patches = pairwise(Euclidean(),
        fault["xyz"][:,1:3], fault["xyz"][:,1:3], dims=1)
    println("Done")

    if options["flags"]["flag_parallel"] > 1
        # Construct a ParforProgressbar object
        progress = ProgressMeter.Progress(
            n_ICs2invert; desc = "Inverting ICs", enabled = true
            )
        
        m = zeros(2*n_patches,n_ICs2invert);
        #Cm = cell(1, n_ICs2invert);
        Cm = [zeros(2*n_patches,2*n_patches) for i=1:n_ICs2invert]
        # Cm = Matrix{Float64}[]
        # for i in 1:n_ICs2invert
        #     push!(Cm, zeros(2*n_patches,2*n_patches))
        # end

        # I don't know why, but it is actually faster if not parallel...
        #Threads.@threads for i in 1:n_ICs2invert
        for i in 1:n_ICs2invert
            if options["inversion"]["verbose"] == true
                println("   IC #", ind_comps[i])
            end
            m0, Cm0 = create_priors(n_patches, options, ind_sigma0_comps[i],
                                    d_patches);
            m[:,i],Cm[i] = invert_comp(G, dobs[:,i], Cdobs[i], m0, Cm0);
            ProgressMeter.next!(progress)
        end
        if options["inversion"]["verbose"] == true
            println("Done")
            println("")
        end
    else

    end

    return m, Cm
end

function create_model(m, Cm, decomp, fault, options, ind_comps)
    println("Creating slip model... ")
    result = Dict()
    result["timeline"] = decomp["timeline"]

    result["obs"],result["cov_obs"],result["var_obs"], result["obs_strike"],
    result["obs_dip"], result["var_obs_strike"], result["var_obs_dip"] = 
        calc_slip_cov_slip_var_slip(m,Cm,
            decomp["S"][ind_comps,ind_comps],decomp["V"][:,ind_comps],
            decomp["var_V"][:,ind_comps],options["inversion"]["calc_slip_cov"]);

    n_patches,n_samples = size(result["obs"]);

    rake_requested = options["inversion"]["rake_pos"];
    rake_lim1 = rake_requested - 90;
    rake_lim2 = rake_requested + 90;
    rake_obs = atand.(result["obs_dip"],result["obs_strike"]);
    ind_flip = rake_obs .< rake_lim1 .|| rake_obs .> rake_lim2;
    result["obs"][ind_flip] = -result["obs"][ind_flip];

    obs_local_vec = zeros(n_samples, n_patches, 3);
    # Construct a ParforProgressbar object
    progress = ProgressMeter.Progress(
        n_samples; desc = "Creating slip local vector", enabled = true
        )
    Threads.@threads for t in 1:n_samples
        for p=1:n_patches
            obs_local_vec[t, p, :] = fault2local(
                [result["obs_strike"][p,t],result["obs_dip"][p,t],0.0],
                fault["strike"][p],fault["dip"][p]);
        end
        ProgressMeter.next!(progress)
    end

    result["obs_strike_local"] = transpose(obs_local_vec[:,:,1]);
    result["obs_dip_local"]    = transpose(obs_local_vec[:,:,2]);
    println("Done")
    
    return result
end

function calc_slip(L,S,V)
    # calc_slip calculates the slip distribution at P patches for T epochs
    # using the decomposition of N components.
    # -------------------------------------------------------------------------
    # INPUT
    # L: Matrix containing the model vectors derived inverting the spatial
    #    vectors of the matrix U of the decomposition. Size: 2PxN
    # S: Diagonal matrix containing the weights associated to the inverted
    #    components. Size: NxN
    # V: Matrix containing the temporal information of the inverted components.
    #    Size: TxN
    # -------------------------------------------------------------------------
    # OUTPUT
    # slip_strike: Matrix containing the slip value for every patch and epoch
    #              in the strike direction. Size: PxT
    # slip_dip: Matrix containing the slip value for every patch and epoch in
    #           the dip direction. Size: PxT
    # slip: Matrix containing the slip value for every patch and epoch. Size:
    #       PxT
    # -------------------------------------------------------------------------
    # 
    # Adriano Gualandi - 20 Oct 2016
    # California Institute of Technology
    # Geological and Planetary Sciences Division
    # Adapted to Julia by Adriano Gualandi - 04 Nov 2024
    # University of Cambridge, UK
    # Department of Earth Sciences, Bullard Laboratories
    
    # Extract the number of patches P
    P = Int(size(L,1)/2);
    # Extract the spatial slip component L relative to the strike and dip
    L_strike = L[1:P,:];
    L_dip    = L[P+1:2*P,:];
    # Calculate the slip along the strike and dip directions (eq. 2) and the
    # final slip (eq. 1). A matrix notation is used to speed up and avoid a
    # loop through the patches.
    slip_strike = L_strike*S*transpose(V);
    slip_dip    = L_dip*S*transpose(V);
    slip = sqrt.(slip_strike.^2 .+ slip_dip.^2);

    return slip_strike,slip_dip,slip
end

function calc_slip_cov_slip_var_slip(L,cov_L,S,V,var_V,flag_compute_cov)
    # calc_slip_cov_slip_var_slip computes the slip and the corresponding
    # covariance matrix on the entire fault for all the epochs taken into
    # account. It reconstructs the slip on the patch p at time t using the
    # following equation:
    # 
    # slip_p(t) = sqrt(slip_strike_p(t)^2 + slip_dip_p^2)               (eq. 1)
    # 
    # where slip_strike_p(t) and slip_dip_p(t) are the slip on patch p at time
    # t along the strike and dip directions, respectively. These quantities are
    # calculated as:
    # 
    # slip_strike_p(t) = L_strike_p*S*V_t                              (eq. 2a)
    # slip_dip_p(t)    = L_dip_p*S*V_t                                 (eq. 2b)
    # 
    # where L_strike and L_dip are the matrices of size PxN containing the
    # result of the inversion of the spatial components U, S is a matrix of
    # size NxN, and V is the matrix of size NxT. P is the total number of
    # patches, N is the total number of inverted components, and T is the total
    # number of epochs. Thus, L_strike_p and L_dip_p are the p-th rows of the
    # L_strike and L_dip matrices, respectively, and V_t is the t-th rows of
    # the V matrix.
    # The function calculates also the covariance matrix of the slip at every
    # single epoch using equation (S16) of the Supplementary Material of
    # Gualandi et al. (2016b):
    # 
    #                                                                   (eq. 3)
    # cov(slip_q(t),slip_r(t)) = sum_{j=1}^{3N} sum_{k=1}^{3N}
    #   {
    #    partiald(slip_q(t)/x_q^j(t))|_E[x_q^j(t)] * 
    #    * cov(x_q^j(t),y_r^k(t)) *
    #    * partiald(slip_r(t)/y_r^k(t))|_E[y_r^k(t)]
    #   } =
    #  = ((nabla_{x_q(t)}(slip_q(t)))|_E[x_q(t)])' * 
    #    * C_{x_q(t),y_r(t)} *
    #    * ((nabla_{y_r(t)}(slip_r(t)))|_E[y_r(t)])
    # 
    # where sum_{j=1}^{3N} and sum_{k=1}^{3N} indicate the sum from j,k=1 to
    # 3N, partiald(f(z)/z) is the partial derivative of f(z) with respect to z,
    # the vectors x_q(t) and y_r(t) are the set of random variables expressed
    # by:
    # 
    # x_q(t) = [L_strike_q', L_dip_q', V_t']'                          (eq. 4a)
    # y_r(t) = [L_strike_r', L_dip_r', V_t']'                          (eq. 4b)
    # 
    # The size of these two vectors is 3Nx1. Indeed, L_strike_p, L_dip_p, and
    # V_t have size 1xN, so that L_strike_p', L_dip_p', and V_t' have size Nx1
    # each. The j-th and k-th elements of these two vectors are written as
    # x_q^j(t) and y_r^k(t).
    # The symbol |_E[] indicates that the partial derivative is calculated in
    # the point given by the expected value of the random variable indicated in
    # the brackets.
    # In compact matrix notation, the double sum and the elementwise partial
    # derivative are written using the nabla operator for the gradient, and the
    # covariance matrix C.
    # Since the slip output is a matrix of size PxT, while the covariance is a
    # matrix PxP for every epoch, the function outputs also the matrix with the
    # variance associated to every single patch at every single epoch:
    # var_slip. This matrix is useful because when we plot the time series
    # relative to a given patch we use this value to estimate the standard
    # deviation on the slip, i.e. its error bar, which neglects the potential
    # covariances with other patches and epochs. The var_slip matrix is
    # calculated using (eq. 3) for q=r=p, with p that goes from 1 to P. Using
    # (eq. 3) and the assumptions of non-correlation between different
    # components as well as spatial and temporal contributions, the variance of
    # the slip_p(t) is given by:
    # 
    #                                                                   (eq. 5)
    # var(slip_p(t)) = sum_{n=1}^{N}
    #   [partiald(slip_p(t)/L_strike_p^n)^2 * var(L_strike_p^n) + 
    #   + partiald(slip_p(t)/L_dip_p^n)^2 * var(L_dip_p^n) +
    #   + 2*partiald(slip_p(t)/L_strike_p^n)*partiald(slip_p(t)/L_dip_p^n)*
    #     *cov(L_strike_p^n,L_dip_p^n) +
    #   + partiald(slip_p(t)/V^n(t))^2 * var(V^n(t))]
    # 
    # where *^n indicates the element corresponding to the n-th inverted
    # component.
    # -------------------------------------------------------------------------
    # INPUT
    # L: Matrix containing the model vectors derived inverting the spatial
    #    vectors of the matrix U of the decomposition. Size: 2PxN
    # cov_L: Cell variable containing the covariance matrix of the model
    #        vectors. The n-th element of the cell contains the covariance
    #        matrix of the n-th model vector. Cell size: Nx1. Each cell
    #        contains a matrix of size: 2Px2P
    # S: Diagonal matrix containing the weights associated to the inverted
    #    components. Size: NxN
    # V: Matrix containing the temporal information of the inverted components.
    #    Size: TxN
    # var_V: Matrix containing the variance associated to the matrix V. The
    #        values are stored in a simple matrix because we neglect the
    #        covariances through time. It is thus easier to save the variances
    #        in a single matrix instead of having several diagonal matrices of
    #        size TxT in a cell variable of size Nx1. Size: TxN.
    # flag_compute_cov: Since the calculation of the cov_slip variable is
    #                   demanding from a computation point of view, the user
    #                   has the option to skip this calculation (flag must be 0
    #                   in this case). If the flag is different from 0, then
    #                   the function calculates the cov_slip variable. Default
    #                   value: 0
    # -------------------------------------------------------------------------
    # OUTPUT
    # slip: Matrix containing the slip value for every patch and epoch. Size:
    #       PxT
    # cov_slip: Cell variable containing the covariance matrix of the slip. It
    #           is composed of T cells, each containing a matrix of size 2Px2P.
    # var_slip: Matrix containing the variance associated to the slip value for
    #           every patch and epoch. Size: PxT
    # -------------------------------------------------------------------------
    # References:
    # Gualandi et al., Tectonophysics, 2016b. Pre- and post-seismic deformation
    # related to the 2015, Mw 7.8 Gorkha Earthquake, Nepal
    # 
    # Adriano Gualandi - 19 Oct 2016
    # California Institute of Technology
    # Geological and Planetary Science Division
    # Adapted to Julia by Adriano Gualandi - 04 Nov 2024
    # University of Cambridge, UK
    # Department of Earth Sciences, Bullard Laboratories

    ########################
    ### USEFUL VARIABLES ###
    ########################
    # Extract the number of patches P
    P = Int(size(L,1)/2);
    # Extract the number of components inverted N and the number of epochs T
    N,T = size(transpose(V));
    # Extract the spatial slip component L relative to the strike and dip
    L_strike = L[1:P,:];
    L_dip    = L[P+1:2*P,:];
    # Initialize the covariances and variances variables to speed up the loop
    cov_L_strike_L_strike = zeros(N,P,P)
    cov_L_strike_L_dip = zeros(N,P,P)
    cov_L_dip_L_dip = zeros(N,P,P)
    var_L_strike = zeros(P,N);
    var_L_dip    = zeros(P,N);
    # For every inverted component...
    for nn=1:N
        # Calculate the covariance matrices between the strike and strike
        # directions, between the strike and dip directions, and between the
        # dip and dip directions
        cov_L_strike_L_strike[nn,:,:] = cov_L[nn][1:P,1:P];
        cov_L_strike_L_dip[nn,:,:]    = cov_L[nn][P+1:2*P,1:P];
        cov_L_dip_L_dip[nn,:,:]       = cov_L[nn][P+1:2*P,P+1:2*P];
        # To speed up future calculation, store the variance in a simple matrix
        # instead of being the diagonal of matrices in different cells
        var_L_strike[:,nn]  = diag(cov_L_strike_L_strike[nn,:,:]);
        var_L_dip[:,nn]     = diag(cov_L_dip_L_dip[nn,:,:]);
    end

    ############
    ### SLIP ###
    ############
    # Calculate the slip along the strike and dip directions (eq. 2) and the
    # final slip (eq. 1). A matrix notation is used to speed up and avoid a
    # loop through the patches.
    slip_strike = L_strike*S*transpose(V);
    slip_dip    = L_dip*S*transpose(V);
    slip = sqrt.(slip_strike.^2 .+ slip_dip.^2);

    #########################
    ### COVARIANCE MATRIX ###
    #########################
    # If flag_compute_cov is set to 0 then set cov_slip to NaN; otherwise
    # compute the cov_slip variable
    if flag_compute_cov == 0
        cov_slip     = NaN;
    else
        # Initialize the variables containing the set of random variables
        # defined by (eq. 4) and the gradient of the slip with the respect to
        # these random variables used in (eq. 3).
        # See also equation (S19) of Gualandi et al. (2016b) for the definition
        # of the gradient of the slip on the patch p at time t. That equation
        # has a misprint: the last N elements should have L_strike_{p1},
        # L_dip_{p1}, ..., L_strike_{pN}, L_dip_{pN} instead of a dependence on
        # time. Furthermore, here I'm using the index n=1,...,N for the number
        # of components, while in Gualandi et al. (2016b) I used r=1,...,R. I
        # prefer now to use N to avoid confusion with the index r in the couple
        # (q,r) that goes from 1 to P.
        set_rvs = [[] for i=1:P, j=1:T]
        #der_slip_set_rvs = [zeros(3*N,N) for i=1:P, j=1:T]
        der_slip_set_rvs = zeros(P,T,3*N)
        # For every epoch...
        for tt=1:T
            # And for every patch...
            for pp=1:P
                # Define the set of random variables of (eq. 4)
                set_rvs[pp,tt] = vcat(L_strike[pp,:], L_dip[pp,:], V[tt,:]);
                # Find the gradient of the slip with respect to those random
                # variables (eq. (S19) of Gualandi et al., 2016b)
                der_slip_set_rvs[pp,tt,:] = vcat(
                    (diag(S).*V[tt,:])*(slip_strike[pp,tt]/slip[pp,tt]),
                    (diag(S).*V[tt,:])*(slip_dip[pp,tt]/slip[pp,tt]),
                    (diag(S)) .* (((L_strike[pp,:]*slip_strike[pp,tt] .+
                    L_dip[pp,:]*slip_dip[pp,tt])/slip[pp,tt])));
            end
        end
        # It remains to find the covariance matrix of the set of random
        # variables defined above. Initialize the output cov_slip variable, and
        # the variable of the covariance matrix of the set of random variables.
        cov_slip = zeros(T,P,P)
        cov_set_rvs = zeros(T,P,P,3*N,3*N)
        # Initialize also a variable to keep track of the progressing time
        elapsed_time = zeros(T);
        # # Construct a ParforProgressbar object
        # progress = ProgressMeter.Progress(
        #     T; desc = "Calculating slip covariance", enabled = true
        #     )
        # For every epoch...
        for tt=1:T
            # Start the clock
            elapsed_time[tt] = @elapsed begin
                # Initialize the matrix cov_slip[tt] for epoch tt
                #cov_slip[tt,:,:] = zeros(P,P);
                # Initialize the matrix cov_set_rvs[tt] for epoch tt
                #cov_set_rvs[tt] = cell(P,P);
                # For every patch...
                for qq=1:P
                    # And again for every patch...
                    for rr=1:P
                        # The covariance matrix that we are building here is the
                        # one expressed by equation (S17) of the Supplementary
                        # Material of Gualandi et al. (2016b). This is a 3Nx3N
                        # matrix that must be created for every couple of
                        # patches q and r, and for every epoch t. We can divide
                        # the matrix in 9 blocks: top-left, top-center,
                        # top-right, middle-left, middle-center, middle-right,
                        # bottom-left, bottom-center, and bottom-right. We
                        # initialize such 3Nx3N matrix to zeros to speed up the
                        # calculations. Initialize to zero the first N elements
                        # of the nn row, in order to create the top-left block
                        # matrix, corresponding to the covariances between
                        # L_strike_q(t) and L_strike_r(t) for different inverted
                        # components nn.
                        
                        #cov_set_rvs[tt][qq,tt] = zeros(3*N,3*N);
                        # TOP BLOCK MATRICES
                        # For every inverted component...
                        for nn=1:N
                            # cov_set_rvs{tt}{qq,rr}(nn,1:N) = zeros(1,N);
                            # Only the diagonal values of the top-left block are
                            # different from zero, and their value is given by the
                            # covariance for the component nn between L_strike_q(t)
                            # and L_strike_r(t).
                            # cov_set_rvs[tt][qq,rr][nn,nn] = 
                            #     cov_L_strike_L_strike[nn][qq,rr];
                            cov_set_rvs[tt,qq,rr,nn,nn] = 
                                cov_L_strike_L_strike[nn,qq,rr];
                            # cov_set_rvs{tt}{qq,rr}(nn,N+1:2*N) = zeros(1,N);
                            # Also the top-central block matrix is diagonal, with
                            # diagonal values given by the covariance for the
                            # component nn between L_strike_q(t) and L_dip_r(t).
                            cov_set_rvs[tt,qq,rr,nn,N+nn] = 
                                cov_L_strike_L_dip[nn,qq,rr];
                            # The remaining top-right matrix is zero because we
                            # assume that there is no correlation between the
                            # spatial and temporal contributions
                            # cov_set_rvs[tt][qq,rr][nn,2*N+1:3*N] = zeros(1,N);
                        end
                        # MIDDLE BLOCK MATRICES
                        # For every inverted component...
                        for nn=1:N
                            # cov_set_rvs[tt][qq,rr][N+nn,1:N] = zeros(1,N);
                            # Only the diagonal values of the central-left block
                            # are different from zero, and their value is given by
                            # the covariance for the component nn between
                            # L_strike_q(t) and L_dip_r(t).
                            #cov_set_rvs[tt][qq,rr][N+nn,nn] = cov_L_strike_L_dip[nn][rr,qq];
                            cov_set_rvs[tt,qq,rr,N+nn,nn] = cov_L_strike_L_dip[nn,rr,qq];
                            # cov_set_rvs[tt][qq,rr][N+nn,N+1:2*N] = zeros(1,N);
                            # Also the middle-central block matrix is diagonal,
                            # with diagonal values given by the covariance for the
                            # component nn between L_dip_q(t) and L_dip_r(t).
                            #cov_set_rvs[tt][qq,rr][N+nn,N+nn] = cov_L_dip_L_dip[nn][qq,rr];
                            cov_set_rvs[tt,qq,rr,N+nn,N+nn] = cov_L_dip_L_dip[nn,qq,rr];
                            # The remaining middle-right matrix is zero because we
                            # assume that there is no correlation between the
                            # spatial and temporal contributions
                            # cov_set_rvs[tt][qq,rr][N+nn,2*N+1:3*N] = zeros(1,N);
                        end
                        # BOTTOM BLOCK MATRICES
                        # For every inverted component...
                        for nn=1:N
                            # The only non zero block among the bottom block
                            # matrices is the bottom-right one because we are
                            # assuming that there is no correlation between the
                            # spatial and temporal contributions. The bottom-right
                            # block is diagonal because the temporal components are
                            # supposedly independent one from the other, and thus
                            # they are uncorrelated. The diagonal values are given
                            # by the variance of V(t) for the nn component.
                            # cov_set_rvs[tt][qq,rr][2*N+nn,1:N] = zeros(1,N);
                            # cov_set_rvs[tt][qq,rr][2*N+nn,N+1:2*N] = zeros(1,N);
                            # cov_set_rvs[tt][qq,rr][2*N+nn,2*N+1:3*N] = zeros(1,N);
                            #cov_set_rvs[tt][qq,rr][2*N+nn,2*N+nn] = var_V[tt,nn];
                            cov_set_rvs[tt,qq,rr,2*N+nn,2*N+nn] = var_V[tt,nn];
                        end
                        # Finally, the element (q,r) of the PxP matrix cov_slip{tt}
                        # is obtained using equation (S16) of the Supplementary
                        # Material of Gualandi et al. (2016b), i.e. (eq. 3) in the
                        # description of this function
                        # cov_slip[tt][qq,rr] = 
                        #     transpose(der_slip_set_rvs[qq,tt])*
                        #     cov_set_rvs[tt][qq,rr]*der_slip_set_rvs[rr,tt];
                        # println(size(transpose(der_slip_set_rvs[qq,tt,:])))
                        # println(size(cov_set_rvs[tt,qq,rr,:,:]))
                        # println(size(der_slip_set_rvs[rr,tt,:]))
                        cov_slip[tt,qq,rr] = 
                            transpose(der_slip_set_rvs[qq,tt,:]) *
                            cov_set_rvs[tt,qq,rr,:,:] * der_slip_set_rvs[rr,tt,:];
                    end
                end
                # After iterating on both qq and rr we have the full covariance
                # matrix between the slip at patch q and the slip at patch r at a
                # given time t.
                # cov_set_rvs_filename = [dirs.dir_scen,'/',dirs.dir_case,'/cov_set_rvs/cov_set_rvs_',num2str(tt),'.mat'];
                # cov_set_rvs_tt = cov_set_rvs{tt};
                # save(cov_set_rvs_filename,'cov_set_rvs_tt');
                # Keep track of the progress
                #elapsed_time(tt) = toc;
            end
            println("tt = "*string(tt)*"/"*string(T)*" - Partial time: "*
                string(elapsed_time[tt])*" s,      Total time: "*
                string(sum(elapsed_time))*" s");
        end
    end

    ################
    ### VARIANCE ###
    ################
    # Using (eq. 5) we can iterate over the number of inverted components N,
    # and using the matrix notation we can avoid the loop over the epochs T.
    var_slip          = zeros(P,T)
    var_slip_strike   = zeros(P,T)
    var_slip_dip      = zeros(P,T)

    # Construct a ParforProgressbar object
    progress = ProgressMeter.Progress(
        N; desc = "Calculating slip variance", enabled = true
        )
    # It seems to be faster if not in parallel...
    #Threads.@threads for nn=1:N
    for nn=1:N
        # Calculate partiald(slip_p(t)/L_strike_p^n) for all the patches and
        # epochs. This is thus a matrix PxT.
        der_slip_L_strike = 
            (slip_strike*(S[nn,nn]).*repeat(transpose(V[:,nn]),P,1))./slip;
        # Calculate partiald(slip_p(t)/L_dip_p^n) for all the patches and
        # epochs. This is thus a matrix PxT.
        der_slip_L_dip = 
            (slip_dip*(S[nn,nn]).*repeat(transpose(V[:,nn]),P,1))./slip;
        # Calculate partiald(slip_p(t)/V^n(t)) for all the patches and epochs.
        # This is thus a matrix PxT.
        der_slip_V = S[nn,nn]*(
            (slip_strike.*repeat(L_strike[:,nn],1,T)) .+ 
            (slip_dip.*repeat(L_dip[:,nn],1,T)))./slip;
        # Calculate the variance related to the component nn. The calculation
        # is performed simultaneously for all the patches and epochs. This is
        # thus a matrix PxT.
        var_slip_comp = 
            (der_slip_L_strike.^2).*repeat(var_L_strike[:,nn],1,T) .+ 
            (der_slip_L_dip.^2).*repeat(var_L_dip[:,nn],1,T) .+
            (der_slip_V.^2).*repeat(transpose(var_V[:,nn]),P,1) .+
            2*((der_slip_L_strike).*(der_slip_L_dip)).*
            repeat(diag(cov_L_strike_L_dip[nn,:,:]),1,T);
        # var_slip_comp{nn}(:,1) = var_slip_comp{nn}(:,2);
        # Calculate the final variance on the slip summing up the contributions
        # from all the inverted components.
        var_slip = var_slip + var_slip_comp;
        var_slip_comp = nothing

        var_slip_strike_comp = S[nn,nn].^2 * (
            var_L_strike[:,nn,:]*var_V[:,nn,:]' +
            var_L_strike[:,nn,:]*(V[:,nn,:].^2)' +
            L_strike[:,nn,:].^2 * var_V[:,nn,:]' );
        var_slip_dip_comp = S[nn,nn].^2 * (
            var_L_dip[:,nn,:]*var_V[:,nn,:]' +
            var_L_dip[:,nn,:]*(V[:,nn,:].^2)' +
            L_dip[:,nn,:].^2 * var_V[:,nn,:]' );
        var_slip_strike = var_slip_strike + var_slip_strike_comp;
        var_slip_dip = var_slip_dip + var_slip_dip_comp;
        var_slip_strike_comp = nothing
        var_slip_dip_comp = nothing

        ProgressMeter.next!(progress)
    end

    return slip, cov_slip, var_slip, slip_strike, slip_dip, var_slip_strike,
        var_slip_dip
end

function calc_slip_potency(slip)
    slip_potency = Dict()
    slip_potency["timeline"] = slip["timeline"]
    slip_potency["obs_strike"] = slip[]
end