

function extract_optimizer_moment_matrix(V::Matrix{Float64},lu0::Int64,basis_sigma0::Matrix{UInt64},n::Int64,l::Int64,opt_val::Float64,f::Polynomial{true},h::Vector{Polynomial{true,Float64}},x::Vector{PolyVar{true}})
    #extraction of optimizers by moment matrix
    
    V=rref_with_pivots!(Matrix(V'),1e-2)
    U=Matrix(V[1]')
    rk=size(U,2)

    println("Rank of moment matrix = ", rk)
    
    
    sol=Vector{Float64}[]
    if rk >0
        # extraction of Henrion and Lasserre

        w=basis_sigma0[:,V[2]]
        N=Vector{Array{Float64}}(undef,n)
        flag=UInt8(1)
 
        
        for i in 1:n
            @inbounds kk=UInt16[]
            for j=1:size(w)[2]
                @inbounds xwj=w[:,j]
                @inbounds xwj[i]+=1
                for t=1:lu0
                    if xwj==basis_sigma0[:,t]
                        if kk==UInt16[]
                            @inbounds kk=t
                        else
                            @inbounds kk=[kk;t]
                        end
                        break
                    end
                end
            end
            
            if rk!=length(kk)
                @inbounds flag=0
                break
            else
                @inbounds N[i]=U[kk,:]
            end
        end

        if flag==1

            # Create random convex combination
            rands = rand(n);rands = rands/sum(rands)
            M = zeros(length(V[2]),rk)
            println(size(M))
            for i=1:n
                println(size(N[i]))
                @inbounds M+=rands[i]*N[i]
            end

            F= schur(M);
            L=F.Z
            # Extract solution
            for i=1:rk
                @inbounds atom=Float64[]
                for j = 1:n
                    @inbounds coordinatej=L[:,i]'*N[j]*L[:,i]
                    @inbounds push!(atom,coordinatej[1])
                end
                println("------------------------------------")
                #println("atom ",i," = ",atom)
                println("atom ",i,":")
                @inbounds flag=1
                
                @inbounds check=abs(polynomial(f)(x => atom)-opt_val)#/minimum(abs.(coefficients(f)))
                
                #println("  check gap of lower bound  = ",check)
                if abs(check)>1e-1
                    @inbounds flag=0
                end


                for i=1:l
                    @inbounds check=abs(polynomial(h[i])(x => atom))#/minimum(abs.(coefficients(h[i])))
                    #println("  check equality constraint ",i," = ",check)
                    if abs(check)>1e-1
                        @inbounds flag=0
                    end
                end

                if flag ==1
                    @inbounds sol=atom
                    println("####################################")
                    #println("Optimal solution: opt_sol = ",atom)
                    println("It is an approximate optimal solution!")
                    println("####################################")
                else
                    println("It is not an approximate optimal solution!")
                end

            end
        end
    end

    return sol

end




function extract_optimizer(Gr::Matrix{Float64},lu0::Int64,basis_sigma0::Matrix{UInt64},n::Int64,l::Int64,opt_val::Float64,f::Polynomial{true},h::Vector{Polynomial{true,Float64}},x::Vector{PolyVar{true}})
    #extraction of optimizers by Gram matrix
    V=nullspace(Gr,atol=1e-5)
    rk=size(V,2)

    println("Dimension of the null space of Gram matrix = ", rk)
    sol=Vector{Float64}[]
    if rk >0
        # extraction of Henrion and Lasserre
        V= rref_with_pivots!(Matrix(V'),1e-2)
        U=Matrix(V[1]')
        w=basis_sigma0[:,V[2]]
        N=Vector{Array{Float64}}(undef,n)
        flag=UInt8(1)
 
        
        for i in 1:n
            @inbounds kk=UInt16[]
            for j=1:size(w)[2]
                @inbounds xwj=w[:,j]
                @inbounds xwj[i]+=1
                for t=1:lu0
                    if xwj==basis_sigma0[:,t]
                        if kk==UInt16[]
                            @inbounds kk=t
                        else
                            @inbounds kk=[kk;t]
                        end
                        break
                    end
                end
            end
            
            if rk!=length(kk)
                @inbounds flag=0
                break
            else
                @inbounds N[i]=U[kk,:]
            end
        end

        if flag==1

            # Create random convex combination
            rands = rand(n);rands = rands/sum(rands)
            M = zeros(length(V[2]),rk)
            for i=1:n
                @inbounds M+=rands[i]*N[i]
            end

            F= schur(M);
            L=F.Z
            # Extract solution
            for i=1:rk
                @inbounds atom=Float64[]
                for j = 1:n
                    @inbounds coordinatej=L[:,i]'*N[j]*L[:,i]
                    @inbounds push!(atom,coordinatej[1])
                end
                println("------------------------------------")
                #println("atom ",i," = ",atom)
                println("atom ",i,":")
                @inbounds flag=1
                
                @inbounds check=abs(polynomial(f)(x => atom)-opt_val)#/minimum(abs.(coefficients(f)))
                
                #println("  check gap of lower bound  = ",check)
                if abs(check)>1e-1
                    @inbounds flag=0
                end


                for i=1:l
                    @inbounds check=abs(polynomial(h[i])(x => atom))#/minimum(abs.(coefficients(h[i])))
                    #println("  check equality constraint ",i," = ",check)
                    if abs(check)>1e-1
                        @inbounds flag=0
                    end
                end

                if flag ==1
                    @inbounds sol=atom
                    println("####################################")
                    #println("Optimal solution: opt_sol = ",atom)
                    println("It is an approximate optimal solution!")
                    println("####################################")
                else
                    println("It is not an approximate optimal solution!")
                end

            end
        end
    end

    return sol

end





function extract_optimizer2(Gr::Matrix{Float64},lu0::Int64,basis_sigma0::Matrix{UInt64},n::Int64,l_g::Int64,l_h::Int64,opt_val::Float64,f::Polynomial{true},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},x::Vector{PolyVar{true}})
    #extraction of optimizers
    V=nullspace(Gr,atol=1e-5)
    rk=size(V,2)

    println("Dimension of the null space of Gram matrix = ", rk)
    sol=Vector{Float64}[]
    if rk >0
        # extraction of Henrion and Lasserre
        V= rref_with_pivots!(Matrix(V'),1e-2)
        U=Matrix(V[1]')

        w=basis_sigma0[:,V[2]]
        N=Vector{Array{Float64}}(undef,n)
        flag=UInt8(1)
 
        
        for i in 1:n
            @inbounds kk=UInt16[]
            for j=1:size(w)[2]
                @inbounds xwj=w[:,j]
                @inbounds xwj[i]+=1
                for t=1:lu0
                    if xwj==basis_sigma0[:,t]
                        if kk==UInt16[]
                            @inbounds kk=t
                        else
                            @inbounds kk=[kk;t]
                        end
                        break
                    end
                end
            end
            
            if rk!=length(kk)
                @inbounds flag=0
                break
            else
                @inbounds N[i]=U[kk,:]
            end
        end

        if flag==1

            # Create random convex combination
            rands = rand(n);rands = rands/sum(rands)
            M = zeros(length(V[2]),rk)
            for i=1:n
                @inbounds M+=rands[i]*N[i]
            end

            F= schur(M);
            L=F.Z
            # Extract solution
            for i=1:rk
                @inbounds atom=Float64[]
                for j = 1:n
                    @inbounds coordinatej=L[:,i]'*N[j]*L[:,i]
                    @inbounds push!(atom,coordinatej[1])
                end
                println("------------------------------------")
                #println("atom ",i," = ",atom)
                println("atom ",i,":")
                @inbounds flag=1
                
                @inbounds check=polynomial(f)(x => atom)-opt_val
                
                #println("  check gap of lower bound  = ",check)
                if abs(check)>1e-1
                    @inbounds flag=0
                end


                for i=1:l_g
                    @inbounds check=polynomial(g[i])(x => atom)
                    #println("  check inequality constraint ",i," = ",check)
                    if check<-1e-1
                        @inbounds flag=0
                    end
                end
                
                for i=1:l_h
                    @inbounds check=polynomial(h[i])(x => atom)
                    #println("  check equality constraint ",i," = ",check)
                    if abs(check)>1e-1
                        @inbounds flag=0
                    end
                end

                if flag ==1
                    @inbounds sol=atom
                    println("####################################")
                    #println("Optimal solution: opt_sol = ",atom)
                    println("It is an approximate optimal solution!")
                    println("####################################")
                else
                    println("It is not an approximate optimal solution!")
                end

            end
        end
    end

    return sol

end
