function ASC_PolySys(x,h,k,L;method="LMBM",EigAlg="Mix",tol=1e-3)
    @time begin
    @polyvar x_slack
    n=length(x)
    l=length(h)
    bar_h=Vector{Polynomial{true,Float64}}(undef,l)
    for j in 1:l
        bar_h[j]=h[j](x => x.*sqrt(L))/L
    end

    sol=SpectralPOP.ASC([x;x_slack],[1.0-sum([x;x_slack].^2);bar_h],k,method=method,tol=tol,EigAlg=EigAlg)
    
    if sol!=[]
        sol=sol[1:n].*sqrt(L)
    end
    end
    return sol
end



function ASC(x,h,k;method="LMBM",EigAlg="Mix",tol=1e-3)
    @time begin
    n=length(x)
    l=length(h)
    
    ch=deepcopy(h)
    
    # Define centers and square of radius
    a = Matrix{Float64}(LinearAlgebra.I, n, n)
    omega = Vector{Float64}(undef,n)
    f=sum(x.^2)+.0
    
    sol=Vector{Float64}[]
    
    t=UInt64(1)
    while t <=n
        println("------------------------------------")
        println("Determine omega",t,":")
        if t>1
            @inbounds push!(ch,omega[t-1]-f)
        end        

        @inbounds f=sum((x-a[:,t]).^2)
        @inbounds omega[t],sol= CTP_POP(x,f,ch,k,1.0;method=method,EigAlg=EigAlg,tol=tol)
        #@inbounds omega[t],sol=nonsmooth_hierarchy_old(x,f,ch,k,method=method,tol=tol)
        println("omega",t," = ", omega[t])
        println("------------------------------------")

        if sol!=[]
            break
        elseif abs(omega[t])<1e-5
                sol=a[:,t]
                break
        end
        t+=1    
    end
    
    if t==n+1
        flag=UInt8(1)
        atom=(-omega./2).+1
        println("atom=",atom)
        for i=1:l
            check=polynomial(h[i])(x => atom)
            println("check equality ",i," = ",check)
            if abs(check)>1e-3
                flag=0
            end
        end
        if flag==1
            sol=atom
        end
    end
    println("------")
    end
    
    return sol
end






#=function ASC_SumOfSquares(x,h,k)
    @time begin
    n=length(x)
    l=length(h)
    
    ch=deepcopy(h)
    
    # Define centers and square of radius
    a = Matrix{Float64}(LinearAlgebra.I, n, n)
    omega = Vector{Float64}(undef,n)
    f=sum(x.^2)+.0
    
    sol=Vector{Float64}[]
    g=Vector{Polynomial{true,Float64}}([])
    
    t=UInt64(1)
    while t <=n
        println("------------------------------------")
        println("Determine omega",t,":")
        if t>1
            @inbounds push!(ch,omega[t-1]-f)
        end        

        @inbounds f=sum((x-a[:,t]).^2)

        @inbounds omega[t],sol=SumofSquares_POP_WithExtraction(x,f,g,ch,k)
        println("omega",t," = ", omega[t])
        println("------------------------------------")

        if sol!=[]
            break
        end
        t+=1    
    end
    
    if t==n+1
        flag=UInt8(1)
        atom=(-omega./2).+1
        println("atom=",atom)
        for i=1:l
            check=polynomial(h[i])(x => atom)
            println("check equality ",i," = ",check)
            if abs(check)>1e-3
                flag=0
            end
        end
        if flag==1
            sol=atom
        end
    end
    end
    return sol
end=#