function get_info(x::Vector{PolyVar{true}},f::Polynomial{true,Float64},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}})
    n=length(x)
    m=length(g)
    l=length(h)
            
   
    
    lmon_g=Vector{UInt64}(undef,m)
    coe_g=Vector{Vector{Float64}}(undef,m)
    
    lmon_h=Vector{UInt64}(undef,l)
    coe_h=Vector{Vector{Float64}}(undef,l)
    
    
    supp_g=Vector{SparseMatrixCSC{UInt64}}(undef,m)
    supp_h=Vector{SparseMatrixCSC{UInt64}}(undef,l)

    lmon_f,supp_f,coe_f=info(f,x,n,sparse=true)
    
    for i in 1:m
        lmon_g[i],supp_g[i],coe_g[i]=info(g[i],x,n,sparse=true)
        
    end
                    
    for i in 1:l
        lmon_h[i],supp_h[i],coe_h[i]=info(h[i],x,n,sparse=true)
    end
    return n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f
end




function save_data(x,f,g,h,R)
    
    n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f=get_info(x,f,g,h)
    
    d=maximum([sum(supp_f[:,i]) for i in 1:lmon_f])
    
    output_file = open("/home/hoanganh/Desktop/math-topics/SpectralPOP/codes/dataPOP/densePOPsphere_deg$(d)_var$(n)_nineq$(m)_neq$(l).jl","w")
    
    write(output_file, "n = ")
    show(output_file, n)
    write(output_file, "; \n \n")
    
    write(output_file, "m = ")
    show(output_file, m)
    write(output_file, "; \n \n")
    
    write(output_file, "l = ")
    show(output_file, l)
    write(output_file, "; \n \n")
    
    write(output_file, "lmon_g = ")
    show(output_file, lmon_g)
    write(output_file, "; \n \n")
    
    I_g=Vector{Vector{UInt64}}(undef,m)
    J_g=Vector{Vector{UInt64}}(undef,m)
    V_g=Vector{Vector{UInt64}}(undef,m)
    
    for i in 1:m
        I_g[i],J_g[i],V_g[i]=findnz(supp_g[i])
    end
    
    write(output_file, "I_g = ")
    show(output_file, I_g)
    write(output_file, "; \n \n")
    
    write(output_file, "J_g = ")
    show(output_file, J_g)
    write(output_file, "; \n \n")
    
    write(output_file, "V_g = ")
    show(output_file, V_g)
    write(output_file, "; \n \n")
        
    
    write(output_file, "supp_g = Vector{SparseMatrixCSC{UInt64}}([sparse(I_g[i],J_g[i],V_g[i],n,lmon_g[i]) for i in 1:m])")
    write(output_file, "; \n \n")
    
    write(output_file, "coe_g = ")
    show(output_file, coe_g)
    write(output_file, "; \n \n")
    
    write(output_file, "lmon_h = ")
    show(output_file, lmon_h)
    write(output_file, "; \n \n")
    
    I_h=Vector{Vector{UInt64}}(undef,l)
    J_h=Vector{Vector{UInt64}}(undef,l)
    V_h=Vector{Vector{UInt64}}(undef,l)
    
    for i in 1:l
        I_h[i],J_h[i],V_h[i]=findnz(supp_h[i])
    end
    
    write(output_file, "I_h = ")
    show(output_file, I_h)
    write(output_file, "; \n \n")
    
    write(output_file, "J_h = ")
    show(output_file, J_h)
    write(output_file, "; \n \n")
    
    write(output_file, "V_h = ")
    show(output_file, V_h)
    write(output_file, "; \n \n")
        
    
    write(output_file, "supp_h = Vector{SparseMatrixCSC{UInt64}}([sparse(I_h[i],J_h[i],V_h[i],n,lmon_h[i]) for i in 1:l])")
    write(output_file, "; \n \n")
    
    write(output_file, "coe_h = ")
    show(output_file, coe_h)
    write(output_file, "; \n \n")
    
    write(output_file, "lmon_f = ")
    show(output_file, lmon_f)
    write(output_file, "; \n \n")
    
    I_f,J_f,V_f=findnz(supp_f)
    
    write(output_file, "I_f = ")
    show(output_file, I_f)
    write(output_file, "; \n \n")
    
    write(output_file, "J_f = ")
    show(output_file, J_f)
    write(output_file, "; \n \n")
    
    write(output_file, "V_f = ")
    show(output_file, V_f)
    write(output_file, "; \n \n")
        
    
    write(output_file, "supp_f = sparse(I_f,J_f,V_f,n,lmon_f)")
    write(output_file, "; \n \n")
    
    write(output_file, "coe_f = ")
    show(output_file, coe_f)
    write(output_file, "; \n \n")
    
    
    write(output_file, "R = ")
    show(output_file, R)
    write(output_file, "; \n \n")
    
    close(output_file)
    
end


function save_deciding_copos(n,A)
    
    
    output_file = open("/home/hoanganh/Desktop/math-topics/SpectralPOP/codes/dataPOP/deciding_copos_size$(n).jl","w")
    
    write(output_file, "n = ")
    show(output_file, n)
    write(output_file, "; \n \n")
    
    write(output_file, "A = ")
    show(output_file, A)
    write(output_file, "; \n \n")
    
    close(output_file)
    
end

function save_deciding_conv(n,d,c)
    
    
    output_file = open("/home/hoanganh/Desktop/math-topics/SpectralPOP/codes/dataPOP/deciding_conv_var$(n)deg$(d).jl","w")
    
    write(output_file, "n = ")
    show(output_file, n)
    write(output_file, "; \n \n")
    
    write(output_file, "d = ")
    show(output_file, d)
    write(output_file, "; \n \n")
    
    write(output_file, "c = ")
    show(output_file, c)
    write(output_file, "; \n \n")
    
    close(output_file)
    
end

function save_deciding_nonneg(n,d,c)
    
    
    output_file = open("/home/hoanganh/Desktop/math-topics/SpectralPOP/codes/dataPOP/deciding_nonneg_var$(n)deg$(d).jl","w")
    
    write(output_file, "n = ")
    show(output_file, n)
    write(output_file, "; \n \n")
    
    write(output_file, "d = ")
    show(output_file, d)
    write(output_file, "; \n \n")
    
    write(output_file, "c = ")
    show(output_file, c)
    write(output_file, "; \n \n")
    
    close(output_file)
    
end

function get_POP(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f)
    
    @polyvar x[1:n]

    f=get_poly(lmon_f, supp_f, coe_f,x,n)
    
    g=Vector{Polynomial{true,Float64}}(undef,m)
    
    for j in 1:m
        g[j]=get_poly(lmon_g[j], supp_g[j], coe_g[j],x,n)
    end
    
    h=Vector{Polynomial{true,Float64}}(undef,l)
    
    for j in 1:l
        h[j]=get_poly(lmon_h[j], supp_h[j], coe_h[j],x,n)
    end
    
    return x,f,g,h
    
end


function get_poly(lmon, supp, coe,x,n)
    f=Polynomial{true,Float64}(1.0*x[1]-x[1])
    for j=1:lmon
        f+=coe[j]*prod(x[i]^supp[i,j] for i=1:n)
    end
    return f
end
    



function info(f,x::Vector{PolyVar{true}},n::Int64;sparse=false)
    
    mon=monomials(f)
    coe=coefficients(f)
    lmon=length(mon)
    if sparse==false
        supp=zeros(UInt64,n,lmon)
    else
        supp=spzeros(UInt64,n,lmon)
    end
    @simd for i in 1:lmon
        @simd for j in 1:n
            @inbounds supp[j,i]=DynamicPolynomials.degree(mon[i],variable(x[j]))
        end
    end
    return lmon, supp, coe
end

function get_basis(n::Int64,d::Int64)
    
    lb=binomial(n+d,d)
    basis=zeros(UInt64,n,lb)
    i=UInt64(0)
    t=UInt64(1)
    while i<d+1
        if basis[n,t]==i
           if i<d
              @inbounds t+=1
              @inbounds basis[1,t]=i+1
              @inbounds i+=1
           else 
                @inbounds i+=1
           end
        else 
            j=UInt64(1)
             while basis[j,t]==0
                   @inbounds j+=1
             end
             if j==1
                @inbounds t+=1
                @inbounds basis[:,t]=basis[:,t-1]
                @inbounds basis[1,t]=basis[1,t]-1
                @inbounds basis[2,t]=basis[2,t]+1
                else t+=1
                  @inbounds basis[:,t]=basis[:,t-1]
                  @inbounds basis[1,t]=basis[j,t]-1
                  @inbounds basis[j,t]=0
                  @inbounds basis[j+1,t]=basis[j+1,t]+1
             end
        end
    end
    return basis
end


function sym(lmon_f::Int64,supp_f::Vector{Int64},coe_f::Vector{Float64},polymat::SparseMatrixCSC{Int64},size_mommat::Int64)
    
    A=spzeros(Float64,size_mommat,size_mommat)
    r=UInt64(1)
    
    for i in 1:size_mommat, j in i:size_mommat
        r=findfirst(x->x==polymat[i,j],supp_f)

        if r!=nothing
            if i==j
                @inbounds A[i,j]=coe_f[r]
            else
                @inbounds A[i,j]=.5*coe_f[r]
                @inbounds A[j,i]=.5*coe_f[r]
            end
            @inbounds lmon_f-=1
            if lmon_f==0
                break
            end
            @inbounds supp_f=supp_f[1:end .!= r]
            @inbounds coe_f=coe_f[1:end .!= r]                
        end
        
    end

return A
end


function semisym(lmon_f::Int64,supp_f::Vector{Int64},coe_f::Vector{Float64},Findinpolymat::Vector{Int64},dim_mom::Int64,svec2::Vector{Float64})
    
    a=spzeros(Float64,dim_mom)
   
    
    @simd for r in 1:lmon_f
        #@inbounds a[Findinpolymat[supp_f[r]]]=-coe_f[r]*svec2[Findinpolymat[supp_f[r]]]
        @inbounds t=Findinpolymat[supp_f[r]]
        @inbounds a[t]=-coe_f[r]*svec2[t]
                    
    end
    return a

end


function semisym_dense(lmon_f::Int64,supp_f::Vector{Int64},coe_f::Vector{Float64},Findinpolymat::Vector{Int64},dim_mom::Int64,svec2::Vector{Float64})
    
    a=zeros(Float64,dim_mom)
   
    
    @simd for r in 1:lmon_f
        
        @inbounds t=Findinpolymat[supp_f[r]]
        @inbounds a[t]=-coe_f[r]*svec2[t]
                    
    end
    return a

end


function semisym_dense2(lmon_f::Int64,supp_f::Vector{Int64},coe_f::Vector{Float64},Findinpolymat::Vector{Int64},dim_mom::Int64)
    
    a=zeros(Float64,dim_mom)
   
    
    @simd for r in 1:lmon_f
        
        @inbounds t=Findinpolymat[supp_f[r]]
        @inbounds a[t]=-coe_f[r]
                    
    end
    return a

end


function OthantBasicMommat2(polymat::SparseMatrixCSC{Int64},leng_v::Int64,dim_mom::Int64,size_mommat::Int64)

            
    #println("Orthogonal moment basis doesn's exist!")
    A=spzeros(Float64,leng_v,dim_mom)
    t=UInt64(1)
    for i in 1:size_mommat, j in i:size_mommat
        if i==j
            @inbounds A[polymat[i,j],t]=1
        else
            @inbounds A[polymat[i,j],t]=sqrt(2)
        end
        @inbounds t+=1

    end
    N=spzeros(Float64,dim_mom,dim_mom-leng_v)
    t=UInt64(1)
    for r in 1:leng_v
        @inbounds I,V=findnz(A[r,:])
        @inbounds leng_I=length(I)
        if leng_I>1           
            for j in 2:leng_I
                @inbounds N[I[1],t]=V[j]
                @inbounds N[I[j],t]=-V[1]
                @inbounds t+=1
            end
        end
    end
    
        
return N
end


function OthantBasicMommat3(Imagepolymat::SparseMatrixCSC{Float64,Int64},leng_v::Int64,dim_mom::Int64,Inde::Vector{Vector{UInt64}})

            
    #println("Orthogonal moment basis doesn's exist!")
    
    N=spzeros(Float64,dim_mom,dim_mom-leng_v)
    t,I,V,T,leng_I=UInt64(1),Int64(0),Float64(0),UInt64(0),Int64(0)
    @simd for r in 1:leng_v
        @inbounds I,V=findnz(Imagepolymat[r,:])
        @inbounds leng_I=length(I)
        if leng_I>1
            @inbounds T=Inde[I[1]]
            @simd for j in 2:leng_I
                
                @inbounds N[I[1],t]=-V[j]/(sqrt(2)^(1-0^(T[2]-T[1])))
                
                @inbounds T=Inde[I[j]]
                
                @inbounds N[I[j],t]=V[1]/(sqrt(2)^(1-0^(T[2]-T[1])))
               
                @inbounds t+=1
            end
        end
    end
    
    
       
return N
end


function OthantBasicMommat(Imagepolymat::SparseMatrixCSC{Float64,Int64},leng_v::Int64,dim_mom::Int64,invP::Vector{Float64},Inde::Vector{Vector{UInt64}})

            
    #println("Orthogonal moment basis doesn's exist!")
    
    N=spzeros(Float64,dim_mom,dim_mom-leng_v)
    t,I,V,T,leng_I=UInt64(1),Int64(0),Float64(0),UInt64(0),Int64(0)
    @simd for r in 1:leng_v
        @inbounds I,V=findnz(Imagepolymat[r,:])
        @inbounds leng_I=length(I)
        if leng_I>1
            @inbounds T=Inde[I[1]]
            @simd for j in 2:leng_I
                
                @inbounds N[I[1],t]=-V[j]*invP[T[1]]*invP[T[2]]/(sqrt(2)^(1-0^(T[2]-T[1])))
                
                @inbounds T=Inde[I[j]]
                
                @inbounds N[I[j],t]=V[1]*invP[T[1]]*invP[T[2]]/(sqrt(2)^(1-0^(T[2]-T[1])))
               
                @inbounds t+=1
            end
        end
    end
    
    
       
return N
end







function makesym(a::Vector{Float64},invInde::Matrix{UInt64},size_mommat::Int64)

    B=zeros(Float64,size_mommat,size_mommat)
    for i in 1:size_mommat, j in i:size_mommat
        @inbounds B[i,j]=a[invInde[i,j]]
        @inbounds B[j,i]= copy(B[i,j])
    end
    
    return B
end

function makesym2(a::SparseArrays.SparseVector{Float64,Int64},Inde::Vector{Vector{UInt64}},size_mommat::Int64)

    I,V=findnz(a)
    length_V=length(V)
    B=spzeros(Float64,size_mommat,size_mommat)
    for t in 1:length_V
        @inbounds i,j=Inde[I[t]]
        if i==j
            @inbounds B[i,j]=V[t]
        else
            @inbounds B[j,i]=V[t]/sqrt(2)
            @inbounds B[i,j]=B[j,i]
        end
    end
    
    return B
end

function makesemisym(B::Matrix{Float64},invInde::Matrix{UInt64},dim_mom::Int64,s::Int64)
    a=zeros(Float64,dim_mom)
    for i in 1:s, j in i:s               
        @inbounds a[invInde[i,j]]=B[i,j]               
    end   
    return a
end

function makesemisym_sparse(B::SparseMatrixCSC{Float64},invInde::Matrix{UInt64},dim_mom::Int64,s::Int64)
    a=zeros(Float64,dim_mom)
    for i in 1:s, j in i:s               
        @inbounds a[invInde[i,j]]=B[i,j]               
    end   
    return a
end

function invdiaP(supp_theta::Vector{Int64}, coe_theta::Vector{Int64},size_mommat::Int64,diagu::Array{Int64,1})
    
    D=zeros(Float64,size_mommat)
    
    
    @simd for j in 1:size_mommat
        @inbounds i=findfirst(x->x==diagu[j],supp_theta)
        if supp_theta[i]==diagu[j]
            @inbounds D[j]=1/sqrt(coe_theta[i])
        end
    end
    
    return D
end

function invdiaP2(supp_theta::Vector{UInt64}, coe_theta::Vector{Int64},size_mommat::Int64,diagu::Array{UInt64,1})
    
    D=zeros(Float64,size_mommat)
    
    
    @simd for j in 1:size_mommat
        @inbounds i=findfirst(x->x==diagu[j],supp_theta)
        if supp_theta[i]==diagu[j]
            @inbounds D[j]=1/sqrt(coe_theta[i])
        end
    end
    
    return D
end

function bfind(A::Matrix{UInt64},l::Int64,a::Vector{UInt64},n::Int64)
    if l==0
        return 0
    end
    low=UInt64(1)
    high=l
    while low<=high
        @inbounds mid=Int(ceil(1/2*(low+high)))
        @inbounds order=comp(A[:,mid],a,n)
        if order==0
           return mid
        elseif order<0
           @inbounds low=mid+1
        else
           @inbounds high=mid-1
        end
    end
    return 0
end

function comp(a::Vector{UInt64},b::Vector{UInt64},n::Int64)
    i=UInt64(1)
    while i<=n
          if a[i]<b[i]
             return -1
          elseif a[i]>b[i]
             return 1
          else
             @inbounds i+=1
          end
    end
    if i==n+1
       return 0
    end
end








#=projection(ui::SparseVector{Float64},gi::SparseVector{Float64})=gi.*(ui'*gi)

function calculateV(ui::SparseVector{Float64},position::Int64,G::SparseMatrixCSC{Float64})

    @fastmath @inbounds @simd for i=1:position-1
        ui-=projection(ui,G[i,:])
    end
    return ui
end

calculateG(vi::SparseVector{Float64})=vi./norm(vi)

function gramschimidt(b::SparseMatrixCSC{Float64})
	#base = float(b)
    l=size(b,1)
	
    
    b[1,:] = calculateG(b[1,:])
	@simd for i=2:l
		@inbounds b[i,:] = calculateG(calculateV(b[i,:],i,b[1:i-1,:]))
	end	
    
	return b
end

function gramschimidt2(V::SparseMatrixCSC{Float64})
    n = size(V,1)

    V[1,:] = V[1,:]./norm(V[1,:])
    @fastmath @inbounds @simd for i = 2:n
      @fastmath @inbounds @simd for j = 1:i-1
        V[i,:] -= V[j,:].*((V[j,:]'*V[i,:])/norm(V[j,:]))
      end
      V[i,:] = V[i,:]./norm(V[i,:]);
    end
	return V
end=#

