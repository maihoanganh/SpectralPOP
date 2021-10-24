function makevec(B::SparseMatrixCSC{Float64},l_vec::UInt64,sqrt2::Float16,IndVec::SparseMatrixCSC{UInt64})
         b=spzeros(Float64,l_vec)
         I,J,V=findnz(B)
         l_I=length(I)
         @simd for i in 1:l_I
            if I[i]==J[i]
                @inbounds b[IndVec[I[i],J[i]]]=V[i]
            else
                @inbounds b[IndVec[I[i],J[i]]]=V[i]*sqrt2
            end
          end
     
         return b
     end  

function makevec2(B::SparseMatrixCSC{Float64},l_vec::UInt64,sqrt2::Float16,IndVec::SparseMatrixCSC{UInt64})
         b=spzeros(Float64,l_vec)
         I,J,V=findnz(B)
         l_I=length(I)
         @simd for i in 1:l_I
            if I[i]==J[i]
                @inbounds b[IndVec[I[i],J[i]]]=V[i]
            elseif I[i]<J[i]
                @inbounds b[IndVec[I[i],J[i]]]=V[i]*sqrt2
            end
          end
     
         return b
end

        
function makevec_dense(B::SparseMatrixCSC{Float64},l_vec::UInt64,sqrt2::Float16,IndVec::SparseMatrixCSC{UInt64})
     b=zeros(Float64,l_vec)
     I,J,V=findnz(B)
     l_I=length(I)
     @simd for i in 1:l_I
            if I[i]==J[i]
                @inbounds b[IndVec[I[i],J[i]]]=V[i]
            else
                @inbounds b[IndVec[I[i],J[i]]]=V[i]*sqrt2
            end
      end

     return b
end         

    
    
  function makemat(b::SparseVector{Float64},l_vec::UInt64,sk_bar::Int64,sqrt2::Float16,reIndVec::Vector{Vector{UInt64}})
        B=spzeros(Float64,sk_bar,sk_bar)
        for i in 1:l_vec
            if reIndVec[i][1]==reIndVec[i][2]
                B[reIndVec[i][1],reIndVec[i][2]]=b[i]
            else
                B[reIndVec[i][1],reIndVec[i][2]]=b[i]/sqrt2
            end
            B[reIndVec[i][2],reIndVec[i][1]]=B[reIndVec[i][1],reIndVec[i][2]]
        end
    
        return B
     end

function makemat_dense(b::Vector{Float64},l_vec::UInt64,sk_bar::Int64,sqrt2::Float16,reIndVec::Vector{Vector{UInt64}})
        B=spzeros(Float64,sk_bar,sk_bar)
        for i in 1:l_vec
            if reIndVec[i][1]==reIndVec[i][2]
                B[reIndVec[i][1],reIndVec[i][2]]=b[i]
            else
                B[reIndVec[i][1],reIndVec[i][2]]=b[i]/sqrt2   
            end
            B[reIndVec[i][2],reIndVec[i][1]]=B[reIndVec[i][1],reIndVec[i][2]]
        end
        return B
     end
   



function BTP_POP(x::Vector{PolyVar{true}},f::Polynomial{true},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64,R::Float64;ak::Float64=0.0,EigAlg="Mix",tol=1e-5,showNormGrad=false,scale=false)

   println("**Bounded trace property (BTP):") 
    
  @time begin
           @time begin
    
    if showNormGrad
        global norm_grad=Vector{Float64}([])
    end
    n=length(x)
    l_g=length(g)
    l_h=length(h)
            
            
            
    v=get_basis(n,2*k)
    s2k=size(v,2)
    sort_v=sortslices(v,dims=2)
    re_ind=Vector{UInt64}(undef,s2k)
    @simd for j in 1:s2k
        @inbounds re_ind[bfind(sort_v,s2k,v[:,j],n)]=j
    end
    sk=binomial(k+n,n)      
    Order(alpha::Vector{UInt64})=re_ind[bfind(sort_v,s2k,alpha,n)]        
            
   
    ceil_g=[@inbounds ceil(Int64,maxdegree(g[i])/2) for i in 1:l_g]
    sk_g=[@inbounds binomial(k-ceil_g[i]+n,n) for i in 1:l_g]
    s2k_g=[@inbounds binomial(2*(k-ceil_g[i])+n,n) for i in 1:l_g]
    s2k_h=[@inbounds binomial(2*(k-ceil(Int64,maxdegree(h[j])/2))+n,n) for j in 1:l_h]
      
   
    
 
    bound_ak=sk
    lmon_gj,supp_gj,coe_gj=Int64(0),zeros(UInt64,1,1),zeros(Float64,1)
    
    
    @simd for j in 1:l_g
        @inbounds lmon_gj,supp_gj,coe_gj=info(g[j],x,n)
        bound_ak+=norm(coe_gj,1)*sk_g[j]
            end      
    bound_ak*=R^k        
            
    println("  Largest upper bound of psd matrix: bound_ak=",bound_ak)
        
   if ak==0 
      ak=bound_ak
   end
                
            
   

    s=sk+sum(sk_g)+1
    numBloc=l_g+2
    println("  Number of blocks: numBloc=",numBloc)
    println("  Size of block-diagonal matrix: s=",s)
    println("  Size of blocks: ",[sk;sk_g;1])
   
    A=SparseMatrixCSC{Float64}[]
    t_A=UInt64(1)
    
    
    M=spzeros(UInt64,sk,sk)
    IndM=[Vector{Int64}[] for j in 1:s2k]
    r=UInt64(0)

    for i in 1:sk, j in i:sk
        @inbounds r=Order(v[:,i]+v[:,j])
        @inbounds M[i,j]=r
        @inbounds append!(IndM[r],[[i,j]]) 
    end
    l_IndM=[length(IndM[r]) for r in 1:s2k]
    
    for r in 1:s2k
        if l_IndM[r]>1
            for i in 2:l_IndM[r]
                matA=spzeros(Float64,s,s)
                if IndM[r][1][1]==IndM[r][1][2]     
                    matA[IndM[r][1][1],IndM[r][1][2]]=1
                else
                    matA[IndM[r][1][1],IndM[r][1][2]]=.5
                end
                    
                if IndM[r][i][1]==IndM[r][i][2]     
                    matA[IndM[r][i][1],IndM[r][i][2]]=-1
                else
                    matA[IndM[r][i][1],IndM[r][i][2]]=-.5
                end
                
                append!(A,[matA])
           end
       end
    end
   
    
    
    
    IndMg=Vector{Vector{Vector{Vector{UInt64}}}}(undef,l_g)
    l_IndMg=Vector{Vector{UInt64}}(undef,l_g)
    t_Blo=sk
    for j in 1:l_g
        IndMg[j]=[[zeros(UInt64,2)] for r in 1:s2k_g[j]]
 
        r=UInt64(0)

        for p in 1:sk_g[j], q in p:sk_g[j]
            @inbounds r=Order(v[:,p]+v[:,q])
            if IndMg[j][r][1]==[0,0]
               IndMg[j][r][1]=[p,q]
            else
                @inbounds append!(IndMg[j][r],[[p,q]]) 
            end
        end
        l_IndMg[j]=[length(IndMg[j][t]) for t in 1:s2k_g[j]]
        for r in 1:s2k_g[j]
            if l_IndMg[j][r]>1
                for i in 2:l_IndMg[j][r]
                    matA=spzeros(Float64,s,s)
                    if IndMg[j][r][1][1]==IndMg[j][r][1][2]     
                        matA[IndMg[j][r][1][1]+t_Blo,IndMg[j][r][1][2]+t_Blo]=1
                    else
                        matA[IndMg[j][r][1][1]+t_Blo,IndMg[j][r][1][2]+t_Blo]=.5
                    end

                    if IndMg[j][r][i][1]==IndMg[j][r][i][2]
                        matA[IndMg[j][r][i][1]+t_Blo,IndMg[j][r][i][2]+t_Blo]=-1
                    else
                        matA[IndMg[j][r][i][1]+t_Blo,IndMg[j][r][i][2]+t_Blo]=-.5
                    end
                    
                    append!(A,[matA])
               end
           end
        end
        t_Blo+=sk_g[j]
     end
            
     
     
            
    
            
            
            
    
    lmon_gj,supp_gj,coe_gj,I=Int64(0),zeros(UInt64,1,1),zeros(Float64,1),zeros(UInt64,2)
    
    t_Blo=sk
    @simd for j in 1:l_g
        @inbounds lmon_gj,supp_gj,coe_gj=info(g[j],x,n)
        @simd for r in 1:s2k_g[j]
                  matA=spzeros(Float64,s,s)
                   @inbounds I=IndMg[j][Order(v[:,r])][1]
                  if I[1]==I[2]
                      @inbounds matA[I[1]+t_Blo,I[2]+t_Blo]=-1
                  else
                      @inbounds matA[I[1]+t_Blo,I[2]+t_Blo]=-.5
                  end
            @simd for p in 1:lmon_gj   
                      @inbounds I=IndM[Order(supp_gj[:,p]+v[:,r])][1]
                      if I[1]==I[2]
                          @inbounds matA[I[1],I[2]]=coe_gj[p]
                      else
                          @inbounds matA[I[1],I[2]]=.5*coe_gj[p]
                      end
                   end
                    @inbounds append!(A,[matA])
               end
        t_Blo+=sk_g[j]            
    end             
            
    
    
  
            
            
    lmon_hj,supp_hj,coe_hj,I=Int64(0),zeros(UInt64,1,1),zeros(Float64,1),zeros(UInt64,2)
    
    
    @simd for j in 1:l_h
        @inbounds lmon_hj,supp_hj,coe_hj=info(h[j],x,n)
        @simd for r in 1:s2k_h[j]
            matA=spzeros(Float64,s,s)
            @simd for p in 1:lmon_hj   
                      @inbounds I=IndM[Order(supp_hj[:,p]+v[:,r])][1]
                      if I[1]==I[2]
                          @inbounds matA[I[1],I[2]]=coe_hj[p]
                      else
                          @inbounds matA[I[1],I[2]]=.5*coe_hj[p]
                      end
                   end
                    @inbounds append!(A,[matA])
               end       
    end 
    
    
    matA=spzeros(Float64,s,s)
    matA[1,1]=1
    append!(A,[matA])
    m=length(A)
            
    println("  Number of trace equality constraints: m=",m) 
           

    l_vec=UInt64(.5*(sk*(sk+1)+sum(sk_g[i]*(sk_g[i]+1) for i in 1:l_g)))+1

            
            
            
    IndVec=spzeros(UInt64,s,s)
    reIndVec=Vector{Vector{UInt64}}(undef,l_vec)
    t_vec=UInt(1)
    for i in 1:sk, j in i:sk
        IndVec[i,j]=t_vec
        reIndVec[t_vec]=[i,j]
        t_vec+=1
    end
    t_Blo=sk
    for r in 1:l_g
        for i in 1:sk_g[r], j in i:sk_g[r]
            IndVec[i+t_Blo,j+t_Blo]=t_vec
            reIndVec[t_vec]=[i+t_Blo,j+t_Blo]
            t_vec+=1
        end
        t_Blo+=sk_g[r]
    end
    
    IndVec[t_Blo+1,t_Blo+1]=t_vec
    reIndVec[t_vec]=[t_Blo+1,t_Blo+1]
            
            
    sqrt2=Float16(sqrt(2)) 
            
    a=spzeros(Float64,l_vec,m)
    @simd for j in 1:m
        @inbounds a[:,j]=makevec(A[j],l_vec,sqrt2,IndVec)
    end    
    
    
    A=[makemat(a[:,j],l_vec,s,sqrt2,reIndVec) for j in 1:m]
            
    C=spzeros(Float64,s,s)
    lmon_f,supp_f,coe_f=info(f,x,n)
    @simd for p in 1:lmon_f   
      @inbounds I=IndM[Order(supp_f[:,p])][1]
      if I[1]==I[2]
          @inbounds C[I[1],I[2]]=-coe_f[p]
      else
          @inbounds C[I[1],I[2]]=-0.5*coe_f[p]  
      end
   end 
            
            
   c=makevec_dense(C,l_vec,sqrt2,IndVec)  
   if scale         
       a,c,norm_c,opnorm_a=rescale_dense(a,c,m)
    else
        norm_c,opnorm_a=1.0,1.0
    end
        end
     

    println("***LMBM solver:***")
        
       
     
    function my_eigs(matA2::SparseMatrixCSC{Float64})
        eig_val=Vector{Float64}(undef,numBloc)
        eig_vec=Vector{SparseVector{Float64}}(undef,numBloc)
        eigval,eigvec=LargEig(Matrix(matA2[1:sk,1:sk]),sk,EigAlg=EigAlg)
        eig_val[1]=eigval
        eig_vec[1]=spzeros(Float64,s)
        eig_vec[1][1:sk]=eigvec
            
        t_Blo=sk
        @simd for j in 1:l_g
            @inbounds eigval,eigvec=LargEig(Matrix(matA2[1+t_Blo:sk_g[j]+t_Blo,1+t_Blo:sk_g[j]+t_Blo]),sk_g[j],EigAlg=EigAlg)
            @inbounds eig_val[j+1]=eigval
            @inbounds eig_vec[j+1]=spzeros(Float64,s)
            @inbounds eig_vec[j+1][1+t_Blo:sk_g[j]+t_Blo]=eigvec
            t_Blo+=sk_g[j]
        end

       eig_val[numBloc]=matA2[t_Blo+1,t_Blo+1]
       eig_vec[numBloc]=spzeros(Float64,s)
       eig_vec[numBloc][s]=1
       
       eig_val_star,i_star=findmax(eig_val)
       return eig_val_star,eig_vec[i_star]
    end

        
        
    function phi(nvar::Cint,xp::Ptr{Cdouble},gp::Ptr{Cdouble})
        zvar=unsafe_wrap(Array, xp, (convert(Int, nvar),))
        grad=unsafe_wrap(Array, gp, (convert(Int, nvar),))  
        eigval,eigvec=my_eigs(makemat_dense(c-a*zvar,l_vec,s,sqrt2,reIndVec))
        grad[:]=-a'*makevec2(eigvec*eigvec',l_vec,sqrt2,IndVec)
        grad[nvar]+=1/opnorm_a/ak
        
        if showNormGrad
            global norm_grad=[norm_grad;norm(grad)]
        end
        return(convert(Cdouble,eigval+zvar[nvar]/opnorm_a/ak))
    end

    opt_val,z= lmbm(phi,zeros(Float64,m);printinfo=true,tol=tol)
    opt_val=-opt_val*norm_c*ak



    println("####################################")
    println("opt_val = ",opt_val)
    println("####################################")
    
        
    G=makemat_dense(c-a*z,l_vec,s,sqrt2,reIndVec)
    Gr=Matrix(G[1:sk,1:sk])
    eigval,eigvec=LargEig(Gr,sk,EigAlg=EigAlg)
    Gr-=Matrix{Int64}(LinearAlgebra.I,sk,sk).*eigval
        
    opt_sol=extract_optimizer2(Gr,sk,v[:,1:sk],n,l_g,l_h,opt_val,f,g,h,x)  
    if showNormGrad
        println("norm_grad=",norm_grad)
    end    
    end
    return opt_val,opt_sol
end