function CTP_POP_on_Ball(x,f,g,h,k,R;EigAlg="Arpack",method="LMBM",tol=1e-5,showNormGrad=false,scale=false)
    @time begin
    l_g=length(g)
    l_h=length(h)
    n=length(x)
    @polyvar x_slack[1:l_g+1]
    R_bar=R
    invb=0.0
    sol=0.0
    for j in 1:l_g
        println("**Computing the upper bound ",j)
        invb,sol = CTP_POP([x;x_slack[1]],-g[j],[h;g[1]-x_slack[1]^2],k,R;method=method,EigAlg=EigAlg,tol=tol,scale=scale)
        R_bar-=invb
        println("--------------------------")
    end

    println("  Radius of big ball: R_bar=",R_bar)

    bar_h=[h;[g[j]-x_slack[j]^2 for j in 1:l_g];R_bar-sum([x;x_slack].^2)]
    
    opt_val,opt_sol = CTP_POP([x;x_slack],f,bar_h,k,R_bar;method=method,tol=tol,EigAlg=EigAlg,showNormGrad=showNormGrad)
    if opt_sol!=[]
        opt_sol=opt_sol[1:n]
    end
    end
    return opt_val,opt_sol
end




function CTP_POP_UniqueBallIneq(x,f,h,k,R;method="LMBM",EigAlg="Arpack",tol=1e-5,showNormGrad=false,scale=false)
    @time begin
    @polyvar x_slack
    opt_val,opt_sol = CTP_POP([x;x_slack],f,[h;R-sum(x.^2)-x_slack^2],k,R;method=method,EigAlg=EigAlg,tol=tol,showNormGrad=showNormGrad,scale=scale)
    if opt_sol!=[]
        opt_sol=opt_sol[1:length(x)]
    end
    end
    return opt_val,opt_sol
end





function CTP_POP(x::Vector{PolyVar{true}},f::Polynomial{true,Float64},h::Vector{Polynomial{true,Float64}},k::Int64,R::Float64;method="LMBM",EigAlg="Mix",tol=1e-5,showEvaluation=false,showNormGrad=false,scale=false)
    println("**Constant trace property (CTP):") 
    @time begin
        if showEvaluation
            global linear_oper=0
            global adjoint_oper=0
            global max_size=0
            global num_eig=0
        end
        if showNormGrad
            global norm_grad=Vector{Float64}([])
        end
        @time n,l,v,s,m,a0,a,Ib,Vb,invInde,invP,norm_a0,opnorm_a = ConvertStandardSDP(x,f,h,k,scale=scale)
        @time opt_val,sol=SpectralSDP(s,m,a0,a,Ib,Vb,invInde,(R+1)^k,norm_a0,opnorm_a,method=method,EigAlg=EigAlg,tol=tol,showNormGrad=showNormGrad,showEvaluation=showEvaluation)
     
        
        println("------------------------------------")
        println("**Numerical result:")
        println("====================================")
        println("opt_val=",opt_val)
        println("====================================")
        if method=="LMBM" || method=="PB"
            @time opt_sol=ExtractionOptSol(n,l,v,s,a0,a,invInde,sol,invP,opt_val,f,h,x,EigAlg=EigAlg,showEvaluation=showEvaluation)
        elseif method=="SketchyCGAL"
            invPmat=diagm(invP)
            @time opt_sol=extract_optimizer_moment_matrix(invPmat*sol*invPmat,s,v,n,l,opt_val,f,h,x)
        else
            opt_sol=Vector{Float64}([])
        end
        
        if showNormGrad
            println("----------------------------")
            println("norm_grad=",norm_grad)
            println("----------------------------")
        end
        if showEvaluation
            println("----------------------------")
            println("linear_oper=",linear_oper)
            println("adjoint_oper=",adjoint_oper)
            println("max_size=",max_size)
            println("num_eig=",num_eig)
            println("----------------------------")
        end
    end
    return opt_val,opt_sol
end

function ExtractionOptSol(n::Int64,l::Int64,v::Matrix{UInt64},s::Int64,a0::Vector{Float64},a::SparseMatrixCSC{Float64},invInde::SparseMatrixCSC{UInt64},z::Vector{Float64},invP::Vector{Float64},opt_val::Float64,f::Polynomial{true,Float64},h::Vector{Polynomial{true,Float64}},x::Vector{PolyVar{true}};EigAlg="Arpack",showEvaluation=false)
    
    P=diagm(invP.^-1)
    Gr=AdjOper(a0+a*z,invInde,s,showEvaluation=showEvaluation)
    eigval,eigvec=LargEig(Gr,s,EigAlg=EigAlg,showEvaluation=showEvaluation)
    for j in 1:s
        Gr[j,j]-=eigval
    end
    Gr=P*Gr*P
    return extract_optimizer(Gr,s,v,n,l,opt_val,f,h,x)
 
end    

function SpectralSDP(s::Int64,m::Int64,a0::Vector{Float64},a::SparseMatrixCSC{Float64},Ib::Vector{UInt64},Vb::Vector{Float64},invInde::SparseMatrixCSC{UInt64,Int64},CT::Float64,norm_a0,opnorm_a;method="LMBM",EigAlg="Arpack",tol=1e-5,showNormGrad=false,showEvaluation=false)
   
    if method=="LMBM"
        println("**LMBM solver:")


        function phi(nvar::Cint,xp::Ptr{Cdouble},gp::Ptr{Cdouble})
            zvar=unsafe_wrap(Array, xp, (convert(Int, nvar),))
            grad=unsafe_wrap(Array, gp, (convert(Int, nvar),))
            eigval,eigvec=LargEig(AdjOper(a0+a*zvar,invInde,s,showEvaluation=showEvaluation),s,EigAlg=EigAlg,showEvaluation=showEvaluation)
            @fastmath grad[:]=LinearOper(s,a,eigvec,showEvaluation=showEvaluation)
            grad[Ib]+=Vb/opnorm_a/CT
            if showNormGrad
                global norm_grad=[norm_grad;norm(grad)]
            end
            return(convert(Cdouble,eigval+(zvar[Ib]'*Vb)[1]/opnorm_a/CT))
        end

        optval, optsol= lmbm(phi,zeros(Float64,m);printinfo=true,tol=tol)
        return -optval*norm_a0*CT, optsol

    elseif method=="PB"
        println("**PB solver:")
        function evaluate_f(z::Vector{Float64})
            eigval,eigvec=LargEig(AdjOper(a0+a*z,invInde,s),s,EigAlg=EigAlg)
            subgrads=LinearOper(s,a,eigvec)
            subgrads[Ib]+=Vb/opnorm_a/CT
            return eigval+(z[Ib]'*Vb)[1]/opnorm_a/CT, subgrads
        end

        bundle = SpectralSOS.Model{ProximalMethod}(m, evaluate_f,tol)

        runb(bundle)

        optsol=getsolution(bundle)
        optval=getobjectivevalue(bundle)
        return -optval*norm_a0*CT, optsol
    elseif method=="SketchyCGAL"
        println("**SketchyCGAL solver:")
        C=AdjOper(-a0,invInde,s)
        Primitive1(x::Vector{Float64})= C*x
        Primitive2(y::Vector{Float64},x::Vector{Float64})= AdjOper(-a*y,invInde,s,showEvaluation=showEvaluation)*x
        Primitive3(x::Vector{Float64})= LinearOper(s,-a,x,showEvaluation=showEvaluation)

        cons = Vector{Float32}([1;1])
        b = zeros(Float64,m)

        b[Ib]=Vb/opnorm_a/CT

        R = UInt16(1) # rank/sketch size parameter
        maxit = UInt64(1e6) # limit on number of iterations

        optval, U, Delt = CGAL(UInt32(s),Primitive1,Primitive2,Primitive3, cons, b, R, maxit,STOPTOL=tol,showEvaluation=showEvaluation,EigAlg=EigAlg)
        optsol=U*Delt*U'
        
        return optval*norm_a0*CT, optsol
    else 
        println("No CTP-SDP method!")
    end
    
end



              

function AdjOper(a::Vector{Float64},invInde::SparseMatrixCSC{UInt64},s::Int64;showEvaluation=false)
    if showEvaluation
        global adjoint_oper+=1
    end
    B=zeros(Float64,s,s)
    for i in 1:s, j in i:s
        @inbounds B[i,j]=a[invInde[i,j]]
        @inbounds B[j,i]= copy(B[i,j])
    end
    
    return B
end
               
function LinearOper(s::Int64,a::SparseMatrixCSC{Float64},u::Vector{Float64};showEvaluation=false)
    if showEvaluation
        global linear_oper+=1
    end
    return a'*[@inbounds u[i]*u[j]*(2-0^(j-i)) for i in 1:s for j in i:s]      
end
    
  
function LargEig(mat::Matrix{Float64},s::Int64;EigAlg="Arpack",showEvaluation=false)
    if showEvaluation
        global num_eig+=1
        global max_size=maximum([max_size,s])
    end
    if EigAlg=="Arpack"
       E=eigs(mat,nev = 1,which=:LR) 
       return E[1][1],E[2][:,1]
    elseif EigAlg=="Normal"
       E=eigen(Symmetric(mat),s:s)
       return E.values[1],E.vectors[:,1]
    elseif EigAlg=="Mix"
       try 
           E=eigs(mat,nev = 1,which=:LR) 
           return E[1][1],E[2][:,1]
       catch
           E=eigen(Symmetric(mat),s:s)
           return E.values[1],E.vectors[:,1]
       end
    else
       println("No eigenvalue algorithm!!!")
    end  
    
end

    
function ConvertStandardSDP(x::Vector{PolyVar{true}},f::Polynomial{true},h::Vector{Polynomial{true,Float64}},k::Int64;scale=false)

    println("**Converting the moment relaxation to the standard SDP:")
    
    n=length(x)
    l_h=length(h)
            
            
            
    v=get_basis(n,2*k)
    s2k=size(v,2)
    sort_v=sortslices(v,dims=2)
    re_ind=Vector{Int64}(undef,s2k)
    @simd for j in 1:s2k
        @inbounds re_ind[bfind(sort_v,s2k,v[:,j],n)]=j
    end
    sk=binomial(k+n,n)      
    println("  Size of psd matrix: sk=",sk)
                    
    Order(alpha::Vector{UInt64})=re_ind[bfind(sort_v,s2k,alpha,n)]        
            
    s2k_h=[@inbounds binomial(2*(k-ceil(Int64,maxdegree(h[j])/2))+n,n) for j in 1:l_h]
    
    d=Int64(0.5*sk*(sk+1))
    m=d-s2k+1+sum(s2k_h)
    println("  Number of trace equality constraints: m=",m)
   
        
    
    Ib=Vector{UInt64}([m])
    Vb=Vector{Float64}([1])
        
    
    IndM=[Vector{Int64}[] for j in 1:s2k]
    diagM=Vector{Int64}(undef,sk)
    invInde=spzeros(UInt64,sk,sk)
    
    r=Int64(0)
    t=UInt64(1)

    for i in 1:sk, j in i:sk
        @inbounds r=Order(v[:,i]+v[:,j])
        @inbounds append!(IndM[r],[[i,j]])
        @inbounds invInde[i,j]=t
        t+=1
        if i==j
            @inbounds diagM[i]=r
        end
    end
    
    lmon,supp,coe=info((1+sum(x.^2))^k,x,n)
    invP=invdiaP([@inbounds Order(supp[:,j]) for j in 1:lmon],coe,sk,diagM)
    l_IndM=[length(IndM[r]) for r in 1:s2k]
    
    a=spzeros(Float64,d,m)
   
    t=UInt64(1)
    I=zeros(UInt64,2)
                
    for r in 1:s2k
        if l_IndM[r]>1
            for i in 2:l_IndM[r]
                I=IndM[r][1]
                if I[1]==I[2]     
                    a[invInde[I[1],I[2]],t]=invP[I[1]]^2
                else
                    a[invInde[I[1],I[2]],t]=0.5*invP[I[1]]*invP[I[2]]
                end
                
                I=IndM[r][i]
                if I[1]==I[2]     
                    a[invInde[I[1],I[2]],t]=-invP[I[1]]^2
                else
                    a[invInde[I[1],I[2]],t]=-0.5*invP[I[1]]*invP[I[2]]
                end
                t+=1
              
           end
           
            IndM[r]=Vector{Int64}[IndM[r][end]]
       end
    end
   
            
            
    I=zeros(UInt64,2)
    
    
    @simd for j in 1:l_h
        @inbounds lmon,supp,coe=info(h[j],x,n)
        @simd for r in 1:s2k_h[j]
            
            @simd for p in 1:lmon 
                @inbounds I=IndM[Order(supp[:,p]+v[:,r])][1]
                      if I[1]==I[2]
                          @inbounds a[invInde[I[1],I[2]],t]=coe[p]*invP[I[1]]^2
                      else
                          @inbounds a[invInde[I[1],I[2]],t]=0.5*coe[p]*invP[I[1]]*invP[I[2]]
                      end
                   end
                t+=1
               end       
    end 

   
    a[1,t]=-1
    t+=1
    
    a0=zeros(Float64,d)
    
    
    lmon,supp,coe=info(f,x,n)
    @simd for p in 1:lmon
        @inbounds I=IndM[Order(supp[:,p])][1]
      if I[1]==I[2]
          @inbounds a0[invInde[I[1],I[2]]]=-coe[p]*invP[I[1]]^2
      else
          @inbounds a0[invInde[I[1],I[2]]]=-0.5*coe[p]*invP[I[1]]*invP[I[2]]
      end
   end
                
#    for j=1:m
#         a[j,:]=a[j,:]./norm(a[j,:])
#     end
   
#    opnorm_a=svds(a, nsv = 1)[1].S[1]
#    a=a./opnorm_a
                
#    norm_a0=norm(a0)
#    a0=a0./norm_a0
               
#    println(size(a))
#    error()
   if  scale             
       a,a0,norm_a0,opnorm_a=rescale_dense(a,a0,m)
   else
       norm_a0,opnorm_a=1.0,1.0
   end
                
   
            
   return n,l_h,v[:,1:sk],sk,m,a0,a,Ib,Vb,invInde,invP,norm_a0,opnorm_a
end

            
function rescale_dense(a::SparseMatrixCSC{Float64},a0::Vector{Float64},m::Int64)
    
     
    
    a_ind1,a_ind2,a_val=findnz(a)
    a_len=length(a_val)

    norm_a=zeros(Float64,m)
    
    @fastmath @inbounds @simd  for t in 1:a_len
        norm_a[a_ind2[t]]+=a_val[t]^2
    end
        
    norm_a=sqrt.(norm_a)
   
    
   
    @fastmath @inbounds @simd  for t in 1:a_len
        a_val[t]/=norm_a[a_ind2[t]]
    end
    
    
    opnorm_a=svds(sparse(a_ind1,a_ind2,a_val), nsv = 1)[1].S[1]
    
    a_val=a_val./opnorm_a

    
    norm_a0=norm(a0)
    a0=a0./norm_a0
    
    
    return sparse(a_ind1,a_ind2,a_val,length(a0),m),a0,norm_a0,opnorm_a
end