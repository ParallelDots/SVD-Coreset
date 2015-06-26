module coreset
using StatsBase
export getCoreset 
function getCoreset(P::Array{Float64,2},beta::Int64,minFunctionSize::Int64,bicriteriaRatio::Float64,sampleSize::Int64)
	n,d = size(P)
	coreset = []
	while n>=minFunctionSize
       ind = rand(1:n,beta)
       ind = unique(ind)
       B = P[ind,:]
       QS = full(qrfact(B')[:Q])
       PQS = P * QS
       OPQS = P - PQS*QS'
       sqdists = sum(OPQS.*OPQS,2)
       sqdists = squeeze(sqdists,2)
       mIndexes = sortperm(sqdists)
       sort!(sqdists)
       m = Int64(ceil(n*bicriteriaRatio)) #No int convert earlier
       PQS = PQS[mIndexes[1:m],:]
       sqdists = sqdists[1:m]
       
       
       # deleteat!(P,mIndexes[1:m]) # Will not work
       donttakeindices = mIndexes[1:m]
       takeindices = []
       for i in 1:n
       if !(i in donttakeindices)
              push!(takeindices,i)
       end
       end
       takeindices = convert(Array{Int64,1},takeindices)
       #P = sub(P,takeindices,[1:d]) # Not taking sub as P created is small, also sub creates problem while matrix multiply, its just a view
       P = P[takeindices,1:d]
       n = size(P)[1]
       wv = WeightVec(sqdists)
       ind = convert(Array{Int64,1},sample(1:m,wv,sampleSize))
       ind = unique(ind)
       mIndexes = mIndexes[ind]
       d_ind = Dict()
       for i in ind
       d_ind[i] = get(d_ind,i,0)+1
       end
       h = map(x->d_ind[x],ind)
       sqrtw = sqrt(sum(sqdists) / sampleSize ./ sqdists[ind] .* h)
       C1 = sqrtw .* OPQS[mIndexes, :]
       U = full(qrfact(PQS)[:Q])
       sqdists = sum(U .* U, 2)
       ind = convert(Array{Int64,1},sample(1:m,wv,sampleSize))
       ind = unique(ind)
       d_ind = Dict()
       for i in ind
       d_ind[i] = get(d_ind,i,0)+1
       end
       h = map(x->d_ind[x],ind)
       sqrtw = sqrt(sum(sqdists) / sampleSize ./ sqdists[ind] .* h)
       L = PQS[ind, :] * QS'
       C2 = sqrtw .* L
       if coreset == []
       	coreset = vcat(C1,C2)
       else 
       	coreset = vcat(coreset,C1,C2)
       end
    end
    return coreset

end

end
#P = rand(10000,1000)
#beta = 5
#minFunctionSize = 10
#bicriteriaRatio = 0.5
#sampleSize = 50
#K = getCoreset(P,beta,minFunctionSize,bicriteriaRatio,sampleSize)
#println(K)
#println(size(K))
