#Womack prior -----
#K: numero de grupos considerados aparte del control (i.e., p-1 en el paper)
#podds: Rho (zeta)
# Entrega un vector con la probabilidad de los modelos de dicho nivel (nivel = predictores activos)
womack = function(K, podds)
{
out = rep( -Inf, K + 1 )
names(out) = 0 : K
out[ K + 1 ] = 0
	for( k in (K-1):0 )
		{
			j = (k+1):K
			bb = out[j+1] + lgamma(j+1) - lgamma(k+1) - lgamma(j+1-k) + log(podds)
			out[k+1] = log( sum( exp(bb - max(bb) )) ) + max(bb)
			out = out - max(out) -log(sum( exp(out - max(out) ) ) )		
		}
		exp( out + lbeta(c(0:K) + 1, K - c(0:K) + 1 ) + log(K+1) )
}

# la función por la cantidad de elementos en cada nivel (binomio de newton) deben sumar 1
#n = 1000
#sum(womack(n,1)*sapply(0:n, function(k) choose(n, k))) 
