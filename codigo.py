mport pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

np.random.seed(123456789)

def corre(data2):
    covarianza_calculada = np.zeros([data2.shape[1],data2.shape[1]])
   
    for i in range(np.shape(covarianza_calculada)[0]):
        covarianza_calculada[i,i] = 1
        a = data2[:,i] -np.mean(data2[:,i])
        for j in range(i+1,np.shape(covarianza_calculada)[1]):
            b = data2[:,j] - np.mean(data2[:,j])
            covarianza_calculada[i,j] = np.sum(a*b)/((data2.shape[0]-1)*np.var(a)**0.5*np.var(b)**0.5)
            covarianza_calculada[j,i] = covarianza_calculada[i,j]
           
    return covarianza_calculada

def filtrar_genes(lista_genes,n_porcentaje,b, correlacion):
    lista_nueva_genes = pd.DataFrame()
    diez_porciento_genes = int(a.shape[0]*n_porcentaje)
   
    min_correlacion = np.argsort(np.triu(correlacion)+2*np.tril(np.ones(np.shape(correlacion))),axis = 0)[:5,:]

       
    for i in range(min_correlacion.shape[0]):
        a1 = b.T[:,i]-np.mean(b.T[:,i])
        for j in min_correlacion[:,i]:
            b1 = b.T[:,j] - np.mean(b.T[:,j])
            dif_genes = np.argsort((a1*b1)/(np.var(a1)**0.5*np.var(b1)**0.5))
            lista_nueva_genes = pd.concat([lista_nueva_genes,lista_genes.iloc[dif_genes[:diez_porciento_genes]]])
            lista_nueva_genes = lista_nueva_genes.drop_duplicates()
           
    return lista_nueva_genes

def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df

a = pd.read_csv('PollenRaw.csv', index_col = 0)
mayor_0 = a>0
a = a[mayor_0.sum(axis =  1)>10]
a = quantileNormalize(a)

b = a.to_numpy() + 0.001
X_filter = b.T
X_filter=stats.zscore(X_filter,axis = 1)


correlacion = np.corrcoef(X_filter)
np.fill_diagonal(correlacion,0)
plt.figure()
plt.imshow(correlacion)

lista_genes = pd.DataFrame(a.index)

lista_nueva_genes = filtrar_genes(lista_genes, 0.8, X_filter, correlacion)              
a_final = a.iloc[lista_nueva_genes.index]
X_filter_final = a_final.to_numpy().T
X_filter_final = stats.zscore(X_filter_final,axis = 1)
correlacion2 = np.corrcoef(X_filter_final)
np.fill_diagonal(correlacion2,0)
plt.figure()
plt.imshow(correlacion2)





       
       
       


#
#correlation = np.corrcoef(b)
#plt.figure()
#plt.imshow(correlation)

#def win(y,wl,ll):
#    
#    r=np.convolve(y,np.ones([2*wl+1]))
#    r=r[wl+1:ll+wl+1]/(2*wl+1);
#    r[0:wl]=r[wl+1]
#    r[ll-wl+1:ll+1]=r[ll-wl]
#
#    return r
#
#def cef(y,ind,ranks,wl,ll):
#    
#    cey=win(y[ind],wl,ll);
#    return cey[ranks]
#
#def ace(X, o_it, i_it, i_crit, o_crit, wl):
#    
#    dim, ll = X.shape
#    dim = dim - 1
#    ind   = np.zeros([dim + 1, ll], dtype=np.int8)
#    ranks = np.zeros([dim + 1, ll], dtype=np.int8)
#    
#    for d in range(dim + 1):
#        ind[d,:] = np.argsort(X[d,:])
#        ranks[d, np.argsort(X[d,:])] = np.arange(ll)
#        
#    phi=(ranks-(ll-1)/2.)/ np.sqrt(ll*(ll-1)/12);
#    
#    i_eps=1
#    o_eps=1
#    o_it1=1
#    o_crit1=1
#    
#    while o_it1<=o_it and o_crit1>o_crit:
#        i_it1=1
#        i_crit1=1
#        while i_it1<=i_it and i_crit1>i_crit:
#            for d in range(1,dim):
#                sum0=0
#                for dd in range(1,dim):
#                    if dd !=d:
#                        sum0=sum0+phi[dd,:]
#                phi[d,:]=cef(phi[0,:]-sum0,ind[d,:],ranks[d,:],wl,ll)
#            
#            i_crit1=i_eps
#            
#            if dim==1:
#                sum0=phi[1,:]
#            else:
#                sum0=np.sum(phi[1:dim,:])
#                
#            i_eps=np.sum((sum0-phi[0,:]**2))/ll
#            i_crit1=abs(i_crit1-i_eps)
#            
#        phi[0,:]=cef(sum0,ind[0,:],ranks[0,:],wl,ll);
#        phi[0,:]=(phi[0,:]-np.mean(phi[0,:]))/np.std(phi[0,:])
#        
#        o_crit1=o_eps
#        o_eps=np.sum((sum0-phi[0,:])**2)/ll
#        o_crit1=abs(o_crit1-o_eps);
#        
#    psi=np.corrcoef(phi[0,:],sum0)
#    psi=psi[0,1];
#    
#    return psi
#    
#
#prueba = np.zeros([b.shape[0],b.shape[0]])
#
#for i in range(b.shape[0]):
#    for j in range(i+1,b.shape[0]):
#        prueba[i,j] = ace(np.array([b[i,:], b[j,:]]),50,10,1e-2,10*2.2204e-16,5)      
#        prueba[j,i] = prueba[i,j]        
#    
#    
#    
