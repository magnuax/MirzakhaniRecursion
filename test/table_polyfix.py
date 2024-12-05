import pickle
import sympy as sp


with open("data/mytable_poly.pkl", "rb") as file:
    V = pickle.load(file)
    

for g_key in V.keys():
    for n_key in V[g_key].keys():
        print(g_key, n_key)
        V[g_key][n_key] = sp.Poly(V[g_key][n_key])


with open("mytable_poly_fix.pkl", "wb") as file:
    pickle.dump(V, file)
    
    
    