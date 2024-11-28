import itertools
import sympy as sp

def unique_partitions(lst):
    """
    Generate unique unordered partitions of a list into two non-empty subsets.
    """
    partitions = []
    for i in range(1, len(lst)):
        for subset in itertools.combinations(lst, i):
            complement = tuple((set(lst) - set(subset)))
            yield (tuple((subset)), complement)
            #partitions.append( tuple([tuple(subset), tuple(complement)]) )
            #eturn partitions
            
def partitions(L_list):
    """
    Generate all possible partitions of L_list into two non-empty subsets.    
    
    Parameters:
    - L_list: List of boundary lengths
    Yields:
    - Tuple of two lists (L_I, L_J)
    """
    indices = list(range(len(L_list)))

    partitions = []
    for k in range(0, len(L_list)):
        for I_indices in itertools.combinations(indices, k):
            J_indices = [idx for idx in indices if idx not in I_indices]
            
            L_I = tuple( [L_list[idx] for idx in I_indices] )
            L_J = tuple( [L_list[idx] for idx in J_indices] )
            
            partitions.append( tuple([L_I, L_J]) )
            #yield (L_I, L_J)
    return tuple(partitions)
            
def m(alpha, L):
    """
    Computes symmetric monomial
    """
    n = len(L)
    k = len(alpha)
    # Pad alpha with zeros to match the length of L
    padded_alpha = alpha + [0] * (n - k)
    
    # Use a set to avoid duplicate permutations
    unique_permutations = set(itertools.permutations(padded_alpha))
    
    result = 0
    for perm in unique_permutations:
        # Compute the product for the current permutation
        product = 1
        for Li, beta_i in zip(L, perm):
            product *= Li**(2 * beta_i)  # L_i^(2*beta_i)
        result += product
    
    return result

def test_m():
    L1 = sp.Symbol("L1")
    L2 = sp.Symbol("L2")
    L3 = sp.Symbol("L3")
    L = [L1, L2, L3]

    alpha_1 = [3, 2, 1]
    alpha_2 = [2]
    alpha_3 = [1, 1, 1]
    alphas = [alpha_1, alpha_2, alpha_3]
    
    expected_1 = (L1**6*L2**4*L3**2 + L1**6*L2**2*L3**4 + L1**4*L2**6*L3**2 
                  + L1**4*L2**2*L3**6 + L1**2*L2**6*L3**4 + L1**2*L2**4*L3**6)
    expected_2 = L1**4 + L2**4 + L3**4
    expected_3 = L1**2*L2**2*L3**2
    
    expected = [expected_1, expected_2, expected_3]
    computed = [m(alpha, L) for alpha in alphas]
    
    for i in range(3):
        msg = f"Test case α={alphas[i]} failed"
        assert computed[i].equals(expected[i]), msg    
    
if __name__ == "__main__":
    test_m()
    
    
    elements = ("L2", "L3", "L4")
    result = partitions(elements)

    print(result)

    for I, J in result:
        print(f"I: {I}, J: {J}")
